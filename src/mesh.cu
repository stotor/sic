// mesh.cu - the initial lattice (SIC) or particles (PIC), the Morton sort, and refinement by edge bisection.
#include "sic.hpp"
#include <cub/util_allocator.cuh>
#include <thrust/copy.h>
#include <thrust/gather.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/scan.h>
#include <thrust/scatter.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <algorithm>
#include <cmath>

namespace sic {
namespace {
using Key = unsigned long long;
// Thrust's temporary storage and the sort's scratch from a cache (a cudaFree per sort or scan would synchronize the
// device every time): blocks of powers of two, kept up to a quarter of the device's memory (larger releases are freed;
// cub frees them all and retries if the device runs out of memory). No synchronization after an algorithm unless its
// result is needed on the host (copy_if's count).
struct Cached {
    using value_type = char;
    char* allocate(std::ptrdiff_t n) {
        void* p;
        check(pool().DeviceAllocate(&p, n));
        return static_cast<char*>(p);
    }
    void deallocate(char* p, std::size_t) { check(pool().DeviceFree(p)); }
    static cub::CachingDeviceAllocator& pool() {
        static cub::CachingDeviceAllocator cache(2u, 1u, cub::CachingDeviceAllocator::INVALID_BIN, [] {
            std::size_t free, total;
            check(cudaMemGetInfo(&free, &total));
            return total / 4;
        }());
        return cache;
    }
};
inline auto fast() { static Cached alloc; return thrust::cuda::par_nosync(alloc); }
inline auto counted() { static Cached alloc; return thrust::cuda::par(alloc); }

// ---- the initial state ------------------------------------------------------------------------------------------
// A Fadeev chain at y = yc: the density factor F = (1 - eps^2) / (cosh(y'/lam) + eps cos(x/lam))^2, which carries the
// chain's current J_z = (B0/lam) F (Fadeev et al. 1965).
__host__ __device__ inline double fadeev_F(double x, double y, double lam, double eps, double yc) {
    const double c = cosh((y - yc) / lam) + eps * cos(x / lam);
    return (1 - eps * eps) / (c * c);
}
// One species: SIC vertices on a lattice of spacing h = 1/vpc cells, at the centres of its cells along a periodic axis
// (the full code's lattice) and on the walls along a walled one; PIC particles on a regular sub-lattice of each cell.
struct Lattice {
    Grid g;
    int size[3], vpc, ppc, side;  // lattice vertices per axis, per cell and axis; PIC: particles per cell, per axis
    Real q, u[3], perturb, h;     // an element's (lattice volume) or a particle's charge before the profile; the flow
    double dx, lambda, eps, chain_y;  // the Fadeev profile (chain_y < 0: uniform), c/w_p
    bool pic;
    __device__ Real offset(int d) const { return g.wall[d] ? Real(0) : h / 2; }
    __device__ double profile(const double* x) const {  // at a point in cells
        return chain_y < 0 ? 1.0 : fadeev_F(x[0] * dx, x[1] * dx, lambda, eps, chain_y);
    }
};
// A cold flow: the drift, plus (tests) u_d += perturb sin(2 pi x_d / L_d) along the active axes; no normal momentum on
// a wall.
template<int D> __device__ Vertex<D> cold(Pos<D> x, const Lattice& s) {
    Vertex<D> v{x, x, {s.u[0], s.u[1], s.u[2]}, 0};
    for (int d = 0; d < D; ++d) {
        const double phase = (double(x.c[d]) + x.f[d]) / s.g.n[d];
        v.u[d] += s.perturb * Real(sin(6.283185307179586 * phase));
        if (s.g.wall[d] && x.f[d] == 0 && (x.c[d] == 0 || x.c[d] == s.g.n[d])) v.u[d] = 0;
    }
    return v;
}
template<int D> __global__ void k_vertices(Vertex<D>* out, int count, Lattice s) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= count) return;
    Pos<D> x;
    if (s.pic) {  // side^D particles on a regular sub-lattice of each cell
        int cell = k / s.ppc, p = k % s.ppc;
        for (int d = 0; d < D; ++d) {
            x.c[d] = cell % s.g.n[d];
            cell /= s.g.n[d];
            x.f[d] = (Real(p % s.side) + Real(.5)) / s.side;
            p /= s.side;
        }
        double at[3] = {0, 0, 0};
        for (int d = 0; d < D; ++d) at[d] = x.c[d] + double(x.f[d]);
        out[k] = cold<D>(x, s);
        out[k].q = Real(s.q * s.profile(at));
        return;
    }
    const int lattice = s.size[0] * s.size[1] * s.size[2];
    const bool midpoint = k >= lattice;  // curved 1D: the element midpoints follow the lattice
    int rest = midpoint ? k - lattice : k;
    for (int d = 0; d < D; ++d) {
        const int i = rest % s.size[d];
        rest /= s.size[d];
        x.c[d] = i / s.vpc;
        x.f[d] = 0;
        move(x, d, Real(i % s.vpc) * s.h + s.offset(d) + (midpoint ? s.h / 2 : Real(0)));
    }
    out[k] = cold<D>(x, s);
}
// The elements of the lattice cell whose lowest corner is lattice vertex lo: a segment (1D), two triangles (2D) or the
// six Kuhn tetrahedra (3D: one per order of the axes along the path from corner 000 to corner 111). A corner past the
// last vertex of a periodic axis is vertex 0 with image 1. A profile weights each element by its mean over the element
// (the degree-2 rule with D + 1 points, as kinetic).
template<int D> __global__ void k_cells(Element<D>* out, int* mid, int cells, Lattice s) {
    const int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= cells) return;
    int lo[3] = {0, 0, 0}, rest = c;
    for (int d = 0; d < D; ++d) {
        const int span = s.size[d] - s.g.wall[d];
        lo[d] = rest % span;
        rest /= span;
    }
    double at[D + 1][3] = {};  // the corners' positions (cells)
    auto corner = [&](Element<D>& e, int k, int bits) {  // bit d: one lattice step along d
        int index = 0;
        for (int d = D - 1; d >= 0; --d) {
            int i = lo[d] + (bits >> d & 1);
            at[k][d] = i * double(s.h) + s.offset(d);
            e.img[k][d] = i == s.size[d];
            if (i == s.size[d]) i = 0;
            index = index * s.size[d] + i;
        }
        e.v[k] = index;
    };
    constexpr int per_cell = D == 1 ? 1 : D == 2 ? 2 : 6;
    const double a = rule_a(D), b = (1 - a) / D;
    for (int t = 0; t < per_cell; ++t) {
        Element<D>& e = out[per_cell * c + t];
        corner(e, 0, 0);
        if constexpr (D == 1) corner(e, 1, 1);
        if constexpr (D == 2) { corner(e, 1, t ? 3 : 1); corner(e, 2, t ? 2 : 3); }
        if constexpr (D == 3) {
            const int order[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
            for (int k = 0, bits = 0; k < 3; ++k) corner(e, k + 1, bits |= 1 << order[t][k]);
        }
        double mean = 1;
        if (s.chain_y >= 0) {
            mean = 0;
            for (int p = 0; p <= D; ++p) {  // point p: weight a on corner p, b on the others
                double x[3] = {0, 0, 0};
                for (int k = 0; k <= D; ++k)
                    for (int d = 0; d < D; ++d) x[d] += (k == p ? a : b) * at[k][d];
                mean += s.profile(x) / (D + 1);
            }
        }
        e.q = Real(s.q * mean);
    }
    if (mid) mid[c] = s.size[0] + c;
}

// What the vote kernel needs to record refinement votes (vote, below).
struct Votes { Real limit; unsigned long long* keys; int *owner, capacity, *count, *status; };
Votes votes(Work& w, Real limit) {
    return {limit, ptr(w.keys), ptr(w.owner), int(w.keys.size()), ptr(w.status) + 1, ptr(w.status)};
}
// ---- refinement votes -----------------------------------------------------------------------------------------------
// An element votes for its longest edge if that is longer than the limit: records (key << 1 | 1, element) for that edge
// and (key << 1, element) for its other long edges (mesh.cu: an edge whose records all vote is bisected). A key holds
// the edge's lower vertex from bit key_lower, the upper one from bit key_upper (vertex indices < 2^28: --max-vertices),
// and below them the upper's image relative to the lower, image + 1 in two bits per axis. The length is the same bits
// from either end (the separation in the other direction is its exact negative). Curved 1D: the curve's length, folds
// included.
constexpr int key_lower = 34, key_upper = 6;
template<int D> __host__ __device__ Real separation(const Pos<D>& a, const Pos<D>& b, int d, int image, const Grid& g) {
    return Real(b.c[d] - a.c[d] + image * g.n[d]) + (b.f[d] - a.f[d]);
}
template<int D, bool CURVED>
__device__ inline void vote(const Vertex<D>* vx, const Element<D>& e, int mid, int element, const Grid& g, Votes v) {
    using Key = unsigned long long;
    constexpr int E = D * (D + 1) / 2;
    Pos<D> x[D + 1];
#pragma unroll
    for (int k = 0; k <= D; ++k) x[k] = vx[e.v[k]].x;
    Key key[E], best_key = 0;
    Real length[E], best = -1;  // the longest edge (ties: the lower key), in scalars: no array indexed by a variable
    bool ok = true;
    int longest = 0;
#pragma unroll
    for (int i = 0, n = 0; i < D; ++i)
#pragma unroll
        for (int j = i + 1; j <= D; ++j, ++n) {
            const bool flip = e.v[i] > e.v[j];
            key[n] = Key(flip ? e.v[j] : e.v[i]) << key_lower | Key(flip ? e.v[i] : e.v[j]) << key_upper;
            length[n] = 0;
#pragma unroll
            for (int d = 0; d < D; ++d) {
                const int image = e.img[j][d] - e.img[i][d];
                ok = ok && image >= -1 && image <= 1;
                key[n] |= Key((flip ? -image : image) + 1) << (2 * d);
                const Real s = separation(x[i], x[j], d, image, g);
                length[n] += s * s;
            }
            length[n] = sqrt(length[n]);
            if constexpr (CURVED) {
                const Real x1 = separation(x[0], x[1], 0, e.img[1][0] - e.img[0][0], g);
                const Real xm = separation(x[0], vx[mid].x, 0, 0, g), s = turning(0, x1, xm);
                const Real xs = quadratic(Real(0), x1, xm, s);
                length[n] = s > 0 && s < 1 ? fabs(xs) + fabs(x1 - xs) : fabs(x1);
            }
            ok = ok && isfinite(length[n]);
            if (length[n] > best || (length[n] == best && key[n] < best_key)) {
                best = length[n];
                best_key = key[n];
                longest = n;
            }
        }
    if (!ok) { atomicOr(v.status, BAD_REFINE); return; }
    if (!(best > v.limit)) return;
#pragma unroll
    for (int k = 0; k < E; ++k) {
        if (!(length[k] > v.limit)) continue;
        const int slot = atomicAdd(v.count, 1);
        if (slot < v.capacity) { v.keys[slot] = key[k] << 1 | Key(k == longest); v.owner[slot] = element; }
    }
}

// ---- Morton order -----------------------------------------------------------------------------------------------
template<int D> __device__ Key morton(const Pos<D>& p, const Grid& g) {
    Key key = 0;
    for (int d = 0; d < D; ++d) {
        const Key x = Key(wrap(p.c[d], g.n[d]));
        for (int b = 0; b < 63 / D; ++b) key |= (x >> b & 1) << (b * D + d);
    }
    return key;
}
template<int D> __global__ void k_vertex_keys(const Vertex<D>* v, int n, Grid g, Key* keys) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n) keys[k] = morton(v[k].x, g);
}
// Renumbers an element's corners after the vertex sort and makes its sort key: the Morton key of its corner 0, with its
// path class inserted above the low 3D bits (blocks of 8^D cells). The class is, per axis, how the element meets the
// window's plane (no plane: 0; else the corners on its upper side, 1..D), at its current positions: elements of one
// class take one path through the deposit, so a warp of them does not diverge (until they move across the planes;
// the six Kuhn tetrahedra of a lattice are six classes, times the sub-cell positions with vpc > 1).
template<int D>
__global__ void k_renumber(Element<D>* e, int* mid, int n, const int* moved_to, const Vertex<D>* v, Grid g, Key* keys) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    for (int c = 0; c <= D; ++c) e[k].v[c] = moved_to[e[k].v[c]];
    if (mid) mid[k] = moved_to[mid[k]];
    const Pos<D> a = v[e[k].v[0]].x;
    const Real sh = g.order % 2 ? Real(0) : Real(.5);  // shift(order): the window planes sit at integers - sh
    Key shape = 0;
    for (int d = 0; d < D; ++d) {
        Real x[D + 1], lo = Real(INFINITY), hi = -Real(INFINITY);
        for (int c = 0; c <= D; ++c) {
            x[c] = a.f[d] + separation(a, v[e[k].v[c]].x, d, e[k].img[c][d] - e[k].img[0][d], g);
            lo = fmin(lo, x[c]);
            hi = fmax(hi, x[c]);
        }
        const Real b = floor(lo + sh);  // the window's lower branch, as frame() finds it
        int path = 0;
        if (hi + sh >= b + 1)
            for (int c = 0; c <= D; ++c) path += x[c] + sh >= b + 1;
        shape = shape * (D + 1) + path;
    }
    const Key m = morton(a, g);
    keys[k] = (m >> 3 * D) << (3 * D + 6) | shape << 3 * D | (m & ((Key(1) << 3 * D) - 1));
}

// ---- refinement -------------------------------------------------------------------------------------------------
// Every element votes for its longest edge if that is longer than the limit, and an edge is bisected when every element
// around it voted for it. The mesh stays conforming, and the globally longest edge is always bisected, so passes repeat
// until no edge is too long. A new vertex interpolates the positions at both time levels, so the density at either
// level is unchanged, and the momentum. An edge is its two vertices (the lower first) and the image of the second
// relative to the first, packed in a key (sic.hpp: vote); a record is (key << 1 | vote, element), and sorting puts an
// edge's non-voting records first.
// Curved 1D: the length is the curve's (folds included), a split adds the points at s = 1/4 and 3/4 as the children's
// midpoints, and the old midpoint becomes their shared corner.
template<int D, bool CURVED>
__global__ void k_votes(const Vertex<D>* vx, const Element<D>* el, const int* mid, const int* list, int listed,
                        Grid g, Votes votes) {
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= listed) return;
    const int i = list ? list[t] : t;  // all elements, or the listed ones
    vote<D, CURVED>(vx, el[i], CURVED ? mid[i] : 0, i, g, votes);
}
__device__ inline int first_vote(const Key* keys, int r) {  // the first of the records equal to record r
    while (r > 0 && keys[r - 1] == keys[r]) --r;
    return r;
}
__device__ inline bool bisected(const Key* keys, int r0) {  // r0: an edge's first vote; no non-voting record before it
    return (keys[r0] & 1) && (r0 == 0 || keys[r0 - 1] >> 1 != keys[r0] >> 1);
}
struct Voted { __device__ bool operator()(Key k) const { return k & 1; } };
// Each element that splits adds one child; each bisected edge adds `created` vertices (low word) - scanned for ranks.
__device__ inline Key flag(const Key* keys, int r, int r0, int created) {
    return bisected(keys, r0) ? Key(1) << 32 | Key(r == r0 ? created : 0) : 0;
}
__global__ void k_flags(const Key* keys, int n, int created, Key* flags) {
    const int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r < n) flags[r] = flag(keys, r, first_vote(keys, r), created);
}
template<int D, bool CURVED>
__global__ void k_split(Vertex<D>* vx, Element<D>* el, int* mid, const Key* keys, const int* owner, const Key* scan,
                        int n, int nv, int ne, Grid g, int* status) {
    const int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n) return;
    constexpr int created = CURVED ? 2 : 1;
    const int r0 = first_vote(keys, r);
    if (!bisected(keys, r0)) return;
    const int vertex = nv + int((scan[r0] - flag(keys, r0, r0, created)) & 0xffffffffu);
    const int e = owner[r], child = ne + int((scan[r] - flag(keys, r, r0, created)) >> 32);
    const Element<D> parent = el[e];
    Element<D> a = parent, b = parent;
    a.q = b.q = parent.q / 2;
    if constexpr (CURVED) {  // children (corner 0, midpoint) and (midpoint, corner 1)
        const Vertex<D> p0 = vx[parent.v[0]], p1 = vx[parent.v[1]], pm = vx[mid[e]];
        const int image = parent.img[1][0] - parent.img[0][0];
        for (int h = 0; h < 2; ++h) {
            const Real s = h ? Real(.75) : Real(.25);
            Vertex<D> v = p0;
            for (int l = 0; l < 2; ++l) {  // both time levels, from corner 0 (see refine)
                const Pos<D>&x0 = l ? p0.old : p0.x, &x1 = l ? p1.old : p1.x, &xm = l ? pm.old : pm.x;
                const Real s1 = separation(x0, x1, 0, image, g), sm = separation(x0, xm, 0, 0, g);
                move(l ? v.old : v.x, 0, quadratic(Real(0), s1, sm, s));
            }
            for (int c = 0; c < 3; ++c) v.u[c] = quadratic(p0.u[c], p1.u[c], pm.u[c], s);
            // (the push reflects vertices that cross a wall during the step, not one created beyond it)
            const bool outside = v.x.c[0] < 0 || v.x.c[0] > g.n[0] || (v.x.c[0] == g.n[0] && v.x.f[0] > 0);
            if (!isfinite(v.x.f[0]) || (g.wall[0] && outside)) atomicOr(status, BAD_REFINE);
            vx[vertex + h] = v;
        }
        a.v[1] = b.v[0] = mid[e];
        a.img[1][0] = b.img[0][0] = parent.img[0][0];
        mid[e] = vertex;
        mid[child] = vertex + 1;
    } else {
        const Key edge = keys[r] >> 1;
        const int lower = int(edge >> key_lower);
        const int upper = int(edge >> key_upper & ((1u << (key_lower - key_upper)) - 1));
        int i = 0, j = 0;  // the corners of the edge in this element
        for (int k = 0; k <= D; ++k) {
            if (parent.v[k] == lower) i = k;
            if (parent.v[k] == upper) j = k;
        }
        if (r == r0) {
            const Vertex<D> A = vx[lower], B = vx[upper];
            Vertex<D> m = A;
            for (int d = 0; d < D; ++d) {
                const int image = int(edge >> (2 * d) & 3) - 1;
                move(m.x, d, separation(A.x, B.x, d, image, g) / 2);
                move(m.old, d, separation(A.old, B.old, d, image, g) / 2);
            }
            for (int c = 0; c < 3; ++c) m.u[c] = (A.u[c] + B.u[c]) / 2;
            vx[vertex] = m;
        }
        a.v[j] = b.v[i] = vertex;  // the midpoint, in the frame of the lower vertex
        for (int d = 0; d < D; ++d) a.img[j][d] = b.img[i][d] = parent.img[i][d];
    }
    el[e] = a;
    el[child] = b;
}
template<class T> void grow(Device<T>& v, std::size_t n) {
    if (v.capacity() < n) v.reserve(n + n / 4);
    v.resize(n);
}
// A temporary device array from the cache (a Device array would cudaMalloc and cudaFree, which synchronize). Every use,
// the copy back in permute included, is on the default stream: it may be released before the GPU is done with it,
// since its next user on that stream runs after. (The destructor ignores DeviceFree's status: it must not throw.)
template<class T> struct Scratch {
    T* p;
    std::size_t n;
    explicit Scratch(std::size_t n) : p(reinterpret_cast<T*>(Cached().allocate(std::max<std::size_t>(n, 1) * sizeof(T)))), n(n) {}
    ~Scratch() { (void)Cached::pool().DeviceFree(p); }
    Scratch(const Scratch&) = delete;
    Scratch& operator=(const Scratch&) = delete;
    thrust::device_ptr<T> begin() const { return thrust::device_pointer_cast(p); }
    thrust::device_ptr<T> end() const { return begin() + n; }
};
void sorting(Scratch<Key>& keys, Scratch<int>& order) {  // the permutation that sorts keys (sorted in place)
    thrust::sequence(fast(), order.begin(), order.end());
    thrust::sort_by_key(fast(), keys.begin(), keys.end(), order.begin());
}
template<class T> void permute(Device<T>& v, const Scratch<int>& order) {  // v[i] <- v[order[i]]
    Scratch<T> out(v.size());
    thrust::gather(fast(), order.begin(), order.end(), v.begin(), out.begin());
    check(cudaMemcpyAsync(ptr(v), out.p, v.size() * sizeof(T), cudaMemcpyDeviceToDevice));
}
}  // namespace

template<int D> Mesh<D> initial_mesh(const Options& o, const Species& sp) {
    const Grid& g = o.g;
    const double h = 1.0 / o.vpc, unit = sp.charge * sp.density * o.dx * o.dx;  // charge per cell volume (inside)
    Lattice s{g, {1, 1, 1}, o.vpc, o.ppc, 1, Real(0), {Real(sp.u[0]), Real(sp.u[1]), Real(sp.u[2])}, Real(sp.perturb),
              Real(h), o.dx, o.fadeev.lambda, o.fadeev.eps, sp.chain_y, o.pic};
    Mesh<D> m;
    m.qm = Real(sp.charge / sp.mass);
    if (o.pic) {
        s.side = int(std::lround(std::pow(double(o.ppc), 1.0 / D)));
        s.q = Real(unit / o.ppc);
        m.vertex.resize(std::size_t(g.cells()) * o.ppc);
    } else {
        s.q = Real(unit * (D == 1 ? h : D == 2 ? h * h / 2 : h * h * h / 6));
        std::size_t lattice = 1, cells = 1;
        for (int d = 0; d < D; ++d) {
            s.size[d] = g.n[d] * o.vpc + g.wall[d];
            lattice *= s.size[d];
            cells *= s.size[d] - g.wall[d];
        }
        m.vertex.resize(lattice + (o.curved ? cells : 0));
        m.element.resize(cells * (D == 1 ? 1 : D == 2 ? 2 : 6));
        if (o.curved) m.mid.resize(cells);
        k_cells<D><<<blocks(cells), threads>>>(ptr(m.element), o.curved ? ptr(m.mid) : nullptr, int(cells), s);
    }
    k_vertices<D><<<blocks(m.vertex.size()), threads>>>(ptr(m.vertex), int(m.vertex.size()), s);
    check(cudaGetLastError());
    return m;
}

// Vertices (or particles) and elements in the Morton order of their cells: neighbours in space, neighbours in memory.
template<int D> void sort_mesh(Mesh<D>& m, const Grid& g) {
    const int nv = int(m.vertex.size()), ne = int(m.element.size());
    m.vertices_at_sort = m.vertex.size();
    Scratch<Key> keys(nv);
    Scratch<int> order(nv);
    k_vertex_keys<D><<<blocks(nv), threads>>>(ptr(m.vertex), nv, g, keys.p);
    sorting(keys, order);
    permute(m.vertex, order);
    if (!ne) return;
    Scratch<int> moved_to(nv);  // the inverse permutation: each vertex's new index
    thrust::scatter(fast(), thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(nv), order.begin(),
                    moved_to.begin());
    Scratch<Key> element_keys(ne);
    Scratch<int> elements(ne);
    k_renumber<D><<<blocks(ne), threads>>>(ptr(m.element), m.mid.empty() ? nullptr : ptr(m.mid), ne, moved_to.p,
                                           ptr(m.vertex), g, element_keys.p);
    sorting(element_keys, elements);
    permute(m.element, elements);
    if (!m.mid.empty()) permute(m.mid, elements);
    check(cudaGetLastError());
}

// Bisects every edge longer than o.refine whose elements all voted for it, until none is left. New vertices interpolate
// both time levels and the momentum, so the subdivision is exact at either level (main.cu refines at the start of a
// step: the push then overwrites the old level).
template<int D> std::size_t refine(Mesh<D>& m, const Options& o, Work& w) {
    if (o.pic || !(o.refine > 0) || m.element.empty()) return 0;
    const std::size_t before = m.vertex.size();
    const bool curved = D == 1 && o.curved;
    int listed = -1;  // the first pass votes every element; a later one only those that can have changed (w.todo)
    for (int pass = 0; pass < 64; ++pass) {
        const int nv = int(m.vertex.size()), ne = int(m.element.size()), voters = listed < 0 ? ne : listed;
        const int* mid = m.mid.empty() ? nullptr : ptr(m.mid);
        // This pass's records; if they did not fit the buffer (the count is complete even then), again with room for all.
        int words[2];  // the status bits and the number of records, read together
        for (;;) {
            const Votes votes = sic::votes(w, o.refine);
            check(cudaMemsetAsync(votes.count, 0, sizeof(int)));
            auto kernel = k_votes<D, false>;
            if constexpr (D == 1) if (curved) kernel = k_votes<D, true>;
            kernel<<<blocks(voters), threads>>>(ptr(m.vertex), ptr(m.element), mid, listed < 0 ? nullptr : ptr(w.todo),
                                                voters, o.g, votes);
            check(cudaMemcpy(words, ptr(w.status), sizeof words, cudaMemcpyDeviceToHost));
            if (words[1] <= votes.capacity) break;  // else again, with more room for the records
            w.keys.resize(words[1] + words[1] / 2);
            w.owner.resize(words[1] + words[1] / 2);
        }
        const int n = words[1];
        check_status(words[0]);
        if (n == 0) return m.vertex.size() - before;
        thrust::sort_by_key(fast(), w.keys.begin(), w.keys.begin() + n, w.owner.begin());
        grow(w.flags, n);
        k_flags<<<blocks(n), threads>>>(ptr(w.keys), n, curved ? 2 : 1, ptr(w.flags));
        thrust::inclusive_scan(fast(), w.flags.begin(), w.flags.begin() + n, w.flags.begin());
        const Key total = w.flags[n - 1];
        const std::size_t added = total & 0xffffffffu, children = total >> 32;
        check_growth(nv + added, ne + children, o.max_vertices);
        grow(m.vertex, nv + added);
        grow(m.element, ne + children);
        if (curved) grow(m.mid, ne + children);
        auto split = k_split<D, false>;
        if constexpr (D == 1) if (curved) split = k_split<D, true>;
        split<<<blocks(n), threads>>>(ptr(m.vertex), ptr(m.element), curved ? ptr(m.mid) : nullptr, ptr(w.keys),
                                      ptr(w.owner), ptr(w.flags), n, nv, ne, o.g, ptr(w.status));
        check(cudaGetLastError());
        // Only these can vote differently now: the elements that voted (split, or waiting for a neighbour; one record
        // each) and the new children. The others' edges are unchanged and were within the limit.
        grow(w.todo, n + children);
        const auto end =
            thrust::copy_if(counted(), w.owner.begin(), w.owner.begin() + n, w.keys.begin(), w.todo.begin(), Voted());
        listed = int(end - w.todo.begin());
        thrust::sequence(fast(), w.todo.begin() + listed, w.todo.begin() + listed + children, ne);
        listed += int(children);
    }
    throw std::runtime_error("refinement did not finish in 64 passes");
}

#define SIC_INSTANTIATE(D)                                       \
    template Mesh<D> initial_mesh<D>(const Options&, const Species&); \
    template void sort_mesh<D>(Mesh<D>&, const Grid&);           \
    template std::size_t refine<D>(Mesh<D>&, const Options&, Work&);
SIC_INSTANTIATE(1)
SIC_INSTANTIATE(2)
SIC_INSTANTIATE(3)
}  // namespace sic
