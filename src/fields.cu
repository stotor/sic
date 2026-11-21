// fields.cu - the vertex push, the Yee field update and the diagnostics.
#include "sic.hpp"
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>
#include <algorithm>
#include <random>
#include <vector>

namespace sic {
namespace {
template<int D, int N>
__global__ void k_push(Vertex<D>* vertices, int count, Grid g, FieldPtr fields, Real qm, int* status) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= count) return;
    Vertex<D> v = vertices[k];
    if (!push<D, N>(v, g, fields, qm)) atomicOr(status, BAD_VALUE);
    vertices[k] = v;
}
// Yee grid: rho on the nodes, J_d and E_d half a cell along d, B_d half a cell along the other two axes. E and B at
// integer times: B -= dt/2 curl E; E += dt (curl B - J); B -= dt/2 curl E. Derivatives along inactive axes vanish.
struct Node {
    int i[3], up[3], down[3];  // grid indices along each axis; the neighbouring nodes (inactive axes: itself)
    __device__ Node(const Grid& g, int k) {
        i[0] = k % g.n[0];
        i[1] = k / g.n[0] % g.n[1];
        i[2] = k / (g.n[0] * g.n[1]);
        for (int d = 0, stride = 1; d < 3; stride *= g.n[d], ++d) {
            up[d] = i[d] + 1 == g.n[d] ? k - i[d] * stride : k + stride;
            down[d] = i[d] == 0 ? k + (g.n[d] - 1) * stride : k - stride;
        }
    }
};
struct FieldRW { Real *e[3], *b[3]; };
__global__ void k_magnetic(FieldRW f, Grid g) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= g.cells()) return;
    const Node n(g, k);
    const Real h = g.dt / 2;
    auto dE = [&](int c, int d) { return f.e[c][n.up[d]] - f.e[c][k]; };
    f.b[0][k] -= h * (dE(2, 1) - dE(1, 2));
    f.b[1][k] -= h * (dE(0, 2) - dE(2, 0));
    f.b[2][k] -= h * (dE(1, 0) - dE(0, 1));
}
// Also applies the conducting walls: zero tangential E on node 0 (normal B there then keeps its initial value: its
// curl E is made of tangential E on the wall).
__global__ void k_electric(FieldRW f, Grid g, double* jx, double* jy, double* jz) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= g.cells()) return;
    const Node n(g, k);
    auto dB = [&](int c, int d) { return f.b[c][k] - f.b[c][n.down[d]]; };
    double* j[3] = {jx, jy, jz};
    const Real curl[3] = {dB(2, 1) - dB(1, 2), dB(0, 2) - dB(2, 0), dB(1, 0) - dB(0, 1)};
    for (int c = 0; c < 3; ++c) f.e[c][k] += g.dt * (curl[c] - Real(j[c][k]));
    for (int d = 0; d < 3; ++d)
        if (g.wall[d] && n.i[d] == 0)
            for (int c = 0; c < 3; ++c)
                if (c != d) f.e[c][k] = 0;
}
// (gamma - 1) for a momentum u
__device__ inline Real excess(const Real* u) {
    const Real u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    return u2 / (sqrt(1 + u2) + 1);
}
}  // namespace

template<int D> void push_vertices(Mesh<D>& m, const Grid& g, const Fields& f, Work& w) {
    const int n = int(m.vertex.size());
    with_order<D>(g.order, [&](auto order) {
        k_push<D, decltype(order)::value><<<blocks(n), threads>>>(ptr(m.vertex), n, g, pointers(f), m.qm, ptr(w.status));
    });
    check(cudaGetLastError());
}
void advance_fields(Fields& f, const Grid& g) {
    FieldRW p{{ptr(f.e[0]), ptr(f.e[1]), ptr(f.e[2])}, {ptr(f.b[0]), ptr(f.b[1]), ptr(f.b[2])}};
    k_magnetic<<<blocks(g.cells()), threads>>>(p, g);
    k_electric<<<blocks(g.cells()), threads>>>(p, g, ptr(f.j[0]), ptr(f.j[1]), ptr(f.j[2]));
    k_magnetic<<<blocks(g.cells()), threads>>>(p, g);
    check(cudaGetLastError());
}
// max |div E - (rho - background)|, without the wall nodes (the conductor's surface charge is not evolved).
double gauss_error(const Fields& f, const Grid& g, const Device<double>& background) {
    const Real *ex = ptr(f.e[0]), *ey = ptr(f.e[1]), *ez = ptr(f.e[2]);
    const double *rho = ptr(f.rho), *bg = ptr(background);
    auto first = thrust::make_counting_iterator(0);
    return thrust::transform_reduce(first, first + g.cells(), [=] __device__(int k) -> double {
        const Node n(g, k);
        for (int d = 0; d < 3; ++d)
            if (g.wall[d] && n.i[d] == 0) return 0;
        const double div =
            double(ex[k] - ex[n.down[0]]) + double(ey[k] - ey[n.down[1]]) + double(ez[k] - ez[n.down[2]]);
        const double r = fabs(div - (rho[k] - bg[k]));
        return isfinite(r) ? r : INFINITY;
    }, 0.0, thrust::maximum<double>());
}
double field_energy(const Fields& f, bool magnetic) {  // sum of E^2 / 2 or B^2 / 2
    const Device<Real>* F = magnetic ? f.b : f.e;
    const Real *x = ptr(F[0]), *y = ptr(F[1]), *z = ptr(F[2]);
    auto first = thrust::make_counting_iterator(0);
    return thrust::transform_reduce(first, first + int(F[0].size()), [=] __device__(int k) -> double {
        const double a[3] = {x[k], y[k], z[k]};
        return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) / 2;
    }, 0.0, thrust::plus<double>());
}
// Kinetic energy, sum of (gamma - 1) m, and momentum, sum of m u, in one pass (they read the same vertices). An element
// averages gamma - 1 over quadrature points of its interpolated momentum (the degree-2 rule with D + 1 points; curved
// 1D: the quadratic through the midpoint): a diagnostic, the deposit is exact. Its momentum is q / qm times its mean u,
// the mean of its corners' (curved: (u0 + u1 + 4 u_mid) / 6, the quadratic's), exact. The energy of an element in
// Real, the sums in double (in double precision, double arithmetic on a GeForce GPU is 1/64 of the float rate: in a
// single-precision run this reduction cost twice the current deposit). The momenta are the leapfrog's, at the half
// step before the positions.
// The momentum at point p of element e's rule: weight a on corner p, b on the others.
template<int D>
__device__ inline void rule_u(const Vertex<D>* v, const Element<D>& e, int mid, bool curved, int p, Real (&u)[3]) {
    const Real a = Real(rule_a(D)), b = (1 - a) / D;
    for (int c = 0; c < 3; ++c) {
        u[c] = 0;
        if (D == 1 && curved) {  // the quadratic at s = the weight of corner 1
            const Real u0 = v[e.v[0]].u[c], u1 = v[e.v[1]].u[c], um = v[mid].u[c];
            u[c] = quadratic(u0, u1, um, p == 1 ? a : b);
        } else {
            for (int k = 0; k <= D; ++k) u[c] += (k == p ? a : b) * v[e.v[k]].u[c];
        }
    }
}
struct AddKinetic {
    __host__ __device__ Kinetic operator()(Kinetic a, const Kinetic& b) const {
        a.energy += b.energy;
        for (int c = 0; c < 3; ++c) a.momentum[c] += b.momentum[c];
        return a;
    }
};
template<int D> Kinetic kinetic(const Mesh<D>& m, bool curved) {
    const Vertex<D>* v = ptr(m.vertex);
    const Element<D>* el = ptr(m.element);
    const int* mid = ptr(m.mid);
    const double qm = m.qm;
    const bool pic = m.element.empty();
    auto first = thrust::make_counting_iterator(0);
    const int count = int(pic ? m.vertex.size() : m.element.size());
    return thrust::transform_reduce(first, first + count, [=] __device__(int i) -> Kinetic {
        Kinetic k{};
        if (pic) {
            const double mass = double(v[i].q) / qm;
            k.energy = mass * double(excess(v[i].u));
            for (int c = 0; c < 3; ++c) k.momentum[c] = mass * double(v[i].u[c]);
            return k;
        }
        const double mass = double(el[i].q) / qm;
        Real sum = 0;
        for (int p = 0; p <= D; ++p) {
            Real u[3];
            rule_u(v, el[i], curved ? mid[i] : 0, curved, p, u);
            sum += excess(u);
        }
        k.energy = mass * double(sum) / (D + 1);
        for (int c = 0; c < 3; ++c) {
            double u = 0;
            if (D == 1 && curved) u = (double(v[el[i].v[0]].u[c]) + v[el[i].v[1]].u[c] + 4.0 * v[mid[i]].u[c]) / 6;
            else {
                for (int j = 0; j <= D; ++j) u += v[el[i].v[j]].u[c];
                u /= D + 1;
            }
            k.momentum[c] = mass * u;
        }
        return k;
    }, Kinetic{}, AddKinetic());
}
// The spectrum: per block a histogram in shared memory, then added to the bins.
__device__ inline int spectrum_bin(Real e) {
    const int k = e > 0 ? int(floor((log10(double(e)) - spectrum_log0) * spectrum_per_decade)) : 0;
    return k < 0 ? 0 : k >= spectrum_bins ? spectrum_bins - 1 : k;
}
template<int D>
__global__ void k_spectrum(const Vertex<D>* v, const Element<D>* el, const int* mid, int count, bool curved,
                           double scale, double* bins) {
    __shared__ double h[spectrum_bins];
    for (int k = threadIdx.x; k < spectrum_bins; k += blockDim.x) h[k] = 0;
    __syncthreads();
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < count; i += gridDim.x * blockDim.x) {
        if (!el) { atomicAdd(h + spectrum_bin(excess(v[i].u)), double(v[i].q) * scale); continue; }  // PIC
        for (int p = 0; p <= D; ++p) {
            Real u[3];
            rule_u(v, el[i], curved ? mid[i] : 0, curved, p, u);
            atomicAdd(h + spectrum_bin(excess(u)), double(el[i].q) * scale / (D + 1));
        }
    }
    __syncthreads();
    for (int k = threadIdx.x; k < spectrum_bins; k += blockDim.x)
        if (h[k] != 0) atomicAdd(bins + k, h[k]);
}
template<int D> void spectrum(const Mesh<D>& m, bool curved, double scale, double* bins) {
    const bool pic = m.element.empty();
    const int count = int(pic ? m.vertex.size() : m.element.size());
    if (!count) return;
    k_spectrum<D><<<std::min(loop_blocks(), blocks(count)), threads>>>(ptr(m.vertex), pic ? nullptr : ptr(m.element),
                                                                       curved ? ptr(m.mid) : nullptr, count, curved,
                                                                       scale, bins);
    check(cudaGetLastError());
}

// ---- the problems' initial magnetic field ---------------------------------------------------------------------------
// Fadeev field: the discrete curl of A_z on the nodes (div B = 0 to round-off; in 3D, where y is periodic, up to the
// ~1e-6 B0 of the five-image sum at the y seam, ly >= 2 pi lambda), plus a uniform guide field B_z. 2D (P22, the
// full code's make_fadeev.py): one chain at Ly/2 between the walls, with the merger kick
// dA_z = kick B0 lam cos(x / 2 lam) cos(pi (y - Ly/2) / Ly), a field-only perturbation (its current is not carried).
// 3D (R3D, the full code's fadeev_Az): mirror chains at Ly/4 (+J_z) and 3 Ly/4 (-J_z), periodic in y.
namespace {
double ln_cosh(double z) { const double a = std::fabs(z); return a + std::log1p(std::exp(-2 * a)) - M_LN2; }
double fadeev_Az(const Options& o, double x, double y) {
    const Fadeev& p = o.fadeev;
    const double lam = p.lambda, B0 = p.b0, Ly = o.L[1], kick = p.kick * B0 * lam * std::cos(x / (2 * lam));
    if (o.g.dim == 2)
        return -B0 * lam * std::log(std::cosh((y - Ly / 2) / lam) + p.eps * std::cos(x / lam)) +
               kick * std::cos(M_PI * (y - Ly / 2) / Ly);
    // the y-only part makes B_x a double Harris profile; the x-dependent tails decay only as e^(-|y'|/lam): their
    // periodic images are summed (else B_y jumps at the y seam)
    const double y1 = Ly / 4, y2 = 3 * Ly / 4;
    auto tail = [&](double yp) { return std::log1p(p.eps * std::cos(x / lam) / std::cosh(yp / lam)); };
    double A = -B0 * lam * (ln_cosh((y - y1) / lam) - ln_cosh((y - y2) / lam)) + B0 * y;
    for (int n = -2; n <= 2; ++n) A -= B0 * lam * (tail(y - y1 - n * Ly) - tail(y - y2 - n * Ly));
    return A + kick * std::sin(2 * M_PI * y / Ly);
}
// Random seed modes, the full code's draws (the same field in every run with the same rng). 1D (sic_seed_bz): B_z
// modes m = 1..M along x with random phases. 2D and 3D (sic_seed_b): M modes with integer wave numbers in [-K, K]
// along the active axes (not all zero), random unit directions e of the vector potential a = A / |k| e cos(k.x + ph),
// so a mode's B = curl a has amplitude A |sin(e, k)|; B is the discrete curl of a on the cell edges (div B = 0).
void seed_modes(const Options& o, std::vector<double> (&b)[3]) {
    const Seed& s = o.seed;
    const int D = o.g.dim, n[3] = {o.g.n[0], o.g.n[1], o.g.n[2]};
    const double A = s.amplitude, h[3] = {o.L[0] / n[0], o.L[1] / n[1], o.L[2] / n[2]};
    if (A == 0 || s.modes <= 0) return;
    std::mt19937_64 rng((unsigned long long)s.rng);
    if (D == 1) {
        std::uniform_real_distribution<double> u(0, 2 * M_PI);
        std::vector<double> phase(s.modes + 1);
        for (int m = 1; m <= s.modes; ++m) phase[m] = u(rng);
        for (int i = 0; i < n[0]; ++i)
            for (int m = 1; m <= s.modes; ++m)
                b[2][i] += A * std::cos(2 * M_PI * m * ((i + .5) * h[0]) / o.L[0] + phase[m]);  // (half a cell up)
        return;
    }
    std::uniform_int_distribution<int> kd(-s.kmax, s.kmax);
    std::uniform_real_distribution<double> u(0, 1);
    struct Mode { double k[3], e[3], a, phase; };
    std::vector<Mode> modes(s.modes);
    for (Mode& m : modes) {
        int k[3] = {0, 0, 0};
        while (k[0] == 0 && k[1] == 0 && k[2] == 0)
            for (int d = 0; d < D; ++d) k[d] = kd(rng);
        double kk = 0;
        for (int d = 0; d < 3; ++d) {
            m.k[d] = 2 * M_PI * k[d] / o.L[d];
            kk += m.k[d] * m.k[d];
        }
        const double cz = 2 * u(rng) - 1, phi = 2 * M_PI * u(rng), sz = std::sqrt(std::max(0.0, 1 - cz * cz));
        m.e[0] = sz * std::cos(phi), m.e[1] = sz * std::sin(phi), m.e[2] = cz;
        m.a = A / std::sqrt(kk);
        m.phase = 2 * M_PI * u(rng);
    }
    const int cells = n[0] * n[1] * n[2];
    std::vector<double> a[3];
    for (auto& v : a) v.assign(cells, 0.0);
    for (int id = 0; id < cells; ++id) {  // a_c on the edge of cell id along c: half a cell along c
        const int q[3] = {id % n[0], id / n[0] % n[1], id / (n[0] * n[1])};
        for (int c = 0; c < 3; ++c) {
            double x[3] = {q[0] * h[0], q[1] * h[1], q[2] * h[2]}, sum = 0;
            x[c] += .5 * h[c];
            for (const Mode& m : modes) sum += m.a * m.e[c] * std::cos(m.k[0] * x[0] + m.k[1] * x[1] + m.k[2] * x[2] + m.phase);
            a[c][id] = sum;
        }
    }
    for (int id = 0; id < cells; ++id) {  // B_x at (i, j+1/2, k+1/2), B_y at (i+1/2, j, k+1/2), B_z at (i+1/2, j+1/2, k)
        const int q[3] = {id % n[0], id / n[0] % n[1], id / (n[0] * n[1])};
        int next[3];
        for (int d = 0, stride = 1; d < 3; stride *= n[d], ++d) next[d] = q[d] + 1 < n[d] ? id + stride : id - q[d] * stride;
        b[0][id] += (a[2][next[1]] - a[2][id]) / h[1] - (a[1][next[2]] - a[1][id]) / h[2];
        b[1][id] += (a[0][next[2]] - a[0][id]) / h[2] - (a[2][next[0]] - a[2][id]) / h[0];
        b[2][id] += (a[1][next[0]] - a[1][id]) / h[0] - (a[0][next[1]] - a[0][id]) / h[1];
    }
}
}  // namespace
void initial_fields(const Options& o, Fields& f) {
    const int nx = o.g.n[0], ny = o.g.n[1], nz = o.g.n[2], cells = o.g.cells();
    const double h = o.dx;
    std::vector<double> b[3];
    for (auto& v : b) v.assign(cells, 0.0);
    if (o.problem == "fadeev")
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {  // B_x at (i, j+1/2), B_y at (i+1/2, j); the y wall is node 0
                    const int id = i + nx * (j + ny * k);
                    const double x = i * h, y = j * h, A = fadeev_Az(o, x, y);
                    b[0][id] = (fadeev_Az(o, x, y + h) - A) / h;
                    b[1][id] = -(fadeev_Az(o, x + h, y) - A) / h;
                    b[2][id] = o.fadeev.guide * o.fadeev.b0;
                }
    seed_modes(o, b);
    for (int c = 0; c < 3; ++c) {
        std::vector<Real> inside(cells);
        for (int id = 0; id < cells; ++id) inside[id] = Real(b[c][id] * h);
        f.b[c] = inside;
    }
}

#define SIC_INSTANTIATE(D)                                                               \
    template void push_vertices<D>(Mesh<D>&, const Grid&, const Fields&, Work&);         \
    template Kinetic kinetic<D>(const Mesh<D>&, bool);                                   \
    template void spectrum<D>(const Mesh<D>&, bool, double, double*);
SIC_INSTANTIATE(1)
SIC_INSTANTIATE(2)
SIC_INSTANTIATE(3)
}  // namespace sic
