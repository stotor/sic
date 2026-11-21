// sic.hpp - the types, shape functions, gather and push shared by every part of the code, and the functions each
// source file provides.
// Units. Input and output: plasma units, c = epsilon_0 = m_e = e = 1, lengths in c/w_p and times in 1/w_p of density
// n0 = 1. Inside the code lengths are in cells and times in dx/c (cell size dx in c/w_p): a species of density n has
// n dx^2 charge per cell volume per unit charge, E and B are dx times their plasma-unit values, J and rho dx^2 times
// (the species' charges: mesh.cu; the initial fields: fields.cu; the time step and every output: main.cu).
#pragma once
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace sic {
#ifdef SIC_SINGLE
using Real = float;
#else
using Real = double;
#endif
template<class T> using Device = thrust::device_vector<T>;
template<class T> T* ptr(Device<T>& a) { return thrust::raw_pointer_cast(a.data()); }
template<class T> const T* ptr(const Device<T>& a) { return thrust::raw_pointer_cast(a.data()); }
inline void check(cudaError_t e) { if (e != cudaSuccess) throw std::runtime_error(cudaGetErrorString(e)); }
constexpr int threads = 128;
inline int blocks(std::size_t n) { return int((n + threads - 1) / threads); }
inline int loop_blocks() {  // a fixed grid for kernels that loop over their work: 4 blocks per multiprocessor
    static int sms = 0;
    if (!sms) check(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, 0));
    return 4 * sms;
}

// Errors found on the GPU set bits of a status word; the host reads it at the diagnostics, at refinement, after the
// images, every 64 steps and at the end of the run.
enum Status : int { BAD_VALUE = 1, BAD_PIECES = 2, BAD_REFINE = 4, BAD_IMAGE = 8 };
inline void check_status(int status) {
    if (status & BAD_VALUE) throw std::runtime_error("a non-finite position or momentum");
    if (status & BAD_PIECES) throw std::runtime_error("an element needs more bisection levels than the deposit allows");
    if (status & BAD_REFINE) throw std::runtime_error("refinement: a non-finite edge or a vertex beyond a wall");
    if (status & BAD_IMAGE) throw std::runtime_error("an image value is not finite");
}

__host__ __device__ inline int wrap(int i, int n) { int r = i % n; return r < 0 ? r + n : r; }

struct Grid {  // what the kernels need to know about the run
    int dim = 1, n[3] = {16, 1, 1}, order = 1;
    bool wall[3] = {};      // conducting walls at node 0 of an axis (the periodic seam)
    bool galerkin = false;  // the gather: Galerkin (order N - 1 along a field's half-cell axes) or the same shape
    Real dt = Real(.1);
    __host__ __device__ int cells() const { return n[0] * n[1] * n[2]; }
    __host__ __device__ int index(int i, int j, int k) const {
        return wrap(i, n[0]) + n[0] * (wrap(j, n[1]) + n[1] * wrap(k, n[2]));
    }
};

// One cold species: a SIC sheet on the initial lattice, or PIC particles. Charge and mass in e and m_e; an infinite
// mass is immobile, part of the fixed neutralizing background (main.cu does not simulate it).
struct Species {
    std::string name;
    double charge = -1, mass = 1, density = 1;  // density n / n0 (a Fadeev sheet: its peak factor times the profile)
    double u[3] = {0, 0, 0};                    // drift momentum, gamma v / c
    double chain_y = -1;                        // Fadeev sheet: the y of its chain (c/w_p); < 0: uniform
    double perturb = 0;                         // (tests) u_d += perturb sin(2 pi x_d / L_d) along the active axes
};
struct Fadeev { double b0 = 1.6, lambda = 1, eps = .4, kick = -.5, guide = 0, bg = 0; int islands = 2; };  // (bg: main.cu)
struct Seed { double amplitude = 0; int modes = 0, kmax = 0; long long rng = 1; };  // random B modes (main.cu)
// A phase-space image's plane: axis[k] 0..2 the positions x, y, z (over the whole box), 3..5 the momenta ux, uy, uz
// (per mass, units of c; over [lo[k], hi[k]]).
struct PhaseView { int axis[2]; double lo[2], hi[2]; };

struct Options {  // the run as given on the command line (main.cu: parse)
    Grid g;                     // inside: cells and dx/c
    std::string problem;        // weibel or fadeev (empty: the species and fields a test sets)
    double dx = 1, L[3] = {16, 1, 1};  // cell size and box, c/w_p
    std::vector<Species> species;
    Fadeev fadeev;
    Seed seed;
    bool pic = false, curved = false, snapshot = true;
    int ppc = 0, vpc = 1, steps = 20, diag_every = 1, fields_every = 0, sort_every = 100;  // ppc 0: 4 (8 in 3D)
    Real refine = 0;                     // cells
    std::size_t max_vertices = 2000000;  // per species
    int render_every = 0, render_substeps = 1;  // density images (2D) and column images (3D): every N steps, K per step
    double render_scale = 10, render_window[4] = {0, 0, 0, 0};  // pixels per cell; u0 u1 v0 v1 (c/w_p; 0 0 0 0: box)
    int render_axis = 2;                        // 3D: the view axis (the image's axes u, v: the other two, in order)
    double render_slab[2] = {0, 0};             // 3D: the slab along the view axis (c/w_p; 0 0: the whole box)
    std::vector<std::string> render_species;   // (empty: every species)
    bool render_smooth = false;                // 2D: linear between the vertices' lumped densities (else exact, flat)
    int spectra_every = 0, dump_every = 0, phase_every = 0;  // energy spectra, mesh dumps, phase-space images
    std::vector<PhaseView> phase;              // the phase-space images' planes
    int phase_size[2] = {1600, 1000};          // their pixels (A, B)
    std::string output = "out";
};
// SIC elements or PIC particles at most (int indices; refinement records, up to six per tetrahedron, stay in an int).
constexpr std::size_t max_elements = 200000000;
inline void check_diagnostics(double gauss, double kinetic, double field, int step) {  // main.cu: record()
    if (!std::isfinite(gauss) || !std::isfinite(kinetic) || !std::isfinite(field))
        throw std::runtime_error("a non-finite diagnostic at step " + std::to_string(step));
}
// Reaching the vertex cap ends a run normally (main.cu): refinement leaves a conforming mesh before it grows past it.
struct VertexCap : std::runtime_error { using std::runtime_error::runtime_error; };
inline void check_growth(std::size_t vertices, std::size_t elements, std::size_t max_vertices) {  // refinement's caps
    if (vertices > max_vertices) throw VertexCap("refinement reached the vertex cap (--max-vertices)");
    if (elements > max_elements) throw std::runtime_error("refinement reached 2e8 elements");
}
// The degree-2 rule on a simplex with D + 1 points: point p puts weight a on corner p and (1 - a) / D on the others
// (kinetic; the profile's element charges in initial_mesh).
__host__ __device__ constexpr double rule_a(int D) {
    return D == 1 ? .788675134594812882 : D == 2 ? 2. / 3 : .585410196624968515;
}
// The mirror across conducting walls at 0 and n of an axis, for the gather and its transpose, the deposit's fold
// (deposit.cu: k_fold): even along a half-cell axis (j -> -j - 1 or 2n - j - 1), odd along a node axis (-j or 2n - j,
// flipping sign); the far wall's node n is node 0, the periodic seam.
template<class T> __host__ __device__ inline int mirror(int j, int n, bool half, T& sign) {
    if (half) return j < 0 ? -j - 1 : j >= n ? 2 * n - j - 1 : j;
    if (j < 0 || j > n) { j = j < 0 ? -j : 2 * n - j; sign = -sign; }
    return j == n ? 0 : j;
}

// A position is an integer cell plus an offset in [0,1), so a small move never meets a large float.
template<int D> struct Pos { int c[D]; Real f[D]; };
template<int D> __host__ __device__ inline void move(Pos<D>& p, int d, Real dx) {
    Real x = p.f[d] + dx;
    if (!(fabs(x) < Real(1e6))) { p.f[d] = Real(NAN); return; }  // a runaway or non-finite move: NaN, flagged later
    int s = int(floor(x));
    p.c[d] += s;
    p.f[d] = x - Real(s);
    if (p.f[d] >= 1) { p.c[d] += 1; p.f[d] = 0; }  // a tiny negative x rounded up to 1
}
template<int D> __host__ __device__ inline Real coordinate(const Pos<D>& p, int d, int origin) {
    return Real(p.c[d] - origin) + p.f[d];
}
template<int D> __host__ __device__ inline void place(Pos<D>& p, int d, double x) {  // x: within the int range
    const double c = floor(x);
    p.c[d] = int(c);
    p.f[d] = Real(x - c);
    if (p.f[d] >= 1) { p.c[d] += 1; p.f[d] = 0; }  // (an offset that rounds to 1)
}

// A vertex: position at n+1 and n, momentum at n+1/2; a PIC particle's charge (SIC: 0, the elements carry it; the
// field sits in the struct's alignment padding). An element: D+1 corner vertices, the periodic image of each corner
// (box lengths added to its position), its charge; a curved 1D element also has a midpoint vertex (Mesh::mid) in the
// frame of corner 0. (16-byte alignment lets the compiler load them 16 bytes per instruction instead of 4.)
template<int D> struct alignas(16) Vertex { Pos<D> x, old; Real u[3], q; };
template<int D> struct alignas(16) Element { int v[D + 1]; signed char img[D + 1][D]; Real q; };
template<int D> struct Mesh {  // one species
    Device<Vertex<D>> vertex;    // SIC vertices or PIC particles
    Device<Element<D>> element;  // SIC elements (none for PIC)
    Device<int> mid;             // curved 1D: each element's midpoint vertex
    Real qm = -1;                // charge-to-mass ratio, e / m_e
    std::size_t vertices_at_sort = 0;  // the vertex count at the last sort_mesh
};

struct Fields {
    Device<Real> e[3], b[3];
    Device<double> j[3];  // the current of the last step, summed in double (step() zeroes it first)
    Device<double> rho;   // charge density (diagnostics)
    explicit Fields(const Grid& g) {
        for (int c = 0; c < 3; ++c) { e[c].resize(g.cells(), 0); b[c].resize(g.cells(), 0); j[c].resize(g.cells(), 0); }
        rho.resize(g.cells(), 0);
    }
};
struct FieldPtr { const Real* e[3]; const Real* b[3]; };
inline FieldPtr pointers(const Fields& f) {
    return {{ptr(f.e[0]), ptr(f.e[1]), ptr(f.e[2])}, {ptr(f.b[0]), ptr(f.b[1]), ptr(f.b[2])}};
}

struct Work {  // GPU scratch, kept between steps
    Device<int> status = Device<int>(2, 0);   // [0] status bits, [1] refinement records (read together)
    Device<int> list;                         // deposit: [0] count, [1] next to take, then the elements for pieces
    Device<unsigned long long> keys, flags;   // refinement records and their ranks
    Device<int> owner, todo;                  // the records' elements; the elements to vote in the next pass
    Device<double> pad[3];                    // walls: the deposit's padded grid (deposit.cu: with_walls)
    int read_status() const { return status[0]; }
};

// ---- B-splines ------------------------------------------------------------------------------------------------------
// Order N (0 NGP .. 4 quartic): on branch b = floor(x + shift(N)), t = x + shift(N) - b, node b - N/2 + k carries
// w_k(t) = sum_m spline(N, k, m) t^m, k = 0..N.
__host__ __device__ constexpr Real shift(int N) { return N % 2 ? Real(0) : Real(.5); }
constexpr double binomial(int n, int k) {
    double b = 1;
    for (int i = 1; i <= k; ++i) b = b * (n + 1 - i) / i;
    return b;
}
// w_k(t) = B(t + N - k), B the cardinal spline (1/N!) sum_j (-1)^j C(N+1, j) (x - j)_+^N. These formulas run only at
// compile time, in make_tables() below.
constexpr double spline(int N, int k, int m) {
    double s = 0, sign = 1, factorial = 1;
    for (int j = 0; j <= N - k; ++j, sign = -sign) {
        double power = 1;
        for (int e = 0; e < N - m; ++e) power *= N - k - j;
        s += sign * binomial(N + 1, j) * power;
    }
    for (int i = 2; i <= N; ++i) factorial *= i;
    return s * binomial(N, m) / factorial;
}
// A window spans two branches b, b+1 per axis: nodes b - N/2 .. b - N/2 + N + 1. With u = x + shift - (b + 1), the
// coordinate from the plane between the branches, node i's weight is S_i(u) = sum_q window(N, i, q) phi_q(u) with
// phi_q = u^q (q <= N) and phi_{N+1} = max(u, 0)^N: the polynomial w_i(u + 1) of branch b, plus the jump of the N-th
// derivative at the plane (B-splines are smooth up to N - 1).
constexpr double window(int N, int i, int q) {
    if (q == N + 1) return (i > 0 ? spline(N, i - 1, N) : 0) - (i <= N ? spline(N, i, N) : 0);
    if (i > N) return 0;
    double s = 0;
    for (int m = q; m <= N; ++m) s += spline(N, i, m) * binomial(m, q);
    return s;
}
// The coefficients in constant memory, as tables built at compile time: device code does not fold the formulas
// above, while a constant-memory entry at a fixed index can be an instruction's operand (no separate load).
// flux[N][i][q] = -sum_{i' <= i} window(N, i', q): the current through face i + 1/2 of a window is minus the prefix
// sum of the weights' change, and summing the rows exactly here gives exact zeros where the cumulative weights are
// constant (no atomics for them).
struct Tables { Real spline[5][5][5], window[5][6][6], flux[5][6][6]; };
constexpr Tables make_tables() {
    Tables t{};
    for (int N = 0; N <= 4; ++N)
        for (int q = 0; q <= N + 1; ++q) {
            double prefix = 0;
            for (int i = 0; i <= N + 1; ++i) {
                if (i <= N && q <= N) t.spline[N][i][q] = Real(spline(N, i, q));
                t.window[N][i][q] = Real(window(N, i, q));
                prefix -= window(N, i, q);
                t.flux[N][i][q] = Real(prefix);
            }
        }
    return t;
}
static __constant__ Tables tables = make_tables();
template<int N> __device__ inline void spline_weights(Real t, Real* w) {
#pragma unroll
    for (int k = 0; k <= N; ++k) {
        Real s = tables.spline[N][k][N];
#pragma unroll
        for (int m = N - 1; m >= 0; --m) s = s * t + tables.spline[N][k][m];
        w[k] = s;
    }
}
template<int N> __device__ inline void basis(Real u, Real* phi) {
    phi[0] = 1;
#pragma unroll
    for (int q = 1; q <= N; ++q) phi[q] = phi[q - 1] * u;
    phi[N + 1] = u >= 0 ? phi[N] : Real(0);
}
// The window's node weights at u (flux = true: minus their prefix sums, the faces' cumulative weights).
template<int N> __device__ inline void window_weights(Real u, Real* S, bool flux = false) {
    Real phi[N + 2];
    basis<N>(u, phi);
#pragma unroll
    for (int i = 0; i < N + 2; ++i) {
        Real s = 0;
#pragma unroll
        for (int q = 0; q < N + 2; ++q) s += (flux ? tables.flux[N][i][q] : tables.window[N][i][q]) * phi[q];
        S[i] = s;
    }
}

// ---- gather and push (SIC vertices and PIC particles) --------------------------------------------------------------
// The order-M weights at coordinate x (cells) into w (M + 1 values); returns the first node.
template<int M> __device__ int stencil(Real x, Real* w) {
    const Real y = x + shift(M);
    const int b = int(floor(y));
    spline_weights<M>(y - b, w);
    return b - M / 2;
}
// Yee staggering: E_d sits half a cell along d, B_d half a cell along the other two axes (H: bit d = half along d).
// The gather (--gather): the same shape, order N along every axis at the field's own position (the full code's
// SAMESHAPE, which every run of the SIC paper used), or Galerkin, order N - 1 along a half-cell axis (the transpose of
// Esirkepov's current there; NGP: 0 everywhere). Conducting walls mirror the stencil: even
// across the wall for a half-cell axis, odd (sign flip) for a node axis.
template<int D, int N, int H> __device__ Real gather(const Real* F, const Grid& g, const Pos<D>& p) {
    constexpr int L = N > 0 ? N - 1 : 0;  // the Galerkin order
    int node[3][N + 1] = {}, count[3] = {1, 1, 1};  // per axis: the stencil's grid offsets (index times stride)
    Real w[3][N + 1] = {{1}, {1}, {1}};
#pragma unroll
    for (int d = 0; d < D; ++d) {
        const bool half = H >> d & 1, galerkin = g.galerkin && half && N > 0;
        const Real x = p.f[d] - (half ? Real(.5) : Real(0));
        const int first = p.c[d] + (galerkin ? stencil<L>(x, w[d]) : stencil<N>(x, w[d]));
        const int n = g.n[d], stride = d == 0 ? 1 : d == 1 ? g.n[0] : g.n[0] * g.n[1];
        count[d] = galerkin ? L + 1 : N + 1;
        int i = wrap(first, n);
#pragma unroll
        for (int a = 0; a < N + 1; ++a, i = i + 1 == n ? 0 : i + 1) {
            int j = i;
            if (g.wall[d]) j = mirror(first + a, n, half, w[d][a]);
            node[d][a] = j * stride;
        }
    }
    Real value = 0;
#pragma unroll
    for (int c = 0; c < (D > 2 ? N + 1 : 1); ++c)
#pragma unroll
        for (int b = 0; b < (D > 1 ? N + 1 : 1); ++b)
#pragma unroll
            for (int a = 0; a < N + 1; ++a)
                if (a < count[0] && b < count[1] && c < count[2])
                    value += w[0][a] * w[1][b] * w[2][c] * F[node[0][a] + node[1][b] + node[2][c]];
    return value;
}
// Relativistic Boris push and move (qm: charge-to-mass ratio); reflecting walls. Returns false for a non-finite result.
template<int D, int N> __device__ bool push(Vertex<D>& v, const Grid& g, const FieldPtr& f, Real qm) {
    v.old = v.x;
    const Real e[3] = {gather<D, N, 1>(f.e[0], g, v.x), gather<D, N, 2>(f.e[1], g, v.x),
                       gather<D, N, 4>(f.e[2], g, v.x)};
    const Real b[3] = {gather<D, N, 6>(f.b[0], g, v.x), gather<D, N, 5>(f.b[1], g, v.x),
                       gather<D, N, 3>(f.b[2], g, v.x)};
    const Real h = qm * g.dt / 2;
    Real* u = v.u;
    for (int d = 0; d < 3; ++d) u[d] += h * e[d];
    Real gamma = sqrt(1 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    Real t[3], t2 = 0;
    for (int d = 0; d < 3; ++d) { t[d] = h * b[d] / gamma; t2 += t[d] * t[d]; }
    const Real p[3] = {u[0] + u[1] * t[2] - u[2] * t[1], u[1] + u[2] * t[0] - u[0] * t[2],
                       u[2] + u[0] * t[1] - u[1] * t[0]};
    const Real s = 2 / (1 + t2);
    u[0] += s * (p[1] * t[2] - p[2] * t[1]) + h * e[0];
    u[1] += s * (p[2] * t[0] - p[0] * t[2]) + h * e[1];
    u[2] += s * (p[0] * t[1] - p[1] * t[0]) + h * e[2];
    gamma = sqrt(1 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    bool ok = isfinite(gamma);
#pragma unroll
    for (int d = 0; d < D; ++d) {
        move(v.x, d, g.dt * u[d] / gamma);
        if (g.wall[d]) {
            int n = g.n[d], c = v.x.c[d];
            Real f = v.x.f[d];
            if (c < 0 || c > n || (c == n && f > 0)) {  // reflect
                v.x.c[d] = c < 0 ? -c : 2 * n - c;
                v.x.f[d] = 0;
                move(v.x, d, -f);
                u[d] = -u[d];
            }
        }
        ok = ok && isfinite(v.x.f[d]);
    }
    return ok;
}

// ---- curved 1D elements ---------------------------------------------------------------------------------------------
// The quadratic through x0 (s = 0), x1 (s = 1) and xm (s = 1/2), its turning point (-1 if none), its root in [lo, hi].
template<class T> __host__ __device__ inline T quadratic(T x0, T x1, T xm, T s) {
    return (1 - s) * (1 - 2 * s) * x0 + s * (2 * s - 1) * x1 + 4 * s * (1 - s) * xm;
}
__host__ __device__ inline Real turning(Real x0, Real x1, Real xm) {
    const Real B = 4 * (xm - x0) - (x1 - x0), C = 2 * ((x1 - x0) - 2 * (xm - x0));
    return C != 0 ? -B / (2 * C) : Real(-1);
}
__host__ __device__ inline Real crossing(Real x0, Real x1, Real xm, Real lo, Real hi) {  // the root of x(s) in [lo, hi]
    const Real B = 4 * (xm - x0) - (x1 - x0), C = 2 * ((x1 - x0) - 2 * (xm - x0));
    if (C == 0) return fmin(hi, fmax(lo, -x0 / B));
    const Real disc = fmax(Real(0), B * B - 4 * C * x0), r = -Real(.5) * (B + copysign(sqrt(disc), B));
    const Real r1 = r / C, r2 = r != 0 ? x0 / r : r1, mid = (lo + hi) / 2;
    return fmin(hi, fmax(lo, fabs(r1 - mid) < fabs(r2 - mid) ? r1 : r2));
}

// ---- the functions of each source file (D = 1, 2, 3) ----------------------------------------------------------------
// mesh.cu
template<int D> Mesh<D> initial_mesh(const Options&, const Species&);
template<int D> void sort_mesh(Mesh<D>&, const Grid&);
template<int D> std::size_t refine(Mesh<D>&, const Options&, Work&);  // returns the vertices added
// deposit.cu
template<int D> void deposit_current(const Mesh<D>&, const Grid&, bool curved, Fields&, Work&);
template<int D>  // (the same into three arrays)
void deposit_current(const Mesh<D>&, const Grid&, bool curved, double* const (&j)[3], Work&);
template<int D> void deposit_density(const Mesh<D>&, const Grid&, bool curved, Device<double>& rho, Work&);  // adds
template<int D> void push_deposit_particles(Mesh<D>&, const Grid&, Fields&, Work&);  // PIC: push and deposit
// 3D column images: the view axis w, the image's axes u, v (the other two, in order) and pixels, the slab along w.
struct View {
    int axis[2], w, n[2];  // the image's axes u, v; the view axis; the image's pixels along u, v
    double scale[2];       // pixels per cell along u, v
    double lo, hi;         // the slab along w (cells, 0 <= lo < hi <= the box)
};
// Adds a species' column of charge per pixel area (charge / pixel^2) along the slab at the sample points, its vertices
// at step - 1 + alpha (a mesh without elements: PIC particles). The sample point of pixel (i, k): (i + sample_u,
// k + sample_v) pixels from the box's corner, off the centre so that the initial lattices' faces miss it.
constexpr double sample_u = .5 + 1 / 7.1, sample_v = .5 + 1 / 13.3;
void column_image(const Mesh<3>&, const Grid&, const View&, double alpha, double* image, Work&);
// Phase-space images: adds a species' charge in the pixels of the plane (view.axis[0], view.axis[1]) to
// image[size[1]][size[0]]: a position axis spans the box in size[k] pixels, a momentum axis [lo[k], hi[k]) (pixel j:
// [lo + j w, lo + (j + 1) w), w = (hi - lo) / size[k]), and only material inside the momentum ranges counts. SIC in 1D
// and 2D (straight elements), PIC in any dimension.
template<int D>
void phase_image(const Mesh<D>&, const Grid&, const PhaseView&, const int (&size)[2], double* image, Work&);
// fields.cu
template<int D> void push_vertices(Mesh<D>&, const Grid&, const Fields&, Work&);
void advance_fields(Fields&, const Grid&);
double gauss_error(const Fields&, const Grid&, const Device<double>& background);
double field_energy(const Fields&, bool magnetic);  // (inside units)
// A species' kinetic energy, sum of (gamma - 1) m, and momentum, sum of m u (inside units; the energy by the degree-2
// rule, the momentum exact).
struct Kinetic { double energy, momentum[3]; };
template<int D> Kinetic kinetic(const Mesh<D>&, bool curved);
// A species' energy spectrum: q x `scale` in each of spectrum_bins log bins of gamma - 1 (spectrum_per_decade a decade,
// from 10^spectrum_log0; below: the first bin, above: the last), each element integrated with kinetic's rule.
constexpr int spectrum_bins = 300, spectrum_per_decade = 30, spectrum_log0 = -6;
template<int D> void spectrum(const Mesh<D>&, bool curved, double scale, double* bins);  // adds (GPU, spectrum_bins)
void initial_fields(const Options&, Fields&);  // the problem's B (Fadeev field, random seed modes), inside units

// Calls f(std::integral_constant<int, N>{}) for the run's shape order (orders 3 and 4 only in 1D): the run-time order
// becomes a compile-time constant, decltype(order)::value inside the (generic) lambda f, instantiated for every order.
template<int D, class F> void with_order(int order, F&& f) {
    if (order < 0 || order > (D == 1 ? 4 : 2)) throw std::runtime_error("unsupported shape order");
    switch (order) {
    case 0: f(std::integral_constant<int, 0>{}); break;
    case 1: f(std::integral_constant<int, 1>{}); break;
    case 2: f(std::integral_constant<int, 2>{}); break;
    case 3: if constexpr (D == 1) f(std::integral_constant<int, 3>{}); break;
    default: if constexpr (D == 1) f(std::integral_constant<int, 4>{});
    }
}

// One time step (main.cu, tests). SIC refines every species first, so that the vertex cap stops a run before anything
// moved; then each species pushes its vertices and deposits its current (PIC: both in one kernel), the fields advance,
// and a mesh is re-sorted every o.sort_every steps or when it grew by an eighth (sort_every 0: never).
template<int D> void step(int n, std::vector<Mesh<D>>& species, const Options& o, Fields& f, Work& w) {
    for (Mesh<D>& m : species) refine(m, o, w);
    for (int c = 0; c < 3; ++c) check(cudaMemsetAsync(ptr(f.j[c]), 0, f.j[c].size() * sizeof(double)));
    for (Mesh<D>& m : species) {
        if (o.pic) push_deposit_particles(m, o.g, f, w);
        else {
            push_vertices(m, o.g, f, w);
            deposit_current(m, o.g, o.curved, f, w);
        }
    }
    advance_fields(f, o.g);
    for (Mesh<D>& m : species) {
        const std::size_t sorted = m.vertices_at_sort;
        if (o.sort_every > 0 && (n % o.sort_every == 0 || m.vertex.size() > sorted + sorted / 8)) sort_mesh(m, o.g);
    }
}
}  // namespace sic
