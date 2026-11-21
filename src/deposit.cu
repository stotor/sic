// deposit.cu - charge-conserving deposition. A SIC element carries its charge uniformly over its material (reference)
// simplex. Its density is the exact integral of the B-spline over the element, and its current is Esirkepov's, also
// integrated exactly, so Gauss's law holds to round-off. An element whose old and new corners lie within two branch
// cells per axis is integrated in registers by signed corners (below); any other is bisected into pieces that do.
// PIC particles are points: the same window and current assembly without the integrals.
// Template parameters: D dimension, N shape order, K corners (D + 1; a curved 1D element: 3), C output components
// (1: density, 3: current), MASK the axes whose window plane has cut a part, A the axis being cut. Loops over these
// are unrolled (#pragma unroll), so arrays indexed by them can stay in registers (3D TSC spills a little);
// `#pragma unroll 1` keeps a loop rolled.
#include "sic.hpp"

namespace sic {
namespace {
// The out-of-plane velocity components, 3 - D, as an array size: at least 1 (C++ has no empty arrays; 3D has none).
constexpr int velocities(int D) { return D < 3 ? 3 - D : 1; }
// The window: N + 2 nodes per axis, T in all; node t has digit d along axis d (base N + 2). No loops or recursion in
// these helpers: the compiler then folds them wherever the arguments are constants (unrolled loops).
template<int D, int N> struct Win {
    static constexpr int W = N + 2, T = D == 1 ? W : D == 2 ? W * W : W * W * W;
    __host__ __device__ static constexpr int stride(int d) { return d == 0 ? 1 : d == 1 ? W : W * W; }
    __host__ __device__ static constexpr int digit(int t, int d) { return t / stride(d) % W; }
};
__host__ __device__ constexpr int bits(int a) { return (a & 1) + (a >> 1 & 1) + (a >> 2 & 1); }  // (3 bits)
// Esirkepov's weights: the time average of a product of m <= 3 shape factors, k of them at the new time level,
// k! (m - k)! / (m + 1)!.
__device__ __forceinline__ Real esirkepov(int m, int k) {
    if (k == 0 || k == m) return m == 0 ? Real(1) : m == 1 ? Real(.5) : m == 2 ? Real(1. / 3) : Real(.25);
    return m == 2 ? Real(1. / 6) : Real(1. / 12);
}
struct Out { double* p[3]; };  // the current components, or the density in p[0]

// ---- an element in local coordinates and its window -------------------------------------------------------------
// x[l][k][d]: corner k's coordinate along axis d at time level l (0 old, 1 new), in cells from the integer cell o;
// v[k][c]: the velocity of out-of-plane component D + c. K corners: D + 1 for a simplex, 3 for a curved 1D element
// (its endpoints, then its midpoint).
template<int D, int K> struct Local { Real x[2][K][D], v[K][velocities(D)], q; int o[D]; };
template<int D, int K> constexpr bool curved = D == 1 && K == 3;

template<int D, int K>
__device__ __forceinline__ Local<D, K> load(const Vertex<D>* vx, const Element<D>& e, int mid, const Grid& g) {
    Local<D, K> L;
    L.q = e.q;
#pragma unroll
    for (int d = 0; d < D; ++d) L.o[d] = vx[e.v[0]].old.c[d] + e.img[0][d] * g.n[d];
#pragma unroll
    for (int k = 0; k < K; ++k) {
        const int c = k <= D ? k : 0;  // the midpoint is in the frame of corner 0
        const Vertex<D>& p = vx[k <= D ? e.v[k] : mid];
#pragma unroll
        for (int d = 0; d < D; ++d) {
            const int origin = L.o[d] - e.img[c][d] * g.n[d];
            L.x[0][k][d] = coordinate(p.old, d, origin);
            L.x[1][k][d] = coordinate(p.x, d, origin);
        }
        const Real inverse_gamma = 1 / sqrt(1 + p.u[0] * p.u[0] + p.u[1] * p.u[1] + p.u[2] * p.u[2]);
#pragma unroll
        for (int c = 0; c < 3 - D; ++c) L.v[k][c] = p.u[D + c] * inverse_gamma;
    }
    return L;
}
// Per axis: the lower branch b of the window (base: the grid index of its first node), the plane between branches b
// and b + 1 in local coordinates, whether the element reaches across it (planes), and whether it reaches beyond b + 1
// (misfit; fits: on no axis). The element's extent: its corners, and a curved element's extreme point.
template<int D> struct Frame { int base[D], planes, misfit; Real plane[D]; bool fits, finite; };
template<int D, int N, int K> __device__ __forceinline__ Frame<D> frame(const Local<D, K>& L, bool current) {
    Frame<D> f;
    f.planes = 0;
    f.misfit = 0;
    f.finite = true;
#pragma unroll
    for (int d = 0; d < D; ++d) {
        Real lo = Real(INFINITY), hi = -Real(INFINITY);
#pragma unroll
        for (int l = 0; l < 2; ++l) {
            if (l == 0 && !current) continue;
#pragma unroll
            for (int k = 0; k < K; ++k) {
                lo = fmin(lo, L.x[l][k][d]);
                hi = fmax(hi, L.x[l][k][d]);
                f.finite = f.finite && isfinite(L.x[l][k][d]);
            }
            if constexpr (curved<D, K>) {  // the extreme of a curved element
                const Real s = turning(L.x[l][0][0], L.x[l][1][0], L.x[l][2][0]);
                if (s > 0 && s < 1) {
                    const Real x = quadratic(L.x[l][0][0], L.x[l][1][0], L.x[l][2][0], s);
                    lo = fmin(lo, x);
                    hi = fmax(hi, x);
                }
            }
        }
        const Real b = floor(lo + shift(N));
        if (!(hi + shift(N) < b + 2)) f.misfit |= 1 << d;
        if (hi + shift(N) >= b + 1) f.planes |= 1 << d;
        f.plane[d] = b + 1 - shift(N);
        f.base[d] = f.finite ? L.o[d] + int(b) - N / 2 : 0;
    }
    f.fits = !f.misfit && f.finite;
    return f;
}

// ---- the integrals over a simplex: signed corners ------------------------------------------------------------------
// In the window, node i's weight along an axis is sum_q window(N, i, q) phi_q(u) (sic.hpp), so the integral of a
// node's weight over a simplex E is a combination of the basis integrals B[q] = int_E prod_d phi_{q_d}(u_d). Where
// q_d = N + 1 (the truncated power max(u_d, 0)^N), the integral runs over the part of E with u_d >= 0. Cutting by one
// plane gives that part as signed simplices: no corner on the positive side: nothing; one: the corner simplex at it;
// all: E; all but one: E minus the corner at the one outside; in 3D two: the wedge, as three tetrahedra. Cutting by
// the next axis splits each of these again; on every part the integrand is a polynomial, integrated exactly (the
// simplex moments below). The polynomial part of the shape is integrated over the whole simplex; only the kinks need
// the parts. Measures are material measures: a folded or collapsed element keeps its charge. A point particle is a
// simplex whose corners coincide.
template<int D> struct Simplex { Real u[D + 1][D], v[D + 1][velocities(D)], mu; };

template<int D> __device__ __forceinline__ void copy_corner(const Simplex<D>& S, int from, Simplex<D>& Q, int to) {
#pragma unroll
    for (int d = 0; d < D; ++d) Q.u[to][d] = S.u[from][d];
#pragma unroll
    for (int c = 0; c < 3 - D; ++c) Q.v[to][c] = S.v[from][c];
}
// Compare-exchange of corners i, j by coordinate A: the larger first. The sorting networks put the corners with
// u_A >= 0 first.
template<int D, int A> __device__ __forceinline__ void order_pair(Simplex<D>& S, int i, int j) {
    const bool swap = S.u[i][A] < S.u[j][A];
#pragma unroll
    for (int d = 0; d < D; ++d) {
        const Real a = S.u[i][d], b = S.u[j][d];
        S.u[i][d] = swap ? b : a;
        S.u[j][d] = swap ? a : b;
    }
#pragma unroll
    for (int c = 0; c < 3 - D; ++c) {
        const Real a = S.v[i][c], b = S.v[j][c];
        S.v[i][c] = swap ? b : a;
        S.v[j][c] = swap ? a : b;
    }
}
template<int D, int A> __device__ __forceinline__ void sort_corners(Simplex<D>& S) {
    if constexpr (D == 1) order_pair<D, A>(S, 0, 1);
    if constexpr (D == 2) { order_pair<D, A>(S, 0, 1); order_pair<D, A>(S, 1, 2); order_pair<D, A>(S, 0, 1); }
    if constexpr (D == 3) {
        order_pair<D, A>(S, 0, 1); order_pair<D, A>(S, 2, 3); order_pair<D, A>(S, 0, 2);
        order_pair<D, A>(S, 1, 3); order_pair<D, A>(S, 1, 2);
    }
}
// The point of edge i-j where u_A = 0 (set exactly) into corner k of Q; returns its fraction from corner i.
template<int D, int A>
__device__ __forceinline__ Real cut(const Simplex<D>& S, int i, int j, Simplex<D>& Q, int k) {
#ifdef SIC_SINGLE  // (the approximate division: 2 ulp, two instructions instead of the IEEE routine's twenty)
    const Real t = __fdividef(S.u[i][A], S.u[i][A] - S.u[j][A]);
#else
    const Real t = S.u[i][A] / (S.u[i][A] - S.u[j][A]);
#endif
#pragma unroll
    for (int d = 0; d < D; ++d) Q.u[k][d] = d == A ? Real(0) : S.u[i][d] + t * (S.u[j][d] - S.u[i][d]);
#pragma unroll
    for (int c = 0; c < 3 - D; ++c) Q.v[k][c] = S.v[i][c] + t * (S.v[j][c] - S.v[i][c]);
    return t;
}
template<int D, int A> __device__ __forceinline__ void corner(const Simplex<D>& S, int i, Simplex<D>& Q) {
    Q.mu = S.mu;
#pragma unroll
    for (int j = 0; j <= D; ++j) {
        if (j == i) copy_corner(S, i, Q, j);
        else Q.mu *= cut<D, A>(S, i, j, Q, j);
    }
}
template<int D> __device__ __forceinline__ int parts(int in) {  // the signed simplices of the part u_A >= 0
    return in == 0 ? 0 : in == 1 || in == D + 1 ? 1 : in == D ? 2 : 3;
}
// Simplex k of the part u_A >= 0 of S (sorted: its 'in' corners first) into Q; returns its sign.
template<int D, int A> __device__ __forceinline__ Real part(const Simplex<D>& S, int in, int k, Simplex<D>& Q) {
    if (in == D + 1 || (D > 1 && in == D && k == 0)) { Q = S; return 1; }
    if (in == 1) { corner<D, A>(S, 0, Q); return 1; }
    if (in == D) { corner<D, A>(S, D, Q); return -1; }
    if constexpr (D == 3) {  // the wedge of corners 0, 1 (in) and 2, 3 (out): a prism, as three tetrahedra
        if (k == 0) {
            copy_corner(S, 0, Q, 0);
            const Real t02 = cut<D, A>(S, 0, 2, Q, 1), t03 = cut<D, A>(S, 0, 3, Q, 2);
            copy_corner(S, 1, Q, 3);
            Q.mu = S.mu * t02 * t03;
        } else if (k == 1) {
            const Real t02 = cut<D, A>(S, 0, 2, Q, 0), t03 = cut<D, A>(S, 0, 3, Q, 1);
            copy_corner(S, 1, Q, 2);
            const Real t12 = cut<D, A>(S, 1, 2, Q, 3);
            Q.mu = S.mu * t03 * t12 * (1 - t02);
        } else {
            const Real t03 = cut<D, A>(S, 0, 3, Q, 0);
            copy_corner(S, 1, Q, 1);
            const Real t12 = cut<D, A>(S, 1, 2, Q, 2), t13 = cut<D, A>(S, 1, 3, Q, 3);
            Q.mu = S.mu * (1 - t03) * t12 * t13;
        }
    }
    return 1;
}

// Simplex moments, two ways. Orders 0 and 1: the mean of a product of up to three affine fields, in closed form.
// Orders >= 2: all the moments at once, from the simplex's generating function (below).
// The mean over a simplex (uniform material measure) of a product of up to three affine fields: factor f is the
// coordinate along axis f (f < D) or the velocity of out-of-plane component f - D; -1: none.
template<int D> __device__ __forceinline__ Real field(const Simplex<D>& P, int k, int f) {
    return f < 0 ? Real(1) : f < D ? P.u[k][f] : P.v[k][f - D];
}
template<int D> __device__ __forceinline__ Real mean(const Simplex<D>& P, int f0, int f1, int f2) {
    Real s0 = 0, s1 = 0, s2 = 0, s01 = 0, s02 = 0, s12 = 0, s012 = 0;
#pragma unroll
    for (int k = 0; k <= D; ++k) {
        const Real a = field(P, k, f0), b = field(P, k, f1), c = field(P, k, f2);
        s0 += a; s1 += b; s2 += c; s01 += a * b; s02 += a * c; s12 += b * c; s012 += a * b * c;
    }
    const int m = (f0 >= 0) + (f1 >= 0) + (f2 >= 0);
    constexpr Real r1 = Real(1. / (D + 1)), r2 = Real(r1 / (D + 2)), r3 = Real(r1 / ((D + 2) * (D + 3)));
    if (m == 0) return 1;
    if (m == 1) return s0 * r1;
    if (m == 2) return (s01 + s0 * s1) * r2;
    return (s0 * s1 * s2 + s0 * s12 + s1 * s02 + s2 * s01 + 2 * s012) * r3;
}
// All the moments at once. Over a simplex with corners y_k, the mean of y^a = prod_v y_v^a_v is D! a! / (D + |a|)!
// times the coefficient of s^a in prod_k 1 / (1 - y_k . s). Multiplying a truncated series G by one factor
// 1 / (1 - y . s) is the recurrence G[a] += sum_v y_v G[a - e_v], over the entries in increasing order (the lower
// entries are then already multiplied). The variables v: the D coordinates (powers up to N) and, for a current, the
// out-of-plane velocities (powers up to 1). The compiler drops the moments a part does not use.
template<int D, int N, int C> struct Series {  // entry r: the powers a_v = the digits of r (radix: the highest + 1)
    __host__ __device__ static constexpr int radix(int v) { return v < D ? N + 1 : C > 1 ? 2 : 1; }
    __host__ __device__ static constexpr int stride(int v) { return v ? radix(0) * (v == 2 ? radix(1) : 1) : 1; }
    __host__ __device__ static constexpr int digit(int r, int v) { return r / stride(v) % radix(v); }
    static constexpr int count = radix(0) * radix(1) * radix(2);
};
__host__ __device__ constexpr double factorial(int n) { return n < 2 ? 1 : n * factorial(n - 1); }
struct Scales { Real scale[3][5][3][3]; };  // [D - 1][a]: D! a! / (D + |a|)!
constexpr Scales make_scales() {
    Scales t{};
    for (int D = 1; D <= 3; ++D)
        for (int a = 0; a < 5; ++a)
            for (int b = 0; b < 3; ++b)
                for (int c = 0; c < 3; ++c)
                    t.scale[D - 1][a][b][c] =
                        Real(factorial(D) * factorial(a) * factorial(b) * factorial(c) / factorial(D + a + b + c));
    return t;
}
__constant__ Scales scales = make_scales();
template<int D, int N, int C>
__device__ __forceinline__ void moments(const Simplex<D>& P, Real (&G)[Series<D, N, C>::count]) {
    using Z = Series<D, N, C>;
    constexpr int R = Z::count;
    G[0] = 1;
#pragma unroll
    for (int r = 1; r < R; ++r) G[r] = 0;
#pragma unroll
    for (int k = 0; k <= D; ++k)  // G <- G / (1 - y_k . s)
#pragma unroll
        for (int r = 1; r < R; ++r)
#pragma unroll
            for (int v = 0; v < 3; ++v)
                if (Z::digit(r, v) > 0) G[r] += field(P, k, v) * G[r - Z::stride(v)];
#pragma unroll
    for (int r = 0; r < R; ++r) G[r] *= scales.scale[D - 1][Z::digit(r, 0)][Z::digit(r, 1)][Z::digit(r, 2)];
}
// Basis integrals summed over the mixed simplices of a step: m[t][q] += coef[t] * sign * int_P prod_d phi_{q_d}(u_d),
// component t < D (or the density, C = 1) unweighted, component t >= D weighted by the velocity of component t - D.
// For the current along axis t, the integrals with q_t = 0 do not involve that axis and cancel exactly between the
// time levels.
template<int D, int N, int C> struct Moments {
    using S = Win<D, N>;
    Real m[C][S::T], coef[C];
    // The entries of subset MASK: q_d = N + 1 on the axes in MASK, 0..N on the others (entry r: the digits of r in base
    // N + 1, on the free axes in order).
    // Component t needs entry q unless t is an axis and q does not involve it (q_t = 0: cancels between the levels).
    __host__ __device__ static constexpr bool used(int t, int q) { return C == 1 || t >= D || S::digit(q, t) != 0; }
    __host__ __device__ static constexpr int entries(int mask) {
        const int n = N + 1, f = D - bits(mask);
        return f == 0 ? 1 : f == 1 ? n : f == 2 ? n * n : n * n * n;
    }
    __host__ __device__ static constexpr int entry(int mask, int r) {
        const int n = N + 1, free0 = !(mask & 1), free1 = !(mask >> 1 & 1);
        const int d0 = free0 ? r % n : N + 1, r1 = free0 ? r / n : r;
        const int d1 = free1 ? r1 % n : N + 1, r2 = free1 ? r1 / n : r1;
        const int d2 = mask >> 2 & 1 ? N + 1 : r2 % n;
        return d0 + (D > 1 ? d1 * S::stride(1) : 0) + (D > 2 ? d2 * S::stride(2) : 0);
    }
    template<int MASK> __device__ __forceinline__ void add(const Simplex<D>& P, Real sign) {
        constexpr int count = entries(MASK);
        if constexpr (N <= 1) {
#pragma unroll
            for (int r = 0; r < count; ++r) {
                const int q = entry(MASK, r);
                int f[3] = {-1, -1, -1}, n = 0;
#pragma unroll
                for (int d = 0; d < D; ++d)
                    if (N == 1 && S::digit(q, d) > 0) f[n++] = d;
                const Real plain = sign * P.mu * mean(P, f[0], f[1], f[2]);
#pragma unroll
                for (int t = 0; t < C; ++t) {
                    if (!used(t, q)) continue;
                    if (C == 1 || t < D) { m[t][q] += coef[t] * plain; continue; }
                    int g[3] = {f[0], f[1], f[2]};
                    g[n] = t;
                    m[t][q] += coef[t] * sign * P.mu * mean(P, g[0], g[1], g[2]);
                }
            }
        } else {
            using Z = Series<D, N, C>;
            Real G[Z::count];
            moments<D, N, C>(P, G);
#pragma unroll
            for (int r = 0; r < count; ++r) {
                const int q = entry(MASK, r);
                int a = 0;  // the moment of q's monomial (u^N on the cut axes)
#pragma unroll
                for (int d = 0; d < D; ++d) a += (S::digit(q, d) > N ? N : S::digit(q, d)) * Z::stride(d);
#pragma unroll
                for (int t = 0; t < C; ++t)
                    if (used(t, q)) m[t][q] += coef[t] * sign * P.mu * G[C == 1 || t < D ? a : a + Z::stride(t)];
            }
        }
    }
};
// Every subset of the axes with a plane, in increasing order: add P's integrals for the axes cut so far (MASK), then
// cut P by each later axis and recurse into the parts. (P is sorted in place: the integrals do not depend on the
// order of its corners.)
template<int D, int N, int MASK, int A, class M>
__device__ __forceinline__ void branch(Simplex<D> P, Real sign, int planes, M& acc);
template<int D, int N, int MASK, int NEXT, class M>
__device__ __forceinline__ void visit(const Simplex<D>& P, Real sign, int planes, M& acc) {
    acc.template add<MASK>(P, sign);  // ("template": add is a member template of a dependent type)
    branch<D, N, MASK, NEXT>(P, sign, planes, acc);
}
template<int D, int N, int MASK, int A, class M>
__device__ __forceinline__ void branch(Simplex<D> P, Real sign, int planes, M& acc) {
    if constexpr (A < D) {
        if (planes >> A & 1) {
            sort_corners<D, A>(P);
            int in = 0;
#pragma unroll
            for (int k = 0; k <= D; ++k) in += P.u[k][A] >= 0;
            const int n = parts<D>(in);
#pragma unroll 1
            for (int k = 0; k < n; ++k) {
                Simplex<D> Q;
                const Real s = part<D, A>(P, in, k, Q);
                visit<D, N, MASK | 1 << A, A + 1>(Q, sign * s, planes, acc);
            }
        }
        branch<D, N, MASK, A + 1>(P, sign, planes, acc);
    }
}

// ---- node weights, currents and the atomic additions ------------------------------------------------------------
// Basis integrals to node weights: the window matrix along every axis; along flux_axis (the current's own axis; -1:
// none), the face matrix instead: then node i of that axis holds the current through face i + 1/2, divided by q / dt.
template<int D, int N> __device__ __forceinline__ void to_nodes(Real (&m)[Win<D, N>::T], int flux_axis) {
    using S = Win<D, N>;
#pragma unroll
    for (int d = 0; d < D; ++d) {
        Real r[S::T];
#pragma unroll
        for (int o = 0; o < S::T; ++o) {
            const int i = S::digit(o, d), first = o - i * S::stride(d);
            Real s = 0;
#pragma unroll
            for (int q = 0; q < S::W; ++q)
                s += (d == flux_axis ? tables.flux[N][i][q] : tables.window[N][i][q]) * m[first + q * S::stride(d)];
            r[o] = s;
        }
#pragma unroll
        for (int o = 0; o < S::T; ++o) m[o] = r[o];
    }
}
// Node weights W (one per component; C = 1: density) to the grid: Esirkepov's current along an axis, q W / dt (face
// values, the last face of the window carries none); out-of-plane components and density, q W.
template<int D, int N, int C>
__device__ __forceinline__ void scatter(const Real (&W)[C][Win<D, N>::T], const int* base, Real q, const Grid& g,
                                        Out out) {
    using S = Win<D, N>;
    int node[3][S::W];  // wrapped grid index of each window node per axis, scaled by the axis stride
#pragma unroll
    for (int d = 0; d < 3; ++d) {
        const int n = g.n[d], stride = d == 0 ? 1 : d == 1 ? g.n[0] : g.n[0] * g.n[1];
        int i = d < D ? wrap(base[d], n) : 0;
#pragma unroll
        for (int k = 0; k < S::W; ++k) {
            node[d][k] = i * stride;
            i = i + 1 == n ? 0 : i + 1;
        }
    }
    auto index = [&](int t) {
        return node[0][S::digit(t, 0)] + (D > 1 ? node[1][S::digit(t, 1)] : 0) + (D > 2 ? node[2][S::digit(t, 2)] : 0);
    };
    auto add = [&](int c, int t, Real value) {
        if (value != 0) atomicAdd(out.p[c] + index(t), double(value));
    };
#pragma unroll
    for (int c = 0; c < C; ++c)
#pragma unroll
        for (int t = 0; t < S::T; ++t) {
            if (C == 1 || c >= D) add(c, t, q * W[c][t]);
            else if (S::digit(t, c) < N + 1) add(c, t, q / g.dt * W[c][t]);
        }
}
// One element (or piece) that fits its window. Current: the 2^D mixed simplices, axis d at time level a_d, weighted
// by Esirkepov; density: the new time level.
template<int D, int N, bool CURRENT>
__device__ __forceinline__ void deposit_simplex(const Local<D, D + 1>& L, const Frame<D>& f, const Grid& g,
                                                Out out) {
    constexpr int C = CURRENT ? 3 : 1;
    Moments<D, N, C> acc;
#pragma unroll
    for (int c = 0; c < C; ++c)
#pragma unroll
        for (int q = 0; q < Win<D, N>::T; ++q) acc.m[c][q] = 0;
#pragma unroll 1
    for (int a = CURRENT ? 0 : (1 << D) - 1; a < 1 << D; ++a) {
        Simplex<D> E;
        E.mu = 1;
#pragma unroll
        for (int k = 0; k <= D; ++k) {
#pragma unroll
            for (int d = 0; d < D; ++d) E.u[k][d] = L.x[a >> d & 1][k][d] - f.plane[d];
#pragma unroll
            for (int c = 0; c < 3 - D; ++c) E.v[k][c] = L.v[k][c];
        }
        const int ones = bits(a);
        if (CURRENT) {
#pragma unroll
            for (int c = 0; c < C; ++c)  // along an axis: the time difference on it, the average on the others
                acc.coef[c] = c < D ? Real(a >> c & 1 ? 1 : -1) * esirkepov(D - 1, ones - (a >> c & 1))
                                    : esirkepov(D, ones);
        } else acc.coef[0] = 1;
        visit<D, N, 0, 0>(E, Real(1), f.planes, acc);
    }
#pragma unroll
    for (int c = 0; c < C; ++c) to_nodes<D, N>(acc.m[c], C > 1 && c < D ? c : -1);
    scatter<D, N, C>(acc.m, f.base, L.q, g, out);
}
// A curved 1D element: its weights by Gauss quadrature in the material coordinate s, cut at the turning point and
// where the curve crosses the window's plane. On each interval the integrand is a polynomial in s of degree at most
// 2N + 2 (a weight of the quadratic x(s), times a quadratic velocity): 6 points are exact to degree 11.
template<int N, bool CURRENT>
__device__ __forceinline__ void deposit_curved(const Local<1, 3>& L, const Frame<1>& f, const Grid& g, Out out) {
    constexpr int C = CURRENT ? 3 : 1, W = N + 2;
    const Real node[6] = {.033765242898423986, .16939530676686774, .38069040695840155,
                          .61930959304159845,  .83060469323313226, .96623475710157601};
    const Real weight[6] = {.08566224618958517, .18038078652406930, .23395696728634552,
                            .23395696728634552, .18038078652406930, .08566224618958517};
    Real m[C][W];
#pragma unroll
    for (int c = 0; c < C; ++c)
#pragma unroll
        for (int i = 0; i < W; ++i) m[c][i] = 0;
#pragma unroll 1
    for (int l = CURRENT ? 0 : 1; l < 2; ++l) {
        const Real x0 = L.x[l][0][0] - f.plane[0], x1 = L.x[l][1][0] - f.plane[0], xm = L.x[l][2][0] - f.plane[0];
        // [0, 1] in monotone parts (split at the turning point), each split where it crosses the plane
        const Real s = turning(x0, x1, xm), bounds[3] = {0, s > 0 && s < 1 ? s : Real(1), 1};
        Real ends[5];
        int count = 0;
        for (int p = 0; p < 2; ++p) {
            if (!(bounds[p] < bounds[p + 1])) continue;
            ends[count++] = bounds[p];
            const Real a = quadratic(x0, x1, xm, bounds[p]), b = quadratic(x0, x1, xm, bounds[p + 1]);
            if ((a >= 0) != (b >= 0)) ends[count++] = crossing(x0, x1, xm, bounds[p], bounds[p + 1]);
        }
        ends[count] = 1;
        for (int p = 0; p < count; ++p) {
            const Real lo = ends[p], h = ends[p + 1] - ends[p];
            if (h <= 0) continue;
#pragma unroll
            for (int k = 0; k < 6; ++k) {
                const Real sk = lo + h * node[k], wk = h * weight[k], u = quadratic(x0, x1, xm, sk);
                Real S[W], F[W];  // node weights; face weights (the x current)
                window_weights<N>(u, S);
                if constexpr (CURRENT) {
                    window_weights<N>(u, F, true);
                    const Real vy = quadratic(L.v[0][0], L.v[1][0], L.v[2][0], sk);
                    const Real vz = quadratic(L.v[0][1], L.v[1][1], L.v[2][1], sk);
#pragma unroll
                    for (int i = 0; i < W; ++i) {
                        m[0][i] += (l == 0 ? -wk : wk) * F[i];
                        m[1][i] += Real(.5) * wk * S[i] * vy;
                        m[2][i] += Real(.5) * wk * S[i] * vz;
                    }
                } else {
#pragma unroll
                    for (int i = 0; i < W; ++i) m[0][i] += wk * S[i];
                }
            }
        }
    }
    scatter<1, N, C>(m, f.base, L.q, g, out);
}
// (Inlined into both kernels: no call frame in local memory.)
template<int D, int N, int K, bool CURRENT>
__device__ __forceinline__ void deposit_local(const Local<D, K>& L, const Frame<D>& f, const Grid& g, Out out) {
    if constexpr (curved<D, K>) deposit_curved<N, CURRENT>(L, f, g, out);
    else deposit_simplex<D, N, CURRENT>(L, f, g, out);
}

// ---- kernels ------------------------------------------------------------------------------------------------------
// Every element that fits its window is deposited here; the others go to a list for k_pieces. Consecutive blocks take
// stretches of the element list 1/64 of it apart: in Morton order the blocks running at once would otherwise all
// deposit into the same few cells, and their atomic adds would queue at the same addresses (2D folded sheets: the
// deposit 1.7x slower). The launch has 64 x stride blocks (spread_blocks); block 64 a + r takes the elements of block
// r stride + a (a transpose: each element once), lanes unchanged.
constexpr int spread = 64;
inline int spread_blocks(int count) { return spread * ((count + spread * threads - 1) / (spread * threads)); }
template<int D, int N, int K, bool CURRENT>
__global__ void k_elements(const Vertex<D>* vx, const Element<D>* el, const int* mid, int count, Grid g, Out out,
                           int* list, int* status) {
    const int stride = gridDim.x / spread, b = blockIdx.x;
    const int e = ((b % spread) * stride + b / spread) * blockDim.x + threadIdx.x;
    if (e >= count) return;
    const Local<D, K> L = load<D, K>(vx, el[e], curved<D, K> ? mid[e] : 0, g);
    const Frame<D> f = frame<D, N>(L, CURRENT);
    if (!f.finite) atomicOr(status, BAD_VALUE);
    else if (!f.fits) list[2 + atomicAdd(list, 1)] = e;
    else deposit_local<D, N, K, CURRENT>(L, f, g, out);
}
// The listed elements, bisected depth first into pieces that fit. A piece is a matrix lam of barycentric rows over
// the element's corners (curved: over its control points; a child is the quadratic on half the parameter interval).
// A piece is bisected along the axes it does not fit. 24 levels take segments 2^24 cells long; a triangle halves its
// extent every two levels, and 40 take it about 10^6 cells across (a phase image's cells are pixels, up to 65536 along
// a momentum); a tetrahedron halves it only every 3-4 levels, and 40 take it to about 700 cells. A deeper piece is dropped with BAD_PIECES (the run then
// stops). With fewer listed elements than threads (a 1D fixed mesh: a few hundred, each cut into tens of pieces), each
// element's first `fan` levels are shared out among 2^fan threads: entry i follows path i mod 2^fan (bit j: the child at
// level j), and a piece that fits before level fan is deposited by the one path whose remaining bits are 0, so that it
// counts once.
template<int D, int N, int K, bool CURRENT>
__global__ void k_pieces(const Vertex<D>* vx, const Element<D>* el, const int* mid, Grid g, Out out, int* list,
                         int* status) {
    constexpr int depth_max = D == 1 ? 24 : 40;
    // Each thread takes the entry of its own index, then the next unclaimed one from a work queue in global memory
    // (list[1]: a thread that finishes takes another, whatever the others still do; nothing listed: no claims).
    const int all = gridDim.x * blockDim.x, listed = list[0];
    int fan = 0;
    while (fan < 10 && (long(listed) << (fan + 1)) <= all) ++fan;
    const int entries = listed << fan;  // (fan > 0: at most all; fan 0: the listed elements)
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < entries; i = all + atomicAdd(list + 1, 1)) {
        const int e = list[2 + (i >> fan)], path = i & ((1 << fan) - 1);
        const Local<D, K> L = load<D, K>(vx, el[e], curved<D, K> ? mid[e] : 0, g);
        Real stack[depth_max + 1][K][K];
        int depth[depth_max + 1], top = 1;
        for (int a = 0; a < K; ++a)
            for (int b = 0; b < K; ++b) stack[0][a][b] = a == b;
        depth[0] = 0;
        while (top > 0) {
            --top;
            Real lam[K][K];
            for (int a = 0; a < K; ++a)
                for (int b = 0; b < K; ++b) lam[a][b] = stack[top][a][b];
            const int level = depth[top];
            Local<D, K> P;
            P.q = ldexp(L.q, -level);
#pragma unroll
            for (int d = 0; d < D; ++d) P.o[d] = L.o[d];
#pragma unroll
            for (int k = 0; k < K; ++k) {
#pragma unroll
                for (int l = 0; l < 2; ++l)
#pragma unroll
                    for (int d = 0; d < D; ++d) {
                        Real s = 0;
#pragma unroll
                        for (int j = 0; j < K; ++j) s += lam[k][j] * L.x[l][j][d];
                        P.x[l][k][d] = s;
                    }
#pragma unroll
                for (int c = 0; c < 3 - D; ++c) {
                    Real s = 0;
#pragma unroll
                    for (int j = 0; j < K; ++j) s += lam[k][j] * L.v[j][c];
                    P.v[k][c] = s;
                }
            }
            const Frame<D> f = frame<D, N>(P, CURRENT);
            if (f.fits) {
                if (level >= fan || (path >> level) == 0) deposit_local<D, N, K, CURRENT>(P, f, g, out);
                continue;
            }
            if (level == depth_max || top + 2 > depth_max + 1) { atomicOr(status, BAD_PIECES); continue; }
            Real child[2][K][K];
            if constexpr (curved<D, K>) {  // [0, 1/2]: (x(0), x(1/2), x(1/4)); [1/2, 1]: (x(1/2), x(1), x(3/4))
                const Real M[2][3][3] = {{{1, 0, 0}, {0, 0, 1}, {.375, -.125, .75}},
                                         {{0, 0, 1}, {0, 1, 0}, {-.125, .375, .75}}};
                for (int h = 0; h < 2; ++h)
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b)
                            child[h][a][b] = M[h][a][0] * lam[0][b] + M[h][a][1] * lam[1][b] + M[h][a][2] * lam[2][b];
            } else {  // bisect the edge with the largest extent along the axes that do not fit (the others stay)
                int ei = 0, ej = 1;
                Real longest = -1;
                for (int a = 0; a < K; ++a)
                    for (int b = a + 1; b < K; ++b) {
                        Real extent = 0;
                        for (int l = 0; l < 2; ++l)
                            for (int d = 0; d < D; ++d)
                                if (f.misfit >> d & 1) extent = fmax(extent, fabs(P.x[l][a][d] - P.x[l][b][d]));
                        if (extent > longest) { longest = extent; ei = a; ej = b; }
                    }
                for (int h = 0; h < 2; ++h)
                    for (int a = 0; a < K; ++a)
                        for (int b = 0; b < K; ++b)
                            child[h][a][b] = a == (h ? ei : ej) ? (lam[ei][b] + lam[ej][b]) / 2 : lam[a][b];
            }
            for (int h = 0; h < 2; ++h) {
                if (level < fan && h != (path >> level & 1)) continue;  // (before level fan: the path's child only)
                for (int a = 0; a < K; ++a)
                    for (int b = 0; b < K; ++b) stack[top][a][b] = child[h][a][b];
                depth[top++] = level + 1;
            }
        }
    }
}
// PIC: push, then deposit the move from the old to the new position (a window always holds it: |dx| < 1). Esirkepov's
// current of a point factorizes: along axis c, the change of the face weights on c times the time average of the
// node weights over the other axes.
template<int D, int N> __device__ __forceinline__ Real averaged(const Real (&weight)[2][D][N + 2], int mask, int t) {
    using S = Win<D, N>;
    Real sum = 0;
#pragma unroll
    for (int a = 0; a < 1 << D; ++a) {  // the time levels of the axes in mask
        if (a & ~mask) continue;
        Real p = esirkepov(bits(mask), bits(a));
#pragma unroll
        for (int d = 0; d < D; ++d)
            if (mask >> d & 1) p *= weight[a >> d & 1][d][S::digit(t, d)];
        sum += p;
    }
    return sum;
}
template<int D, int N>
__global__ void k_particles(Vertex<D>* particles, int count, Grid g, Grid deposit_grid, FieldPtr fields, Real qm,
                            Out out, int* status) {
    using S = Win<D, N>;
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= count) return;
    Vertex<D> p = particles[k];
    const bool ok = push<D, N>(p, g, fields, qm);
    particles[k].x = p.x;  // (a particle's old position stays in registers)
#pragma unroll
    for (int c = 0; c < 3; ++c) particles[k].u[c] = p.u[c];
    if (!ok) { atomicOr(status, BAD_VALUE); return; }
    Real weight[2][D][S::W], face[2][D][S::W];  // node and face weights per time level and axis
    int base[D];
#pragma unroll
    for (int d = 0; d < D; ++d) {
        const Real x0 = p.old.f[d] + shift(N), x1 = coordinate(p.x, d, p.old.c[d]) + shift(N), b = floor(fmin(x0, x1));
        if (fmax(x0, x1) >= b + 2) { atomicOr(status, BAD_VALUE); return; }
        base[d] = p.old.c[d] + int(b) - N / 2;
#pragma unroll
        for (int l = 0; l < 2; ++l) {
            window_weights<N>((l ? x1 : x0) - (b + 1), weight[l][d]);
            window_weights<N>((l ? x1 : x0) - (b + 1), face[l][d], true);
        }
    }
    const Real inverse_gamma = 1 / sqrt(1 + p.u[0] * p.u[0] + p.u[1] * p.u[1] + p.u[2] * p.u[2]);
    Real J[3][S::T];
#pragma unroll
    for (int c = 0; c < D; ++c)  // lines along c: the transverse average once per line
#pragma unroll
        for (int line = 0; line < S::T / S::W; ++line) {
            const int t = line % S::stride(c) + line / S::stride(c) * S::stride(c) * S::W;  // (digit c set to 0)
            const Real transverse = averaged<D, N>(weight, ((1 << D) - 1) & ~(1 << c), t);
#pragma unroll
            for (int i = 0; i < S::W; ++i) J[c][t + i * S::stride(c)] = (face[1][c][i] - face[0][c][i]) * transverse;
        }
#pragma unroll
    for (int c = D; c < 3; ++c)
#pragma unroll
        for (int t = 0; t < S::T; ++t) J[c][t] = p.u[c] * inverse_gamma * averaged<D, N>(weight, (1 << D) - 1, t);
    scatter<D, N, 3>(J, base, p.q, deposit_grid, out);
}
template<int D, int N>
__global__ void k_particle_density(const Vertex<D>* particles, int count, Grid g, Out out) {
    using S = Win<D, N>;
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= count) return;
    const Pos<D> x = particles[k].x;
    const Real q = particles[k].q;
    Real weight[D][S::W], W[1][S::T];
    int base[D];
#pragma unroll
    for (int d = 0; d < D; ++d) {
        const Real b = floor(x.f[d] + shift(N));
        base[d] = x.c[d] + int(b) - N / 2;
        window_weights<N>(x.f[d] + shift(N) - (b + 1), weight[d]);
    }
#pragma unroll
    for (int t = 0; t < S::T; ++t) {
        W[0][t] = 1;
#pragma unroll
        for (int d = 0; d < D; ++d) W[0][t] *= weight[d][S::digit(t, d)];
    }
    scatter<D, N, 1>(W, base, q, g, out);
}

// ---- conducting walls --------------------------------------------------------------------------------------------
// A window reaches up to N/2 + 1 nodes beyond a wall. With walls the deposit writes into a padded grid, G = N/2 + 2
// extra nodes past each wall (an unwrapped node index j < 0 kept at n' + j), and k_fold adds those back through the
// gather's mirror (sic.hpp: gather): odd across a wall along a node axis, even along a half-cell axis (a current along
// its own axis): the fold is the transpose of the gather's mirror, and charge is conserved off the wall nodes.
Grid padded(const Grid& g) {
    Grid p = g;
    for (int d = 0; d < g.dim; ++d)
        if (g.wall[d]) p.n[d] = g.n[d] + 2 * (g.order / 2 + 2);
    return p;
}
__global__ void k_fold(const double* pad, Grid g, Grid p, int half_axis, double* out) {  // half_axis -1: density
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= p.cells()) return;
    const double value = pad[k];
    if (value == 0) return;
    int q[3] = {k % p.n[0], k / p.n[0] % p.n[1], k / (p.n[0] * p.n[1])};
    double sign = 1;
    for (int d = 0; d < 3; ++d) {
        if (!g.wall[d]) continue;
        const int n = g.n[d];
        const int j = q[d] < n + (p.n[d] - n) / 2 ? q[d] : q[d] - p.n[d];  // the unwrapped node index
        q[d] = mirror(j, n, d == half_axis, sign);
    }
    atomicAdd(out + g.index(q[0], q[1], q[2]), sign * value);
}
// Runs deposit(Out, deposit grid) into the padded scratch if the grid has walls, and folds it into `into` (component
// count C: the density's 1 or the current's 3).
template<class F> void with_walls(const Grid& g, Work& w, double* const (&into)[3], int C, F&& deposit) {
    if (!g.wall[0] && !g.wall[1] && !g.wall[2]) { deposit(Out{{into[0], into[1], into[2]}}, g); return; }
    const Grid p = padded(g);
    for (int c = 0; c < C; ++c) {
        w.pad[c].resize(p.cells());
        check(cudaMemsetAsync(ptr(w.pad[c]), 0, p.cells() * sizeof(double)));
    }
    deposit(Out{{ptr(w.pad[0]), C > 1 ? ptr(w.pad[1]) : nullptr, C > 1 ? ptr(w.pad[2]) : nullptr}}, p);
    for (int c = 0; c < C; ++c) k_fold<<<blocks(p.cells()), threads>>>(ptr(w.pad[c]), g, p, C > 1 ? c : -1, into[c]);
    check(cudaGetLastError());
}

template<int D, bool CURRENT>
void deposit(const Mesh<D>& m, const Grid& g, bool curved, Out out, Work& w) {
    const int ne = int(m.element.size());
    if (!ne) return;
    if (int(w.list.size()) < ne + 2) w.list.resize(ne + 2);
    check(cudaMemsetAsync(ptr(w.list), 0, 2 * sizeof(int)));
    with_order<D>(g.order, [&](auto order) {
        constexpr int N = decltype(order)::value;
        auto launch = [&](auto corners) {
            constexpr int K = decltype(corners)::value;
            k_elements<D, N, K, CURRENT><<<spread_blocks(ne), threads>>>(ptr(m.vertex), ptr(m.element), ptr(m.mid), ne, g, out,
                                                                   ptr(w.list), ptr(w.status));
            k_pieces<D, N, K, CURRENT><<<loop_blocks(), threads>>>(ptr(m.vertex), ptr(m.element), ptr(m.mid), g, out,
                                                                    ptr(w.list), ptr(w.status));
        };
        if constexpr (D == 1 && N >= 1) {
            if (curved) { launch(std::integral_constant<int, 3>{}); return; }
        }
        launch(std::integral_constant<int, D + 1>{});
    });
    check(cudaGetLastError());
}
}  // namespace

// (The deposit kernels take the deposit grid: along a walled axis its node count differs, which changes nothing else
// they compute, since no element has a periodic image along a walled axis.)
template<int D> void deposit_current(const Mesh<D>& m, const Grid& g, bool curved, double* const (&j)[3], Work& w) {
    with_walls(g, w, j, 3, [&](Out out, const Grid& dg) { deposit<D, true>(m, dg, curved, out, w); });
}
template<int D> void deposit_current(const Mesh<D>& m, const Grid& g, bool curved, Fields& f, Work& w) {
    double* const j[3] = {ptr(f.j[0]), ptr(f.j[1]), ptr(f.j[2])};
    deposit_current(m, g, curved, j, w);
}
template<int D> void deposit_density(const Mesh<D>& m, const Grid& g, bool curved, Device<double>& rho, Work& w) {
    double* const into[3] = {ptr(rho), nullptr, nullptr};
    with_walls(g, w, into, 1, [&](Out out, const Grid& dg) {
        if (!m.element.empty()) { deposit<D, false>(m, dg, curved, out, w); return; }
        const int np = int(m.vertex.size());
        with_order<D>(g.order, [&](auto order) {
            k_particle_density<D, decltype(order)::value><<<blocks(np), threads>>>(ptr(m.vertex), np, dg, out);
        });
        check(cudaGetLastError());
    });
}
template<int D> void push_deposit_particles(Mesh<D>& m, const Grid& g, Fields& f, Work& w) {
    const int np = int(m.vertex.size());
    double* const j[3] = {ptr(f.j[0]), ptr(f.j[1]), ptr(f.j[2])};
    with_walls(g, w, j, 3, [&](Out out, const Grid& dg) {
        with_order<D>(g.order, [&](auto order) {
            k_particles<D, decltype(order)::value>
                <<<blocks(np), threads>>>(ptr(m.vertex), np, g, dg, pointers(f), m.qm, out, ptr(w.status));
        });
        check(cudaGetLastError());
    });
}
#define SIC_INSTANTIATE(D)                                                                         \
    template void deposit_current<D>(const Mesh<D>&, const Grid&, bool, Fields&, Work&);           \
    template void deposit_current<D>(const Mesh<D>&, const Grid&, bool, double* const (&)[3], Work&); \
    template void deposit_density<D>(const Mesh<D>&, const Grid&, bool, Device<double>&, Work&);   \
    template void push_deposit_particles<D>(Mesh<D>&, const Grid&, Fields&, Work&);
#if !defined(SIC_DIM) || SIC_DIM == 1  // the Makefile builds one object per dimension, in parallel
SIC_INSTANTIATE(1)
#endif
#if !defined(SIC_DIM) || SIC_DIM == 2
SIC_INSTANTIATE(2)
#endif
#if !defined(SIC_DIM) || SIC_DIM == 3
SIC_INSTANTIATE(3)

// ---- 3D column images -----------------------------------------------------------------------------------------------
// The column of charge along w through the slab at a sample point: each tetrahedron adds its density (charge / volume)
// times the length of the sample's ray inside it and inside the slab's periodic copies (the full code's 3D render). Its
// set-up is in double precision; the samples in single precision unless the tetrahedron is thin (Column::thin: single
// precision's 6e-8 times (w extent x shadow area) / volume, the relative error of its column, would exceed 6e-4, under
// the frames' float16; a thin one's column stays finite as its volume goes to zero, where single precision loses it,
// and double precision runs at 1/64 of the single rate on GeForce GPUs: all-double frames took 2.8x as long).
// The corners are unwrapped by the element's images and moved by whole periods so that corner 0 lies in the box (the
// pixel indices stay small; the image is periodic in u and v). A ray on a face parallel to it belongs to the
// tetrahedron on the side of a symbolic displacement of the sample, along +u first, then +v, so that the two sharing the
// face do not both count it. The sample points sit off the pixel centres (sic.hpp: sample_u, sample_v), so that the initial
// lattices' faces miss them even where the arithmetic is inexact. A tetrahedron of zero volume adds nothing (its
// charge lies in a plane); one whose shadow holds more than big_shadow sample points is listed for a block of its own.
namespace {
constexpr long big_shadow = 4096;
template<class T> struct Column {  // one tetrahedron, ready to draw (T: the samples' precision)
    double u0, v0, w0;    // corner 0: u, v in pixels, w in cells
    T m[4][3];            // the gradients of its barycentric coordinates in (u, v, w)
    T inverse_b[4];       // 1 / m[r][2] (where it is not 0)
    double density;       // charge per pixel^2 per cell
    int i0, i1, k0, k1;   // its shadow's sample points: columns i0..i1, rows k0..k1
    int j0, j1;           // the slab's copies it can reach: [lo, hi) + j n_w
    bool thin;            // its samples need double precision: (w extent x shadow area) / volume above 1e4
    __device__ long size() const { return long(i1 - i0 + 1) * (k1 - k0 + 1); }
    template<class U> __device__ Column<U> as() const {
        Column<U> c{u0, v0, w0, {}, {}, density, i0, i1, k0, k1, j0, j1, thin};
        for (int r = 0; r < 4; ++r) {
            for (int a = 0; a < 3; ++a) c.m[r][a] = U(m[r][a]);
            c.inverse_b[r] = U(inverse_b[r]);
        }
        return c;
    }
};
// false: nothing to draw (no volume, or not finite: then the status bit)
__device__ bool column_setup(const Vertex<3>* vx, const Element<3>& el, const Grid& g, double alpha, const View& view,
                             int* status, Column<double>& c) {
    const int axis[3] = {view.axis[0], view.axis[1], view.w};
    double p[4][3];
    for (int k = 0; k < 4; ++k)
        for (int a = 0; a < 3; ++a) {
            const Pos<3>& x0 = vx[el.v[k]].old;
            const Pos<3>& x1 = vx[el.v[k]].x;
            const int d = axis[a];
            const double y0 = x0.c[d] + double(x0.f[d]), y1 = x1.c[d] + double(x1.f[d]);
            p[k][a] = y0 + alpha * (y1 - y0) + double(el.img[k][d]) * g.n[d];
        }
    for (int a = 0; a < 3; ++a) {  // whole periods: corner 0 into the box; u, v in pixels
        const double n = g.n[axis[a]], shift = floor(p[0][a] / n) * n;
        for (int k = 0; k < 4; ++k) p[k][a] = (p[k][a] - shift) * (a < 2 ? view.scale[a] : 1.0);
    }
    double t[3][3];  // T = [p1 - p0 | p2 - p0 | p3 - p0]
    for (int a = 0; a < 3; ++a)
        for (int k = 0; k < 3; ++k) t[a][k] = p[k + 1][a] - p[0][a];
    const double det = t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1]) -
                       t[0][1] * (t[1][0] * t[2][2] - t[1][2] * t[2][0]) +
                       t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0]);
    if (!isfinite(det) || !isfinite(double(el.q))) { atomicOr(status, BAD_IMAGE); return false; }
    if (det == 0) return false;
    // Row k of T^-1 is the gradient of barycentric coordinate k + 1: entry (k, a) = cofactor of T[a][k] / det.
    for (int k = 0; k < 3; ++k)
        for (int a = 0; a < 3; ++a) {
            const int a1 = (a + 1) % 3, a2 = (a + 2) % 3, k1 = (k + 1) % 3, k2 = (k + 2) % 3;
            c.m[k + 1][a] = (t[a1][k1] * t[a2][k2] - t[a1][k2] * t[a2][k1]) / det;
        }
    for (int a = 0; a < 3; ++a) c.m[0][a] = -(c.m[1][a] + c.m[2][a] + c.m[3][a]);
    for (int r = 0; r < 4; ++r) c.inverse_b[r] = c.m[r][2] != 0 ? 1 / c.m[r][2] : 0;
    c.u0 = p[0][0], c.v0 = p[0][1], c.w0 = p[0][2];
    c.density = double(el.q) / (fabs(det) / 6);
    double lo[3], hi[3];
    for (int a = 0; a < 3; ++a) {
        lo[a] = hi[a] = p[0][a];
        for (int k = 1; k < 4; ++k) lo[a] = fmin(lo[a], p[k][a]), hi[a] = fmax(hi[a], p[k][a]);
    }
    c.i0 = int(ceil(lo[0] - sample_u)), c.i1 = int(floor(hi[0] - sample_u));
    c.k0 = int(ceil(lo[1] - sample_v)), c.k1 = int(floor(hi[1] - sample_v));
    const double nw = g.n[view.w];
    c.j0 = int(floor((lo[2] - view.hi) / nw)) + 1, c.j1 = int(floor((hi[2] - view.lo) / nw));
    c.thin = (hi[2] - lo[2]) * (hi[0] - lo[0]) * (hi[1] - lo[1]) > 1e4 * fabs(det) / 6;
    return c.i0 <= c.i1 && c.k0 <= c.k1 && c.j0 <= c.j1;  // (j0 > j1: outside every copy of the slab)
}
template<class T> __device__ void column_sample(const Column<T>& c, int i, int k, const View& view, int nw, double* image) {
    const T du = T(i + sample_u - c.u0), dv = T(k + sample_v - c.v0);
    T lo = -INFINITY, hi = INFINITY;  // the ray's w - w0 inside the tetrahedron
    for (int r = 0; r < 4; ++r) {
        const T a = T(r == 0) + c.m[r][0] * du + c.m[r][1] * dv;  // coordinate r = a + m[r][2] (w - w0)
        if (c.m[r][2] > 0) lo = max(lo, -a * c.inverse_b[r]);
        else if (c.m[r][2] < 0) hi = min(hi, -a * c.inverse_b[r]);
        else if (a < 0 || (a == 0 && (c.m[r][0] < 0 || (c.m[r][0] == 0 && c.m[r][1] < 0)))) return;  // (ties: +u, +v)
    }
    if (!(hi > lo)) return;
    double length = 0;
    for (int j = c.j0; j <= c.j1; ++j)
        length += fmax(0.0, fmin(c.w0 + hi, view.hi + j * nw) - fmax(c.w0 + lo, view.lo + j * nw));
    if (length > 0) atomicAdd(image + std::size_t(wrap(k, view.n[1])) * view.n[0] + wrap(i, view.n[0]), c.density * length);
}
__global__ void k_columns(const Vertex<3>* vx, const Element<3>* el, int count, Grid g, double alpha, View view,
                          double* image, int* list, int* status) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= count) return;
    Column<double> c;
    if (!column_setup(vx, el[e], g, alpha, view, status, c)) return;
    if (c.size() > big_shadow) { list[2 + atomicAdd(list, 1)] = e; return; }
    auto draw = [&](const auto& t) {
        for (int k = t.k0; k <= t.k1; ++k)
            for (int i = t.i0; i <= t.i1; ++i) column_sample(t, i, k, view, g.n[view.w], image);
    };
    if (c.thin) draw(c);
    else draw(c.as<float>());
}
// The listed tetrahedra: a block each, its threads striding over the shadow.
__global__ void k_columns_big(const Vertex<3>* vx, const Element<3>* el, Grid g, double alpha, View view, double* image,
                              const int* list, int* status) {
    for (int t = blockIdx.x; t < list[0]; t += gridDim.x) {
        Column<double> c;
        if (!column_setup(vx, el[list[2 + t]], g, alpha, view, status, c)) continue;
        auto draw = [&](const auto& t) {
            const int width = t.i1 - t.i0 + 1;
            for (long s = threadIdx.x; s < t.size(); s += blockDim.x)
                column_sample(t, t.i0 + int(s % width), t.k0 + int(s / width), view, g.n[view.w], image);
        };
        if (c.thin) draw(c);
        else draw(c.as<float>());
    }
}
// PIC: a particle's charge into the pixel column it is in, if it lies in the slab.
__global__ void k_column_particles(const Vertex<3>* vx, int count, Grid g, View view, double* image) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= count) return;
    const Pos<3>& x = vx[k].x;
    const int nw = g.n[view.w];
    const double w = wrap(x.c[view.w], nw) + double(x.f[view.w]);  // (in [0, nw + 1))
    if (!(w - view.lo - floor((w - view.lo) / nw) * nw < view.hi - view.lo)) return;
    int at[2];
    for (int a = 0; a < 2; ++a) {
        const int d = view.axis[a];
        at[a] = wrap(int(floor((wrap(x.c[d], g.n[d]) + double(x.f[d])) * view.scale[a])), view.n[a]);
    }
    atomicAdd(image + std::size_t(at[1]) * view.n[0] + at[0], double(vx[k].q));
}
}  // namespace

void column_image(const Mesh<3>& m, const Grid& g, const View& view, double alpha, double* image, Work& w) {
    const int ne = int(m.element.size()), nv = int(m.vertex.size());
    if (!ne) {
        if (nv) k_column_particles<<<blocks(nv), threads>>>(ptr(m.vertex), nv, g, view, image);
        check(cudaGetLastError());
        return;
    }
    if (int(w.list.size()) < ne + 2) w.list.resize(ne + 2);
    check(cudaMemsetAsync(ptr(w.list), 0, 2 * sizeof(int)));
    k_columns<<<blocks(ne), threads>>>(ptr(m.vertex), ptr(m.element), ne, g, alpha, view, image, ptr(w.list), ptr(w.status));
    k_columns_big<<<loop_blocks(), 256>>>(ptr(m.vertex), ptr(m.element), g, alpha, view, image, ptr(w.list),
                                            ptr(w.status));
    check(cudaGetLastError());
}
#endif

// ---- phase-space images ---------------------------------------------------------------------------------------------
// Coordinates in double, in pixels: positions over the whole box, an element's corners unwrapped by its images and
// moved by whole periods so that corner 0 lies in the box (as the column images); momenta from lo. Along a momentum
// axis an element is cut exactly at the image's edges in its own material coordinates, the part inside keeping its
// share of the charge, so nothing outside the range reaches the image however far an element extends. PIC: each
// particle in its pixel. 1D: a segment's charge is uniform along it; each pixel takes the fraction of the segment
// inside it (a walk through the pixels between the ends of the part inside: no bisection depth, for a fixed mesh's
// segments thousands of pixels long). 2D: an element maps affinely onto a triangle of the plane, uniform in charge
// (collapsed ones too: material measures); a cut one becomes the up to five triangles of its polygon; the order-0
// deposit integrates them over each pixel exactly (as the density images), in batches, on a grid with a guard pixel
// beyond each edge that is not periodic (round-off at the cuts). A point on an upper wall is in the last pixel. (3D
// SIC: a tetrahedron's shadow is not uniform; not drawn.)
namespace {
struct PhaseAxis {
    int axis, n;          // 0..2 a position, 3..5 a momentum; pixels
    double lo, scale;     // the pixel coordinate: position (cells) x scale, or (u - lo) x scale
    bool periodic, wall;  // a position axis: periodic, or between walls
};
PhaseAxis phase_axis(const Grid& g, const PhaseView& view, const int (&size)[2], int k) {
    PhaseAxis a;
    a.axis = view.axis[k], a.n = size[k];
    const bool momentum = a.axis >= 3;
    a.lo = momentum ? view.lo[k] : 0;
    a.scale = momentum ? a.n / (view.hi[k] - view.lo[k]) : double(a.n) / g.n[a.axis];
    a.wall = !momentum && g.wall[a.axis];
    a.periodic = !momentum && !a.wall;
    return a;
}
// Corner k of element e along axis a, in pixels (an upper wall's point: just below it, in the last pixel).
template<int D>
__device__ double corner(const Vertex<D>* v, const Element<D>& e, int k, const PhaseAxis& a, const Grid& g) {
    if (a.axis >= 3) return (double(v[e.v[k]].u[a.axis - 3]) - a.lo) * a.scale;
    const int d = a.axis;
    const double n = g.n[d];
    auto at = [&](int j) { const Pos<D>& x = v[e.v[j]].x; return x.c[d] + double(x.f[d]) + double(e.img[j][d]) * n; };
    const double y = (at(k) - (a.periodic ? floor(at(0) / n) * n : 0)) * a.scale;
    return a.wall ? fmin(y, a.n - 1e-6) : y;
}
template<int D>
__global__ void k_phase_points(const Vertex<D>* v, int count, PhaseAxis a, PhaseAxis b, Grid g, double* image,
                               int* status) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= count) return;
    const PhaseAxis ax[2] = {a, b};
    int at[2];
    for (int j = 0; j < 2; ++j) {
        const int d = ax[j].axis;
        const double y = d >= 3 ? (double(v[k].u[d - 3]) - ax[j].lo) * ax[j].scale
                                : ((ax[j].periodic ? wrap(v[k].x.c[d], g.n[d]) : v[k].x.c[d]) + double(v[k].x.f[d])) *
                                      ax[j].scale;
        if (!isfinite(y)) { atomicOr(status, BAD_IMAGE); return; }
        const double i = floor(y);
        if (d < 3) at[j] = int(fmin(fmax(i, 0.0), ax[j].n - 1.0));  // (an upper wall's point, or round-off: the last)
        else if (i >= 0 && i < ax[j].n) at[j] = int(i);
        else return;  // outside the momentum range
    }
    atomicAdd(image + std::size_t(at[1]) * a.n + at[0], double(v[k].q));
}
// 1D: segment e from corner 0 (t = 0) to corner 1 (t = 1); its part inside the momentum ranges, t0..t1, walked through
// the pixels (next: the value of t at the next pixel edge along each axis; per: t per pixel).
template<int D>  // (D = 1)
__global__ void k_phase_segments(const Vertex<D>* v, const Element<D>* el, int count, PhaseAxis a, PhaseAxis b, Grid g,
                                 double* image, int* status) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= count) return;
    const PhaseAxis ax[2] = {a, b};
    double p[2], d[2], t0 = 0, t1 = 1;  // per axis: corner 0, corner 1 - corner 0
    for (int j = 0; j < 2; ++j) {
        p[j] = corner(v, el[e], 0, ax[j], g);
        d[j] = corner(v, el[e], 1, ax[j], g) - p[j];
        if (!isfinite(p[j]) || !isfinite(d[j])) { atomicOr(status, BAD_IMAGE); return; }
        if (ax[j].axis < 3) continue;
        if (d[j] == 0) {
            if (!(p[j] >= 0 && p[j] < ax[j].n)) return;
        } else {
            const double s0 = -p[j] / d[j], s1 = (ax[j].n - p[j]) / d[j];
            t0 = fmax(t0, fmin(s0, s1)), t1 = fmin(t1, fmax(s0, s1));
        }
    }
    if (!(t0 < t1)) return;
    int cell[2], step[2];
    double next[2], per[2];
    for (int j = 0; j < 2; ++j) {
        const double y = p[j] + t0 * d[j];
        cell[j] = int(floor(y));
        step[j] = d[j] > 0 ? 1 : -1;
        per[j] = d[j] != 0 ? 1 / fabs(d[j]) : INFINITY;
        next[j] = d[j] > 0 ? t0 + (cell[j] + 1 - y) / d[j] : d[j] < 0 ? t0 + (y - cell[j]) / -d[j] : INFINITY;
    }
    const double q = el[e].q;
    for (double t = t0;;) {
        const double s = fmin(fmin(next[0], next[1]), t1);
        int at[2];
        bool inside = s > t;
        for (int j = 0; j < 2; ++j) {
            at[j] = ax[j].periodic ? wrap(cell[j], ax[j].n) : cell[j];
            inside = inside && at[j] >= 0 && at[j] < ax[j].n;  // (round-off at the cuts and walls)
        }
        if (inside) atomicAdd(image + std::size_t(at[1]) * a.n + at[0], q * (s - t));
        if (s >= t1) break;
        t = s;
        const int j = next[0] <= next[1] ? 0 : 1;
        cell[j] += step[j];
        next[j] += per[j];
    }
}
// 2D: element first + i as triangles of the plane, each with three vertices of its own (pixel coordinates, the guard
// added, order 0's half-pixel shift taken off), appended from made[0]: the element itself, the polygon left of it by the
// cuts at the momentum edges (as a fan, a triangle's charge its share of the material: area in (s, t), the element's
// own coordinates of corners 1 and 2), or nothing (wholly outside).
template<int D>  // (D = 2)
__global__ void k_phase_triangles(const Vertex<D>* v, const Element<D>* el, int first, int count, PhaseAxis a,
                                  PhaseAxis b, Grid g, Vertex<2>* out, Element<2>* tri, int* made, int* status) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count) return;
    const Element<D>& e = el[first + i];
    const PhaseAxis ax[2] = {a, b};
    double P[8][4];  // the polygon's vertices: s, t, then the pixel coordinates along A and B
    int m = 3;
    for (int k = 0; k < 3; ++k) {
        P[k][0] = k == 1, P[k][1] = k == 2;
        for (int j = 0; j < 2; ++j) P[k][2 + j] = corner(v, e, k, ax[j], g);
    }
    for (int j = 0; j < 2; ++j) {
        const double lo = fmin(fmin(P[0][2 + j], P[1][2 + j]), P[2][2 + j]);
        const double hi = fmax(fmax(P[0][2 + j], P[1][2 + j]), P[2][2 + j]);
        if (!isfinite(lo) || !isfinite(hi)) { atomicOr(status, BAD_IMAGE); return; }
        if (ax[j].axis >= 3 && (hi < 0 || lo >= ax[j].n)) return;  // wholly outside
    }
    for (int j = 0; j < 2; ++j) {
        if (ax[j].axis < 3) continue;
        for (int side = 0; side < 2; ++side) {  // keep y >= 0, then y <= n (Sutherland-Hodgman)
            double Q[8][4];
            int r = 0;
            for (int k = 0; k < m; ++k) {
                const double *A = P[k], *B = P[(k + 1) % m];
                const double da = side ? ax[j].n - A[2 + j] : A[2 + j], db = side ? ax[j].n - B[2 + j] : B[2 + j];
                if (da >= 0) {
                    for (int c = 0; c < 4; ++c) Q[r][c] = A[c];
                    ++r;
                }
                if ((da >= 0) != (db >= 0)) {
                    const double f = da / (da - db);
                    for (int c = 0; c < 4; ++c) Q[r][c] = A[c] + f * (B[c] - A[c]);
                    ++r;
                }
            }
            m = r;
            for (int k = 0; k < m; ++k)
                for (int c = 0; c < 4; ++c) P[k][c] = Q[k][c];
        }
    }
    if (m < 3) return;
    const int t = atomicAdd(made, m - 2);
    for (int f = 0; f < m - 2; ++f) {
        const double* T[3] = {P[0], P[f + 1], P[f + 2]};
        Element<2> r{};
        for (int k = 0; k < 3; ++k) {
            Vertex<2> p{};
            for (int j = 0; j < 2; ++j) place(p.x, j, T[k][2 + j] + (ax[j].periodic ? 0 : 1) - .5);
            p.old = p.x;
            r.v[k] = 3 * (t + f) + k;
            out[r.v[k]] = p;
        }
        const double share = fabs((T[1][0] - T[0][0]) * (T[2][1] - T[0][1]) - (T[2][0] - T[0][0]) * (T[1][1] - T[0][1]));
        r.q = Real(double(e.q) * share);
        tri[t + f] = r;
    }
}
__global__ void k_phase_crop(const double* plane, int width, int na, int nb, int ga, int gb, double* image) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < na * nb) image[k] += plane[std::size_t(k / na + gb) * width + k % na + ga];
}
}  // namespace

template<int D>
void phase_image(const Mesh<D>& m, const Grid& g, const PhaseView& view, const int (&size)[2], double* image, Work& w) {
    const PhaseAxis a = phase_axis(g, view, size, 0), b = phase_axis(g, view, size, 1);
    const int nv = int(m.vertex.size()), ne = int(m.element.size());
    if (!ne) {
        if (nv) k_phase_points<D><<<blocks(nv), threads>>>(ptr(m.vertex), nv, a, b, g, image, ptr(w.status));
    } else if constexpr (D == 1) {
        k_phase_segments<D><<<blocks(ne), threads>>>(ptr(m.vertex), ptr(m.element), ne, a, b, g, image, ptr(w.status));
    } else if constexpr (D == 2) {
        const int ga = a.periodic ? 0 : 1, gb = b.periodic ? 0 : 1;  // the guard pixels
        Grid p;
        p.dim = 2, p.order = 0, p.n[0] = a.n + 2 * ga, p.n[1] = b.n + 2 * gb, p.n[2] = 1;
        Device<double> plane(std::size_t(p.n[0]) * p.n[1], 0.0);
        constexpr int batch = 1 << 18;  // elements per batch (at most 5 triangles each)
        Mesh<2> pieces;
        pieces.vertex.resize(15 * std::size_t(std::min(ne, batch)));
        Device<int> made(1);
        for (int first = 0; first < ne; first += batch) {
            const int count = std::min(batch, ne - first);
            pieces.element.resize(5 * std::size_t(count));
            check(cudaMemsetAsync(ptr(made), 0, sizeof(int)));
            k_phase_triangles<D><<<blocks(count), threads>>>(ptr(m.vertex), ptr(m.element), first, count, a, b, g,
                                                             ptr(pieces.vertex), ptr(pieces.element), ptr(made),
                                                             ptr(w.status));
            check(cudaGetLastError());
            const int triangles = made[0];
            if (!triangles) continue;  // (an empty element list would be PIC particles)
            pieces.element.resize(triangles);
            deposit_density(pieces, p, false, plane, w);
        }
        k_phase_crop<<<blocks(std::size_t(a.n) * b.n), threads>>>(ptr(plane), p.n[0], a.n, b.n, ga, gb, image);
    } else {
        throw std::runtime_error("phase images: SIC in 1D and 2D");
    }
    check(cudaGetLastError());
}
#define SIC_PHASE(D) \
    template void phase_image<D>(const Mesh<D>&, const Grid&, const PhaseView&, const int (&)[2], double*, Work&);
#if !defined(SIC_DIM) || SIC_DIM == 1
SIC_PHASE(1)
#endif
#if !defined(SIC_DIM) || SIC_DIM == 2
SIC_PHASE(2)
#endif
#if !defined(SIC_DIM) || SIC_DIM == 3
SIC_PHASE(3)
#endif
}  // namespace sic
