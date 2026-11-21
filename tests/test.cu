// tests/test.cu - checks of the deposit, refinement, field solver and short runs against independent references.
// Host code here is fixtures and oracles only. Run on a GPU: build/test (double), build/test-single.
#include "sic.hpp"
#include <cuda_fp16.h>
#include <thrust/fill.h>
#include <thrust/host_vector.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <random>
#include <sstream>
#include <vector>
#include <unistd.h>

using namespace sic;
namespace {
int failures = 0;
void expect(bool ok, const char* what, double value = 0) {
    std::printf("%s %-58s %.3g\n", ok ? "PASS" : "FAIL", what, value);
    failures += !ok;
}
const double tol = sizeof(Real) == 8 ? 1e-10 : 2e-4;  // relative to the element charge
// The larger error; NaN wins (std::max would drop it), so a NaN anywhere fails the test.
double worse(double a, double b) { return std::isnan(a) || std::isnan(b) ? NAN : std::max(a, b); }
double flagged(const Work& w) { return w.read_status() ? NAN : 0; }  // a status bit set on the GPU: a failure
std::mt19937_64 rng(12345);
double uniform(double a, double b) { return std::uniform_real_distribution<double>(a, b)(rng); }

template<int D> Pos<D> at(const double* x) {
    Pos<D> p;
    for (int d = 0; d < D; ++d) {
        p.c[d] = int(std::floor(x[d]));
        p.f[d] = Real(x[d] - p.c[d]);
    }
    return p;
}
template<int D> double coord(const Pos<D>& p, int d) { return p.c[d] + double(p.f[d]); }
Grid grid(int D, int n, int order) {
    Grid g;
    g.dim = D;
    g.order = order;
    for (int d = 0; d < 3; ++d) g.n[d] = d < D ? n : 1;
    g.dt = Real(.5 / std::sqrt(double(D)));
    return g;
}
template<int D> Mesh<D> mesh(const std::vector<Vertex<D>>& v, const std::vector<Element<D>>& e) {
    Mesh<D> m;
    m.vertex = thrust::host_vector<Vertex<D>>(v.begin(), v.end());
    m.element = thrust::host_vector<Element<D>>(e.begin(), e.end());
    return m;
}
std::vector<double> host(const Device<double>& a) { thrust::host_vector<double> h = a; return {h.begin(), h.end()}; }
template<int D> void deposit_density(const Mesh<D>& m, const Grid& g, bool curved, Fields& f, Work& w) {  // f.rho only
    thrust::fill(f.rho.begin(), f.rho.end(), 0.0);
    sic::deposit_density(m, g, curved, f.rho, w);
}

// A random element near the middle of a 32-cell box: corners spread by `size` cells, moved by up to 0.4 cells.
template<int D>
void random_element(double size, std::vector<Vertex<D>>& v, std::vector<Element<D>>& e, double q = -0.7) {
    Element<D> el{};
    el.q = Real(q);
    double centre[3];
    for (int d = 0; d < D; ++d) centre[d] = uniform(14, 18);
    for (int k = 0; k <= D; ++k) {
        double x[3], y[3];
        for (int d = 0; d < D; ++d) {
            x[d] = centre[d] + uniform(-size / 2, size / 2);
            y[d] = x[d] + uniform(-.4, .4);
        }
        Vertex<D> p{at<D>(y), at<D>(x), {Real(uniform(-.5, .5)), Real(uniform(-.5, .5)), Real(uniform(-.5, .5))}, 0};
        el.v[k] = int(v.size());
        v.push_back(p);
    }
    e.push_back(el);
}
// E[prod of affine fields] over a simplex under the uniform measure, by expanding in barycentric coordinates:
// E[lambda^a] = D! a! / (D + |a|)!. fields[j][k]: field j at corner k.
double simplex_mean(const std::vector<std::vector<double>>& fields, int D) {
    const int m = int(fields.size());
    double total = 0;
    std::vector<int> pick(m, 0);
    for (;;) {
        int count[4] = {0, 0, 0, 0};
        double product = 1;
        for (int j = 0; j < m; ++j) { product *= fields[j][pick[j]]; ++count[pick[j]]; }
        double weight = 1;
        for (int k = 0; k <= D; ++k) for (int c = 2; c <= count[k]; ++c) weight *= c;
        for (int c = D + 1; c <= D + m; ++c) weight /= c;
        total += product * weight;
        int j = 0;
        while (j < m && ++pick[j] > D) pick[j++] = 0;
        if (j == m) return total;
    }
}

// 1. Density: sum_nodes p(node) W = q E[E_eps p(x + eps)] for every tensor monomial of degree <= N per axis (B-splines
// reproduce them; eps has the spline's moments), on elements from 0.05 to 4 cells (the large ones take pieces).
template<int D, int N> void density_moments() {
    const Grid g = grid(D, 32, N);
    double worst = 0;
    for (double size : {.05, .3, .9, 1.7, 4.0})
        for (int trial = 0; trial < 6; ++trial) {
            std::vector<Vertex<D>> v;
            std::vector<Element<D>> e;
            random_element<D>(size, v, e);
            Mesh<D> m = mesh(v, e);
            Fields f(g);
            Work w;
            deposit_density(m, g, false, f, w);
            worst = worse(worst, flagged(w));
            const std::vector<double> rho = host(f.rho);
            const double var = (N + 1) / 12.0, mu4 = -(N + 1) / 120.0 + 3 * var * var;
            const double moment[5][5] = {{1}, {0, 1}, {var, 0, 1}, {0, 3 * var, 0, 1}, {mu4, 0, 6 * var, 0, 1}};
            int k[3] = {0, 0, 0};
            for (;;) {
                double lhs = 0;
                for (int n = 0; n < g.cells(); ++n) {
                    const int i[3] = {n % g.n[0], n / g.n[0] % g.n[1], n / (g.n[0] * g.n[1])};
                    double p = rho[n];
                    for (int d = 0; d < D; ++d) p *= std::pow(double(i[d]), k[d]);
                    lhs += p;
                }
                double rhs = 0;  // q E[prod_d m_{k_d}(x_d)], expanded into monomials
                int j[3] = {0, 0, 0};
                for (;;) {
                    double c = 1;
                    std::vector<std::vector<double>> fields;
                    for (int d = 0; d < D; ++d) {
                        c *= moment[k[d]][j[d]];
                        for (int r = 0; r < j[d]; ++r) {
                            std::vector<double> x;
                            for (int corner = 0; corner <= D; ++corner) x.push_back(coord(v[e[0].v[corner]].x, d));
                            fields.push_back(x);
                        }
                    }
                    if (c != 0) rhs += c * simplex_mean(fields, D);
                    int d = 0;
                    while (d < D && ++j[d] > k[d]) j[d++] = 0;
                    if (d == D) break;
                }
                rhs *= e[0].q;
                worst = worse(worst, std::fabs(lhs - rhs) / std::pow(20.0, k[0] + k[1] + k[2]));
                int d = 0;
                while (d < D && ++k[d] > N) k[d++] = 0;
                if (d == D) break;
            }
        }
    char what[96];
    std::snprintf(what, sizeof what, "%dD order %d: density moments", D, N);
    expect(worst < tol, what, worst);
}
// 2. Additivity: an element deposits what its 2^L material sub-elements deposit together (density and current): checks
// every cut configuration of the signed corners against pieces that need fewer cuts.
template<int D> Vertex<D> between(const Vertex<D>& a, const Vertex<D>& b) {  // positions and velocities (not momenta)
    double x[3], y[3], v[3], v2 = 0, ga = 1, gb = 1;
    for (int d = 0; d < D; ++d) {
        x[d] = (coord(a.x, d) + coord(b.x, d)) / 2;
        y[d] = (coord(a.old, d) + coord(b.old, d)) / 2;
    }
    for (int c = 0; c < 3; ++c) { ga += double(a.u[c]) * a.u[c]; gb += double(b.u[c]) * b.u[c]; }
    for (int c = 0; c < 3; ++c) { v[c] = (a.u[c] / std::sqrt(ga) + b.u[c] / std::sqrt(gb)) / 2; v2 += v[c] * v[c]; }
    Vertex<D> m{at<D>(x), at<D>(y), {}, 0};
    for (int c = 0; c < 3; ++c) m.u[c] = Real(v[c] / std::sqrt(1 - v2));
    return m;
}
template<int D>
void subdivide(std::vector<Vertex<D>>& v, const Element<D>& e, int levels, std::vector<Element<D>>& out) {
    if (!levels) { out.push_back(e); return; }
    int a = 0, b = 1;  // bisect the longest edge
    double longest = -1;
    for (int i = 0; i <= D; ++i)
        for (int j = i + 1; j <= D; ++j) {
            double s = 0;
            for (int d = 0; d < D; ++d) s += std::pow(coord(v[e.v[i]].x, d) - coord(v[e.v[j]].x, d), 2);
            if (s > longest) { longest = s; a = i; b = j; }
        }
    v.push_back(between(v[e.v[a]], v[e.v[b]]));
    Element<D> c0 = e, c1 = e;
    c0.q = c1.q = e.q / 2;
    c0.v[b] = c1.v[a] = int(v.size()) - 1;
    subdivide(v, c0, levels - 1, out);
    subdivide(v, c1, levels - 1, out);
}
template<int D, int N> void additivity() {
    const Grid g = grid(D, 32, N);
    double worst = 0;
    for (double size : {.4, 1.2, 3.0})
        for (int trial = 0; trial < 6; ++trial) {
            std::vector<Vertex<D>> v;
            std::vector<Element<D>> one, many;
            random_element<D>(size, v, one);
            subdivide(v, one[0], 3 * D, many);
            std::vector<double> result[2][4];
            for (int which = 0; which < 2; ++which) {
                Mesh<D> m = mesh(v, which ? many : one);
                Fields f(g);
                Work w;
                deposit_current(m, g, false, f, w);
                deposit_density(m, g, false, f, w);
                for (int c = 0; c < 3; ++c) result[which][c] = host(f.j[c]);
                result[which][3] = host(f.rho);
                worst = worse(worst, flagged(w));
            }
            for (int c = 0; c < 4; ++c)
                for (std::size_t n = 0; n < result[0][c].size(); ++n)
                    worst = worse(worst, std::fabs(result[0][c][n] - result[1][c][n]) * (c < 3 ? g.dt : 1));
        }
    char what[96];
    std::snprintf(what, sizeof what, "%dD order %d: element = sum of its sub-elements", D, N);
    expect(worst < tol, what, worst);
}
// 3. Continuity: rho_new - rho_old + dt div J = 0 for random elements of every size.
template<int D, int N> void continuity(bool curved = false) {
    const Grid g = grid(D, 32, N);
    double worst = 0;
    for (double size : {.1, .8, 2.5}) {
        const double before = worst;
        for (int trial = 0; trial < 6; ++trial) {
            std::vector<Vertex<D>> v;
            std::vector<Element<D>> e;
            random_element<D>(size, v, e);
            Mesh<D> m = mesh(v, e);
            if (curved) {  // a midpoint off the chord makes the element curved; it moves less than the deposit allows
                double x[3], y[3];
                for (int d = 0; d < D; ++d) {
                    x[d] = (coord(v[0].x, d) + coord(v[1].x, d)) / 2 + uniform(-.3, .3) * size;
                    y[d] = x[d] + uniform(-.4, .4);
                }
                v.push_back({at<D>(x), at<D>(y), {.1, .2, -.3}, 0});
                m = mesh(v, e);
                m.mid = thrust::host_vector<int>(1, int(v.size()) - 1);
            }
            Fields f(g);
            Work w;
            deposit_current(m, g, curved, f, w);
            deposit_density(m, g, curved, f, w);
            const std::vector<double> jx = host(f.j[0]), jy = host(f.j[1]), jz = host(f.j[2]), rho1 = host(f.rho);
            thrust::host_vector<Vertex<D>> old = m.vertex;
            for (auto& p : old) p.x = p.old;
            m.vertex = old;
            deposit_density(m, g, curved, f, w);
            const std::vector<double> rho0 = host(f.rho);
            worst = worse(worst, flagged(w));
            for (int n = 0; n < g.cells(); ++n) {
                const int i[3] = {n % g.n[0], n / g.n[0] % g.n[1], n / (g.n[0] * g.n[1])};
                double div = jx[n] - jx[g.index(i[0] - 1, i[1], i[2])];
                if (D > 1) div += jy[n] - jy[g.index(i[0], i[1] - 1, i[2])];
                if (D > 2) div += jz[n] - jz[g.index(i[0], i[1], i[2] - 1)];
                worst = worse(worst, std::fabs(rho1[n] - rho0[n] + g.dt * div));
            }
        }
        if (worst > before && worst > tol) std::printf("     (size %.2f: %.3g)\n", size, worst);
    }
    char what[96];
    std::snprintf(what, sizeof what, "%dD order %d%s: continuity of charge", D, N, curved ? " curved" : "");
    expect(worst < tol, what, worst);
}
// 4. A PIC particle deposits what a collapsed SIC element (all corners on it) deposits, after one free push.
template<int D, int N> void point_is_collapsed_simplex() {
    const Grid g = grid(D, 16, N);
    double worst = 0;
    for (int trial = 0; trial < 10; ++trial) {
        double x[3];
        for (int d = 0; d < D; ++d) x[d] = uniform(6, 10);
        const Vertex<D> p{at<D>(x), at<D>(x), {Real(uniform(-2, 2)), Real(uniform(-2, 2)), Real(uniform(-2, 2))}, 0};
        Element<D> el{};
        el.q = Real(-.3);
        for (int k = 0; k <= D; ++k) el.v[k] = k;
        Vertex<D> particle = p;
        particle.q = el.q;
        Mesh<D> sic = mesh<D>(std::vector<Vertex<D>>(D + 1, p), {el}), pic = mesh<D>({particle}, {});
        Fields a(g), b(g);
        Work w;
        push_vertices(sic, g, a, w);
        deposit_current(sic, g, false, a, w);
        push_deposit_particles(pic, g, b, w);
        for (int c = 0; c < 3; ++c) {
            const std::vector<double> ja = host(a.j[c]), jb = host(b.j[c]);
            for (std::size_t n = 0; n < ja.size(); ++n) worst = worse(worst, std::fabs(ja[n] - jb[n]) * g.dt);
        }
        worst = worse(worst, flagged(w));
    }
    char what[96];
    std::snprintf(what, sizeof what, "%dD order %d: PIC point = collapsed SIC element", D, N);
    expect(worst < tol, what, worst);
}
// 5. Refinement of a deformed lattice: all edges within the limit, every face shared by exactly two elements (periodic
// box), and the same density at both time levels as before (the subdivision is exact).
template<int D> void refinement(int order) {
    Options o;
    o.g = grid(D, 8, order);
    o.refine = Real(.6);
    Mesh<D> m = initial_mesh<D>(o, Species{});
    thrust::host_vector<Vertex<D>> v = m.vertex;
    for (auto& p : v)
        for (int d = 0; d < D; ++d) {  // a smooth deformation, different at the two time levels
            const double x = coord(p.x, d);
            move(p.x, d, Real(.3 * std::sin(2 * M_PI * x / 8 + d)));
            move(p.old, d, Real(.2 * std::cos(2 * M_PI * x / 8 - d)));
            p.u[d] = Real(uniform(-.5, .5));
        }
    m.vertex = v;
    // Momentum: sum_e q_e (mean of its corners' u). Bisecting an edge at the mean of its ends' u leaves it unchanged.
    auto momentum = [](const Mesh<D>& mesh, double& scale) {
        const thrust::host_vector<Vertex<D>> p = mesh.vertex;
        const thrust::host_vector<Element<D>> el = mesh.element;
        std::vector<double> sum(3, 0);
        scale = 0;
        for (const auto& e : el) {
            for (int k = 0; k <= D; ++k)
                for (int c = 0; c < 3; ++c) sum[c] += double(e.q) * p[e.v[k]].u[c] / (D + 1);
            scale += std::fabs(e.q);
        }
        return sum;
    };
    double scale = 0;
    const std::vector<double> p0 = momentum(m, scale);
    Fields f(o.g);
    Work w;
    std::vector<double> before[2], after[2];
    for (int pass = 0; pass < 2; ++pass) {
        for (int level = 0; level < 2; ++level) {
            Mesh<D> copy = m;
            thrust::host_vector<Vertex<D>> c = copy.vertex;
            if (level == 0) for (auto& p : c) p.x = p.old;
            copy.vertex = c;
            deposit_density(copy, o.g, false, f, w);
            (pass ? after : before)[level] = host(f.rho);
        }
        if (!pass) refine(m, o, w);
    }
    double change = flagged(w);
    for (int level = 0; level < 2; ++level)
        for (std::size_t n = 0; n < before[level].size(); ++n)
            change = worse(change, std::fabs(before[level][n] - after[level][n]));
    // edges and faces on the host
    v = m.vertex;
    const thrust::host_vector<Element<D>> el = m.element;
    double longest = 0;
    std::map<std::vector<int>, int> faces;
    for (const auto& e : el) {
        for (int i = 0; i <= D; ++i)
            for (int j = i + 1; j <= D; ++j) {
                double s = 0;
                for (int d = 0; d < D; ++d) {
                    const double dx = coord(v[e.v[j]].x, d) + 8 * e.img[j][d] - coord(v[e.v[i]].x, d) - 8 * e.img[i][d];
                    s += dx * dx;
                }
                longest = worse(longest, std::sqrt(s));
            }
        for (int skip = 0; skip <= D; ++skip) {  // a face: its vertices with images relative to the lowest one
            std::vector<std::pair<int, std::vector<int>>> corners;
            for (int k = 0; k <= D; ++k)
                if (k != skip) corners.push_back({e.v[k], std::vector<int>(e.img[k], e.img[k] + D)});
            std::sort(corners.begin(), corners.end());
            std::vector<int> key;
            for (auto& c : corners) {
                key.push_back(c.first);
                for (int d = 0; d < D; ++d) key.push_back(c.second[d] - corners[0].second[d]);
            }
            ++faces[key];
        }
    }
    bool conforming = true;
    for (auto& face : faces) conforming = conforming && face.second == 2;
    char what[96];
    std::snprintf(what, sizeof what, "%dD refinement: %zu elements, edges <= limit", D, el.size());
    expect(longest <= o.refine * (1 + 1e-6), what, longest);
    std::snprintf(what, sizeof what, "%dD refinement: every face shared by two elements", D);
    expect(conforming, what);
    std::snprintf(what, sizeof what, "%dD refinement: density unchanged at both time levels", D);
    expect(change < tol, what, change);
    const std::vector<double> p1 = momentum(m, scale);
    double drift = 0;
    for (int c = 0; c < 3; ++c) drift = worse(drift, std::fabs(p1[c] - p0[c]) / scale);
    std::snprintf(what, sizeof what, "%dD refinement: momentum unchanged (midpoints interpolate it)", D);
    expect(drift < tol, what, drift);
}
// 6. A vacuum standing wave E_y = sin(k x): exactly cos(omega n dt) sin(k x), Yee: sin(omega dt / 2) = dt sin(k/2).
template<int D> void vacuum_wave() {
    Grid g = grid(D, 16, 1);
    Fields f(g);
    thrust::host_vector<Real> ey(g.cells());
    const double k = 2 * M_PI / 16;
    for (int n = 0; n < g.cells(); ++n) ey[n] = Real(std::sin(k * (n % 16)));
    f.e[1] = ey;
    const int steps = 200;
    for (int s = 0; s < steps; ++s) advance_fields(f, g);
    const double omega = 2 / double(g.dt) * std::asin(double(g.dt) * std::sin(k / 2));
    ey = f.e[1];
    double worst = 0;
    for (int n = 0; n < g.cells(); ++n)
        worst = worse(worst, std::fabs(ey[n] - std::cos(omega * steps * double(g.dt)) * std::sin(k * (n % 16))));
    char what[96];
    std::snprintf(what, sizeof what, "%dD vacuum wave after %d steps (Yee dispersion)", D, steps);
    expect(worst < (sizeof(Real) == 8 ? 1e-11 : 2e-5), what, worst);
}
// 7. Langmuir oscillation of a cold species: the field energy peaks at t = pi / (2 omega), omega^2 = n Z^2 / M in plasma
// units, with cells of dx c/w_p (density n dx^2 per cell inside); Gauss's law holds. main.cu's step, re-sorted every
// 25 steps.
template<int D> void langmuir(bool pic, int order, bool curved = false, bool refined = false, bool wall = false,
                              double charge = -1, double mass = 1, double dx = 1) {
    Options o;
    o.g = grid(D, 16, order);
    o.g.dt = Real(.05);
    o.dx = dx;
    for (int d = 0; d < D; ++d) o.g.wall[d] = wall;  // walls on every active axis
    o.pic = pic;
    o.curved = curved;
    o.ppc = D == 1 ? 8 : D == 2 ? 4 : 8;
    o.refine = refined ? Real(.7) : Real(0);
    o.sort_every = 25;
    Species s{"s", charge, mass};
    s.perturb = .01;
    const double omega = std::fabs(charge) / std::sqrt(mass), t_peak = M_PI / (2 * omega);  // plasma units
    std::vector<Mesh<D>> m;
    m.push_back(initial_mesh<D>(o, s));
    Fields f(o.g);
    Work w;
    sort_mesh(m[0], o.g);
    deposit_density(m[0], o.g, curved, f, w);
    const Device<double> background = f.rho;
    double best = 0, peak = 0, gauss = 0;
    bool finite = true;
    const int steps = int(1.6 * t_peak / (double(o.g.dt) * dx)) + 1;
    // (the kinetic energy after step n is at t - dt/2, the fields at t: their sum at t takes the mean with step n + 1's)
    std::vector<double> kinetic = {sic::kinetic(m[0], curved).energy}, field = {0};
    for (int n = 1; n <= steps; ++n) {
        step(n, m, o, f, w);
        const double energy = field_energy(f, false) + field_energy(f, true);
        kinetic.push_back(sic::kinetic(m[0], curved).energy);
        field.push_back(energy);
        finite = finite && std::isfinite(energy);
        if (energy > best) { best = energy; peak = n * double(o.g.dt) * dx; }
        deposit_density(m[0], o.g, curved, f, w);
        gauss = worse(gauss, gauss_error(f, o.g, background) / (dx * dx));
    }
    char what[96];
    std::string extra = charge != -1 || mass != 1 ? " Z " + std::to_string(int(charge)) + " M " + std::to_string(int(mass)) : "";
    if (dx != 1) extra += " dx " + std::to_string(dx).substr(0, 4);
    std::snprintf(what, sizeof what, "%dD %s order %d%s%s%s%s: Langmuir peak at pi/2w", D, pic ? "PIC" : "SIC", order,
                  curved ? " curved" : "", refined ? " refined" : "", wall ? " wall" : "", extra.c_str());
    // (walls on every axis shift the peak by a step: the earlier compact code gives the same history to round-off)
    const double step_time = double(o.g.dt) * dx;
    expect(std::fabs(peak - t_peak) < (wall ? .1 : .06) * t_peak / (M_PI / 2) + step_time && finite &&
           w.read_status() == 0, what, peak / t_peak);
    std::snprintf(what, sizeof what, "%dD %s order %d%s%s%s%s: Gauss's law", D, pic ? "PIC" : "SIC", order,
                  curved ? " curved" : "", refined ? " refined" : "", wall ? " wall" : "", extra.c_str());
    expect(gauss < (sizeof(Real) == 8 ? 1e-12 : 1e-5) * std::fabs(charge), what, gauss);
    if (charge != -1 || mass != 1 || dx != 1) {  // (the kinetic energy's 1/qm against the field energy, in these units)
        const std::size_t last = field.size() - 2;
        const double begin = (kinetic[0] + kinetic[1]) / 2 + field[0];
        const double end = (kinetic[last] + kinetic[last + 1]) / 2 + field[last];
        std::snprintf(what, sizeof what, "%dD %s order %d%s: total energy conserved", D, pic ? "PIC" : "SIC", order,
                      extra.c_str());
        expect(std::fabs(end / begin - 1) < .01, what, end / begin - 1);
    }
}
// 8. One 1D element, straight or curved, against an independent oracle: its density, J_x and both transverse currents
// node by node, from an explicit B-spline formula integrated over the material parameter s (split where either time
// level crosses a knot; 6-point Gauss-Legendre on each piece, exact there), with unequal control-point velocities.
double bspline(int N, double x) {  // the cardinal B-spline of degree N, centred at 0
    double s = 0, c = 1, factorial = 1;  // c: (-1)^j binomial(N + 1, j)
    for (int j = 0; j <= N + 1; c = -c * (N + 1 - j) / (j + 1), ++j) {
        const double y = x + (N + 1) / 2.0 - j;
        if (y > 0) s += c * std::pow(y, N);
    }
    for (int i = 2; i <= N; ++i) factorial *= i;
    return s / factorial;
}
double curve(const double (&p)[3], double s) {  // the quadratic through p[0] (s = 0), p[1] (s = 1), p[2] (s = 1/2)
    return p[0] + s * (-3 * p[0] - p[1] + 4 * p[2]) + s * s * (2 * p[0] + 2 * p[1] - 4 * p[2]);
}
void roots(const double (&p)[3], double c, std::vector<double>& out) {  // the s in (0, 1) where curve(p, s) = c
    const double a = 2 * p[0] + 2 * p[1] - 4 * p[2], b = -3 * p[0] - p[1] + 4 * p[2], r = p[0] - c;
    std::vector<double> t;
    if (std::fabs(a) < 1e-13 * (std::fabs(b) + 1)) {
        if (b != 0) t.push_back(-r / b);
    } else if (b * b - 4 * a * r >= 0) {
        const double h = -(b + std::copysign(std::sqrt(b * b - 4 * a * r), b)) / 2;
        t.push_back(h / a);
        if (h != 0) t.push_back(r / h);
    }
    for (double x : t)
        if (x > 0 && x < 1) out.push_back(x);
}
template<int N> void element_oracle(bool curved) {
    const Grid g = grid(1, 32, N);
    const double q = Real(-0.7), node[6] = {.033765242898423986, .16939530676686774, .38069040695840155,
                                      .61930959304159845, .83060469323313226, .96623475710157601};
    const double weight[6] = {.08566224618958517, .18038078652406930, .23395696728634552,
                              .23395696728634552, .18038078652406930, .08566224618958517};
    double worst = 0;
    for (int trial = 0; trial < 8; ++trial) {
        double x[2][3], u[3][3];  // [old, new][end 0, end 1, midpoint]; momenta of the three control points
        const double size = trial < 4 ? .7 : 2.5, x0 = uniform(14, 17), x1 = x0 + (trial % 2 ? -1 : 1) * size;
        x[1][0] = x0;
        x[1][1] = x1;
        x[1][2] = !curved ? (x0 + x1) / 2 : trial == 0 ? x1 + (x1 - x0) / 2 : (x0 + x1) / 2 + uniform(-.3, .3) * size;
        for (int k = 0; k < 3; ++k) x[0][k] = x[1][k] + uniform(-.4, .4);  // (trial 0, curved: a fold)
        if (!curved) x[0][2] = (x[0][0] + x[0][1]) / 2;
        for (int k = 0; k < 3; ++k)
            for (int c = 0; c < 3; ++c) u[k][c] = uniform(-.5, .5);
        std::vector<Vertex<1>> v;
        for (int k = 0; k < (curved ? 3 : 2); ++k) {
            Vertex<1> p{at<1>(&x[1][k]), at<1>(&x[0][k]), {}, 0};
            for (int c = 0; c < 3; ++c) p.u[c] = Real(u[k][c]);
            v.push_back(p);
        }
        Element<1> el{};
        el.v[0] = 0;
        el.v[1] = 1;
        el.q = Real(q);
        Mesh<1> m = mesh<1>(v, {el});
        if (curved) m.mid = thrust::host_vector<int>(1, 2);
        Fields f(g);
        Work w;
        deposit_current(m, g, curved, f, w);
        deposit_density(m, g, curved, f, w);
        const std::vector<double> result[4] = {host(f.rho), host(f.j[0]), host(f.j[1]), host(f.j[2])};
        worst = worse(worst, flagged(w));
        // The oracle: velocities per control point (u / gamma), quadratic in s like the positions (straight: linear).
        double vel[2][3], pos[2][3];  // [y, z] and [old, new] at the control points, as the code holds them
        for (int k = 0; k < int(v.size()); ++k) {
            const double a = v[k].u[0], b = v[k].u[1], c = v[k].u[2], gamma = std::sqrt(1 + a * a + b * b + c * c);
            vel[0][k] = b / gamma;
            vel[1][k] = c / gamma;
            pos[0][k] = coord(v[k].old, 0);
            pos[1][k] = coord(v[k].x, 0);
        }
        if (!curved)  // a straight element: the quadratics through the chord's midpoint are its linear interpolants
            for (int l = 0; l < 2; ++l) {
                vel[l][2] = (vel[l][0] + vel[l][1]) / 2;
                pos[l][2] = (pos[l][0] + pos[l][1]) / 2;
            }
        std::vector<double> cut = {0, 1};
        const double shift = N % 2 ? 0 : .5;  // knots at integers (odd N) or half-integers
        for (int l = 0; l < 2; ++l)
            for (int kn = 0; kn < 40; ++kn) roots(pos[l], kn + shift, cut);
        std::sort(cut.begin(), cut.end());
        std::vector<double> expected[4];
        for (auto& e : expected) e.assign(g.cells(), 0);
        for (std::size_t p = 0; p + 1 < cut.size(); ++p)
            for (int k = 0; k < 6; ++k) {
                const double s = cut[p] + (cut[p + 1] - cut[p]) * node[k], ws = (cut[p + 1] - cut[p]) * weight[k];
                const double xo = curve(pos[0], s), xn = curve(pos[1], s), vy = curve(vel[0], s), vz = curve(vel[1], s);
                double flux = 0;  // minus the change of the weights left of face i + 1/2
                for (int i = 5; i < 28; ++i) {
                    const double so = bspline(N, xo - i), sn = bspline(N, xn - i);
                    flux -= sn - so;
                    expected[0][i] += q * ws * sn;
                    expected[1][i] += q / g.dt * ws * flux;
                    expected[2][i] += q * ws * (so + sn) / 2 * vy;
                    expected[3][i] += q * ws * (so + sn) / 2 * vz;
                }
            }
        for (int c = 0; c < 4; ++c)
            for (int i = 0; i < g.cells(); ++i) worst = worse(worst, std::fabs(result[c][i] - expected[c][i]));
    }
    char what[96];
    std::snprintf(what, sizeof what, "1D order %d%s: density and currents = independent oracle", N,
                  curved ? " curved" : "");
    expect(worst < (sizeof(Real) == 8 ? 1e-10 : 1e-5), what, worst);
}
// 9. The harness itself: a NaN error, a status bit set on the GPU and a NaN field each fail a check; the growth caps.
void harness() {
    Work w;
    w.status[0] = BAD_PIECES;
    Fields f(grid(1, 8, 1));
    f.b[0][3] = Real(NAN);
    const bool caught = !(worse(0, NAN) < tol) && !(worse(NAN, 0) < tol) && !(flagged(w) < tol) &&
                        !std::isfinite(field_energy(f, true));
    expect(caught, "test harness: a NaN or a status bit fails a check");
    bool caps = true;
    auto throws = [](auto&& f) { try { f(); } catch (const std::runtime_error&) { return true; } return false; };
    caps = caps && !throws([] { check_growth(10, max_elements, 10); });
    caps = caps && throws([] { check_growth(10, max_elements + 1, 10); });
    caps = caps && throws([] { check_growth(11, 10, 10); });
    expect(caps, "refinement caps: 2e8 elements, --max-vertices");
    // record()'s check of the diagnostics (main.cu): each of the three, NaN or infinite, at any step, stops the run
    bool diagnostics = !throws([] { check_diagnostics(1e-14, 2.5, 3.5, 50); });
    for (double bad : {double(NAN), double(INFINITY), -double(INFINITY)}) {
        diagnostics = diagnostics && throws([=] { check_diagnostics(bad, 2.5, 3.5, 50); });
        diagnostics = diagnostics && throws([=] { check_diagnostics(1e-14, bad, 3.5, 50); });
        diagnostics = diagnostics && throws([=] { check_diagnostics(1e-14, 2.5, bad, 50); });
    }
    expect(diagnostics, "diagnostics check: a NaN or infinite energy or Gauss residual stops the run");
}
// 10. Weibel (filamentation) growth: cold electron beams u = (0, +-0.2, 0) of density 1/2 each, seeded with opposite
// u_x perturbations at k = 4 (the full code's weibel_1d.in; a B_z seed would go almost entirely into the light wave),
// grow at the cold-fluid rate Gamma^2 (Gamma^2 + k^2 + 1/g^3) = beta^2 k^2 / g (omega_p = 1 for the total): the magnetic
// energy as e^(2 Gamma t), measured from t = 25 to 35 (after the decaying mode, before saturation; the local rate is
// within 0.3% of the theory there).
void weibel_growth(bool pic, int order) {
    Options o;
    o.problem = "weibel";
    o.pic = pic;
    o.g = grid(1, 64, order);
    o.L[0] = M_PI / 2;
    o.dx = o.L[0] / 64;
    o.L[1] = o.L[2] = o.dx;
    o.g.dt = Real(.9);
    o.vpc = 4;
    o.ppc = 16;
    Species beam1{"beam1", -1, 1, .5}, beam2{"beam2", -1, 1, .5};
    beam1.u[1] = .2;
    beam2.u[1] = -.2;
    beam1.perturb = 1e-5;
    beam2.perturb = -1e-5;
    std::vector<Mesh<1>> m;
    m.push_back(initial_mesh<1>(o, beam1));
    m.push_back(initial_mesh<1>(o, beam2));
    Fields f(o.g);
    Work w;
    initial_fields(o, f);
    for (Mesh<1>& x : m) sort_mesh(x, o.g);
    const double step_time = double(o.g.dt) * o.dx;
    const int n1 = int(std::lround(25 / step_time)), n2 = int(std::lround(35 / step_time));
    double energy[2] = {0, 0};
    for (int n = 1; n <= n2; ++n) {
        step(n, m, o, f, w);
        if (n == n1 || n == n2) energy[n == n2] = field_energy(f, true);
    }
    const double k = 4, u = .2, g = std::sqrt(1 + u * u), beta = u / g, a = k * k + 1 / (g * g * g);
    const double theory = std::sqrt((-a + std::sqrt(a * a + 4 * beta * beta * k * k / g)) / 2);
    const double measured = std::log(energy[1] / energy[0]) / (2 * (n2 - n1) * step_time);
    char what[96];
    std::snprintf(what, sizeof what, "1D %s order %d: Weibel growth = cold-fluid rate %.4f", pic ? "PIC" : "SIC", order,
                  theory);
    expect(std::fabs(measured / theory - 1) < .005 && w.read_status() == 0, what, measured / theory);
}
// 11b. The problems' initial magnetic energies against the full code's (its energy.txt at step 0, plasma units): W1D's
// seed modes, P22's Fadeev field between walls, R3D's mirror chains and seed: the same formulas and random draws.
void initial_energies() {
    // (full: the full code's value, from a double-precision run for W1D, single for P22 and R3D: tol, its precision)
    struct Case { const char* name; const char* problem; int dim, n[3]; double L[3]; bool wall_y; Seed seed; double full, tol; };
    const double pi = M_PI;
    const Case cases[3] = {
        {"W1D", "weibel", 1, {128, 1, 1}, {12.8, .1, .1}, false, {1e-6, 64, 0, 1}, 2.07909609614e-10, 1e-10},
        {"P22", "fadeev", 2, {192, 192, 1}, {4 * pi, 4 * pi, 4 * pi / 192}, true, {0, 0, 0, 1}, 179.50638496, 3e-7},
        {"R3D", "fadeev", 3, {64, 64, 16}, {4 * pi, 4 * pi, pi}, false, {1e-2, 32, 3, 1}, 488.168990959, 3e-7}};
    for (const Case& c : cases) {
        Options o;
        o.problem = c.problem;
        o.g.dim = c.dim;
        for (int d = 0; d < 3; ++d) { o.g.n[d] = c.n[d]; o.L[d] = c.L[d]; }
        o.dx = c.L[0] / c.n[0];
        o.g.wall[1] = c.wall_y;
        o.seed = c.seed;
        Fields f(o.g);
        initial_fields(o, f);
        const double energy = field_energy(f, true) * std::pow(o.dx, c.dim - 2);
        char what[96];
        std::snprintf(what, sizeof what, "%s initial magnetic energy = the full code's %.12g", c.name, c.full);
        expect(std::fabs(energy / c.full - 1) < std::max(c.tol, sizeof(Real) == 8 ? 0.0 : 2e-6), what, energy / c.full - 1);
    }
}
// 11. The problems' initial fields: div B = 0 (the random seed modes in 2D and 3D, the Fadeev field between walls in 2D)
// and the Fadeev sheets' current = curl B (kick 0; the two sheet species' density deposit times their drift velocity
// against the discrete curl at the E_z nodes away from the walls, as large as the profile's O(dx^2) error).
double divergence(const Fields& f, const Grid& g) {  // max |div B| over the cells
    thrust::host_vector<Real> b[3] = {f.b[0], f.b[1], f.b[2]};
    double worst = 0;
    for (int n = 0; n < g.cells(); ++n) {
        const int i[3] = {n % g.n[0], n / g.n[0] % g.n[1], n / (g.n[0] * g.n[1])};
        const double div = double(b[0][g.index(i[0] + 1, i[1], i[2])]) - b[0][n] +
                           double(b[1][g.index(i[0], i[1] + 1, i[2])]) - b[1][n] +
                           double(b[2][g.index(i[0], i[1], i[2] + 1)]) - b[2][n];
        worst = worse(worst, std::fabs(div));
    }
    return worst;
}
template<int D> void seed_divergence() {
    Options o;
    o.problem = "weibel";
    o.g = grid(D, 16, 1);
    o.dx = .2;
    for (int d = 0; d < 3; ++d) o.L[d] = d < D ? 16 * o.dx : o.dx;
    o.seed = Seed{1e-3, 32, 4, 1};
    Fields f(o.g);
    initial_fields(o, f);
    const double field = std::sqrt(2 * field_energy(f, true) / o.g.cells());
    char what[96];
    std::snprintf(what, sizeof what, "%dD random B seed: div B = 0 (relative to the rms B)", D);
    expect(divergence(f, o.g) / field < (sizeof(Real) == 8 ? 1e-12 : 1e-5) && field > 0, what,
           divergence(f, o.g) / field);
}
void fadeev_state(bool pic) {
    Options o;
    o.problem = "fadeev";
    o.pic = pic;
    o.fadeev.kick = 0;
    o.g = grid(2, 96, 1);
    o.g.wall[1] = true;
    o.L[0] = o.L[1] = 4 * M_PI;
    o.dx = o.L[0] / 96;
    o.L[2] = o.dx;
    o.vpc = 2;
    o.ppc = 16;
    Fields f(o.g);
    initial_fields(o, f);
    std::vector<double> jz(o.g.cells(), 0.0);
    Work w;
    for (double charge : {-1., 1.}) {  // each species carries half: v = B0 / (2 lam), along z for the positrons
        Species s{"sheet", charge};
        s.chain_y = o.L[1] / 2;
        const double v = o.fadeev.b0 / (2 * o.fadeev.lambda);
        s.u[2] = charge * v / std::sqrt(1 - v * v);
        const Mesh<2> m = initial_mesh<2>(o, s);
        deposit_density(m, o.g, false, f, w);
        const std::vector<double> rho = host(f.rho);
        for (int n = 0; n < o.g.cells(); ++n) jz[n] += rho[n] * charge * v;
    }
    const thrust::host_vector<Real> bx = f.b[0], by = f.b[1];
    double worst = 0, peak = 0;
    for (int n = 0; n < o.g.cells(); ++n) {
        const int i = n % 96, j = n / 96;
        if (j < 2 || j > 94) continue;  // (the wall rows)
        const double curl = double(by[n]) - by[o.g.index(i - 1, j, 0)] - (double(bx[n]) - bx[o.g.index(i, j - 1, 0)]);
        worst = worse(worst, std::fabs(curl - jz[n]));
        peak = std::max(peak, std::fabs(jz[n]));
    }
    char what[96];
    std::snprintf(what, sizeof what, "2D %s Fadeev chain: div B = 0", pic ? "PIC" : "SIC");
    expect(divergence(f, o.g) < (sizeof(Real) == 8 ? 1e-12 : 1e-6) * peak, what, divergence(f, o.g) / peak);
    std::snprintf(what, sizeof what, "2D %s Fadeev chain: sheet current = curl B (to the profile's dx^2)", pic ? "PIC" : "SIC");
    expect(worst / peak < (pic ? .01 : .004) && w.read_status() == 0, what, worst / peak);
}
// 12. Conducting walls (y): an element or particle near a wall deposits (density and current) what it and its mirror
// image (charge -q, y -> -y at both time levels, u_y -> -u_y) deposit in a periodic box twice as tall (the image of y
// is 32 - y there, for both walls), on every node off the wall rows: the fold is the gather's mirror, transposed.
template<int N> void wall_fold(bool pic) {
    Grid walled = grid(2, 16, N), doubled = grid(2, 16, N);
    walled.wall[1] = true;
    doubled.n[1] = 32;
    double worst = 0;
    for (int trial = 0; trial < 8; ++trial) {
        // (a particle starts .3 cells or more from a wall: its push, at most .18 cells, then reflects at none)
        const double y0 = trial % 2 ? uniform(pic ? .3 : .15, .9) : uniform(15.1, pic ? 15.7 : 15.85), x0 = uniform(5, 11);
        std::vector<Vertex<2>> v;
        for (int k = 0; k < (pic ? 1 : 3); ++k) {
            double x[2], y[2];  // new, old (inside the box at both levels)
            x[0] = x0 + (pic ? 0 : uniform(-.6, .6));
            x[1] = std::min(15.95, std::max(.05, y0 + (pic ? 0 : uniform(-.6, .6))));
            y[0] = x[0] + uniform(-.3, .3);
            y[1] = std::min(15.95, std::max(.05, x[1] + uniform(-.3, .3)));
            Vertex<2> p{at<2>(x), at<2>(y), {Real(uniform(-.5, .5)), Real(uniform(-.5, .5)), Real(uniform(-.5, .5))}, 0};
            v.push_back(p);
        }
        auto mirrored = [](Vertex<2> p) {
            const double x[2] = {coord(p.x, 0), 32 - coord(p.x, 1)}, y[2] = {coord(p.old, 0), 32 - coord(p.old, 1)};
            p.x = at<2>(x);
            p.old = at<2>(y);
            p.u[1] = -p.u[1];
            return p;
        };
        const Real q = Real(-.7);
        std::vector<Vertex<2>> both = v;
        for (const auto& p : v) both.push_back(mirrored(p));
        Mesh<2> one, two;
        if (pic) {
            both[0].q = q;
            both[1].q = -q;
            one = mesh<2>({both[0]}, {});
            two = mesh<2>(both, {});
        } else {
            Element<2> e{}, image{};
            e.q = q;
            image.q = -q;
            for (int k = 0; k < 3; ++k) { e.v[k] = k; image.v[k] = 3 + k; }
            one = mesh<2>(v, {e});
            two = mesh<2>(both, {e, image});
        }
        Fields a(walled), b(doubled);
        Work w;
        if (pic) {  // (a free push from the positions given: no fields)
            push_deposit_particles(one, walled, a, w);
            push_deposit_particles(two, doubled, b, w);
        } else {
            deposit_current(one, walled, false, a, w);
            deposit_current(two, doubled, false, b, w);
        }
        deposit_density(one, walled, false, a, w);
        deposit_density(two, doubled, false, b, w);
        worst = worse(worst, flagged(w));
        const std::vector<double> ra = host(a.rho), rb = host(b.rho);
        for (int c = 0; c < 4; ++c) {
            const std::vector<double> fa = c < 3 ? host(a.j[c]) : ra, fb = c < 3 ? host(b.j[c]) : rb;
            for (int j = 0; j < 16; ++j)
                for (int i = 0; i < 16; ++i) {
                    if (c != 1 && j == 0) continue;  // (the wall row; J_y sits half a cell up: all rows are inside)
                    worst = worse(worst, std::fabs(fa[walled.index(i, j, 0)] - fb[doubled.index(i, j, 0)]) *
                                             (c < 3 ? double(walled.dt) : 1));
                }
        }
    }
    char what[96];
    std::snprintf(what, sizeof what, "2D %s order %d walls: deposit = element + mirror image", pic ? "PIC" : "SIC", N);
    expect(worst < tol, what, worst);
}
// 12b. 3D column images on built meshes (column_image; 16^3 cells, pixels of a cell, along z): a thin tetrahedron's
// column stays finite and exact as its volume goes to zero (corners 0, (2, 0, 2), (0, 2, 0), (0, 0, eps), charge 1/6:
// the column at a sample point (u, v) inside its shadow is R / 4, R = 1 - (u + v) / 2), inverted and negative too, to
// 1e-4 at eps = 2^-40 (the ray's ends, about 0.6 cells, differ by 1e-12: double precision's limit); two
// Kuhn cubes of unit density side by side, at u in [s, s + 1] and [s + 1, s + 2] with s = sample_u, so that their shared
// face holds sample points (exactly, in double): those rays count once (1, not 2; in double the left face's rays count
// too, by the tie rule); a tetrahedron moved by 1.6e8 cells along u and back along w (beyond 2^31 pixels at scale 16)
// draws the same image to its corners' rounding; PIC particles moved so draw exactly the same, a slab keeping only
// those inside it.
void columns() {
    Work w;
    View view;
    view.axis[0] = 0, view.axis[1] = 1, view.w = 2, view.n[0] = view.n[1] = 16, view.scale[0] = view.scale[1] = 1;
    view.lo = 0, view.hi = 16;
    const Grid g = grid(3, 16, 1);
    Device<double> image(16 * 16);
    auto draw = [&](const std::vector<Vertex<3>>& v, const std::vector<Element<3>>& e) {
        thrust::fill(image.begin(), image.end(), 0.0);
        column_image(mesh(v, e), g, view, 1.0, ptr(image), w);
        return host(image);
    };
    auto vertex = [](double x, double y, double z, long shift = 0, double q = 0) {  // (shift: whole cells along u, -w)
        const double p[3] = {x, y, z};
        Vertex<3> v{};
        v.x = at<3>(p);
        v.x.c[0] += int(shift), v.x.c[2] -= int(shift);
        v.old = v.x;
        v.q = Real(q);
        return v;
    };
    auto element = [](int a, int b, int c, int d, double q) {
        Element<3> e{};
        e.v[0] = a, e.v[1] = b, e.v[2] = c, e.v[3] = d;
        e.q = Real(q);
        return e;
    };
    double thin = 0;
    const double R = 1 - (sample_u + sample_v) / 2;
    for (double eps : {std::ldexp(1.0, -20), std::ldexp(1.0, -26), std::ldexp(1.0, -40)})
        for (int flip = 0; flip < 2; ++flip) {
            const std::vector<Vertex<3>> v = {vertex(0, 0, 0), vertex(2, 0, 2), vertex(0, 2, 0), vertex(0, 0, eps)};
            const double q = flip ? -1.0 / 6 : 1.0 / 6;
            const std::vector<double> a = draw(v, {flip ? element(0, 2, 1, 3, q) : element(0, 1, 2, 3, q)});
            thin = worse(thin, std::fabs(a[0] / (q > 0 ? R / 4 : -R / 4) - 1));
        }
    expect(thin < 1e-4, "3D images: a thin tetrahedron's column, finite and exact as its volume goes to 0", thin);
    std::vector<Vertex<3>> v;  // two cubes of Kuhn tetrahedra (the paths 0, e_a, e_a + e_b, 1 over the axes' orders)
    std::vector<Element<3>> e;
    const int order[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int cube = 0; cube < 2; ++cube) {
        const int first = int(v.size());
        for (int corner = 0; corner < 8; ++corner)
            v.push_back(vertex(sample_u + cube + (corner & 1), corner >> 1 & 1, corner >> 2 & 1));
        for (const auto& o : order) {
            int path[4] = {0, 0, 0, 0};
            for (int k = 0; k < 3; ++k) path[k + 1] = path[k] | 1 << o[k];
            e.push_back(element(first + path[0], first + path[1], first + path[2], first + path[3], 1.0 / 6));
        }
    }
    const std::vector<double> a = draw(v, e);
    const double faces = std::fabs(a[1] - 1) + (sizeof(Real) == 8 ? std::fabs(a[0] - 1) : 0);
    expect(faces < 1e-6, "3D images: rays on a shared face parallel to them count once", faces);  // (single: q's rounding)
    double moved = 0, total = 0;  // a tetrahedron, then moved by whole periods
    std::vector<double> image0;
    for (long shift : {0L, 160000000L}) {
        const std::vector<Vertex<3>> t = {vertex(.2, .3, .1, shift), vertex(2.9, .4, .6, shift),
                                          vertex(.5, 3.1, .35, shift), vertex(1.1, 1.3, 2.7, shift)};
        const std::vector<double> c = draw(t, {element(0, 1, 2, 3, 1)});
        if (!shift) image0 = c;
        for (std::size_t k = 0; k < c.size(); ++k) moved += std::fabs(c[k] - image0[k]), total += std::fabs(image0[k]);
    }
    moved = total > 0 ? moved / total : NAN;
    // PIC: particles at (3.5, 4.5, 1.5) and (5.5, 6.5, 2.5), slab w in [0, 2): only the first, in pixel (3, 4).
    view.hi = 2;
    double pic = 0;
    for (long shift : {0L, 160000000L}) {
        const std::vector<double> c = draw({vertex(3.5, 4.5, 1.5, shift, 1), vertex(5.5, 6.5, 2.5, shift, 1)}, {});
        for (std::size_t k = 0; k < c.size(); ++k) pic = worse(pic, std::fabs(c[k] - (k == 4 * 16 + 3 ? 1 : 0)));
    }
    expect(moved < 1e-6 && pic == 0 && !w.read_status(), "3D images: whole periods move nothing (SIC beyond 2^31 pixels, PIC)",
           worse(moved, pic));
}
// 12c. Phase-space images on built meshes (phase_image), against references that cut each element with each pixel: a
// 1D segment puts q times the fraction of it inside a pixel there (Liang-Barsky), a 2D element q times the fraction of
// its plane triangle's area (the triangle clipped to the pixel), summed over the periodic copies along a position axis.
// Segments: tilted across the seam (an image +1), pixels smaller than the segment; the reviewers' segment ux 0..4 in
// [0, 1) x [-1, 1) (1/4 of it, no wrapped tail); both tails cut (1/7 inside); an extreme zoom (1e12 pixels, 1e-9 of it
// inside, no hang); moved by 1.6e8 cells (beyond 2^31 pixels) the same image; wholly outside and far away nothing; up to
// an upper wall (into the last pixel). PIC: the momentum range's lower edge in, its upper edge out, just inside, just
// outside, far outside, moved by 1.6e8 cells, a point on an upper wall in the last pixel. Triangles: across one momentum
// edge, across both momentum axes' edges, collapsed onto a segment of the plane (the material density along it 2 (1 - s):
// 3/4 inside), 30000 pixels long (the deposit's 40 levels in 2D), moved by 1.6e8 cells, touching an upper wall. (Walls:
// a point on the upper wall is moved 1e-6 pixels down, into the last pixel: tolerance 1e-6.)
namespace phase_ref {
using P2 = std::array<double, 2>;
std::array<double, 2> interval(const P2& a, const P2& b, double x0, double y0) {  // of t: a + t (b - a) in the pixel
    double t0 = 0, t1 = 1;
    for (int j = 0; j < 2; ++j) {
        const double lo = j ? y0 : x0, d = b[j] - a[j];
        if (d == 0) {
            if (!(a[j] >= lo && a[j] < lo + 1)) return {0, 0};
            continue;
        }
        const double s0 = (lo - a[j]) / d, s1 = (lo + 1 - a[j]) / d;
        t0 = std::max(t0, std::min(s0, s1)), t1 = std::min(t1, std::max(s0, s1));
    }
    return {t0, std::max(t0, t1)};
}
double area(const std::vector<P2>& p) {
    double s = 0;
    for (std::size_t k = 0; k < p.size(); ++k) s += p[k][0] * p[(k + 1) % p.size()][1] - p[(k + 1) % p.size()][0] * p[k][1];
    return std::fabs(s) / 2;
}
double area_in(std::vector<P2> p, double x0, double y0) {  // of polygon p in [x0, x0 + 1] x [y0, y0 + 1]
    for (int side = 0; side < 4 && !p.empty(); ++side) {
        const int j = side / 2;
        const double edge = (j ? y0 : x0) + side % 2;
        auto in = [&](const P2& v) { return side % 2 ? v[j] <= edge : v[j] >= edge; };
        std::vector<P2> q;
        for (std::size_t k = 0; k < p.size(); ++k) {
            const P2 &a = p[k], &b = p[(k + 1) % p.size()];
            if (in(a)) q.push_back(a);
            if (in(a) != in(b)) {
                const double f = (edge - a[j]) / (b[j] - a[j]);
                q.push_back({a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])});
            }
        }
        p = q;
    }
    return p.size() < 3 ? 0 : area(p);
}
// Images of segments and triangles (corners in pixels, unwrapped; periodic: copies along axis 0) in na x nb pixels.
std::vector<double> segments(const std::vector<std::array<P2, 2>>& s, const std::vector<double>& q, int na, int nb,
                             bool periodic) {
    std::vector<double> r(std::size_t(na) * nb, 0.0);
    for (std::size_t e = 0; e < s.size(); ++e) {
        const double lo = std::min(s[e][0][0], s[e][1][0]), hi = std::max(s[e][0][0], s[e][1][0]);
        for (int k = periodic ? int(std::floor(lo / na)) : 0; k <= (periodic ? int(std::floor(hi / na)) : 0); ++k)
            for (int j = 0; j < nb; ++j)
                for (int i = 0; i < na; ++i) {
                    const auto t = interval({s[e][0][0] - k * na, s[e][0][1]}, {s[e][1][0] - k * na, s[e][1][1]}, i, j);
                    r[std::size_t(j) * na + i] += q[e] * (t[1] - t[0]);
                }
    }
    return r;
}
std::vector<double> triangles(const std::vector<std::array<P2, 3>>& t, const std::vector<double>& q, int na, int nb,
                              bool periodic) {
    std::vector<double> r(std::size_t(na) * nb, 0.0);
    for (std::size_t e = 0; e < t.size(); ++e) {
        const std::vector<P2> p(t[e].begin(), t[e].end());
        double lo[2] = {INFINITY, INFINITY}, hi[2] = {-INFINITY, -INFINITY};
        for (const P2& c : p)
            for (int j = 0; j < 2; ++j) lo[j] = std::min(lo[j], c[j]), hi[j] = std::max(hi[j], c[j]);
        for (int k = periodic ? int(std::floor(lo[0] / na)) : 0; k <= (periodic ? int(std::floor(hi[0] / na)) : 0); ++k) {
            std::vector<P2> s = p;
            for (P2& c : s) c[0] -= k * na;
            const int i0 = std::max(0, int(std::floor(lo[0] - k * na))), i1 = std::min(na - 1, int(std::floor(hi[0] - k * na)));
            const int j0 = std::max(0, int(std::floor(lo[1]))), j1 = std::min(nb - 1, int(std::floor(hi[1])));
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) r[std::size_t(j) * na + i] += q[e] * area_in(s, i, j) / area(p);
        }
    }
    return r;
}
}  // namespace phase_ref
void phase_images() {
    using phase_ref::P2;
    Work w;
    const double ptol = sizeof(Real) == 8 ? 1e-9 : 3e-5;  // relative to an element's charge
    PhaseView view{};
    int size[2] = {0, 0};
    auto set = [&](int a, int b, double lo0, double hi0, double lo1, double hi1, int na, int nb) {
        view.axis[0] = a, view.axis[1] = b, view.lo[0] = lo0, view.hi[0] = hi0, view.lo[1] = lo1, view.hi[1] = hi1;
        size[0] = na, size[1] = nb;
    };
    auto draw = [&](const auto& m, const Grid& g) {
        Device<double> image(std::size_t(size[0]) * size[1], 0.0);
        phase_image(m, g, view, size, ptr(image), w);
        return host(image);
    };
    auto differ = [](const std::vector<double>& a, const std::vector<double>& b, double scale) {
        double worst = a.size() == b.size() ? 0 : NAN;
        for (std::size_t k = 0; k < a.size() && k < b.size(); ++k) worst = worse(worst, std::fabs(a[k] - b[k]) / scale);
        return worst;
    };
    // the pixel coordinate of a built vertex along the view's axis k (positions in cells, unwrapped; momenta from lo)
    auto pixel = [&](const auto& v, int k, const Grid& g, int image = 0) {
        const int a = view.axis[k];
        if (a >= 3) return (double(v.u[a - 3]) - view.lo[k]) * size[k] / (view.hi[k] - view.lo[k]);
        return (coord(v.x, a) + double(image) * g.n[a]) * size[k] / g.n[a];
    };
    auto v1 = [](double x, double ux, double uy = 0, long shift = 0, double q = 0) {
        const double p[1] = {x};
        Vertex<1> v{};
        v.x = at<1>(p);
        v.x.c[0] += int(shift);
        v.old = v.x;
        v.u[0] = Real(ux), v.u[1] = Real(uy), v.q = Real(q);
        return v;
    };
    auto e1 = [](int image, double q) {  // corners 0 and 1 (corner 1's image along x)
        Element<1> e{};
        e.v[0] = 0, e.v[1] = 1, e.img[1][0] = (signed char)image, e.q = Real(q);
        return e;
    };
    const Grid g1 = grid(1, 16, 1);
    Grid w1 = g1;
    w1.wall[0] = true;
    double seg = 0;
    auto segment = [&](const std::vector<Vertex<1>>& v, const Element<1>& e, const Grid& g, double scale = 0) {
        const std::vector<double> got = draw(mesh(v, {e}), g);
        const std::array<P2, 2> s = {P2{pixel(v[0], 0, g), pixel(v[0], 1, g)}, P2{pixel(v[1], 0, g, e.img[1][0]), pixel(v[1], 1, g)}};
        return differ(got, phase_ref::segments({s}, {double(e.q)}, size[0], size[1], view.axis[0] < 3 && !g.wall[0]),
                      scale ? scale : std::fabs(double(e.q)));
    };
    set(0, 3, 0, 0, 0, .5, 64, 50);
    const std::vector<Vertex<1>> tilted = {v1(13.3, .1), v1(1.9, .35)};
    seg = worse(seg, segment(tilted, e1(1, -.8), g1));
    const std::vector<double> before = draw(mesh(tilted, {e1(1, -.8)}), g1);
    set(3, 4, 0, 1, -1, 1, 16, 16);
    const std::vector<double> codex = draw(mesh<1>({v1(5.5, 0), v1(6.5, 4)}, {e1(0, 1)}), g1);
    double in = 0;
    for (double x : codex) in += x;
    seg = worse(seg, worse(segment({v1(5.5, 0), v1(6.5, 4)}, e1(0, 1), g1), std::fabs(in - .25)));
    set(0, 3, 0, 0, 0, 1, 64, 50);
    seg = worse(seg, segment({v1(4.2, -3), v1(9.7, 4)}, e1(0, 1), g1));
    set(0, 3, 0, 0, 0, 1e-6, 64, 1000);
    seg = worse(seg, segment({v1(2.5, 0), v1(3.5, 1e3)}, e1(0, 1), g1, 1e-12));
    const std::vector<double> far = draw(mesh<1>({v1(2.5, 1e6), v1(3.5, 2e6)}, {e1(0, 1)}), g1);
    for (double x : far) seg = worse(seg, std::fabs(x));
    set(0, 3, 0, 0, 0, .5, 64, 50);
    const std::vector<double> moved = draw(mesh<1>({v1(13.3, .1, 0, 160000000), v1(1.9, .35, 0, 160000000)}, {e1(1, -.8)}), g1);
    const double shifted = differ(moved, before, .8);
    const double wall = segment({v1(12, .25), v1(16, .25)}, e1(0, .5), w1);
    expect(seg < ptol && shifted < 1e-6 && wall < 1e-6 && !w.read_status(),
           "phase images: 1D segments cut exactly at the momentum range", worse(seg, worse(shifted, wall)));
    // PIC (1D): expected pixels (13, 0) 1, (13, 49) 2, (13, 25) 4; walls: (63, 25) 8
    set(0, 3, 0, 0, 0, 1, 64, 50);
    const std::vector<double> points = draw(mesh<1>({v1(3.3, 0, 0, 0, 1), v1(3.3, 1, 0, 0, 10), v1(3.3, 1 - 1e-6, 0, 0, 2),
                                                     v1(3.3, -1e-6, 0, 0, 20), v1(3.3, .5, 0, 160000000, 4),
                                                     v1(3.3, 1e30, 0, 0, 40)}, {}), g1);
    const std::vector<double> walled = draw(mesh<1>({v1(16, .5, 0, 0, 8)}, {}), w1);
    double pic = 0;
    for (std::size_t k = 0; k < points.size(); ++k) {
        const double want = k == 13 ? 1 : k == 49 * 64 + 13 ? 2 : k == 25 * 64 + 13 ? 4 : 0;
        pic = worse(pic, std::fabs(points[k] - want) + std::fabs(walled[k] - (k == 25 * 64 + 63 ? 8 : 0)));
    }
    expect(pic == 0 && !w.read_status(), "phase images: PIC particles in their pixels, the range's edges, walls", pic);
    // triangles (2D)
    auto v2 = [](double x, double y, double ux, double uy, double uz, long shift = 0) {
        const double p[2] = {x, y};
        Vertex<2> v{};
        v.x = at<2>(p);
        v.x.c[1] += int(shift);
        v.old = v.x;
        v.u[0] = Real(ux), v.u[1] = Real(uy), v.u[2] = Real(uz);
        return v;
    };
    Element<2> tri{};
    tri.v[0] = 0, tri.v[1] = 1, tri.v[2] = 2, tri.q = Real(.9);
    const Grid g2 = grid(2, 16, 1);
    Grid w2 = g2;
    w2.wall[1] = true;
    auto triangle = [&](const std::vector<Vertex<2>>& v, const Grid& g, double tolerance) {
        const std::vector<double> got = draw(mesh(v, {tri}), g);
        std::array<P2, 3> t;
        for (int k = 0; k < 3; ++k) t[k] = {pixel(v[k], 0, g), pixel(v[k], 1, g)};
        return differ(got, phase_ref::triangles({t}, {double(tri.q)}, size[0], size[1], false), double(tri.q)) / tolerance;
    };
    double cut = 0;
    set(1, 5, 0, 0, 0, 1, 64, 40);
    const std::vector<Vertex<2>> across = {v2(3, 2.2, 0, 0, .2), v2(9, 6.1, 0, 0, .7), v2(5, 11.7, 0, 0, 1.6)};
    cut = worse(cut, triangle(across, g2, ptol));
    const std::vector<double> before2 = draw(mesh(across, {tri}), g2);
    const std::vector<double> moved2 = draw(mesh<2>({v2(3, 2.2, 0, 0, .2, 160000000), v2(9, 6.1, 0, 0, .7, 160000000),
                                                  v2(5, 11.7, 0, 0, 1.6, 160000000)}, {tri}), g2);
    cut = worse(cut, differ(moved2, before2, double(tri.q)) / 1e-6);
    set(3, 4, 0, 1, 0, 1, 40, 30);
    cut = worse(cut, triangle({v2(3, 3, -.3, .5, 0), v2(4, 3, .8, -.2, 0), v2(3, 4, 1.4, 1.3, 0)}, g2, ptol));
    set(3, 1, 0, 1, 0, 0, 65536, 16);
    cut = worse(cut, triangle({v2(1, 3.1, .01, 0, 0), v2(2, 4.2, .47, 0, 0), v2(3, 3.6, .2, 0, 0)}, g2, ptol));
    set(1, 3, 0, 0, 0, 1, 64, 20);
    cut = worse(cut, triangle({v2(3, 13.5, .3, 0, 0), v2(5, 16, .5, 0, 0), v2(7, 15.2, .7, 0, 0)}, w2, 1e-6));
    // collapsed: corners 0 and 2 meet in the plane (x, ux); along the segment s from (8, 10) to (24, 50) pixels the
    // material density is 2 (1 - s); the range ux < 0.3 keeps s < 1/2
    set(0, 3, 0, 0, 0, .3, 64, 30);
    const std::vector<double> flat = draw(mesh<2>({v2(2, 2, .1, 0, 0), v2(6, 2, .5, 0, 0), v2(2, 9, .1, 0, 0)}, {tri}), g2);
    std::vector<double> line(flat.size(), 0.0);
    for (int j = 0; j < 30; ++j)
        for (int i = 0; i < 64; ++i) {
            const auto t = phase_ref::interval({8, 10}, {24, 50}, i, j);
            line[std::size_t(j) * 64 + i] = double(tri.q) * (2 * (t[1] - t[0]) - (t[1] * t[1] - t[0] * t[0]));
        }
    double total = 0;
    for (double x : flat) total += x;
    cut = worse(cut, worse(differ(flat, line, double(tri.q)), std::fabs(total / double(tri.q) - .75)) / ptol);
    expect(cut < 1 && !w.read_status(), "phase images: 2D elements cut exactly at the momentum ranges", cut);
}
// 12d. Momentum and spectra on built meshes (kinetic, spectrum): a triangle and a tetrahedron whose corners move
// differently, a curved 1D element whose midpoint's momentum is off the chord's, PIC particles of unequal charge (one at
// rest: the first bin; one with gamma - 1 > 1e4: the last). Momentum: q / qm times the corners' mean (curved:
// (u0 + u1 + 4 u_mid) / 6); spectrum: each point of the rule (weight a on its corner, (1 - a) / D on the others; curved:
// the quadratic at s = a and 1 - a) in its bin, with q / (D + 1).
void diagnostic_sums() {
    auto bin = [](const double (&u)[3]) {
        const double u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2], e = u2 / (std::sqrt(1 + u2) + 1);
        const int k = e > 0 ? int(std::floor((std::log10(e) - spectrum_log0) * spectrum_per_decade)) : 0;
        return std::min(std::max(k, 0), spectrum_bins - 1);
    };
    double worst = 0;
    auto compare = [&](const auto& m, bool curved, const double (&p)[3], const std::vector<double>& bins) {
        const Kinetic k = kinetic(m, curved);
        const double* got = k.momentum;
        Device<double> h(spectrum_bins, 0.0);
        spectrum(m, curved, 1.0, ptr(h));
        const std::vector<double> b = host(h);
        double scale = 0;
        for (double x : bins) scale += std::fabs(x);
        for (int c = 0; c < 3; ++c) worst = worse(worst, std::fabs(got[c] - p[c]) / (1 + std::fabs(p[c])));
        for (int k = 0; k < spectrum_bins; ++k) worst = worse(worst, std::fabs(b[k] - bins[k]) / scale);
    };
    auto vertex = [](int D, const double (&u)[3], double q = 0) {
        Vertex<3> v{};
        for (int d = 0; d < D; ++d) v.x.c[d] = 3 + d, v.x.f[d] = Real(.25);
        v.old = v.x;
        for (int c = 0; c < 3; ++c) v.u[c] = Real(u[c]);
        v.q = Real(q);
        return v;
    };
    auto real = [](const Vertex<3>& v, int c) { return double(v.u[c]); };
    const double U[4][3] = {{.1, -.4, .3}, {1.2, .5, -.2}, {-.6, 2, .8}, {.4, .4, .4}};
    for (int D = 2; D <= 3; ++D) {  // a triangle (q -0.7, electrons) and a tetrahedron (q 0.3, positrons)
        const double q = D == 2 ? -.7 : .3, qm = D == 2 ? -1 : 1, a = rule_a(D), b = (1 - a) / D;
        std::vector<Vertex<3>> v;
        for (int k = 0; k <= D; ++k) v.push_back(vertex(D, U[k]));
        double p[3] = {0, 0, 0};
        std::vector<double> bins(spectrum_bins, 0.0);
        for (int r = 0; r <= D; ++r) {
            double u[3] = {0, 0, 0};
            for (int c = 0; c < 3; ++c) {
                for (int k = 0; k <= D; ++k) u[c] += (k == r ? a : b) * real(v[k], c);
                p[c] += q / qm * real(v[r], c) / (D + 1);
            }
            bins[bin(u)] += double(Real(q)) / (D + 1);
        }
        if (D == 2) {
            std::vector<Vertex<2>> v2(3);
            for (int k = 0; k < 3; ++k) {
                v2[k].x.c[0] = v[k].x.c[0], v2[k].x.c[1] = v[k].x.c[1], v2[k].x.f[0] = v2[k].x.f[1] = Real(.25);
                v2[k].old = v2[k].x;
                for (int c = 0; c < 3; ++c) v2[k].u[c] = v[k].u[c];
            }
            Element<2> e{};
            e.v[0] = 0, e.v[1] = 1, e.v[2] = 2, e.q = Real(q);
            Mesh<2> m = mesh(v2, {e});
            m.qm = Real(qm);
            compare(m, false, p, bins);
        } else {
            Element<3> e{};
            e.v[0] = 0, e.v[1] = 1, e.v[2] = 2, e.v[3] = 3, e.q = Real(q);
            Mesh<3> m = mesh(v, {e});
            m.qm = Real(qm);
            compare(m, false, p, bins);
        }
    }
    {  // curved 1D: corners 0, 1, the midpoint 2; q -0.5
        const double u[3][3] = {{.1, .2, 0}, {.5, -.3, .1}, {1.5, .9, -.4}}, q = -.5, a = rule_a(1);
        std::vector<Vertex<1>> v(3);
        for (int k = 0; k < 3; ++k) {
            v[k].x.c[0] = 3 + k, v[k].old = v[k].x;
            for (int c = 0; c < 3; ++c) v[k].u[c] = Real(u[k][c]);
        }
        Element<1> e{};
        e.v[0] = 0, e.v[1] = 1, e.q = Real(q);
        Mesh<1> m = mesh(v, {e});
        m.mid = thrust::host_vector<int>(1, 2);
        double p[3];
        std::vector<double> bins(spectrum_bins, 0.0);
        for (int c = 0; c < 3; ++c) p[c] = q / -1 * (double(v[0].u[c]) + double(v[1].u[c]) + 4 * double(v[2].u[c])) / 6;
        for (double s : {a, 1 - a}) {
            double w[3];
            for (int c = 0; c < 3; ++c) w[c] = quadratic(double(v[0].u[c]), double(v[1].u[c]), double(v[2].u[c]), s);
            bins[bin(w)] += double(Real(q)) / 2;
        }
        compare(m, true, p, bins);
    }
    {  // PIC 1D: charges -1, -2 (at rest), -0.5 (gamma - 1 about 2e4), -1.5
        const double u[4][3] = {{.3, 0, 0}, {0, 0, 0}, {2e4, 0, 0}, {-.7, 1.1, .2}}, q[4] = {-1, -2, -.5, -1.5};
        std::vector<Vertex<1>> v(4);
        double p[3] = {0, 0, 0};
        std::vector<double> bins(spectrum_bins, 0.0);
        for (int k = 0; k < 4; ++k) {
            v[k].x.c[0] = 3 + k, v[k].old = v[k].x, v[k].q = Real(q[k]);
            for (int c = 0; c < 3; ++c) v[k].u[c] = Real(u[k][c]), p[c] += q[k] / -1 * double(v[k].u[c]);
            bins[bin(u[k])] += q[k];
        }
        if (bin(u[1]) != 0 || bin(u[2]) != spectrum_bins - 1) worst = NAN;  // (the fixture: the first and the last bin)
        compare(mesh<1>(v, {}), false, p, bins);
    }
    expect(worst < (sizeof(Real) == 8 ? 1e-12 : 1e-6), "diagnostics: momentum and spectrum of built meshes", worst);
}
// 12e. The deposit's shared bisection (k_pieces' fan): M copies of one long element deposit M times what one does, for
// M = T >> f, f = 0..10 (T: k_pieces' threads; then the fan is f) and M = T + 7 (fan 0, the work queue), density and
// current: a 1D segment of 100 cells (pieces that fit before level fan and after it), one of 3000 cells (12 levels), a
// curved one, a triangle 60 cells across. The depth limit: a segment 2^23 cells long deposits its charge; one 2^25 cells
// long sets BAD_PIECES.
void fan_out() {
    Work w;
    const int T = loop_blocks() * threads;
    double worst = 0;
    auto copies = [&](auto m, int M, const Grid& g, bool curved) {  // rho, then J (current: the element moved)
        using E = typename decltype(m.element)::value_type;
        const thrust::host_vector<E> one = m.element;
        m.element = thrust::host_vector<E>(M, one[0]);
        if (curved) m.mid = thrust::host_vector<int>(M, 2);
        Fields f(g);
        sic::deposit_density(m, g, curved, f.rho, w);
        deposit_current(m, g, curved, f, w);
        std::vector<double> out = host(f.rho);
        for (int c = 0; c < 3; ++c) for (double x : host(f.j[c])) out.push_back(x);
        return out;
    };
    auto check = [&](const auto& m, const Grid& g, bool curved, const std::vector<int>& counts) {
        const std::vector<double> one = copies(m, 1, g, curved);
        double scale = 0;
        for (double x : one) scale = std::max(scale, std::fabs(x));
        for (int M : counts) {
            const std::vector<double> got = copies(m, M, g, curved);
            for (std::size_t k = 0; k < one.size(); ++k) worst = worse(worst, std::fabs(got[k] - M * one[k]) / (M * scale));
        }
    };
    auto segment = [](double length, bool curved) {
        std::vector<Vertex<1>> v(curved ? 3 : 2);
        const double x[3] = {3.2, 3.2 + length, 3.2 + length / 2 + 2.5};
        const double u[3][3] = {{.1, .3, -.2}, {-.2, .1, .4}, {.3, -.2, .1}};
        for (int k = 0; k < int(v.size()); ++k) {
            const double now[1] = {x[k]}, then[1] = {x[k] - .3};
            v[k].x = at<1>(now), v[k].old = at<1>(then);
            for (int c = 0; c < 3; ++c) v[k].u[c] = Real(u[k][c]);
        }
        Element<1> e{};
        e.v[0] = 0, e.v[1] = 1, e.q = Real(-.6);
        return mesh(v, {e});
    };
    std::vector<int> all, deep = {T >> 10, T >> 7};
    for (int f = 0; f <= 10; ++f) all.push_back(T >> f);
    all.push_back(T + 7);
    const Grid g1 = grid(1, 64, 1), g2 = grid(2, 32, 1);
    check(segment(100, false), g1, false, all);
    check(segment(3000, false), g1, false, deep);
    check(segment(100, true), g1, true, {T >> 3});
    std::vector<Vertex<2>> v(3);
    const double x[3][2] = {{2.3, 3.1}, {62.6, 9.4}, {20.1, 41.7}};
    for (int k = 0; k < 3; ++k) {
        const double then[2] = {x[k][0] - .2, x[k][1] + .1};
        v[k].x = at<2>(x[k]), v[k].old = at<2>(then);
        v[k].u[2] = Real(.3 - .2 * k);
    }
    Element<2> e{};
    e.v[0] = 0, e.v[1] = 1, e.v[2] = 2, e.q = Real(.8);
    check(mesh(v, {e}), g2, false, {T >> 9, T >> 5});
    expect(worst < tol && !w.read_status(), "deposit: M copies of a long element deposit M times one (every fan)", worst);
    // the depth limit (1D: 24 levels)
    Fields f(g1);
    double charge[2];
    for (int k = 0; k < 2; ++k) {
        std::vector<Vertex<1>> s(2);
        const double ends[2] = {.5, .5 + std::ldexp(1.0, k ? 25 : 23)};
        for (int j = 0; j < 2; ++j) s[j].x = s[j].old = at<1>(ends + j);
        Element<1> el{};
        el.v[0] = 0, el.v[1] = 1, el.q = 1;
        sic::deposit_density(mesh(s, {el}), g1, false, f.rho, w);
        charge[k] = 0;
        for (double x : host(f.rho)) charge[k] += x;
        thrust::fill(f.rho.begin(), f.rho.end(), 0.0);
        if (!k) charge[0] = w.read_status() ? NAN : std::fabs(charge[0] - 1);
    }
    const bool limit = w.read_status() & BAD_PIECES;
    w.status[0] = 0;
    expect(charge[0] < 1e-6 && limit, "deposit: 2^23 cells deposit their charge, 2^25 set BAD_PIECES", charge[0]);
}
// 13. Density images (the program's renders): a uniform species draws its density in every pixel, flat and smooth,
// for electrons (a signed charge) as for positive species, in sub-step frames too: Weibel beam1 (density 1, charge -1)
// at steps 0 and 1 (frames at 1/2 and 1), and the Fadeev background (0.2) between the walls, edge rows included.
// 3D column images, the mean density along the view axis at the pixel centres: the same for the Weibel beams along
// z, along x, through a slab, and as PIC particles (no centre on a lattice face is counted twice or missed); the R3D
// sheet's image means its total charge over the box; and two slabs that split the box, weighted by their lengths, add up
// to the whole column (the slab's periodic copies).
// The values of a float16 .npy (count -1: not one); shape: its two dimensions (rows, columns), if it has two.
std::vector<float> npy_float16(const std::string& path, long& count, long* shape = nullptr) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    std::vector<float> v;
    count = -1;
    if (!f) return v;
    unsigned char head[10];
    if (std::fread(head, 1, 10, f) == 10) {
        std::string dict(head[8] + 256 * head[9], ' ');
        if (std::fread(&dict[0], 1, dict.size(), f) != dict.size() || dict.find("'<f2'") == std::string::npos) {
            std::fclose(f);
            return v;
        }
        if (shape && std::sscanf(dict.c_str() + dict.find("'shape': (") + 10, "%ld, %ld", shape, shape + 1) != 2)
            shape[0] = shape[1] = -1;
        for (unsigned short bits; std::fread(&bits, 2, 1, f) == 1;) {
            __half_raw raw;
            raw.x = bits;
            v.push_back(__half2float(__half(raw)));
        }
        count = long(v.size());
    }
    std::fclose(f);
    return v;
}
void renders(const std::string& program) {
    const std::string dir = "/tmp/sic_test_render_" + std::to_string(getpid());
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    struct Case { std::string args, frame; double density; };
    const std::string w = "--problem weibel --dim 2 --nx 16 --vpc 1 --steps 1 --render-every 1 --render-substeps 2 "
                          "--render-scale 4 --render beam1 --no-snapshot";
    const std::string w3 = "--problem weibel --dim 3 --nx 8 --vpc 1 --steps 1 --render-every 1 --render-substeps 2 "
                           "--render-scale 2 --render beam1 --no-snapshot";
    const std::vector<Case> cases = {
        {w, "000001/beam1_n_s1of2.npy", 1}, {w, "000001/beam1_n.npy", 1}, {w + " --render-smooth", "000001/beam1_n.npy", 1},
        {"--problem fadeev --dim 2 --nx 32 --steps 0 --render-every 1 --render-scale 4 --render bg_positrons --no-snapshot",
         "000000/bg_positrons_n.npy", .2},
        {w3, "000001/beam1_n_s1of2.npy", 1}, {w3, "000001/beam1_n.npy", 1}, {w3 + " --render-axis x", "000001/beam1_n.npy", 1},
        {w3 + " --render-axis y --render-slab 1.3,2.9", "000001/beam1_n.npy", 1},
        {w3 + " --render-axis z --render-slab 0,1.4", "000001/beam1_n.npy", 1},  // (the copy past the box's top)
        {"--problem weibel --dim 3 --nx 16 --engine pic --ppc 8 --steps 0 --render-every 1 --render-scale 1 --render beam1 "
         "--no-snapshot", "000000/beam1_n.npy", 1}};
    double worst = 0;
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const std::string out = dir + "/run" + std::to_string(i);
        const int status = std::system((program + " " + cases[i].args + " --output " + out + " > " + out + ".log 2>&1").c_str());
        long count = 0;
        const std::vector<float> v = npy_float16(out + "/render/" + cases[i].frame, count);
        if (status != 0 || count < 16 * 16) { worst = NAN; std::printf("     (%s: exit %d, %ld pixels)\n", cases[i].args.c_str(), status, count); continue; }
        for (float x : v) worst = worse(worst, std::fabs(x / cases[i].density - 1));
    }
    expect(worst < 2e-3, "images: a uniform species draws its density in every pixel", worst);
    // A non-uniform image: the Fadeev sheet at step 0 is symmetric under a half turn of the box (its profile and its
    // triangulation), which pins the pixel convention ([k, k + 1) of the box, no wrap); a window is the matching block.
    const double pi = M_PI;
    const std::string f = "--problem fadeev --dim 2 --nx 32 --steps 0 --render-every 1 --render-scale 4 "
                          "--render sheet_positrons --no-snapshot";
    const std::string window = " --render-window " + std::to_string(pi) + "," + std::to_string(3 * pi) + "," +
                               std::to_string(pi) + "," + std::to_string(2 * pi);
    std::vector<float> image[2];
    long count[2] = {0, 0};
    for (int k = 0; k < 2; ++k) {
        const std::string out = dir + "/fadeev" + std::to_string(k);
        std::system((program + " " + f + (k ? window : "") + " --output " + out + " > " + out + ".log 2>&1").c_str());
        image[k] = npy_float16(out + "/render/000000/sheet_positrons_n.npy", count[k]);
    }
    std::filesystem::remove_all(dir);
    double turn = count[0] == 128 * 128 ? 0 : NAN, block = count[1] == 64 * 32 ? 0 : NAN;
    const std::vector<float>& a = image[0];
    for (std::size_t k = 0; !std::isnan(turn) && k < a.size(); ++k)  // (pixel (i, j) against (127 - i, 127 - j))
        turn = worse(turn, std::fabs(a[k] - a[a.size() - 1 - k]) / std::max(1e-3f, std::fabs(a[k])));
    for (int j = 0; !std::isnan(block) && j < 32; ++j)
        for (int i = 0; i < 64; ++i) block = worse(block, std::fabs(image[1][j * 64 + i] - a[(32 + j) * 128 + 32 + i]));
    expect(turn < 2e-3, "images: the Fadeev sheet is symmetric under a half turn", turn);
    expect(block == 0, "images: a window is the matching block of the whole image", block);
    // 3D: R3D sheet 1 (positrons) at step 0 (box 4 pi x 4 pi x pi, 32 x 32 x 8 cells, 2 pixels a cell), along y through
    // both chains: the whole box, the slabs [0, 5) and [5, Ly); along z, the whole image and a window. Every run must end
    // well and its images have the expected shape; run.txt records the view.
    const std::string r = "--problem fadeev --dim 3 --nx 32 --steps 0 --render-every 1 --render-scale 2 --render sheet1_positrons";
    const double ly = 4 * pi, slab[3] = {0, 5, ly};
    std::filesystem::create_directories(dir);
    std::vector<float> column[5];
    long pixels[5], shape[5][2], good = 0;
    const long expected[5][2] = {{16, 64}, {16, 64}, {16, 64}, {64, 64}, {16, 32}};  // (rows: v; columns: u)
    for (int k = 0; k < 5; ++k) {
        const std::string out = dir + "/r3d" + std::to_string(k);
        char view[192] = " --render-axis y";  // (every digit: 4 pi as the program's own box length)
        if (k == 1 || k == 2) std::snprintf(view, sizeof view, " --render-axis y --render-slab %.17g,%.17g", slab[k - 1], slab[k]);
        if (k == 3) std::snprintf(view, sizeof view, " --render-axis z");
        if (k == 4) std::snprintf(view, sizeof view, " --render-axis z --render-window %.17g,%.17g,%.17g,%.17g", pi, 3 * pi, pi, 2 * pi);
        const int status = std::system((program + " " + r + view + " --output " + out + " > " + out + ".log 2>&1").c_str());
        column[k] = npy_float16(out + "/render/000000/sheet1_positrons_n.npy", pixels[k], shape[k]);
        good += status == 0 && shape[k][0] == expected[k][0] && shape[k][1] == expected[k][1];
    }
    double total = 0;
    std::ifstream elements(dir + "/r3d0/elements_sheet1_positrons.csv");
    std::string line;
    std::getline(elements, line);  // (the header)
    while (std::getline(elements, line)) total += std::stod(line.substr(line.find(',') + 1));
    std::ifstream parameters(dir + "/r3d1/run.txt");
    std::string text((std::istreambuf_iterator<char>(parameters)), std::istreambuf_iterator<char>());
    const bool recorded = text.find("render_axis=y\n") != std::string::npos && text.find("render_slab=0 5\n") != std::string::npos;
    std::filesystem::remove_all(dir);
    expect(good == 5 && recorded && total > 0, "3D images: the runs end well, images of the expected shapes, run.txt has the view",
           double(good));
    if (good < 5) return;
    double mean = 0, split = 0, cropped = 0;
    for (float x : column[0]) mean += x;
    mean = mean / pixels[0] / (total / (4 * pi * ly * pi)) - 1;
    for (std::size_t k = 0; k < column[0].size(); ++k)
        split = worse(split, std::fabs((column[1][k] * 5 + column[2][k] * (ly - 5)) / ly - column[0][k]) /
                                 std::max(1e-2f, column[0][k]));
    // Along z: chain 1 lies at y = Ly/4 (row 16 of 64) with its islands' O-points at x = pi, 3 pi (columns 16, 48);
    // the window [pi, 3 pi] x [pi, 2 pi] is the block of columns 16..47, rows 16..31.
    const std::vector<float>& z = column[3];
    int row = 0, left = 0, right = 32;
    for (int j = 0; j < 64; ++j)
        if (z[j * 64 + 20] > z[row * 64 + 20]) row = j;
    for (int i = 0; i < 32; ++i) {
        if (z[row * 64 + i] > z[row * 64 + left]) left = i;
        if (z[row * 64 + 32 + i] > z[row * 64 + right]) right = 32 + i;
    }
    for (int j = 0; j < 16; ++j)
        for (int i = 0; i < 32; ++i) cropped = worse(cropped, std::fabs(column[4][j * 32 + i] - z[(16 + j) * 64 + 16 + i]));
    const bool oriented = std::abs(row - 16) <= 1 && std::abs(left - 16) <= 1 && std::abs(right - 48) <= 1;
    expect(std::fabs(mean) < 1e-3, "3D images: the mean column is the total charge over the box", mean);
    expect(oriented, "3D images: chain 1 and its islands where they are (rows: y, columns: x)", row * 10000 + left * 100 + right);
    expect(cropped == 0, "3D images: a window is the matching block of the whole image", cropped);
    expect(split < 4e-3, "3D images: two slabs that split the box add up to the whole column", split);
}
// 13b. Diagnostics (the program): W1D at step 0 (two cold beams, density 1 each, u = +-(0.1, 0.5, 0), box 12.8): each
// beam's momentum is 12.8 u; its spectrum holds its 12.8 particles in the bin of gamma - 1 = sqrt(1.26) - 1; its (x, ux)
// image over ux in [-0.25, 0.25] (64 x 99 pixels) is one row, 69 or 29, of f = 1 / (0.5 / 99) = 198 in every pixel, and
// integrates to 12.8; its dump has 256 vertices (x, y, z, ux, uy, uz), 256 elements (vertex, image x, y, z per corner)
// and charge -12.8. Weibel 2D (box 6.4^2): the (ux, uy) images integrate to 40.96 and peak at the beams' momenta. P22:
// the (y, uz) image of the sheet positrons, across the walls, integrates to their dumped charge. 3D PIC (box 6.4^3): the
// (x, uz) image and the spectrum hold 262.144 particles.
template<class T> std::vector<T> npy(const std::string& path, const char* descr, std::vector<long>& shape) {
    std::vector<T> v;
    shape.clear();
    std::ifstream f(path, std::ios::binary);
    char head[10];
    if (!f.read(head, 10)) return v;
    std::string dict(std::size_t((unsigned char)head[8]) + 256 * (unsigned char)head[9], ' ');
    if (!f.read(&dict[0], long(dict.size())) || dict.find(std::string("'") + descr + "'") == std::string::npos) return v;
    std::stringstream dims(dict.substr(dict.find("'shape': (") + 10));
    long count = 1;
    for (long n; dims >> n;) {
        shape.push_back(n), count *= n;
        if (dims.peek() == ',') dims.ignore();
    }
    v.resize(std::size_t(count));
    if (!f.read(reinterpret_cast<char*>(v.data()), long(v.size() * sizeof(T)))) v.clear(), shape.clear();
    return v;
}
std::vector<std::vector<double>> csv(const std::string& path, std::vector<std::string>& names) {  // rows after the header
    std::ifstream f(path);
    std::vector<std::vector<double>> rows;
    std::string line;
    names.clear();
    if (std::getline(f, line)) {
        std::stringstream in(line);
        for (std::string t; std::getline(in, t, ',');) names.push_back(t);
    }
    while (std::getline(f, line)) {
        rows.emplace_back();
        std::stringstream in(line);
        for (std::string t; std::getline(in, t, ',');) rows.back().push_back(std::stod(t));
    }
    return rows;
}
void diagnostics(const std::string& program) {
    const std::string dir = "/tmp/sic_test_diag_" + std::to_string(getpid());
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    auto run = [&](const std::string& name, const std::string& args) {
        return std::system((program + " " + args + " --output " + dir + "/" + name + " > " + dir + "/" + name + ".log 2>&1").c_str());
    };
    const double tol = 1e-5;  // (single precision: the element sums' rounding)
    int ok = run("w1", "--problem weibel --dim 1 --steps 0 --spectra-every 1 --dump-every 1 --phase-every 1 "
                       "--phase x,ux,-0.25,0.25 --phase-size 64,99") == 0;
    std::vector<std::string> names;
    const auto history = csv(dir + "/w1/history.csv", names);
    auto column = [&](const std::string& name) {
        for (std::size_t k = 0; k < names.size(); ++k)
            if (names[k] == name && !history.empty() && k < history[0].size()) return history[0][k];
        return double(NAN);
    };
    double worst = 0;
    const double want[4] = {1.28, 6.4, -1.28, -6.4};
    const char* cols[4] = {"px_beam1", "py_beam1", "px_beam2", "py_beam2"};
    for (int k = 0; k < 4; ++k) worst = worse(worst, std::fabs(column(cols[k]) / want[k] - 1));
    expect(ok && worst < tol, "diagnostics: momentum m u per beam (W1D)", worst);
    const auto spectra = csv(dir + "/w1/spectra_beam1.csv", names);
    double total = 0, misplaced = spectra.size() == 1 && names.size() == 2 + spectrum_bins &&
                                  spectra[0].size() == names.size() ? 0 : NAN;
    const double e = std::sqrt(1.26) - 1;
    for (int k = 0; !std::isnan(misplaced) && k < spectrum_bins; ++k) {
        const double count = spectra[0][2 + k], lo = std::stod(names[2 + k]);
        const double hi = k + 1 < spectrum_bins ? std::stod(names[3 + k]) : INFINITY;
        total += count;
        if (!(lo <= e && e < hi)) misplaced += count;
    }
    expect(std::fabs(total / 12.8 - 1) < tol && misplaced == 0, "diagnostics: the spectrum's particles in their bin (W1D)",
           std::fabs(total / 12.8 - 1) + misplaced);
    std::vector<long> shape;
    double image = 0;
    for (int b = 0; b < 2; ++b) {
        const std::vector<float> f = npy<float>(dir + "/w1/phase0/000000/beam" + std::to_string(b + 1) + "_f.npy", "<f4", shape);
        if (shape != std::vector<long>{99, 64}) { image = NAN; break; }
        double sum = 0;
        for (int j = 0; j < 99; ++j)
            for (int i = 0; i < 64; ++i) {
                const double x = f[j * 64 + i];
                sum += x;
                if (j == (b ? 29 : 69) ? std::fabs(x / 198 - 1) > tol : x != 0) image = NAN;
            }
        image = worse(image, std::fabs(sum * (12.8 / 64) * (0.5 / 99) / 12.8 - 1));
    }
    expect(image < tol, "diagnostics: a cold beam's (x, ux) image is its line, f = n / du (W1D)", image);
    std::vector<long> vs, es, cs;
    const std::vector<double> vertices = npy<double>(dir + "/w1/mesh/000000/beam1_vertices.npy", "<f8", vs);
    const std::vector<int> corners = npy<int>(dir + "/w1/mesh/000000/beam1_elements.npy", "<i4", es);
    const std::vector<double> charges = npy<double>(dir + "/w1/mesh/000000/beam1_charges.npy", "<f8", cs);
    double charge = 0;
    for (double q : charges) charge += q;
    const bool dumped = vs == std::vector<long>{256, 6} && es == std::vector<long>{256, 8} && cs == std::vector<long>{256} &&
                        std::fabs(charge / -12.8 - 1) < tol && vertices[3] == Real(.1) && corners[4] == 1;
    expect(dumped, "diagnostics: the mesh dump's shapes and charge (W1D)", charge);
    // Weibel 2D: (ux, uy) over [-1, 1]^2 at 1601 x 1001 (the momenta mid-pixel): totals, and the peaks' pixels.
    ok = run("w2", "--problem weibel --dim 2 --steps 0 --no-snapshot --phase-every 1 --phase ux,uy,-1,1,-1,1 "
                   "--phase-size 1601,1001") == 0;
    double plane = ok ? 0 : NAN;
    for (int b = 0; b < 2 && !std::isnan(plane); ++b) {
        const std::vector<float> f = npy<float>(dir + "/w2/phase0/000000/beam" + std::to_string(b + 1) + "_f.npy", "<f4", shape);
        if (shape != std::vector<long>{1001, 1601}) { plane = NAN; break; }
        std::size_t top = 0;
        double sum = 0;
        for (std::size_t k = 0; k < f.size(); ++k) sum += f[k], top = f[k] > f[top] ? k : top;
        const double ux = -1 + (top % 1601 + .5) * 2 / 1601, uy = -1 + (top / 1601 + .5) * 2 / 1001, sign = b ? -1 : 1;
        plane = worse(plane, std::fabs(sum * (2. / 1601) * (2. / 1001) / 40.96 - 1));
        if (std::fabs(ux - sign * .1) > 1. / 1601 || std::fabs(uy - sign * .5) > 1. / 1001) plane = NAN;
    }
    expect(plane < tol, "diagnostics: (ux, uy) images: the beams' particles at their momenta (Weibel 2D)", plane);
    // P22 (walls in y): the (y, uz) image of the sheet positrons holds their dumped charge; the uniform background's
    // fills every y pixel alike, the first and the last (at the walls) too.
    ok = run("p22", "--problem fadeev --dim 2 --nx 64 --steps 0 --no-snapshot --dump-every 1 --phase-every 1 "
                    "--phase y,uz,-2,2 --phase-size 128,100 --render sheet_positrons,bg_positrons") == 0;
    const std::vector<float> fy = npy<float>(dir + "/p22/phase0/000000/sheet_positrons_f.npy", "<f4", shape);
    const std::vector<double> qy = npy<double>(dir + "/p22/mesh/000000/sheet_positrons_charges.npy", "<f8", cs);
    double inside = 0, dumped_charge = 0;
    for (float x : fy) inside += x;
    for (double q : qy) dumped_charge += q;
    double walls = ok && shape == std::vector<long>{100, 128} ? std::fabs(inside * (4 * M_PI / 128) * (4. / 100) / dumped_charge - 1) : NAN;
    const std::vector<float> fb = npy<float>(dir + "/p22/phase0/000000/bg_positrons_f.npy", "<f4", shape);
    if (shape != std::vector<long>{100, 128}) walls = NAN;
    for (int i = 0; i < 128 && !std::isnan(walls); ++i) {
        double column = 0, first = 0;
        for (int j = 0; j < 100; ++j) column += fb[j * 128 + i], first += fb[j * 128];
        walls = worse(walls, std::fabs(column / first - 1));
    }
    expect(walls < tol, "diagnostics: (y, uz) images across the walls: the species, every row (P22)", walls);
    // 3D PIC: the (x, uz) image and the spectrum.
    ok = run("pic3", "--problem weibel --dim 3 --nx 8 --engine pic --steps 0 --no-snapshot --phase-every 1 "
                     "--phase x,uz,-1,1 --phase-size 32,32 --spectra-every 1") == 0;
    const std::vector<float> fp = npy<float>(dir + "/pic3/phase0/000000/beam1_f.npy", "<f4", shape);
    double particles = 0;
    for (float x : fp) particles += x * (6.4 / 32) * (2. / 32);
    const auto spectrum3 = csv(dir + "/pic3/spectra_beam1.csv", names);
    double counted = 0;
    for (std::size_t k = 2; !spectrum3.empty() && k < spectrum3[0].size(); ++k) counted += spectrum3[0][k];
    const double pic = ok && shape == std::vector<long>{32, 32} ? worse(std::fabs(particles / 262.144 - 1), std::fabs(counted / 262.144 - 1)) : NAN;
    expect(pic < tol, "diagnostics: 3D PIC phase image and spectrum hold the particles", pic);
    // The dumps of a run's last step against its final snapshot (CSV): every row of W1D, Weibel 2D, curved 1D, 1D PIC.
    auto same = [&](const std::string& name, const std::string& species, const std::string& step, bool pic, int D, bool curved) {
        const std::string m = dir + "/" + name + "/mesh/" + step + "/" + species;
        std::vector<long> vs, es, cs;
        const std::vector<double> v = npy<double>(m + "_vertices.npy", "<f8", vs);
        std::vector<std::string> names;
        const auto vc = csv(dir + "/" + name + "/vertices_" + species + ".csv", names);
        const long columns = pic ? 7 : 6;
        double worst = vs == std::vector<long>{long(vc.size()), columns} && !vc.empty() ? 0 : NAN;
        for (std::size_t i = 0; i < vc.size() && !std::isnan(worst); ++i)
            for (long c = 0; c < columns; ++c) worst = worse(worst, std::fabs(v[i * columns + c] - vc[i][1 + c]) / (1 + std::fabs(vc[i][1 + c])));
        if (pic) return worst;
        const std::vector<int> e = npy<int>(m + "_elements.npy", "<i4", es);
        const std::vector<double> q = npy<double>(m + "_charges.npy", "<f8", cs);
        const auto ec = csv(dir + "/" + name + "/elements_" + species + ".csv", names);
        const long width = 4 * (D + 1) + curved;
        if (es != std::vector<long>{long(ec.size()), width} || cs != std::vector<long>{long(ec.size())} || ec.empty()) return double(NAN);
        for (std::size_t i = 0; i < ec.size(); ++i) {
            for (long c = 0; c < width; ++c) if (e[i * width + c] != ec[i][2 + c]) return double(NAN);
            if (std::fabs(q[i] / ec[i][1] - 1) > (sizeof(Real) == 8 ? 1e-15 : 1e-8)) return double(NAN);  // (CSV: Real's digits)
        }
        return worst;
    };
    ok = run("d1", "--problem weibel --dim 1 --steps 3 --dump-every 3") == 0 &&
         run("d2", "--problem weibel --dim 2 --nx 16 --steps 2 --dump-every 2") == 0 &&
         run("dc", "--problem weibel --dim 1 --curved --order 2 --dt .076 --steps 2 --dump-every 2") == 0 &&
         run("dp", "--problem weibel --dim 1 --engine pic --steps 2 --dump-every 1") == 0;
    double rows = ok ? 0 : NAN;
    rows = worse(rows, same("d1", "beam1", "000003", false, 1, false));
    rows = worse(rows, same("d2", "beam2", "000002", false, 2, false));
    rows = worse(rows, same("dc", "beam1", "000002", false, 1, true));
    rows = worse(rows, same("dp", "beam2", "000002", true, 1, false));
    expect(rows < 1e-15, "diagnostics: mesh dumps, every row as the snapshot's (1D, 2D, curved, PIC)", rows);
    // Cadence: every 2 steps of 3, two planes: steps 0 and 2 only.
    ok = run("every", "--problem weibel --dim 1 --steps 3 --no-snapshot --phase-every 2 --phase x,ux,-1,1 "
                      "--phase ux,uz,-1,1,-1,1 --spectra-every 2") == 0;
    bool cadence = ok;
    for (int k = 0; k < 2; ++k)
        for (int n = 0; n < 4; ++n)
            cadence = cadence && std::filesystem::exists(dir + "/every/phase" + std::to_string(k) + "/00000" +
                                                         std::to_string(n) + "/beam2_f.npy") == (n % 2 == 0);
    const auto rows2 = csv(dir + "/every/spectra_beam2.csv", names);
    cadence = cadence && rows2.size() == 2 && rows2[0][0] == 0 && rows2[1][0] == 2;
    expect(cadence, "diagnostics: phase images and spectra every 2 steps, two planes");
    std::filesystem::remove_all(dir);
}
// 14. The program (next to this test binary): every case must end as expected and say why (its log holds the expected
// text; a success, and only a success, prints the summary " s -> "). Covered stages: option checks (bounds, numbers
// beyond a double or a float, species names, the problems' ranges), the initial diagnostics (beam momenta whose squares
// overflow), a GPU status bit after a step (single precision: a seed field so strong that the Boris rotation's t^2
// overflows a float, while every diagnostic, summed in double, stays finite; in double the same run is fine), the
// vertex cap (the run ends normally there) and the initial mesh above it, and a write that fails when a file is
// flushed at its close: a 4 KB file-size limit (SIGXFSZ ignored) against a snapshot file between 4 and 8 KB, the
// stream buffer, so written whole at its close (without the explicit close the run would succeed).
void command_line(const std::string& path) {
    std::string program = "'";  // the path in single quotes for sh (a quote in it: '\'')
    for (char c : path) program += c == '\'' ? std::string("'\\''") : std::string(1, c);
    program += "'";
    const std::string dir = "/tmp/sic_test_cli_" + std::to_string(getpid());
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    const bool single = sizeof(Real) == 4;
    const std::string big = single ? "3e38" : "1e308";
    struct Case { std::string args, expect; bool success, limit; };  // expect: a text of the log
    const std::string w1 = "--problem weibel --dim 1 ", w2 = "--problem weibel --dim 2 ", f2 = "--problem fadeev --dim 2 ";
    const std::vector<Case> cases = {
        {w1 + "--steps 2 --no-snapshot", " s -> ", true, false},
        {w1 + "--curved --dt .095 --steps 1", "curved: 0.76", false, false},
        {w1 + "--curved --order 2 --dt .076 --steps 20 --no-snapshot", " s -> ", true, false},
        {w1 + "--dt 1e-39 --steps 1", "1e-6 <= dt", false, false},
        {w1 + "--dt 2e-7 --steps 1 --no-snapshot", " s -> ", true, false},
        {w1 + "--u 1e400,0,0 --steps 0", "expected a finite number", false, false},
        {w1 + "--u 1e40,0,0 --steps 0 --no-snapshot", single ? "out of range" : " s -> ", !single, false},
        {w1 + "--u " + big + "," + big + ",0 --steps 0", "non-finite diagnostic at step 0", false, false},
        {w1 + "--seed 1e30 --steps 1 --no-snapshot", single ? "non-finite position or momentum" : " s -> ", !single,
         false},
        {w2 + "--order x", "expected a finite number", false, false},
        {"--problem bogus", "weibel or fadeev", false, false},
        {"--problem fadeev --dim 1", "dim 2 or 3", false, false},
        {f2 + "--islands 3", "even, with a kick", false, false},
        {w1 + "--mass nobody=2 --steps 1", "no species name contains nobody", false, false},
        {w1 + "--mass ions=1836 --steps 2 --no-snapshot", "species ions: charge 1, mass 1836", true, false},
        {w1 + "--mass beam1=1e308 --steps 0 --no-snapshot", single ? "out of range" : "a non-finite momentum of beam1",
         false, false},  // (double: the energy 1.6e308 is finite, p_y 6.4e308 is not)
        {f2 + "--nx 32 --mass positrons=25 --steps 2 --no-snapshot", "bg_positrons: charge 1, mass 25", true, false},
        {f2 + "--nx 32 --b0 3 --steps 1", "below c", false, false},
        {w1 + "--render-every 1 --steps 1", "2D and 3D", false, false},
        {"--problem weibel --dim 3 --nx 8 --render-every 1 --render-scale 2 --steps 1 --no-snapshot", " s -> ", true, false},
        {"--problem weibel --dim 3 --nx 8 --render-every 1 --render-scale 2 --render-smooth --steps 0", "2D only", false, false},
        {"--problem weibel --dim 3 --nx 8 --render-axis q --steps 0", "x, y or z", false, false},
        {"--problem weibel --dim 3 --nx 8 --render-every 1 --render-scale 2 --render-slab 1,99 --steps 0", "render-slab", false,
         false},
        {w2 + "--render-axis x --steps 0", "3D only", false, false},
        {w1 + "--phase x,ux,-1,1 --steps 0", "both or neither", false, false},
        {"--problem weibel --dim 3 --nx 8 --phase-every 1 --phase x,ux,-1,1 --steps 0", "SIC in 1D and 2D", false, false},
        {w1 + "--phase-every 1 --phase y,ux,-1,1 --steps 0", "a position axis of the box", false, false},
        {w1 + "--phase-every 1 --phase x,ux --steps 0", "a range lo,hi for each momentum axis", false, false},
        {w1 + "--phase-every 1 --phase x,ux,-1,1,3 --steps 0", "no more", false, false},
        {w2 + "--phase-every 1 --phase x,y --steps 0", "a momentum axis", false, false},
        {w1 + "--phase-every 1 --phase x,ux,1,-1 --steps 0", "lo < hi", false, false},
        {w1 + "--curved --dt .076 --phase-every 1 --phase x,ux,-1,1 --steps 0", "straight elements", false, false},
        {w1 + "--phase-size 8,100 --steps 0", "phase-size", false, false},
        {w1 + "--phase-every 1 --phase ux,ux,0,1,0,1 --steps 0", "two different axes", false, false},
        {w1 + "--phase-every 1 --phase x,ux,0,1e-307 --steps 0", "too wide or too narrow", false, false},
        {w1 + "--phase-every 1 --phase x,ux,-1e308,1e308 --steps 0", single ? "out of range" : "too wide or too narrow",
         false, false},  // (double: the width overflows; single: the ends do)
        {w1 + "--spectra-every 1 --dump-every 1 --phase-every 1 --phase ux,uz,-1,1,-1,1 --steps 1 --no-snapshot", " s -> ",
         true, false},
        {w2 + "--render-every 1 --render-scale .1 --steps 1", "pixels per axis", false, false},
        {w2 + "--engine pic --render-every 1 --render-substeps 2 --steps 1", "PIC: 1", false, false},
        {w2 + "--lx -3 --steps 1", "box lengths", false, false},
        {"--problem fadeev --dim 3 --nx 32 --steps 0 --no-snapshot --charge sheet2_electrons=-2",  // (chain 2: -J_z;
         "species sheet2_electrons: charge -2, mass 1, density 1 x Fadeev profile, u (0, 0, 0.4364", true, false},  // v .4)
        {"--problem fadeev --dim 3 --nx 32 --ly 3 --steps 0", "ly >= 2 pi lambda", false, false},
        {w1 + "--lx 1e-170 --steps 0", "out of range", false, false},
        {w1 + "--engine pic --nx 4 --lx 4e-160 --ppc 4096 --steps 0", "out of range", false, false},  // (q / ppc -> 0)
        {w2 + "--nx 4 --lx 4 --ly 4.0000001 --render-every 1 --render-scale 4 --render-window 0,4,4,4.0000001 --steps 0",
         "render-window", false, false},  // (the box is exactly its cells: 4 x 4)
        {w2 + "--vpc 1 --refine .5 --max-vertices 3000 --steps 3 --no-snapshot", "vertex cap", true, false},
        {w2 + "--vpc 2 --refine 2 --max-vertices 100 --steps 1", "initial mesh", false, false},
        {w1 + "--steps 2 --nx 40", "iostream", false, true},  // a failed close (fields.csv: 6.2 KB, single 4.3)
    };
    int failed = 0;
    auto read = [](const std::string& path) {
        std::FILE* f = std::fopen(path.c_str(), "r");
        std::string text;
        for (int c; f && (c = std::fgetc(f)) != EOF;) text += char(c);
        if (f) std::fclose(f);
        return text;
    };
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const std::string out = dir + "/run" + std::to_string(i), log = out + ".log";
        const std::string limit = cases[i].limit ? "trap '' XFSZ; ulimit -f 8; " : "";
        const int status = std::system((limit + program + " " + cases[i].args + " --output " + out + " > " +
                                        log + " 2>&1").c_str());
        const std::string text = read(log);
        const bool success = cases[i].success;
        const bool ok = (status == 0) == success && text.find(cases[i].expect) != std::string::npos &&
                        (text.find(" s -> ") != std::string::npos) == success;
        if (ok) continue;
        std::printf("     (%s: exit %d, %s)\n", cases[i].args.c_str(), status, text.empty() ? "no log" : text.c_str());
        ++failed;
    }
    const std::string again = program + " --problem weibel --dim 1 --steps 1 --output " + dir + " > " + dir + ".log 2>&1";
    const int status = std::system(again.c_str());
    const std::string text = read(dir + ".log");
    const bool exists = status != 0 && text.find("output directory exists") != std::string::npos &&
                        text.find(" s -> ") == std::string::npos;
    std::filesystem::remove_all(dir);
    std::filesystem::remove(dir + ".log");
    expect(failed == 0 && exists, "command line: every case ends as expected, with its reason", failed);
}
}  // namespace

int main(int, char** argv) {
    std::printf("%s precision\n", sizeof(Real) == 8 ? "double" : "single");
    harness();
    density_moments<1, 0>(); density_moments<1, 1>(); density_moments<1, 2>(); density_moments<1, 3>();
    density_moments<1, 4>(); density_moments<2, 0>(); density_moments<2, 1>(); density_moments<2, 2>();
    density_moments<3, 0>(); density_moments<3, 1>(); density_moments<3, 2>();
    additivity<1, 2>(); additivity<2, 0>(); additivity<2, 1>(); additivity<2, 2>();
    additivity<3, 0>(); additivity<3, 1>(); additivity<3, 2>();
    continuity<1, 1>(); continuity<1, 4>(); continuity<2, 1>(); continuity<2, 2>();
    continuity<3, 1>(); continuity<3, 2>();
    continuity<1, 1>(true); continuity<1, 3>(true);
    point_is_collapsed_simplex<1, 2>(); point_is_collapsed_simplex<2, 1>();
    point_is_collapsed_simplex<3, 1>(); point_is_collapsed_simplex<3, 2>();
    refinement<2>(1); refinement<3>(1); refinement<2>(2); refinement<3>(2);
    vacuum_wave<1>(); vacuum_wave<2>(); vacuum_wave<3>();
    langmuir<1>(false, 1); langmuir<1>(false, 3, true, true); langmuir<1>(true, 2);
    langmuir<1>(false, 2, false, false, true);
    langmuir<2>(false, 1); langmuir<2>(false, 2, false, true); langmuir<2>(true, 1);
    langmuir<3>(false, 1); langmuir<3>(false, 2); langmuir<3>(false, 1, false, true); langmuir<3>(true, 2);
    langmuir<2>(false, 2, false, true, true); langmuir<2>(true, 1, false, false, true);
    langmuir<3>(false, 1, false, true, true); langmuir<3>(true, 2, false, false, true);
    element_oracle<0>(false); element_oracle<1>(false); element_oracle<2>(false); element_oracle<3>(false);
    element_oracle<4>(false); element_oracle<1>(true); element_oracle<2>(true); element_oracle<3>(true);
    element_oracle<4>(true);
    langmuir<1>(false, 1, false, false, false, 1, 4); langmuir<2>(false, 2, false, false, false, -2, 1);
    langmuir<1>(true, 2, false, false, false, 1, 4); langmuir<2>(false, 1, false, true, false, -1, 1, .25);
    langmuir<3>(true, 1, false, false, false, -1, 1, .25);
    weibel_growth(false, 1); weibel_growth(false, 2); weibel_growth(true, 2);
    seed_divergence<2>(); seed_divergence<3>(); fadeev_state(false); fadeev_state(true); initial_energies();
    columns();
    phase_images();
    diagnostic_sums();
    fan_out();
    wall_fold<1>(false); wall_fold<2>(false); wall_fold<2>(true); wall_fold<1>(true);
    const std::string self = argv[0], dir = self.substr(0, self.rfind('/') + 1);
    command_line(dir + (sizeof(Real) == 8 ? "sic" : "sic-single"));
    renders("'" + dir + (sizeof(Real) == 8 ? "sic" : "sic-single") + "'");
    diagnostics("'" + dir + (sizeof(Real) == 8 ? "sic" : "sic-single") + "'");
    std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
    return failures != 0;
}
