// main.cu - the options and the two problems, the time loop and the output (history, field dumps, density images).
#include "sic.hpp"
#include <cuda_fp16.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>
#include <thrust/host_vector.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

using namespace sic;
namespace {
const char* usage =
    "SIC / PIC on the GPU in 1D, 2D, 3D: cold species, plasma units (c/w_p and 1/w_p of density 1; e, m_e).\n"
    "  --problem weibel|fadeev  --dim 1|2|3 (fadeev: 2, 3)\n"
    "  --nx N                         cells along x (the other axes: their box length over the same cell size)\n"
    "  --lx L --ly L --lz L           box, c/w_p (weibel: 12.8 in 1D, 6.4 per axis in 2D and 3D; fadeev: 2 pi lam\n"
    "                                 per island along x, --ly 4 pi lam, --lz pi lam)\n"
    "  --steps N  --dt DT             dt in 1/w_p (default 0.9 dx / sqrt(dim)), dt / dx <= 0.95 / sqrt(dim)\n"
    "  --engine sic|pic  --order N    SIC: 0..4 (1D), 0..2 (2D, 3D); PIC: 1..2\n"
    "  --gather sameshape|galerkin    field interpolation: order N on every axis (default, the paper's runs) or\n"
    "                                 order N - 1 along a field's half-cell axes\n"
    "  --vpc N  --ppc N  --curved     SIC vertices per cell and axis; PIC particles per cell (a square, a cube);\n"
    "                                 quadratic 1D SIC elements\n"
    "  --refine L  --max-vertices N   split SIC edges longer than L cells; the run ends when a species reaches N\n"
    "  --charge NAME=Z  --mass NAME=M charge (e), mass (m_e; inf: immobile) of each species whose name contains NAME\n"
    "  weibel: --u UX,UY,UZ           beam momenta +-u\n"
    "  fadeev: --b0 B --lambda L --eps E --kick K --guide G --bg N --islands N\n"
    "  --seed A  --seed-rng R         random B seed amplitude and generator (weibel 1e-6, fadeev 3D 1e-2; 2D fadeev:\n"
    "                                 none, with --seed 32 modes up to |k| = 3 per axis)\n"
    "  --diag-every N --fields-every N --sort-every N --no-snapshot --output NEW_DIRECTORY\n"
    "  --render-every N --render-substeps K --render-scale S --render-window U0,U1,V0,V1 --render NAME,...\n"
    "  --render-smooth                density images every N steps, K per step, S pixels per cell: 2D exact pixel\n"
    "                                 means (smooth: each element linear between its corners' lumped densities)\n"
    "  --render-axis x|y|z --render-slab W0,W1\n"
    "                                 3D: the mean density along the view axis (default z) through the slab W0..W1\n"
    "                                 (default the box) at the pixel centres; the image's axes: the other two\n"
    "  --phase-every N --phase A,B[,LO,HI...] --phase-size W,H\n"
    "                                 phase-space images every N steps (species: --render), A and B from x y z ux uy uz\n"
    "                                 (positions: the box; a range LO,HI for each momentum); --phase repeats\n"
    "  --spectra-every N --dump-every N  energy spectra (gamma - 1, log bins); vertices and elements (.npy)\n";

struct Override { std::string name; double value; bool mass; };

// The problems; the defaults are the paper's named runs (below). Weibel (W1D, W3D): beams beam1 and beam2 of electrons (density 1
// each, momenta +-u) and immobile ions; a random B seed. Fadeev: cold pair sheets with density n F and a uniform drift
// along z, each species carrying half of its chain's current J_z = (B0/lam) F; 2D (P22): sheet_electrons and
// sheet_positrons in the chain at Ly/2 between walls in y, with a uniform background bg_electrons, bg_positrons;
// 3D (R3D): sheet1_* (y = Ly/4, +J_z) and sheet2_* (3 Ly/4, -J_z), periodic, a random seed, no background.
void set_problem(Options& o, const double* L, int& nx, int& vpc, int& steps, const double* u, double bg,
                 const bool* given) {
    const int D = o.g.dim;
    if (o.problem == "weibel") {
        const double side = D == 1 ? 12.8 : 6.4;
        for (int d = 0; d < D; ++d) o.L[d] = L[d] > 0 ? L[d] : side;
        nx = nx ? nx : D == 1 ? 128 : 32;
        vpc = vpc ? vpc : 2;
        steps = steps >= 0 ? steps : D == 1 ? 500 : 600;
        o.seed = D == 1 ? Seed{1e-6, 64, 0, 1} : Seed{1e-6, 32, 4, 1};
        Species beam1{"beam1"}, beam2{"beam2"}, ions{"ions", 1, INFINITY, 2};
        for (int c = 0; c < 3; ++c) beam1.u[c] = u[c], beam2.u[c] = -u[c];
        o.species = {beam1, beam2, ions};
        return;
    }
    if (o.problem != "fadeev") throw std::runtime_error("problem: weibel or fadeev");
    Fadeev& p = o.fadeev;
    if (D < 2) throw std::runtime_error("fadeev: dim 2 or 3");
    if (!(p.lambda > 0) || !(p.eps >= 0 && p.eps < 1) || p.islands < 1 || p.islands > 1000 ||
        (p.kick != 0 && p.islands % 2))
        throw std::runtime_error("fadeev: lambda > 0, 0 <= eps < 1, islands 1..1000 (even, with a kick)");
    if (L[0] > 0) throw std::runtime_error("fadeev: the box along x is --islands chain periods 2 pi lam");
    o.L[0] = 2 * M_PI * p.lambda * p.islands;
    o.L[1] = L[1] > 0 ? L[1] : 4 * M_PI * p.lambda;
    if (D == 3 && o.L[1] < 2 * M_PI * p.lambda)  // (the mirror chains' tails are summed over 5 images: periodic to 1e-6)
        throw std::runtime_error("fadeev 3D: ly >= 2 pi lambda");
    if (D == 3) o.L[2] = L[2] > 0 ? L[2] : M_PI * p.lambda;
    nx = nx ? nx : (D == 2 ? 96 : 32) * p.islands;
    if (!given[0] && !o.pic) o.refine = Real(D == 2 ? 1 : 2);  // (the decks' refinement and caps)
    if (!given[1]) o.max_vertices = D == 2 ? 20000000 : 5500000;
    vpc = vpc ? vpc : 1;
    steps = steps >= 0 ? steps : D == 2 ? 346 : 600;
    p.bg = std::isnan(bg) ? (D == 2 ? .2 : 0) : bg;
    if (!(p.bg >= 0)) throw std::runtime_error("fadeev: bg >= 0");
    o.g.wall[1] = D == 2;
    if (D == 3) o.seed = Seed{1e-2, 32, 3, 1};
    auto sheet = [&](const std::string& name, double charge, double y) {
        Species s{name, charge};
        s.chain_y = y;
        return s;
    };
    if (D == 2) o.species = {sheet("sheet_electrons", -1, o.L[1] / 2), sheet("sheet_positrons", 1, o.L[1] / 2)};
    else
        o.species = {sheet("sheet1_electrons", -1, o.L[1] / 4), sheet("sheet1_positrons", 1, o.L[1] / 4),
                     sheet("sheet2_electrons", -1, 3 * o.L[1] / 4), sheet("sheet2_positrons", 1, 3 * o.L[1] / 4)};
    if (p.bg > 0)
        for (double charge : {-1., 1.}) o.species.push_back({charge < 0 ? "bg_electrons" : "bg_positrons", charge, 1, p.bg});
}
// After the charge and mass overrides: each sheet species carries half of its chain's current with a uniform drift.
void set_drifts(Options& o) {
    for (Species& s : o.species) {
        if (s.chain_y < 0) continue;
        const double sign = (s.chain_y > o.L[1] / 2 ? -1 : 1) * (s.charge > 0 ? 1 : -1);  // chain 2 (3 Ly/4): -J_z
        const double v = o.fadeev.b0 / (2 * o.fadeev.lambda * std::fabs(s.charge) * s.density);
        if (!(std::fabs(v) < 1)) throw std::runtime_error("fadeev: the sheet drift B0 / (2 lam |q| n) must be below c");
        s.u[2] = sign * v / std::sqrt(1 - v * v);
    }
}

// The image of the whole box and the window's pixels [i0, i1) x [j0, j1) (parse checks them, the Renderer uses them):
// pixel k is [k, k + 1) pixels from the box's corner. The image's axes u, v: the two besides the view axis, in order
// (2D: x, y; the view axis is z).
struct Pixels { int axis[2], n[2], i0, i1, j0, j1; double scale[2]; };  // (scale: pixels per cell)
void image_axes(int view, int (&axis)[2]) { axis[0] = view == 0 ? 1 : 0, axis[1] = view == 2 ? 1 : 2; }
Pixels pixels(const Options& o) {
    Pixels p;
    image_axes(o.render_axis, p.axis);
    for (int d = 0; d < 2; ++d) {
        p.n[d] = int(std::lround(o.g.n[p.axis[d]] * o.render_scale));
        p.scale[d] = double(p.n[d]) / o.g.n[p.axis[d]];
    }
    const double* w = o.render_window;
    auto at = [&](double x, int d) { return int(std::lround(x / o.dx * p.scale[d])); };
    p.i0 = at(w[0], 0), p.i1 = at(w[1], 0), p.j0 = at(w[2], 1), p.j1 = at(w[3], 1);
    return p;
}

Options parse(int argc, char** argv) {
    Options o;
    int nx = 0, vpc = 0, steps = -1;  // (0, -1: the problem's)
    double L[3] = {0, 0, 0}, dt = 0, u[3] = {.1, .5, 0}, bg = NAN, seed = NAN;
    long long seed_rng = -1;
    bool given[2] = {false, false};  // --refine, --max-vertices (else the problem's)
    bool axis_given = false, slab_given = false;  // (3D only)
    std::vector<Override> overrides;
    auto number = [](const std::string& s) {
        std::size_t used = 0;
        double x = NAN;
        try { x = std::stod(s, &used); } catch (const std::exception&) {}  // (not a number, or beyond a double)
        if (used != s.size() || !std::isfinite(x)) throw std::runtime_error("expected a finite number: " + s);
        return x;
    };
    auto real = [&](const std::string& s) {  // (single precision: also finite as a float)
        const double x = number(s);
        if (!std::isfinite(Real(x))) throw std::runtime_error("out of range: " + s);
        return x;
    };
    auto integer = [&](const std::string& s) {
        const double x = number(s);
        if (x != std::floor(x) || x < 0 || x > 2e9) throw std::runtime_error("expected a non-negative integer: " + s);
        return (long long)x;
    };
    auto list = [&](const std::string& s, int count) {
        std::vector<double> v;
        std::stringstream in(s);
        for (std::string item; std::getline(in, item, ',');) v.push_back(real(item));
        if (int(v.size()) != count) throw std::runtime_error("expected " + std::to_string(count) + " numbers: " + s);
        return v;
    };
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") { std::cout << usage; std::exit(0); }
        if (arg == "--curved") { o.curved = true; continue; }
        if (arg == "--no-snapshot") { o.snapshot = false; continue; }
        if (arg == "--render-smooth") { o.render_smooth = true; continue; }
        if (i + 1 == argc) throw std::runtime_error("missing value for " + arg);
        const std::string v = argv[++i];
        if (arg == "--problem") o.problem = v;
        else if (arg == "--dim") o.g.dim = int(integer(v));
        else if (arg == "--nx") nx = int(integer(v));
        else if (arg == "--lx" || arg == "--ly" || arg == "--lz") {
            L[arg[3] - 'x'] = real(v);
            if (!(L[arg[3] - 'x'] > 0)) throw std::runtime_error("box lengths: > 0");
        }
        else if (arg == "--steps") steps = int(integer(v));
        else if (arg == "--dt") dt = real(v);
        else if (arg == "--engine") {
            if (v != "sic" && v != "pic") throw std::runtime_error("engine is sic or pic");
            o.pic = v == "pic";
        } else if (arg == "--order") o.g.order = int(integer(v));
        else if (arg == "--gather") {
            if (v != "sameshape" && v != "galerkin") throw std::runtime_error("gather is sameshape or galerkin");
            o.g.galerkin = v == "galerkin";
        }
        else if (arg == "--vpc") vpc = int(integer(v));
        else if (arg == "--ppc") o.ppc = int(integer(v));
        else if (arg == "--refine") o.refine = Real(real(v)), given[0] = true;
        else if (arg == "--max-vertices") o.max_vertices = std::size_t(integer(v)), given[1] = true;
        else if (arg == "--charge" || arg == "--mass") {
            const std::size_t eq = v.find('=');
            if (eq == std::string::npos || eq == 0) throw std::runtime_error(arg + ": NAME=VALUE");
            const std::string value = v.substr(eq + 1);
            const bool mass = arg == "--mass";
            overrides.push_back({v.substr(0, eq), mass && value == "inf" ? INFINITY : real(value), mass});
        } else if (arg == "--u") { const auto x = list(v, 3); for (int c = 0; c < 3; ++c) u[c] = x[c]; }
        else if (arg == "--b0") o.fadeev.b0 = real(v);
        else if (arg == "--lambda") o.fadeev.lambda = real(v);
        else if (arg == "--eps") o.fadeev.eps = real(v);
        else if (arg == "--kick") o.fadeev.kick = real(v);
        else if (arg == "--guide") o.fadeev.guide = real(v);
        else if (arg == "--bg") bg = real(v);
        else if (arg == "--islands") o.fadeev.islands = int(integer(v));
        else if (arg == "--seed") seed = real(v);
        else if (arg == "--seed-rng") seed_rng = integer(v);
        else if (arg == "--diag-every") o.diag_every = int(integer(v));
        else if (arg == "--fields-every") o.fields_every = int(integer(v));
        else if (arg == "--sort-every") o.sort_every = int(integer(v));
        else if (arg == "--output") o.output = v;
        else if (arg == "--render-every") o.render_every = int(integer(v));
        else if (arg == "--render-substeps") o.render_substeps = int(integer(v));
        else if (arg == "--render-scale") o.render_scale = real(v);
        else if (arg == "--render-window") { const auto x = list(v, 4); for (int k = 0; k < 4; ++k) o.render_window[k] = x[k]; }
        else if (arg == "--render-axis") {
            if (v != "x" && v != "y" && v != "z") throw std::runtime_error("render-axis: x, y or z");
            o.render_axis = v[0] - 'x';
            axis_given = true;
        }
        else if (arg == "--render-slab") { const auto x = list(v, 2); o.render_slab[0] = x[0], o.render_slab[1] = x[1]; slab_given = true; }
        else if (arg == "--render") {
            std::stringstream in(v);
            for (std::string name; std::getline(in, name, ',');) o.render_species.push_back(name);
        } else if (arg == "--spectra-every") o.spectra_every = int(integer(v));
        else if (arg == "--dump-every") o.dump_every = int(integer(v));
        else if (arg == "--phase-every") o.phase_every = int(integer(v));
        else if (arg == "--phase-size") {
            const auto x = list(v, 2);
            for (int k = 0; k < 2; ++k) {
                if (x[k] != std::floor(x[k]) || !(x[k] >= 16 && x[k] <= 65536)) throw std::runtime_error("phase-size: W,H integers 16..65536");
                o.phase_size[k] = int(x[k]);
            }
        } else if (arg == "--phase") {  // A,B and a range for each momentum axis: x,ux,-1,1 or ux,uy,-1,1,-2,2
            static const std::string names[6] = {"x", "y", "z", "ux", "uy", "uz"};
            std::stringstream in(v);
            std::vector<std::string> item;
            for (std::string t; std::getline(in, t, ',');) item.push_back(t);
            PhaseView view{};
            std::size_t next = 2;
            for (int k = 0; k < 2; ++k) {
                view.axis[k] = -1;
                for (int a = 0; a < 6; ++a)
                    if (item.size() > std::size_t(k) && item[k] == names[a]) view.axis[k] = a;
                if (view.axis[k] < 0) throw std::runtime_error("phase: A,B from x y z ux uy uz, then their momentum ranges");
                if (view.axis[k] < 3) continue;
                if (item.size() < next + 2) throw std::runtime_error("phase: a range lo,hi for each momentum axis");
                view.lo[k] = real(item[next]), view.hi[k] = real(item[next + 1]);
                next += 2;
                if (!(view.lo[k] < view.hi[k])) throw std::runtime_error("phase: lo < hi");
            }
            if (item.size() != next) throw std::runtime_error("phase: A,B and a range for each momentum axis, no more");
            if (view.axis[0] < 3 && view.axis[1] < 3) throw std::runtime_error("phase: a momentum axis (else: --render)");
            if (view.axis[0] == view.axis[1]) throw std::runtime_error("phase: two different axes");
            o.phase.push_back(view);
        } else throw std::runtime_error("unknown option " + arg + " (--help)");
    }
    const int D = o.g.dim, N = o.g.order;
    if (D < 1 || D > 3) throw std::runtime_error("dim is 1, 2 or 3");
    set_problem(o, L, nx, vpc, steps, u, bg, given);
    o.vpc = vpc;
    o.steps = steps;
    if (!std::isnan(seed)) {
        if (!o.seed.modes) o.seed = Seed{0, 32, 3, 1};
        o.seed.amplitude = seed;
    }
    if (seed_rng >= 0) o.seed.rng = seed_rng;
    for (const Override& r : overrides) {
        bool found = false;
        for (Species& s : o.species)
            if (s.name.find(r.name) != std::string::npos) {
                (r.mass ? s.mass : s.charge) = r.value;
                found = true;
            }
        if (!found) throw std::runtime_error("no species name contains " + r.name);
    }
    for (const Species& s : o.species)
        if (!(s.charge != 0) || !(s.mass > 0)) throw std::runtime_error("species " + s.name + ": charge != 0, mass > 0");
    set_drifts(o);
    std::vector<Species> moving;  // (an infinite mass: part of the fixed background)
    for (const Species& s : o.species)
        if (std::isfinite(s.mass)) moving.push_back(s);
    o.species = moving;
    if (o.species.empty()) throw std::runtime_error("no mobile species");
    if (o.ppc == 0) o.ppc = D == 3 ? 8 : 4;
    // cubic cells: dx from the box along x, the other axes' cell counts from their lengths
    if (nx < 4 || nx > 65536) throw std::runtime_error("4..65536 cells per active axis");
    o.dx = o.L[0] / nx;
    double cells = 1;
    for (int d = 0; d < 3; ++d) {
        if (d >= D) { o.g.n[d] = 1; o.L[d] = o.dx; continue; }
        const double n = o.L[d] / o.dx;
        o.g.n[d] = int(std::lround(n));
        if (!(std::fabs(n - o.g.n[d]) < 1e-6 * n) || o.g.n[d] < 4 || o.g.n[d] > 65536)
            throw std::runtime_error("4..65536 cells per active axis, cubic cells (box length / cell size an integer)");
        o.L[d] = o.g.n[d] * o.dx;  // (the box, exactly the cells)
        cells *= o.g.n[d];
    }
    // A point of a curved element moves up to 5/4 of its control points' largest move; the deposit takes moves below
    // a cell. (Tiny steps would overflow q / dt in single precision.)
    if (dt == 0) dt = .9 * o.dx / std::sqrt(double(D));
    o.g.dt = Real(dt / o.dx);
    const Real dt_max = Real(.95) / std::sqrt(Real(D)) / (o.curved ? Real(1.25) : Real(1));
    if (!(o.g.dt >= Real(1e-6) && o.g.dt <= dt_max))
        throw std::runtime_error("1e-6 <= dt / dx <= 0.95/sqrt(dim); curved: 0.76");
    if (o.pic ? N < 1 || N > 2 : N < 0 || N > (D == 1 ? 4 : 2))
        throw std::runtime_error("orders: SIC 0..4 in 1D, 0..2 in 2D and 3D; PIC 1..2");
    if (o.curved && (o.pic || D != 1 || N == 0)) throw std::runtime_error("curved elements: 1D SIC, orders 1..4");
    if (o.vpc < 1 || o.vpc > 64) throw std::runtime_error("vpc is 1..64");
    const int side = int(std::lround(std::pow(double(o.ppc), 1.0 / D)));
    if (o.pic && (std::pow(side, D) != o.ppc || o.ppc < 1 || cells * o.ppc > double(max_elements)))
        throw std::runtime_error("ppc: a square in 2D, a cube in 3D; at most 2e8 particles");
    if (!o.pic && cells * std::pow(o.vpc, D) * (D == 3 ? 6 : D) > double(max_elements))
        throw std::runtime_error("at most 2e8 elements");
    if (!(o.refine >= 0) || (o.pic && o.refine > 0)) throw std::runtime_error("refine: a length >= 0, SIC only");
    for (const Species& s : o.species) {  // the charge initial_mesh stores (an element's or a particle's), and q/m
        const double share = o.pic ? 1.0 / o.ppc : std::pow(1.0 / o.vpc, D) / (D == 1 ? 1 : D == 2 ? 2 : 6);
        const Real q = Real(s.charge * s.density * o.dx * o.dx * share), qm = Real(s.charge / s.mass);
        if (!(std::isfinite(q) && q != 0 && std::isfinite(qm) && qm != 0))
            throw std::runtime_error("species " + s.name + ": charge, density, cell size or mass out of range");
    }
    if (o.max_vertices > 100000000) throw std::runtime_error("max-vertices: at most 1e8");
    if ((axis_given || slab_given) && D != 3) throw std::runtime_error("render-axis, render-slab: 3D only");
    if (o.render_every > 0) {
        double* w = o.render_window;
        if (D == 1) throw std::runtime_error("render: 2D and 3D");
        // (PIC keeps no old positions in memory: one frame per step. The box's image: 16 pixels per axis or more,
        // 2^28 in all, so that indices stay ints and the image fits)
        if (o.render_substeps < 1 || o.render_substeps > 16 || (o.pic && (o.render_substeps > 1 || o.render_smooth)))
            throw std::runtime_error("render: substeps 1..16 (PIC: 1, not smooth)");
        if (D == 3 && o.render_smooth) throw std::runtime_error("render-smooth: 2D only");
        int axis[2];
        image_axes(o.render_axis, axis);
        const int a = axis[0], b = axis[1], c = o.render_axis;
        const double px = o.g.n[a] * o.render_scale, py = o.g.n[b] * o.render_scale;  // (before any rounding)
        if (!(px >= 16 && py >= 16 && px * py <= double(1 << 28)))
            throw std::runtime_error("render: the box's image takes 16 pixels per axis or more, 2^28 in all (--render-scale)");
        if (w[0] == 0 && w[1] == 0 && w[2] == 0 && w[3] == 0) w[1] = o.L[a], w[3] = o.L[b];
        const Pixels p = pixels(o);
        if (!(w[0] >= 0 && w[1] > w[0] && w[1] <= o.L[a] && w[2] >= 0 && w[3] > w[2] && w[3] <= o.L[b]) ||
            !(p.i0 >= 0 && p.i0 < p.i1 && p.i1 <= p.n[0] && p.j0 >= 0 && p.j0 < p.j1 && p.j1 <= p.n[1]))
            throw std::runtime_error("render-window: u0 < u1, v0 < v1 inside the box, a pixel or more");
        if (D == 3) {  // (the slab's ends in cells must differ too)
            double* slab = o.render_slab;
            if (slab[0] == 0 && slab[1] == 0) slab[1] = o.L[c];
            if (!(slab[0] >= 0 && slab[1] > slab[0] && slab[1] <= o.L[c] && slab[1] / o.dx > slab[0] / o.dx))
                throw std::runtime_error("render-slab: w0 < w1 inside the box");
        }
    }
    for (const std::string& name : o.render_species) {  // (also the species of the phase-space images)
        bool found = false;
        for (const Species& s : o.species) found = found || s.name == name;
        if (!found) throw std::runtime_error("render: no species " + name);
    }
    if ((o.phase_every > 0) != !o.phase.empty()) throw std::runtime_error("phase-every and phase: both or neither");
    if (o.phase_every > 0) {  // (a 3D element's shadow is not uniform: SIC phase images in 1D and 2D only)
        if (D == 3 && !o.pic) throw std::runtime_error("phase: SIC in 1D and 2D, PIC in any dimension");
        if (o.curved) throw std::runtime_error("phase: straight elements (not --curved)");
        if (double(o.phase_size[0]) * o.phase_size[1] > double(1 << 26)) throw std::runtime_error("phase-size: 2^26 pixels in all");
        for (const PhaseView& view : o.phase)
            for (int k = 0; k < 2; ++k) {
                if (view.axis[k] < 3 && view.axis[k] >= D) throw std::runtime_error("phase: a position axis of the box");
                const double width = view.hi[k] - view.lo[k];  // (a momentum's pixels per unit: finite, as its width)
                if (view.axis[k] >= 3 && !(std::isfinite(width) && std::isfinite(o.phase_size[k] / width)))
                    throw std::runtime_error("phase: a momentum range too wide or too narrow for its pixels");
            }
    }
    return o;
}

template<class T> std::ofstream open(const std::string& path) {
    std::ofstream f;
    f.exceptions(std::ios::failbit | std::ios::badbit);
    f.open(path, std::ios::binary);
    f << std::setprecision(std::numeric_limits<T>::max_digits10);
    return f;
}
// A .npy file (numpy.load): descr '<f8', '<f4', '<f2' (float64, 32, 16) or '<i4' (int32), C order.
std::ofstream npy_file(const std::string& path, const char* descr, const std::vector<long>& shape) {  // its header
    std::string header = std::string("{'descr': '") + descr + "', 'fortran_order': False, 'shape': (";
    for (long n : shape) header += std::to_string(n) + ", ";
    header += "), }";
    while ((10 + header.size() + 1) % 64) header += ' ';
    header += '\n';
    auto f = open<float>(path);
    const std::uint16_t length = std::uint16_t(header.size());
    f.write("\x93NUMPY\x01\x00", 8);
    f.write(reinterpret_cast<const char*>(&length), 2);
    f.write(header.data(), std::streamsize(header.size()));
    return f;
}
void write_npy(const std::string& path, const char* descr, const std::vector<long>& shape, const void* data,
               std::size_t bytes) {
    auto f = npy_file(path, descr, shape);
    f.write(static_cast<const char*>(data), std::streamsize(bytes));
    f.close();
}
std::string step_name(int n) {
    std::ostringstream s;
    s << std::setw(6) << std::setfill('0') << n;
    return s.str();
}

// ---- 2D density images ------------------------------------------------------------------------------------------
// A species' exact density on a grid of pixels: the order-0 deposit (the cell-average of each element) of its elements
// with their vertices at pixel coordinates, at the time step - 1 + j/K (vertices on straight paths between the two time
// levels). Frames as the full code writes them (its studies/render/make_composite.py reads them):
// render/<step>/<species>_n.npy and, for j < K, <species>_n_s<j>of<K>.npy; float16 number density (n0), row 0 at the
// lowest y, the window cropped from the whole box's image.
// The vertices at step - 1 + alpha in pixel coordinates (pixel k: [k, k + 1) of the box), at rest; for smooth shading
// also with u_z of a speed proportional to their lumped density (below).
__global__ void k_scaled(const Vertex<2>* in, Vertex<2>* out, int count, double alpha, double sx, double sy,
                         double top_x, double top_y, const float* lumped_charge, const float* lumped_area, double speed) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= count) return;
    Vertex<2> v = in[k];
    const double s[2] = {sx, sy}, top[2] = {top_x, top_y};  // (top: a point on an upper wall stays in the last pixel)
    for (int d = 0; d < 2; ++d) {
        const double x0 = v.old.c[d] + double(v.old.f[d]), x1 = v.x.c[d] + double(v.x.f[d]);
        place(v.x, d, fmin((x0 + alpha * (x1 - x0)) * s[d] - .5, top[d]));
    }
    v.old = v.x;
    if (lumped_charge) {
        const double z = lumped_area[k] > 0 ? speed * lumped_charge[k] / lumped_area[k] : 0;
        v.u[0] = v.u[1] = 0;
        v.u[2] = Real(z / sqrt(1 - z * z));
    }
    out[k] = v;
}
// Smooth shading (--render-smooth, the full code's): a vertex takes the density of its adjacent elements (a third of
// each one's charge over a third of its area), and an element is drawn with the linear interpolant of its corners'
// densities. That is the order-0 deposit of the out-of-plane current J_z of a copy of the mesh at rest whose elements
// have unit density (charge = area) and whose vertices move along z at speeds proportional to their densities.
__global__ void k_lumped(const Vertex<2>* v, const Element<2>* el, int count, Grid g, double alpha, double sx,
                         double sy, float* lumped_charge, float* lumped_area, Element<2>* unit) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= count) return;
    double x[3][2];
    for (int k = 0; k < 3; ++k)
        for (int d = 0; d < 2; ++d) {
            const Vertex<2>& p = v[el[e].v[k]];
            const double x0 = p.old.c[d] + double(p.old.f[d]), x1 = p.x.c[d] + double(p.x.f[d]);
            x[k][d] = x0 + alpha * (x1 - x0) + double(el[e].img[k][d]) * g.n[d];
        }
    const double a = .5 * fabs((x[1][0] - x[0][0]) * (x[2][1] - x[0][1]) - (x[2][0] - x[0][0]) * (x[1][1] - x[0][1]));
    for (int k = 0; k < 3; ++k) {
        atomicAdd(lumped_charge + el[e].v[k], float(fabs(double(el[e].q)) / 3));
        atomicAdd(lumped_area + el[e].v[k], float(a / 3));
    }
    unit[e] = el[e];
    unit[e].q = Real(a * sx * sy);
}
__global__ void k_frame(const double* image, int nx, int i0, int j0, int w, int h, double scale, unsigned short* out,
                        int* status) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= w * h) return;
    const int i = i0 + k % w, j = j0 + k / w;  // (the window lies inside the image)
    const double v = image[std::size_t(j) * nx + i] * scale;
    if (!isfinite(v)) atomicOr(status, BAD_IMAGE);
    out[k] = __half_as_ushort(__float2half(float(fmin(v, 65504.0))));  // (finite values above float16's largest: it)
}
struct Renderer {
    View view;        // the image's axes, pixels and pixels per cell (2D: the view axis is z), the slab (3D)
    Grid g;           // 2D: the pixel grid of the whole box
    double top[2];    // 2D: the highest pixel coordinate (a walled axis: below the last edge)
    int i0, j0, w, h;
    Device<double> image, unused_current;  // (smooth: the deposit's in-plane currents, zero: its copy is at rest)
    Device<unsigned short> frame;
    Device<Vertex<2>> scratch;             // the vertices in pixel coordinates
    Device<float> lumped_charge, lumped_area;  // smooth: per vertex
    Mesh<2> unit;                              // smooth: the copy at rest, unit-density elements
    explicit Renderer(const Options& o) {
        const Pixels p = pixels(o);  // (parse checked the window and the slab)
        g.dim = 2;
        g.order = 0;
        for (int d = 0; d < 2; ++d) {
            g.n[d] = view.n[d] = p.n[d];
            view.scale[d] = p.scale[d];
            view.axis[d] = p.axis[d];
            top[d] = o.g.wall[p.axis[d]] ? p.n[d] - .5 - 1e-6 : INFINITY;
        }
        view.w = o.render_axis;
        view.lo = o.render_slab[0] / o.dx, view.hi = o.render_slab[1] / o.dx;
        i0 = p.i0, j0 = p.j0, w = p.i1 - p.i0, h = p.j1 - p.j0;
        image.resize(std::size_t(g.n[0]) * g.n[1]);
        frame.resize(std::size_t(w) * h);
    }
    // 3D: the column of charge per pixel area along the slab; the factor to n / n0 divides by the slab's length.
    double columns(const Options& o, const Mesh<3>& m, const Species& s, Work& work, double alpha) {
        column_image(m, o.g, view, alpha, ptr(image), work);
        return view.scale[0] * view.scale[1] / ((view.hi - view.lo) * s.charge * o.dx * o.dx);
    }
    // Each fills image and returns the factor from its values to n / n0. Flat: the species' charge per pixel.
    double flat(const Options& o, Mesh<2>& m, const Species& s, Work& work, double alpha) {
        const int nv = int(m.vertex.size());
        scratch.resize(nv);
        k_scaled<<<blocks(nv), threads>>>(ptr(m.vertex), ptr(scratch), nv, alpha, view.scale[0], view.scale[1], top[0], top[1],
                                          nullptr, nullptr, 0);
        m.vertex.swap(scratch);  // (the deposit reads the mesh's vertices)
        deposit_density(m, g, false, image, work);
        m.vertex.swap(scratch);
        return view.scale[0] * view.scale[1] / (s.charge * o.dx * o.dx);
    }
    // Smooth: J_z of the unit-density copy, the integral of speed x (the linear interpolant of the lumped densities).
    double smooth(const Options& o, const Mesh<2>& m, const Species& s, Work& work, double alpha) {
        const int nv = int(m.vertex.size()), ne = int(m.element.size());
        lumped_charge.assign(nv, 0.f);
        lumped_area.assign(nv, 0.f);
        unit.element.resize(ne);
        k_lumped<<<blocks(ne), threads>>>(ptr(m.vertex), ptr(m.element), ne, o.g, alpha, view.scale[0], view.scale[1],
                                          ptr(lumped_charge), ptr(lumped_area), ptr(unit.element));
        const float *charge = ptr(lumped_charge), *area = ptr(lumped_area);
        auto first = thrust::make_counting_iterator(0);
        const double densest = thrust::transform_reduce(first, first + nv, [=] __device__(int k) -> double {
            return area[k] > 0 ? double(charge[k]) / area[k] : 0.0;
        }, 0.0, thrust::maximum<double>());
        const double speed = densest > 0 ? .5 / densest : 0;  // (speeds up to c/2)
        unit.vertex.resize(nv);
        k_scaled<<<blocks(nv), threads>>>(ptr(m.vertex), ptr(unit.vertex), nv, alpha, view.scale[0], view.scale[1], top[0],
                                          top[1], charge, area, speed);
        unused_current.resize(image.size());
        double* const j[3] = {ptr(unused_current), ptr(unused_current), ptr(image)};
        deposit_current(unit, g, false, j, work);
        return speed > 0 ? 1 / (speed * std::fabs(s.charge) * o.dx * o.dx) : 0;
    }
    template<int D> void write(const Options& o, Mesh<D>& m, const Species& s, Work& work, int step) {
        const int K = o.render_substeps;
        const std::string dir = o.output + "/render/" + step_name(step);
        std::filesystem::create_directories(dir);
        for (int j = 1; j <= K; ++j) {
            thrust::fill(image.begin(), image.end(), 0.0);
            double density;
            if constexpr (D == 3) density = columns(o, m, s, work, double(j) / K);
            else density = o.render_smooth ? smooth(o, m, s, work, double(j) / K) : flat(o, m, s, work, double(j) / K);
            k_frame<<<blocks(frame.size()), threads>>>(ptr(image), g.n[0], i0, j0, w, h, density, ptr(frame),
                                                       ptr(work.status));
            check(cudaGetLastError());
            const thrust::host_vector<unsigned short> host = frame;
            const std::string name = dir + "/" + s.name + "_n" +
                                     (j < K ? "_s" + std::to_string(j) + "of" + std::to_string(K) : "") + ".npy";
            write_npy(name, "<f2", {h, w}, host.data(), host.size() * sizeof(unsigned short));
        }
    }
};

// ---- phase-space images (deposit.cu: phase_image) -------------------------------------------------------------------
// phase<k>/<step>/<species>_f.npy: f(A, B) in particles per unit A per unit B (per unit length along the 3 - D axes the
// box lacks), float32, row 0 at the lowest B; positions over the whole box, momenta over [lo, hi), the momenta at the
// half step before the positions (the leapfrog's).
template<int D> void write_phase(const Options& o, const Mesh<D>& m, const Species& s, Work& w, int step) {
    const int na = o.phase_size[0], nb = o.phase_size[1];
    Device<double> image;
    for (std::size_t k = 0; k < o.phase.size(); ++k) {
        const PhaseView& view = o.phase[k];
        image.assign(std::size_t(na) * nb, 0.0);
        phase_image(m, o.g, view, o.phase_size, ptr(image), w);
        double width[2];  // a pixel, plasma units
        for (int j = 0; j < 2; ++j)
            width[j] = (view.axis[j] < 3 ? o.dx * o.g.n[view.axis[j]] : view.hi[j] - view.lo[j]) / o.phase_size[j];
        const double scale = std::pow(o.dx, D - 2) / s.charge / (width[0] * width[1]);
        const thrust::host_vector<double> h = image;
        std::vector<float> f(h.size());
        for (std::size_t p = 0; p < h.size(); ++p)
            if (!std::isfinite(f[p] = float(h[p] * scale))) throw std::runtime_error("an image value is not finite");
        const std::string dir = o.output + "/phase" + std::to_string(k) + "/" + step_name(step);
        std::filesystem::create_directories(dir);
        write_npy(dir + "/" + s.name + "_f.npy", "<f4", {nb, na}, f.data(), f.size() * sizeof(float));
    }
}

// ---- mesh dumps (--dump-every) --------------------------------------------------------------------------------------
// mesh/<step>/<species>_vertices.npy: float64 rows x, y, z, ux, uy, uz (PIC: and the charge), plasma units, positions
// unwrapped; SIC also <species>_elements.npy: int32 rows of corners (vertex, image x, y, z each; curved: then the
// midpoint), and <species>_charges.npy: float64, an element's charge. As the final snapshot's CSV files. Written in
// blocks of rows (the host holds one block at a time).
template<class T, class Row>
void dump(const std::string& path, const char* descr, const Device<T>& from, int columns, Row row) {
    constexpr std::size_t block = 1 << 20;
    const std::size_t n = from.size();
    std::vector<long> shape = {long(n)};
    if (columns > 1) shape.push_back(columns);
    auto f = npy_file(path, descr, shape);
    thrust::host_vector<T> h;
    std::vector<decltype(row(T{}, std::size_t(0), 0))> rows;
    for (std::size_t first = 0; first < n; first += block) {
        const std::size_t count = std::min(block, n - first);
        h.resize(count);
        thrust::copy(from.begin() + first, from.begin() + first + count, h.begin());
        rows.resize(count * columns);
        for (std::size_t i = 0; i < count; ++i)
            for (int c = 0; c < columns; ++c) rows[i * columns + c] = row(h[i], first + i, c);
        f.write(reinterpret_cast<const char*>(rows.data()), std::streamsize(rows.size() * sizeof(rows[0])));
    }
    f.close();
}
template<int D> void write_dump(const Options& o, const std::vector<Mesh<D>>& meshes, int step) {
    const std::string dir = o.output + "/mesh/" + step_name(step);
    std::filesystem::create_directories(dir);
    const double charge = std::pow(o.dx, D - 2);
    for (std::size_t s = 0; s < meshes.size(); ++s) {
        const Mesh<D>& m = meshes[s];
        const std::string name = dir + "/" + o.species[s].name;
        dump(name + "_vertices.npy", "<f8", m.vertex, o.pic ? 7 : 6, [&](const Vertex<D>& v, std::size_t, int c) {
            return c < 3 ? (c < D ? (double(v.x.c[c]) + double(v.x.f[c])) * o.dx : 0) : c < 6 ? double(v.u[c - 3]) : v.q * charge;
        });
        if (o.pic) continue;
        const thrust::host_vector<int> mid = m.mid;
        dump(name + "_elements.npy", "<i4", m.element, 4 * (D + 1) + (o.curved ? 1 : 0),
             [&](const Element<D>& e, std::size_t i, int c) {
                 const int k = c / 4, d = c % 4 - 1;
                 return k > D ? mid[i] : d < 0 ? e.v[k] : d < D ? int(e.img[k][d]) : 0;
             });
        dump(name + "_charges.npy", "<f8", m.element, 1, [&](const Element<D>& e, std::size_t, int) { return e.q * charge; });
    }
}

// ---- fields in plasma units ---------------------------------------------------------------------------------------
// rho (net: minus the background), J, E, B; inside units times 1/dx^2 (rho, J) or 1/dx (E, B).
struct HostFields { thrust::host_vector<double> rho, j[3]; thrust::host_vector<Real> e[3], b[3]; };
HostFields to_host(const Fields& f, const Device<double>& background) {
    HostFields h;
    h.rho = f.rho;
    const thrust::host_vector<double> bg = background;
    for (std::size_t k = 0; k < h.rho.size(); ++k) h.rho[k] -= bg[k];
    for (int c = 0; c < 3; ++c) { h.j[c] = f.j[c]; h.e[c] = f.e[c]; h.b[c] = f.b[c]; }
    return h;
}
// fields/<step>.npy: float32 [10][nz][ny][nx]: rho, Jx, Jy, Jz, Ex, Ey, Ez, Bx, By, Bz (J of the step before)
void write_fields(const Options& o, const Fields& f, const Device<double>& background, int step) {
    const HostFields h = to_host(f, background);
    const std::size_t cells = h.rho.size();
    const double a = 1 / (o.dx * o.dx), b = 1 / o.dx;
    std::vector<float> out(10 * cells);
    for (std::size_t k = 0; k < cells; ++k) {
        out[k] = float(h.rho[k] * a);
        for (int c = 0; c < 3; ++c) {
            out[(1 + c) * cells + k] = float(h.j[c][k] * a);
            out[(4 + c) * cells + k] = float(h.e[c][k] * b);
            out[(7 + c) * cells + k] = float(h.b[c][k] * b);
        }
    }
    for (float x : out)
        if (!std::isfinite(x)) throw std::runtime_error("a field value beyond float32 in the dump of step " + std::to_string(step));
    std::filesystem::create_directories(o.output + "/fields");
    write_npy(o.output + "/fields/" + step_name(step) + ".npy", "<f4", {10, o.g.n[2], o.g.n[1], o.g.n[0]}, out.data(),
              out.size() * sizeof(float));
}
// The final state: fields.csv, and per species vertices_<name>.csv (SIC vertices or PIC particles: position, momentum)
// and elements_<name>.csv (SIC: charge, corners with their periodic images; curved 1D: the midpoint last), in plasma
// units. Every file is closed explicitly: a failure of the last flush then throws (a destructor would swallow it).
template<int D>
void write_snapshot(const Options& o, const std::vector<Mesh<D>>& meshes, const Fields& f,
                    const Device<double>& background) {
    const HostFields h = to_host(f, background);
    const double a = 1 / (o.dx * o.dx), b = 1 / o.dx, charge = std::pow(o.dx, D - 2);
    auto fields = open<Real>(o.output + "/fields.csv");
    fields << "i,j,k,rho,Jx,Jy,Jz,Ex,Ey,Ez,Bx,By,Bz\n";
    for (int k = 0; k < o.g.n[2]; ++k)
        for (int j = 0; j < o.g.n[1]; ++j)
            for (int i = 0; i < o.g.n[0]; ++i) {
                const int n = o.g.index(i, j, k);
                fields << i << ',' << j << ',' << k << ',' << Real(h.rho[n] * a);
                for (int c = 0; c < 3; ++c) fields << ',' << Real(h.j[c][n] * a);
                for (int c = 0; c < 3; ++c) fields << ',' << Real(h.e[c][n] * b);
                for (int c = 0; c < 3; ++c) fields << ',' << Real(h.b[c][n] * b);
                fields << '\n';
            }
    fields.close();
    for (std::size_t s = 0; s < meshes.size(); ++s) {
        const Mesh<D>& m = meshes[s];
        const thrust::host_vector<Vertex<D>> v = m.vertex;
        const thrust::host_vector<Element<D>> el = m.element;
        const thrust::host_vector<int> mid = m.mid;
        auto vertices = open<double>(o.output + "/vertices_" + o.species[s].name + ".csv");
        vertices << "vertex,x,y,z,ux,uy,uz" << (o.pic ? ",charge\n" : "\n");
        for (std::size_t i = 0; i < v.size(); ++i) {
            vertices << i;
            for (int d = 0; d < 3; ++d) vertices << ',' << (d < D ? (double(v[i].x.c[d]) + double(v[i].x.f[d])) * o.dx : 0.0);
            for (int c = 0; c < 3; ++c) vertices << ',' << v[i].u[c];
            if (o.pic) vertices << ',' << v[i].q * charge;
            vertices << '\n';
        }
        vertices.close();
        if (o.pic) continue;
        auto elements = open<Real>(o.output + "/elements_" + o.species[s].name + ".csv");
        elements << "element,charge";
        for (int k = 0; k <= D; ++k) elements << ",v" << k << ",image_x" << k << ",image_y" << k << ",image_z" << k;
        elements << (o.curved ? ",mid\n" : "\n");
        for (std::size_t i = 0; i < el.size(); ++i) {
            elements << i << ',' << el[i].q * charge;
            for (int k = 0; k <= D; ++k) {
                elements << ',' << el[i].v[k];
                for (int d = 0; d < 3; ++d) elements << ',' << (d < D ? int(el[i].img[k][d]) : 0);
            }
            if (o.curved) elements << ',' << mid[i];
            elements << '\n';
        }
        elements.close();
    }
}

void write_parameters(const Options& o) {
    auto p = open<double>(o.output + "/run.txt");
    const Grid& g = o.g;
    p << "problem=" << o.problem << "\nengine=" << (o.pic ? "pic" : "sic") << "\ndim=" << g.dim << "\norder=" << g.order
      << "\nnx=" << g.n[0] << "\nny=" << g.n[1] << "\nnz=" << g.n[2] << "\nlx=" << o.L[0] << "\nly=" << o.L[1]
      << "\nlz=" << o.L[2] << "\ndx=" << o.dx << "\ndt=" << g.dt * o.dx << "\nsteps=" << o.steps << "\nvpc=" << o.vpc
      << "\nppc=" << o.ppc << "\ncurved=" << o.curved << "\nrefine=" << o.refine << "\nmax_vertices=" << o.max_vertices
      << "\ngather=" << (g.galerkin ? "galerkin" : "sameshape") << "\nwall_x=" << g.wall[0] << "\nwall_y=" << g.wall[1]
      << "\nwall_z=" << g.wall[2]
      << "\nseed=" << o.seed.amplitude << ' ' << o.seed.modes << ' ' << o.seed.kmax << ' ' << o.seed.rng
      << "\ndiag_every=" << o.diag_every << "\nfields_every=" << o.fields_every << "\nsort_every=" << o.sort_every
      << "\nprecision=" << 8 * sizeof(Real) << "\nrender_every=" << o.render_every << "\nrender_substeps="
      << o.render_substeps << "\nrender_scale=" << o.render_scale << "\nrender_smooth=" << o.render_smooth
      << "\nrender_window=" << o.render_window[0] << ' ' << o.render_window[1] << ' ' << o.render_window[2] << ' '
      << o.render_window[3] << "\nrender_axis=" << "xyz"[o.render_axis] << "\nrender_slab=" << o.render_slab[0] << ' '
      << o.render_slab[1] << "\nrender=";
    for (const std::string& r : o.render_species) p << r << ' ';
    p << "\nspectra_every=" << o.spectra_every << "\ndump_every=" << o.dump_every << "\nphase_every=" << o.phase_every
      << "\nphase_size=" << o.phase_size[0] << ' ' << o.phase_size[1] << '\n';
    static const char* axes[6] = {"x", "y", "z", "ux", "uy", "uz"};
    for (const PhaseView& view : o.phase) {
        p << "phase=" << axes[view.axis[0]] << ' ' << axes[view.axis[1]];
        for (int k = 0; k < 2; ++k)
            if (view.axis[k] >= 3) p << ' ' << view.lo[k] << ' ' << view.hi[k];
        p << '\n';
    }
    if (o.problem == "fadeev")
        p << "fadeev=b0 " << o.fadeev.b0 << " lambda " << o.fadeev.lambda << " eps " << o.fadeev.eps << " kick "
          << o.fadeev.kick << " guide " << o.fadeev.guide << " bg " << o.fadeev.bg << " islands " << o.fadeev.islands
          << '\n';
    for (const Species& s : o.species)
        p << "species=" << s.name << " charge " << s.charge << " mass " << s.mass << " density " << s.density << " u "
          << s.u[0] << ' ' << s.u[1] << ' ' << s.u[2] << " chain_y " << s.chain_y << '\n';
    p.close();
}

template<int D> void run(const Options& o) {
    const Grid& g = o.g;
    cudaDeviceProp device{};
    check(cudaGetDeviceProperties(&device, 0));
    std::cout << "GPU: " << device.name << '\n';
    std::vector<Mesh<D>> meshes;
    meshes.reserve(o.species.size());
    for (const Species& s : o.species) {
        meshes.push_back(initial_mesh<D>(o, s));
        std::cout << "species " << s.name << ": charge " << s.charge << ", mass " << s.mass << ", density "
                  << s.density << (s.chain_y >= 0 ? " x Fadeev profile" : "") << ", u (" << s.u[0] << ", " << s.u[1]
                  << ", " << s.u[2] << "), " << meshes.back().vertex.size() << (o.pic ? " particles\n" : " vertices\n");
        if (o.refine > 0 && meshes.back().vertex.size() > o.max_vertices)
            throw std::runtime_error("the initial mesh has more vertices than --max-vertices");
    }
    Fields f(g);
    Work w;
    initial_fields(o, f);
    for (Mesh<D>& m : meshes) sort_mesh(m, g);
    auto density = [&] {
        thrust::fill(f.rho.begin(), f.rho.end(), 0.0);
        for (const Mesh<D>& m : meshes) deposit_density(m, g, o.curved, f.rho, w);
    };
    density();
    const Device<double> background = f.rho;  // the fixed neutralizing background: the initial charge
    if (!std::filesystem::create_directory(o.output)) throw std::runtime_error("output directory exists: " + o.output);
    write_parameters(o);
    auto history = open<double>(o.output + "/history.csv");
    history << "step,time,vertices,elements,kinetic,electric,magnetic,gauss_max";
    for (const Species& s : o.species) history << ",kinetic_" << s.name << ",vertices_" << s.name;
    for (const Species& s : o.species) history << ",px_" << s.name << ",py_" << s.name << ",pz_" << s.name;
    history << '\n';
    const double energy = std::pow(o.dx, D - 2);  // energies in plasma units: inside units times dx^(D - 2)
    auto verify = [&] { check_status(w.read_status()); };
    int recorded = -1;
    auto record = [&](int n) {  // (each row flushed: a cancelled run keeps its history)
        density();
        verify();  // (after the density deposit: a failure there is reported before its row)
        const double gauss = gauss_error(f, g, background) / (o.dx * o.dx);
        const double electric = field_energy(f, false) * energy, magnetic = field_energy(f, true) * energy;
        std::size_t vertices = 0, elements = 0;
        std::vector<Kinetic> k(meshes.size());  // (in plasma units: energies m_e c^2, momenta m_e c)
        double total = 0;
        for (std::size_t s = 0; s < meshes.size(); ++s) {
            k[s] = kinetic(meshes[s], o.curved);
            k[s].energy *= energy;
            for (double& p : k[s].momentum) p *= energy;
            total += k[s].energy;
            vertices += meshes[s].vertex.size();
            elements += meshes[s].element.size();
        }
        history << n << ',' << n * double(g.dt) * o.dx << ',' << vertices << ',' << elements << ',' << total << ','
                << electric << ',' << magnetic << ',' << gauss;
        for (std::size_t s = 0; s < meshes.size(); ++s) history << ',' << k[s].energy << ',' << meshes[s].vertex.size();
        for (std::size_t s = 0; s < meshes.size(); ++s)
            for (double p : k[s].momentum) history << ',' << p;
        history << std::endl;
        check_diagnostics(gauss, total, electric + magnetic, n);
        for (std::size_t s = 0; s < meshes.size(); ++s)
            for (double p : k[s].momentum)
                if (!std::isfinite(p))
                    throw std::runtime_error("a non-finite momentum of " + o.species[s].name + " at step " + std::to_string(n));
        recorded = n;
    };
    std::unique_ptr<Renderer> renderer;
    auto render = [&](int n) {
        if constexpr (D >= 2) {
            for (std::size_t s = 0; s < meshes.size(); ++s) {
                const std::vector<std::string>& r = o.render_species;
                if (r.empty() || std::find(r.begin(), r.end(), o.species[s].name) != r.end())
                    renderer->write(o, meshes[s], o.species[s], w, n);
            }
            try {  // (a failure in a render's deposit, e.g. an element too long for the bisection at pixel scale)
                verify();
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("render: ") + e.what());
            }
        }
    };
    if (o.render_every > 0) renderer = std::make_unique<Renderer>(o);
    // Energy spectra: spectra_<species>.csv, a row per output step: step, time, then the particles (per unit length
    // along the box's missing axes) in each bin of gamma - 1; the header holds the bins' lower edges.
    std::vector<std::ofstream> spectra;
    if (o.spectra_every > 0)
        for (const Species& s : o.species) {
            spectra.push_back(open<double>(o.output + "/spectra_" + s.name + ".csv"));
            spectra.back() << "step,time";
            for (int k = 0; k < spectrum_bins; ++k) spectra.back() << ',' << std::pow(10.0, spectrum_log0 + double(k) / spectrum_per_decade);
            spectra.back() << '\n';
        }
    Device<double> bins(spectrum_bins);
    auto diagnose = [&](int n) {  // spectra, mesh dumps, phase-space images when due
        auto due = [&](int every) { return every > 0 && n % every == 0; };
        if (due(o.spectra_every))
            for (std::size_t s = 0; s < meshes.size(); ++s) {
                thrust::fill(bins.begin(), bins.end(), 0.0);
                spectrum(meshes[s], o.curved, std::pow(o.dx, D - 2) / o.species[s].charge, ptr(bins));
                const thrust::host_vector<double> h = bins;
                spectra[s] << n << ',' << n * double(g.dt) * o.dx;
                for (double x : h) spectra[s] << ',' << x;
                spectra[s] << std::endl;
            }
        if (due(o.dump_every)) write_dump(o, meshes, n);
        if (due(o.phase_every)) {
            const std::vector<std::string>& r = o.render_species;
            try {  // (a failure in an image: an element too long in pixels, a value beyond float32, a write)
                for (std::size_t s = 0; s < meshes.size(); ++s)
                    if (r.empty() || std::find(r.begin(), r.end(), o.species[s].name) != r.end())
                        write_phase(o, meshes[s], o.species[s], w, n);
                verify();
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("phase: ") + e.what());
            }
        }
    };
    record(0);
    if (o.fields_every > 0) write_fields(o, f, background, 0);
    if (renderer) render(0);
    diagnose(0);
    check(cudaDeviceSynchronize());
    const auto start = std::chrono::steady_clock::now();
    int last = 0;
    std::string stopped = "steps";
    bool capped = false;
    for (int n = 1; n <= o.steps; ++n) {
        try {
            step(n, meshes, o, f, w);
        } catch (const VertexCap& e) {
            std::cout << "step " << n << ": " << e.what() << ": the run ends at step " << last << '\n';
            stopped = "vertex cap at step " + std::to_string(n);
            capped = true;  // (refinement may have split edges before reaching the cap: the last row is written again)
            break;
        }
        last = n;
        if (n == o.steps || (o.diag_every > 0 && n % o.diag_every == 0)) record(n);
        else if (n % 64 == 0) verify();
        if (o.fields_every > 0 && n % o.fields_every == 0) {
            if (recorded != n) density();
            write_fields(o, f, background, n);
        }
        if (renderer && n % o.render_every == 0) render(n);
        diagnose(n);
    }
    if (recorded != last || capped) record(last);
    check(cudaDeviceSynchronize());
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    verify();
    history.close();
    for (std::ofstream& file : spectra) file.close();
    auto timing = open<double>(o.output + "/timing.txt");
    timing << "loop_seconds=" << seconds << "\nsteps=" << last << "\nstopped=" << stopped << '\n';
    timing.close();
    if (o.snapshot) write_snapshot(o, meshes, f, background);
    std::size_t vertices = 0, elements = 0;
    for (const Mesh<D>& m : meshes) vertices += m.vertex.size(), elements += m.element.size();
    std::cout << D << "D " << (o.pic ? "PIC" : "SIC") << " " << o.problem << ", order " << g.order << ", "
              << 8 * sizeof(Real) << "-bit: " << last << " steps (t = " << last * double(g.dt) * o.dx << "), "
              << vertices << " vertices, " << elements << " elements, " << seconds << " s -> " << o.output << '\n';
}
}  // namespace

int main(int argc, char** argv) {
    try {
        const Options o = parse(argc, argv);
        if (o.g.dim == 1) run<1>(o);
        else if (o.g.dim == 2) run<2>(o);
        else run<3>(o);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "sic: " << e.what() << '\n';
        return 1;
    }
}
