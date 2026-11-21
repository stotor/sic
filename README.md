# sic

A compact CUDA code for **SIC (simplex-in-cell) and PIC in 1D, 2D and 3D**: cold species, each with its own charge
and mass, exact charge-conserving element deposits, length-based refinement, curved 1D elements, conducting walls. It
sets up two problems, **Weibel** (filamentation of counter-streaming beams) and **Fadeev** (island coalescence), with
the SIC paper's runs (named below) as defaults. Five source files (about 3,350 lines), no `.inc` fragments, no switches between
alternative implementations. Everything runs on the GPU; the host parses options, launches kernels and writes output.

It is a rewrite of two earlier codes by the same group: a large research code (the "full code" below: ~29,500 lines,
92 files, with tuned SIC and PIC engines) and a compact 1D/2D code. It has the compact code's features plus 3D
(tetrahedra on a Kuhn lattice), species and the two problems, and most of the full code's speed (see Performance).

## Build and run

An NVIDIA GPU (the default target is sm_86, RTX 30 series; `make ARCH=sm_XX` for another), CUDA 12 (tested with nvcc
12.9) and a C++17 host compiler (tested with gcc 14).

```sh
make -j6                                     # build/sic (double precision), build/sic-single (single)
./build/sic --problem weibel --dim 1 --output w1d          # the paper's W1D deck: 128 cells, 500 steps, < 1 s
./build/sic-single --problem weibel --dim 3 --output w3d   # W3D: 32^3, beams +-(0.1, 0.5, 0), 600 steps
./build/sic-single --problem fadeev --dim 2 --output p22          # P22 (refine 1, 20M vertices per species): to the cap
./build/sic-single --problem fadeev --dim 2 --mass positrons=25 --output p22_ions          # positive species m 25
./build/sic --help
```

On a shared machine run it through the scheduler (e.g. `srun --gres=gpu:1 ./build/sic-single ...`).

## Problems

**Units.** Plasma units: c = epsilon_0 = m_e = e = 1, lengths in c/w_p and times in 1/w_p of density n0 = 1 (so a
species of density n, charge Z and mass M has w_p^2 = n Z^2 / M). Momenta are u = gamma v / c. The box is set in
these units and `--nx` sets the resolution; cells are cubes (dx = Lx / nx; the other axes need a whole number of
cells). The default time step is 0.9 dx / sqrt(dim) (`--dt` in 1/w_p, at most 0.95 dx / sqrt(dim)).

**Weibel** (`--problem weibel`, 1D, 2D, 3D; the paper's W1D and W3D): two cold electron beams `beam1`, `beam2` of
density 1 each with momenta +-u (`--u`, default (0.1, 0.5, 0)) through immobile ions (`ions`, charge 1, density 2;
`--mass ions=M` makes them a mobile species at rest). The instability grows from a random magnetic seed, the same in
every run with the same `--seed-rng` (the full code's draws): in 1D 64 B_z modes of amplitude `--seed` (default 1e-6)
with random phases; in 2D and 3D 32 divergence-free modes with integer wave numbers up to 4 per axis. Defaults: 1D
128 cells over 12.8 (dx 0.1, dt 0.09, 500 steps); 2D and 3D 32 cells per axis over 6.4 (600 steps); 2 vertices per
cell and axis.

**Fadeev** (`--problem fadeev`, 2D, 3D; the paper's P22 and R3D): a chain of magnetic islands, A_z = -B0 lam ln[cosh(y'
/ lam) + eps cos(x / lam)] (Fadeev et al. 1965), carried by cold pair sheets of density n F, F = (1 - eps^2) / (cosh +
eps cos)^2, each sheet species with a uniform drift along z that carries half of the current J_z = (B0 / lam) F (v =
B0 / (2 lam |Z| n), below c). A field-only kick dA_z = kick B0 lam cos(x / 2 lam) (...) pushes the islands together.
- 2D (P22, the geometry of Stanier et al. 2017): one chain at Ly/2 between conducting walls in y, periodic in x;
  species `sheet_electrons`, `sheet_positrons` and a uniform background at rest `bg_electrons`, `bg_positrons`
  (density `--bg`, default 0.2). Box: `--islands` (default 2) chain periods 2 pi lam along x, `--ly` (default 4 pi lam;
  Stanier's box: 2 pi lam). Defaults 96 cells per island (dx 0.065), 346 steps, kick -0.5 (merger across the central
  X-point), refinement 1 cell with 20M vertices per species: the full code's 2026-11-06 P22 movie run (`final_4k`;
  its fields from `make_fadeev.py 192 192 4pi 4pi 1.6 1 0.4 --pert -0.5`).
- 3D (R3D): two mirror chains at Ly/4 (+J_z) and 3 Ly/4 (-J_z) in a periodic box 4 pi lam x 4 pi lam x pi lam, species
  `sheet1_*` and `sheet2_*` (no background unless `--bg`), a 1% random divergence-free seed (`--seed`), 32 cells per
  island along x, refinement 2 cells with 5.5M vertices per species (the R3D deck's length limit).
- Parameters: `--b0` (1.6), `--lambda` (1), `--eps` (0.4), `--kick` (-0.5), `--guide` (a uniform B_z, in B0; 0),
  `--bg`, `--islands`.
- The sheets' density enters as each element's charge: the profile's mean over the element (a degree-2 rule); PIC
  particles carry the profile at their position.

**Species: charge and mass.** `--charge NAME=Z` and `--mass NAME=M` set the charge (in e) and mass (in m_e) of every
species whose name contains NAME, after the problem has made them: `--mass positrons=25` turns both Fadeev positron
species into ions of mass 25, `--mass ions=1836` makes Weibel's ions mobile. A mass `inf` is immobile: such a species
is not simulated, the fixed neutralizing background stands for it. The background is the initial charge of all
species together, so every run starts neutral. A Fadeev sheet's drift follows its charge and density.

## Options

| Option | Meaning |
|---|---|
| `--problem weibel\|fadeev`, `--dim 1\|2\|3` | the problem (Fadeev: 2D or 3D) |
| `--nx N`, `--lx --ly --lz` | cells along x (4..65536 per axis); box lengths in c/w_p (Fadeev: x from `--islands`) |
| `--steps`, `--dt` | steps; time step in 1/w_p (default 0.9 dx / sqrt(dim); curved elements: dt / dx <= 0.76) |
| `--engine sic\|pic`, `--order N` | SIC orders 0..4 in 1D, 0..2 in 2D and 3D; PIC orders 1..2 |
| `--gather sameshape\|galerkin` | field interpolation: order N along every axis at the field's own position (default; the full code's SAMESHAPE, every run of the SIC paper), or order N - 1 along a field's half-cell axes (Galerkin) |
| `--vpc N`, `--ppc N`, `--curved` | SIC vertices per cell and axis; PIC particles per cell (a square in 2D, a cube in 3D; default 4, 3D 8); 1D SIC quadratic elements (an extra midpoint vertex per element), orders 1..4 |
| `--refine L`, `--max-vertices N` | split SIC edges longer than L cells (0: off); the run ends normally when a species reaches N vertices (default 2M, Fadeev: the decks', see above; the initial mesh above it is an error) |
| `--charge NAME=Z`, `--mass NAME=M` | species charge and mass (above) |
| `--u UX,UY,UZ` | Weibel beam momenta |
| `--b0 --lambda --eps --kick --guide --bg --islands` | Fadeev parameters |
| `--seed A`, `--seed-rng R` | random B seed amplitude and generator |
| `--diag-every N`, `--fields-every N`, `--sort-every N`, `--no-snapshot` | history sampling (0: final only; a sample costs a density deposit and a few reductions); field dumps; Morton re-sort interval (0: at the start only); skip the final CSVs |
| `--render-every N`, `--render-substeps K`, `--render-scale S`, `--render-window U0,U1,V0,V1`, `--render NAME,...`, `--render-smooth` | density images: 2D exact, 3D columns (below) |
| `--render-axis x\|y\|z`, `--render-slab W0,W1` | 3D: the view axis (default z) and the slab along it (default the box) |
| `--phase-every N`, `--phase A,B[,LO,HI...]`, `--phase-size W,H` | phase-space images (below); `--phase` repeats, one image per plane |
| `--spectra-every N`, `--dump-every N` | energy spectra; vertex and element dumps (below) |
| `--output DIR` | a new directory |

## Output

In `--output`: `run.txt` (parameters and species), `history.csv`, `timing.txt` (loop seconds, renders and dumps
included; why the run stopped), and, unless `--no-snapshot`, the final `fields.csv` (rho, J, E, B per node),
`vertices_<species>.csv` and `elements_<species>.csv` (SIC: charge, corners with their images; curved: the midpoint
last). Everything in plasma units. Every file is closed explicitly: a failed write stops the run with an error.
- `history.csv`: step, time, vertices, elements, kinetic, electric and magnetic energy, the largest Gauss residual,
  each species' kinetic energy and vertex count, then each species' momentum px, py, pz (sum of m u, m_e c: exact,
  an element's mean u is its corners' mean; one pass with the energy); flushed every row (a cancelled run keeps it). A
  run ended by the vertex cap writes its last step again with the mesh refined so far (the spectra, images and dumps of
  that step are of the mesh before: the same state). The kinetic energy (and the momentum)
  is from the momenta at t - dt/2: averaged with the next row's it is at t (the total then conserves to about a
  percent where the raw sum swings by +-dt/2 of the energy exchanged). A non-finite energy or momentum stops the run.
- `spectra_<species>.csv` (`--spectra-every`): a row per output step: step, time, then the particles (per unit length
  along the box's missing axes) in 300 log bins of gamma - 1, 30 a decade from 1e-6 (below it: the first bin; above
  1e4: the last); the header holds the bins' lower edges. An element is integrated with the kinetic energy's rule
  (D + 1 points, exact for a quadratic), the full code's 2D spectra in every dimension; PIC: each particle. (An energy
  within round-off of a bin's edge may land in the bin beside it.)
- `phase<k>/<step>/<species>_f.npy` (`--phase-every`, the species of `--render`): float32 [H][W], row 0 at the lowest B:
  f(A, B), particles per unit A per unit B (per unit length along the missing axes), A and B from x, y, z (the whole
  box) and ux, uy, uz (over LO..HI, one range per momentum axis in order: `--phase x,ux,-1,1`, `--phase
  ux,uy,-1,1,-2,2`; pixel k of a momentum axis: [LO + k w, LO + (k + 1) w)); default 1600 x 1000 pixels. Along a
  momentum axis an element is cut exactly at the image's edges in its own material coordinates, the part inside
  keeping its share of the charge: nothing outside the range reaches the image, however far an element extends. Each
  element maps affinely onto a triangle of the plane and is integrated over each pixel exactly (the 2D images'
  deposit, in batches of 2^18 elements; collapsed ones keep their charge), so a cold beam is a sharp line or point; a
  1D segment's charge, uniform along its line, goes to the pixels it crosses by the fraction of its length in each (a
  walk through the pixels: a fixed mesh's segments can be thousands of pixels long); a PIC particle to its pixel. The
  geometry is in double, an element's corners moved by whole periods so that corner 0 lies in the box. Momenta are the
  leapfrog's (t - dt/2), one frame per step; SIC in 1D and 2D (a tetrahedron's shadow is not uniform), PIC in any
  dimension. A point on an upper wall is in the last pixel (moved 1e-6 pixels down, as in the density images). A
  non-finite value stops the run ("phase: ...").
- `mesh/<step>/<species>_vertices.npy` (`--dump-every`): float64 rows x, y, z, ux, uy, uz (PIC: and the charge), positions
  unwrapped; SIC also `<species>_elements.npy` (int32 rows: vertex, image x, y, z for each corner; curved: then the
  midpoint) and `<species>_charges.npy` (float64): the final snapshot's columns, in plasma units (written in blocks
  of 2^20 rows: the host holds one block at a time).
- `fields/<step>.npy` (`--fields-every`): float32 [10][nz][ny][nx]: rho (net), Jx, Jy, Jz (the last step's), Ex, Ey,
  Ez, Bx, By, Bz on their Yee positions (`numpy.load`). With walls, node row 0 is both walls (the periodic seam): its
  rho and tangential J hold both sides.
- `render/<step>/<species>_n.npy` (`--render-every`): the species' number density on a grid of pixels, float16,
  pixel (i, j) the square [i, i + 1) x [j, j + 1) pixels from the box's corner, row 0 at the lowest v, cropped to
  `--render-window` (u0, u1, v0, v1; in 2D u, v = x, y; in 3D the two axes besides the view axis, in order); with `--render-substeps K` (SIC), frames at step - 1 + j/K as `<species>_n_s<j>of<K>.npy`
  (vertices on straight paths between the time levels). Densities above 65504 (float16) are written as 65504. PIC
  images are the particles' charge per pixel (NGP), one frame per step, never smooth. Exact by default: each element's
  charge integrated over each pixel (the order-0 deposit at pixel scale, so anti-aliased). `--render-smooth` draws each
  element linear between its corners' lumped densities (the full code's smooth shading; the total is the same): it
  hides the element mosaic that stretched elements show in the exact density, at 2-3 times the cost per frame (it
  leaves out the charge of a vertex whose adjacent elements all have zero area, a collapsed point). In 3D a frame is
  the mean density along the view axis through the slab at one sample point per pixel (the full code's 3D render), at
  (i + 0.6408, j + 0.5752) pixels: off the centre, so that the initial lattices' faces miss them. Each tetrahedron adds
  its density times the length of the sample's ray inside it and the slab, so folds and caustic surfaces show as bright
  sheets; a tetrahedron of zero volume adds nothing (its charge lies in a plane). PIC: the particles' charge per pixel
  column in the slab. This is the full code's
  frame layout: its `studies/render/make_composite.py RUN OUT.mp4 SPECIES_A SPECIES_B` makes a two-colour movie of
  them.

## Code map and where to start reading

| File | Lines | Contents |
|---|---|---|
| `src/sic.hpp` | 433 | units, types (positions as integer cell + offset, vertices, elements, a species' mesh, fields, options and species), B-splines and the window basis (compile-time tables in constant memory), gather and Boris push, the refinement vote, the functions of each file, the time step |
| `src/deposit.cu` | 1182 | the exact element deposit (signed corners, simplex moments, bisection into pieces; consecutive blocks spread over the element list), curved 1D elements, PIC particles (fused push and deposit), the fold at walls, 3D column images, phase-space images |
| `src/mesh.cu` | 498 | a species' initial lattice or particles (with the Fadeev profile), the Morton sort (by path class within blocks), refinement |
| `src/fields.cu` | 313 | the vertex push, the Yee update, the diagnostics (energies, momentum, spectra), the problems' initial magnetic fields |
| `src/main.cu` | 916 | options and the two problems, the loop, output (history, field and mesh dumps, density images, phase-space frames, spectra, snapshot) |
| `tests/test.cu` | 1766 | tests against independent references (below) |

A reading order:
1. **`sic.hpp`, top to bottom.** The data first (`Pos`, `Vertex`, `Element`, `Mesh`, `Fields`, `Species`,
   `Options`, `Work`), then the shape functions and their constant-memory tables, the gather and push, and at the end
   `step()`, the time step: refine every species, push and deposit each, advance the fields, re-sort if due.
2. **`main.cu`**: the problems (`set_problem`), the options, the loop around `step()`, the outputs.
3. **`fields.cu`**: short kernels (push, Yee update, reductions), then the initial fields.
4. **`deposit.cu`, the heart.** Read it in file order:
   - `Local`, `load` and `frame`: an element in local coordinates, and its window of nodes.
   - `cut`, `corner` and `part`: the signed sub-simplices on one side of a plane.
   - `mean` (orders 0 and 1) and `moments` (orders >= 2): integrals of polynomials over a simplex.
   - `Moments`, `visit` and `branch`: the recursion over the planes.
   - `to_nodes` and `scatter`: node weights onto the grid.
   - `deposit_simplex`: the 2^D mixed simplices of Esirkepov's current.
   - `deposit_curved`, then the kernels `k_elements`, `k_pieces` and the PIC kernels, then the walls' fold.
5. **`mesh.cu`**: the initial lattice, the Morton sort and its path classes, then refinement (votes, keys, splits).
6. **`tests/test.cu`**: each test says in its comment what it checks against what.

C++/CUDA idioms used throughout:
- **Templates.** Kernels and device functions take `D` (dimension), `N` (shape order), `K` (corners) and `C`
  (output components; 1 density, 3 current) as template parameters.
- **Unrolling.** The loops over these are `#pragma unroll`, so the compiler turns the small arrays they index into
  registers; `#pragma unroll 1` keeps a loop rolled.
- **`with_order` (sic.hpp)** turns the run-time order into a compile-time constant through a generic lambda.
- **`if constexpr`** drops the code a dimension or order does not need.
- **Coefficient tables** are built by `constexpr` functions at compile time and stored in `__constant__` memory: device
  code does not fold such functions itself.
- **Sorting, scans and reductions** use thrust (and cub's caching allocator for its temporaries).
- **Errors.** A status word on the GPU collects them (`Status`), and the host reads it every 64 steps, at the
  diagnostics and at the end, so a step needs no host round trip.

## How it works

**Units inside.** Lengths are in cells and times in dx/c: a species of density n has n dx^2 charge per cell volume
per unit charge, E and B are dx times their plasma-unit values, J and rho dx^2 times; energies are dx^(D-2) times.
`main.cu` converts at the input and the output; kernels only see cells.

**Species.** Each species is its own mesh (`Mesh`: vertices, elements, its charge-to-mass ratio), on the same
initial lattice: vertices at the centres of a lattice of spacing 1/vpc cells along a periodic axis (the full code's
lattice) and on the walls along a walled one. Elements carry their charge (PIC particles: each particle). A step
refines every species first (the vertex cap then stops a run before anything has moved), then pushes and deposits
each into the same current, then advances the fields once.

**Positions** are an integer cell plus an offset in [0, 1) per axis, unwrapped (vertices are never wrapped into the
box); an element's corners carry their periodic image shifts. Every geometric computation works in a local frame
around the element, so single precision keeps sub-cell resolution anywhere in the box.

**The element deposit** (`deposit.cu`). An element carries its charge uniformly over its material simplex, so its
density is the integral of the B-spline over the simplex (degenerate and folded elements are fine) and its current is
Esirkepov's, integrated the same way over the 2^D mixed simplices (axis d at the old or the new time level, weighted
k!(m-k)!/(m+1)!). Gauss's law holds to round-off. When an element's old and new corners lie within two branch cells
per axis, every node weight of its window is `sum_q A[i][q] phi_q(u)` per axis, with `phi_q = u^q` (q <= N) and
`phi_(N+1) = max(u, 0)^N`, u measured from the plane between the two branches. The integral of a product of these over
the simplex needs the simplex's part on the positive side of the planes, which is a signed sum of simplices: cutting by
one plane leaves nothing, the corner at the one vertex inside, the whole simplex, the whole minus the corner at the one
vertex outside, or (3D, two and two) the wedge as three tetrahedra. Cutting by the next axis splits each part again.
On every part the integrand is a polynomial, integrated exactly. Orders 0 and 1: closed-form means of products of up
to three affine fields. Higher orders: all moments at once from the simplex's generating function (the mean of y^a
over a simplex with corners y_k is D! a! / (D + |a|)! times the coefficient of s^a in prod_k 1 / (1 - y_k . s); one
factor multiplies a truncated series by the recurrence G[a] += sum_v y_v G[a - e_v]). All of this is straight-line
code in registers, one element per thread. Elements that do not fit the window go to a list; a second kernel bisects
them in material space, depth first and along the axes they do not fit, into pieces that fit (its threads take the
listed elements from a work queue). The same window and current assembly serve PIC particles.

**Walls.** The gather mirrors the stencil across a wall (even along a half-cell axis, odd along a node axis). The
deposit's fold is the transpose of that mirror: with walls the deposit writes into a grid padded by N/2 + 2 nodes past
each wall and adds the padding back through the same mirror, so charge is conserved off the wall nodes at every order.

**Accumulation.** Every element adds its window to the grid with double-precision atomics (native on sm_86, unlike
shared-memory float atomics). A sum of single-precision contributions is then exact, so independent of the order of
the additions, when the sum of their absolute values, in units of the smallest one's last bit, stays below 2^53 (a
sufficient rough bound: (largest / smallest) x their number below 2^29); otherwise its last bit can depend on the
order. A cold drifting lattice stays uniform to the last bit. Float sums were tried and are not faster.

**Refinement** (`mesh.cu`): every element votes for its longest edge if that is longer than the limit, and an edge is
bisected when every element around it voted for it, so the mesh stays conforming in 1D, 2D and 3D alike (the globally
longest edge always qualifies; passes repeat until no edge is too long, a later pass voting only the elements that
can have changed). Refinement runs at the start of a step. New vertices interpolate the positions at both time levels (the density
at either level is unchanged) and the momentum. Curved 1D elements are split at their parameter midpoint.

**Order in memory**: vertices and elements are sorted by the Morton order of their cells at the start, every
`--sort-every` steps, and whenever refinement has grown the mesh by an eighth (`--sort-every 0`: neither). Within
blocks of 8^D cells, elements are grouped by their path class (per axis: whether they meet the window's plane, and
how many corners lie above it), so that a warp's elements take one path through the deposit.

**Initial fields** (`fields.cu`): B from the discrete curl of a potential on the Yee grid, so div B = 0 to round-off (3D Fadeev: up to the ~1e-6 B0 of its five-image sum at the periodic y seam):
the Fadeev A_z on the nodes, the seed's vector potential on the cell edges. E starts at zero, consistent with the
neutral start.

**Images** (`main.cu`): in 2D the vertices are mapped to pixel coordinates (between the two time levels for sub-step
frames) and the order-0 deposit integrates each element over each pixel of a grid covering the box; smooth shading
lumps densities on the vertices and deposits the out-of-plane current of a unit-density copy at rest whose vertices
move along z in proportion to them, which is the integral of the linear interpolant. In 3D one thread per tetrahedron
(a block for one whose shadow holds over 4096 samples) inverts its corner matrix once in double precision
(barycentric coordinates linear in the sample and in w) and, for each sample in its shadow, finds the interval of the
ray where all four are >= 0, clipped to the slab's periodic copies. The samples run in single precision relative to
corner 0 unless the tetrahedron is thin enough for that to cost more than 6e-4 of its column (then double: a thin
tetrahedron's column stays finite as its volume goes to zero; all-single frames were off by up to 30% in a few pixels
of a folded sheet). A ray on a face parallel to it belongs to one of the two tetrahedra sharing the face: the one on
the side of a symbolic displacement of the sample along +u, then +v. The corners move by whole periods so that
corner 0 lies in the box (the image is periodic in u and v). A 4K frame of a 9M-tetrahedron species takes about 1.2 s.

## Conventions

All three momentum and field components exist in every dimension. E and B at integer times, momentum at half steps:
`B -= dt/2 curl E; E += dt (curl B - J); B -= dt/2 curl E` on the Yee grid (E_d and J_d half a cell along d, B_d half
a cell along the other axes). The gather interpolates each component at its own position with the order-N B-spline
along every axis (`--gather galerkin`: order N-1 along half-cell axes, NGP at order 0 and 1 there). Walls sit at node 0 of an axis
(the periodic seam: the walls at 0 and L share it): the gather mirrors across them, tangential E vanishes there and
normal B keeps its initial value, particles and vertices reflect, and the Gauss diagnostic skips the wall nodes (the
conductor's surface charge is not evolved).

## Known behaviour and limits

- **Energy at low orders.** SIC gains energy in a cold Langmuir oscillation driven to wave breaking (1D, 16 cells, dt
  .15, momentum amplitude .8): after 1000 steps x12 at CIC and x16 at TSC with the default gather, x65 at NGP; the
  Galerkin gather makes CIC much worse (x56: its NGP interpolation along a field's own axis makes the force jump at cell
  boundaries).
- **Lattice placement against the shape's knots.** The B-spline knots sit at integers for odd orders (CIC) and at
  half-integers for even orders (NGP, TSC). The lattice sits at cell centres (the full code's): between the knots for
  CIC, on them for TSC, where any expansion makes an element span three cells and the bisection kernel works harder.
- **Bisection depth.** At most 24 levels (segments up to 2^24 cells, triangles up to about 4000 cells across) and 40
  in 3D (tetrahedra up to about 700 cells across); a deeper piece stops the run.
- **Caps.** `--max-vertices` per species (default 2M) and 2e8 elements or particles.
- **PIC** is a plain global-atomics PIC, kept for comparisons and tests; the full code's tiled PIC is 2-2.5x faster.
- **Diagnostics** cost a density deposit and a few reductions per sample.

## Tests

`make build/test build/test-single`, then run both on a GPU (seconds each; they also run the program built next to
them). A status bit set on the GPU or a NaN error fails a test. 149 checks in each precision:
- **density moments** of random elements of every size against an exact oracle (B-splines reproduce polynomials up
  to their order; the simplex means come from barycentric expansion, independent of the deposit code);
- **additivity**: an element deposits exactly what its 2^(3D) material sub-elements deposit, cut configurations
  against pieces with fewer cuts;
- **charge continuity** for straight and curved elements, and a PIC particle against a collapsed SIC element;
- **refinement**: edge lengths, conformity (every face shared by two elements), density and momentum unchanged;
- **the field solver**: the vacuum Yee dispersion, and Langmuir oscillations with Gauss's law, for SIC and PIC,
  refined, curved and walled runs, stepped by `step()` with re-sorts; the plasma frequency w^2 = n Z^2 / M for other
  charges and masses and with cells of dx = 0.25 c/w_p;
- **Weibel growth** of cold beams (1D, k = 4, u = +-0.2, a velocity seed): 0.999 of the cold-fluid rate for SIC CIC
  and TSC and PIC;
- **the problems' initial fields**: div B = 0 for the random seeds and the Fadeev field; the initial magnetic energies
  of W1D, P22 and R3D equal the full code's; the Fadeev sheets' current = curl B to the profile's dx^2 (0.25% at dx
  0.13);
- **walls**: an element or particle near a wall deposits what it and its mirror image deposit in a periodic box twice
  as tall;
- **images**: a uniform species draws its density in every pixel (electrons and positive species, flat and smooth,
  sub-step frames, between walls; 3D columns along z and x, through a slab, PIC); the Fadeev sheet's image is
  symmetric under a half turn of the box (the pixel convention), and a window is the matching block of the whole
  image; in 3D: the Fadeev sheet's mean column is its total charge over the box, two slabs that split the box add up
  to the whole column, its chain and islands lie in the right rows and columns, a window is the matching block, the
  runs record their view; on built meshes, a thin tetrahedron's column stays exact as its volume goes to zero (to
  2^-40), two cubes' shared face (holding samples) counts once, and moving a mesh by whole periods (beyond 2^31
  pixels) or PIC particles changes nothing;
- **phase images on built meshes** against references that cut each element with each pixel (Liang-Barsky for
  segments, a polygon clip for triangles): tilted segments across the seam with pixels smaller than them, segments
  reaching far past the momentum range on either side (1/4 and 1/7 inside, no wrapped tail), an extreme zoom (10^12
  pixels of extent), whole-period moves beyond 2^31 pixels, a wall; triangles across one and both momentum axes'
  edges, collapsed onto a segment of the plane (its material density 2 (1 - s)), 30000 pixels long, at a wall; PIC
  particles at and beside the range's edges, far outside, moved, on a wall;
- **momentum and spectra** of built meshes whose corners move differently (a triangle, a tetrahedron, a curved
  segment with its midpoint off the chord, PIC particles of unequal charge, the first and the last bin);
- **the deposit's shared bisection**: M copies of a long element deposit M times what one does, for every fan
  (M = threads / 2^f, f = 0..10, and the queue past them), density and current, straight, curved, 2D; the depth
  limit (2^23 cells deposit, 2^25 set BAD_PIECES);
- **diagnostics** (the program): W1D's beams: momentum 12.8 u each, their 12.8 particles in the spectrum's bin of
  gamma - 1, their (x, ux) images one row of f = n / du in every pixel; Weibel 2D's (ux, uy) images at the beams'
  momenta with all their particles; P22's (y, uz) image, across the walls, holds the dumped charge, the uniform
  background's every y pixel alike; 3D PIC's phase image and spectrum hold the particles; mesh dumps equal to the
  final snapshot row for row (1D, 2D, curved, PIC); phase images and spectra at their cadence, two planes;
- **1D elements against an independent B-spline oracle**, straight and curved, node by node (density and all three
  currents);
- **the harness itself**, and **the program's command line** (bounds, overflowing inputs, species names, the
  problems' ranges, caps, a failed write).

Against the full code on the paper's decks (its SAMESHAPE binaries; the default gather here):
- W1D (double, fixed meshes): the same magnetic energy history to round-off, 1e-11 through t = 27 and 1e-7 at t = 72
  (orders 1, 2 and 3), and within 2e-4 at t = 90, deep in the chaotic phase.
- P22 (192^2, CIC, refine 1, single): the same initial magnetic energy (179.50638); through t = 8 the field energies
  agree to 1e-4 and the kinetic energy to 7e-4 (the sheets' profile enters as element charges here, as a density map
  there), and the meshes grow alike (0.98-1.00 of the full code's vertices; at 96^2, t = 10, 0.96).
- R3D (3D, single): the same initial magnetic energy as the full code's quiet PIC with the same seed (488.16899), the
  magnetic energy within 0.4% of it through t = 7.
The earlier test history: compute-sanitizer memcheck found no errors; unrefined 1D/2D runs matched the earlier compact
code's histories to 1e-13.

## Performance

Single precision, RTX 3090 Ti at 450 W, ms per step, measured at commit c492fbd with the command line of then (one
electron species: `--drift .1 --perturb 0 --dt .1 --offset .5/vpc`, vertices at cell centres as now; 4096 cells (1D),
256^2 (2D), 64^3 (3D); refinement limits .75 (1D), 1.25 (2D, 3D); perturbed cases with `--perturb .1`; PIC 1024^2 x 16
ppc, 128^3 x 8), alternated with the full code's matching runs. The kernels are the same now; a second species costs
its own deposit.

| case | this code | full code | full / new |
|---|---:|---:|---:|
| 1d-cic | 0.0171 | 0.020-0.030 | 1.2-1.7 |
| 1d-cic-ref | 0.0292 | 0.040 | 1.37 |
| 2d-cic | 0.055-0.081 | 0.080 | 1.0-1.45 |
| 2d-cic-ref | 0.120 | 0.140 | 1.16 |
| 2d-cic-1024 | 0.871 | 0.750 | 0.86 |
| 2d-tsc | 0.100 | 0.130 | 1.30 |
| 2d-cic-pert | 0.0674 | 0.090 | 1.34 |
| 2d-tsc-pert | 0.164 | 0.170 | 1.04 |
| 3d-cic | 3.30 | 1.63 | 0.49 |
| 3d-cic-vpc2 | 10.9 | 7.44 | 0.68 |
| 3d-cic-ref | 7.12 | 12.2 | 1.72 |
| 3d-tsc | 8.55 | 14.5 | 1.70 |
| 3d-cic-pert | 3.44 | 2.78 | 0.81 |
| 3d-tsc-pert | 19.7 | 11.6 | 0.59 |
| 3d-cic-ref-pert | 9.14 | 8.06 | 0.88 |
| 2d-pic | 2.41 | 0.99 | 0.41 |
| 3d-pic | 3.37 | 1.70 | 0.51 |
| 3d-pic-tsc | 7.29 | 2.81 | 0.39 |

P22 at 192^2, CIC, refine 1, to t = 12 (288 steps, 20.5M vertices at the end): 2.36 s, the full code 2.49 s on the
same deck (its meshes within 1.5%; one job, alternated); Weibel 2D (256^2, u_z +-0.5, refine 1, 440 steps) 3.3 s;
W1D (order 1, 8 vertices per cell, 1000 steps, double) 0.28 s, the full code 2.8 s.
A movie run to the 20M-vertex cap per species takes minutes, almost all of it in its 4K images.

Where the time goes:
- **2D folded sheets** (P22 to t = 12): deposit 1.05 s, refinement 0.80 s (votes, splits and a host read per pass),
  Morton sorts 0.28 s, push 0.12 s. In Morton order the blocks running at once deposited into the same few cells and
  their atomic adds queued (the deposit took 1.83 s): consecutive blocks now take stretches 1/64 of the element
  list apart.
- **Refined 3D** (R3D, refine 2, 70 steps: 20 s): 75% in the bisection kernel for elements longer than a cell: half
  of the elements, 3.6 pieces each, each piece as costly as a whole element in the element kernel (7.3 vs 6.5 ns): the
  integrations, not the bisection. (A clipping integrator, the full code's scheme, in one thread per element spilled:
  13x slower.)
- **A few long elements** (W1D's fixed mesh at late times: some 450 per deposit, 20 pieces each): one thread per element
  left the GPU idle, 307 ns a piece. Each element's first levels are now shared out among up to 2^10 threads (when fewer
  elements are listed than threads): 5.9 -> 0.28 s.
- **3D CIC**: the element kernel's integrations, eight mixed simplices of about 24 signed parts each, 77% of it. It
  is not the atomics: without them the kernel is 2-3% faster. Nor is it occupancy: halving the resident warps costs 16%.
- **Perturbed 3D TSC**: the bisection kernel, with vertices on the knots (see above).
- **PIC**: the full code's tiles, as above.

## Where to go next

1. **A cheaper CIC integration in 3D**: the full code integrates cell-clipped pieces of each mixed simplex on separate
   lanes; a factor ~2 on 3D CIC is at stake.
2. **Wider windows (three cells per axis) for bisected pieces**: TSC on the cell-centred lattice.
3. **A tiled PIC** (the full code's design) if PIC speed matters.
