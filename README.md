# PANTHER

**P**hase-separation **AN**alysis **T**ool for **H**eterogeneous polym**ER**s (PANTHER)
is a [MOOSE](https://mooseframework.inl.gov/)-based finite-element application for
studying phase separation in block copolymers and crosslinked polymer networks.

## Overview

PANTHER solves the Cahn–Hilliard equation (and, for void-bearing cases, a coupled
Allen–Cahn equation) to model spinodal decomposition, dissolution, and void
evolution in two- and three-component polymer systems. Simulations run in 2D and
3D and support coupled chemo-mechanical pull tests on the resulting morphologies.

### About this branch (`polyDeg-approx-local`)

This branch concerns **polymer degradation** (`polyDeg`), modeled with a **local
approximation** of the mixing free energy. Instead of the logarithmic
Flory–Huggins form, the local mixing term is represented by a polynomial
(Taylor-expanded) approximation. The classical Flory–Huggins free energy contains
`ln(c)` terms that are expensive to differentiate and ill-behaved near `c = 0` and
`c = 1`; degradation pushes compositions toward those limits, so a smooth local
approximation keeps the solve well-conditioned. The local term is replaced by a
Taylor expansion about a reference composition `(c1_0, c2_0)`:

- **Two-phase cases** use a single-variable polynomial in `c`.
- **Three-phase cases** use a two-variable polynomial in `c1` and `c2`, with the
  third component recovered as `c3 = 1 - c1 - c2`.

The expansion coefficients (`A00`, `A10`, … `A66`) are precomputed and supplied as
constants in each input file. A small stability term `beta*(1/c1 + 1/c2 + …)`
is retained to keep compositions inside the physical range. The original
Flory–Huggins parameters (`chi`, degree of polymerization `N`, `R`, `T`) are kept
as commented-out references next to the polynomial coefficients in each input file.

The most recent commit on this branch also adds explanatory header comments to
every input file and post-processing script so each case is self-documenting.

## Physics and numerics

- **Cahn–Hilliard transport** is solved in the *split* form, introducing the
  chemical potential `w` as an auxiliary variable to keep the discretization
  second-order. MOOSE's `SplitCHWRes`, `SplitCHParsed`, and `CoupledTimeDerivative`
  kernels are used together with PANTHER's custom kernels.
- **Free energy** is assembled from a local mixing term (Taylor polynomial or
  Flory–Huggins) plus a gradient (interfacial) energy `kappa*|grad c|^2`.
  The Cahn number `Cn` sets the interface width via `kappa = Cn^2`.
- **Void evolution** (`*_void_*` cases) adds an Allen–Cahn order parameter `eta`
  that distinguishes polymer from void.
- **Time integration** uses BDF2 with `IterationAdaptiveDT` adaptive stepping;
  the nonlinear system is solved with Newton and an ASM/ILU or LU preconditioner.

## Custom MOOSE objects

Source for the application-specific objects lives under [src/](src/) with headers
in [include/](include/):

| Object | Type | Purpose |
| --- | --- | --- |
| [SplitCHPhaseSep](src/auxkernels/SplitCHPhaseSep.C) | Kernel | Split Cahn–Hilliard kernel for the chemical potential variable with an added nonlocal term. |
| [CHEAux](src/auxkernels/CHEAux.C) | AuxKernel | Rescales a coupled variable to `(value + 1)/2`. |
| [RandomConstraintIC](src/ics/RandomConstraint.C) | InitialCondition | Random IC drawn from a uniform or user-defined distribution, constrained against a coupled variable. |
| [ScaledSumIC](src/ics/ScaledSumIC.C) | InitialCondition | Initializes a variable as a (optionally prefactored) sum of other variables. |
| [DerivativeSpline2Material](src/materials/DerivativeSpline2Material.C) | Material | Evaluates a piecewise-linear spline over triangular regions with automatic first/second derivatives. |

## Repository layout

```
src/, include/      Custom MOOSE kernels, aux kernels, ICs, and materials
problems/           Input files (.i), parameter sweeps, and post-processing
  template/         Templated input for DBC.sh sweeps
  ic_2p/, ic_2pv/,
  ic_3pv/           Saved Exodus initial-condition meshes
  Post_scripts/     Python post-processing and plotting scripts
test/, unit/        Regression and unit tests
doc/                MooseDocs configuration
Makefile, run_tests Build and test entry points
```

## Compilation

1. Install MOOSE following the official
   [installation instructions](https://mooseframework.inl.gov/getting_started/installation/index.html).
2. Clone this repository:
   ```
   git clone https://github.com/baskargroup/panther
   ```
3. Build the application (the Makefile finds MOOSE via the `moose` submodule or
   the `MOOSE_DIR` environment variable):
   ```
   cd panther
   make -j 6
   ```
4. Verify the build:
   ```
   ./run_tests
   ```

This produces the `panther-opt` executable.

## Running simulations

Run any input file with the optimized executable:

```
./panther-opt -i problems/3phase.i
```

or in parallel:

```
mpiexec -n 4 ./panther-opt -i problems/3p_dis_void_ch.i
```

Results are written to `output/` as Exodus (`.e`) and CSV files.

### Problem catalog

Each input file in [problems/](problems/) carries a header comment describing its
case. By family:

**Two-phase**
- [2phase.i](problems/2phase.i) — spinodal decomposition benchmark (Taylor energy).
- [2phase_3d.i](problems/2phase_3d.i) — 3D spinodal decomposition (Taylor energy).
- [2phase_fh.i](problems/2phase_fh.i) — spinodal decomposition with Flory–Huggins energy.
- [2p_void.i](problems/2p_void.i) / [2p_void_3d.i](problems/2p_void_3d.i) — dissolution with fixed circular/spherical voids.
- [2p_void_fh.i](problems/2p_void_fh.i) — void dissolution with Flory–Huggins energy.
- [2p_void_ic.i](problems/2p_void_ic.i) — builds a void initial condition from saved solution data.
- [2pv_3d.i](problems/2pv_3d.i) — 3D polymer/void case with `eta`-defined voids.

**Three-phase**
- [3phase.i](problems/3phase.i) / [3phase_3d.i](problems/3phase_3d.i) — ternary spinodal decomposition in 2D/3D.
- [3p_dis.i](problems/3p_dis.i) — ternary dissolution without explicit void mechanics.
- [3p_dis_void_ch.i](problems/3p_dis_void_ch.i) — coupled dissolution + Allen–Cahn void evolution.
- [3p_dis_void_ch_am.i](problems/3p_dis_void_ch_am.i) — adaptive-mesh variant; `icfile`/`outfile` select input mesh and output prefix.
- [3p_dis_void_ch_valid.i](problems/3p_dis_void_ch_valid.i) — validation setup with diagnostic aux variables.
- [3p_dis_void_ic.i](problems/3p_dis_void_ic.i) / [3p_dis_void_ic_new.i](problems/3p_dis_void_ic_new.i) — three-phase void initial-condition generators.
- [3pvt_dis.i](problems/3pvt_dis.i) — three-phase void transport/dissolution with adaptive refinement.
- [3p_dis_mech.i](problems/3p_dis_mech.i) / [3p_dis_void_mech.i](problems/3p_dis_void_mech.i) — mechanical pull tests driven by a saved composition field.

### Initial conditions

The `ic_2p/`, `ic_2pv/`, and `ic_3pv/` directories hold pre-generated Exodus
meshes used to seed dissolution and void runs. File names encode the parameters,
e.g. `2pv_0.3_ic_0.05_0.2.e` is a two-phase/void case at composition `0.3` with
void radius `0.05` and spacing `0.2`. The `*_ic*.i` input files regenerate these.

### Parameter sweeps

[problems/DBC.sh](problems/DBC.sh) launches a batch of simulations from the
template [problems/template/DBC_split_nd.i](problems/template/DBC_split_nd.i),
substituting composition, `chi`, degree of polymerization `N`, and random seed
into per-case output directories:

```
cd problems
./DBC.sh <number_of_processes>
```

## Post-processing

Python scripts in [problems/Post_scripts/](problems/Post_scripts/) read the Exodus
output and produce plots, images, and animations. Highlights:

- `process_exodus.py`, `py_post_img.py`, `py_post_separate.py` — render fields to images.
- `py_post_gif*.py`, `c_avg_gif.py` — build animations of evolving morphologies.
- `energy_plot.py` — plot total and interfacial free energy over time.
- `c_avg.py`, `minmax.py` — composition statistics.
- `void_dist*.py` — generate void distributions for initial conditions.
- `maxwell_triangle_legend.py` — ternary-composition colour legend.
- `LSA.py` — linear stability analysis.
- [problems/ts_num.py](problems/ts_num.py) — report timestep counts and time ranges of Exodus files.

These scripts require `netCDF4`, `numpy`, and `matplotlib`.

## Testing

Regression tests live in [test/](test/) and unit tests in [unit/](unit/). Run the
full suite with `./run_tests`.

## License

PANTHER is distributed under the GNU LGPL 2.1. See [LICENSE](LICENSE) for details.
