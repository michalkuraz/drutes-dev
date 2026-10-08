# ADEnc SUPG stabilization

Optional residual shock capturing is now available independently; see
[ADENC_SHOCK_CAPTURING.md](ADENC_SHOCK_CAPTURING.md). SUPG itself is unchanged.

Implemented/tested 2026-10-07. This is residual-based streamline-upwind
Petrov–Galerkin stabilization for the single-equation, 2D P1 ADEnc model.
It is not a positivity-preserving limiter or a calibration of river transport.

## Configuration and scope

Optional `drutes.conf/netcdf/supg.conf`:

```
# enable [y/n]
y
# nonnegative, finite multiplier
1.0
```

The root default is enabled with multiplier 2 (changed at user request after
the factor-1 tests below). Missing file, `n`, or multiplier
zero disables the correction. Existing scenario folders and run copies were
not edited; to test SUPG in a new staged scenario, copy this settings file and
the rebuilt executable. Their NetCDF data can remain symlinked. Never rerun
DRUtES in an existing results directory: main deletes `out/*`.

Implementation is `src/models/fluxLS/ncsupg.f90`, bound by `ncpointers.f90`.
The shared `PDE_str` has an optional null-by-default element-stabilization
callback. Shared FEM dispatch calls it after capacity assembly, before adding
capacity to stiffness. Other models have no callback and receive no correction.
The two built Schwarz assemblers have the same one-line hook, but no numerical
Schwarz/subcycling validation was performed; do not infer support from linking.

## Equation, signs and time discretization

Using the existing storage depth H, depth-integrated water flux q and physical
dispersion tensor K, the broken element/quadrature residual is

`R = H (C_new-C_old)/dt + q.grad(C_new) - r C_new - f`.

Here r and f use the solver's reaction/zero-order callback signs; both are zero
in the current ADEnc setup. For P1 concentration, the Hessian is zero. Within
each hydrological cell the coefficients used by ADEnc are piecewise constant,
so the local strong diffusive residual is zero. No distributional jump term at
hydrological-cell interfaces or element interfaces is added. An element may
cross a hydrological-cell boundary; the existing quadrature/mapping is not a
conservative interface remapping. SUPG does not repair that limitation or the
observation-output discrepancy already identified in point 1.

The added weak residual is `integral tau (u.grad(N_i)) R`, with `u=q/H`.
All local contributions carry the existing negative equation sign in DRUtES:

- capacity: `-integral tau (u.grad(N_i)) H N_j`;
- spatial matrix: `-dt integral tau (u.grad(N_i)) (q.grad(N_j)-r N_j)`;
- RHS: capacity correction times the existing `elnode_prev`, plus
  `-dt integral tau (u.grad(N_i)) f`.

The consistent SUPG temporal correction is added AFTER base capacity lumping,
not diagonalized. Both implicit-Euler choices (lumped or consistent base
capacity) and steady-state assembly are handled. The previous nodal vector is
the one supplied by the existing solver; existing Dirichlet-history handling
was not changed. Physical dispersion, convection, storage, GMRES and boundary
conditions are not modified by SUPG.

With P1 gradients, streamline length is
`h_streamline=2/sum_i |d.grad(N_i)|`, `d=u/|u|`. Let
`D_parallel=d^T (K/H) d`. The time scale is

`tau = factor / sqrt((2/dt)^2 + (2|u|/h_streamline)^2
                     + (4 D_parallel/h_streamline^2)^2 + (r/H)^2)`.

The temporal contribution is omitted for steady state. Zero convection gives
zero SUPG. Tau has units of time; depth-integrated q must not be mistaken for
velocity. Directional dispersion is retained, without adding crosswind
stabilization. The combined time/advection/diffusion scale follows the form
used in [transient stabilized finite-element formulations](https://doi.org/10.1140/epjp/s13360-024-05481-9);
this reference supports the scale form, not a validation of ADEnc.

## Verification and observed limitation

Final verification: `python -m pytest -q` passed all 35 Python-level tests,
including the embedded Fortran checks and strip solves (2026-10-07).
Optimized `make build_target` succeeded with NetCDF. A separate isolated
`HAVE_NETCDF=no` build also succeeded, checking that the shared optional hook
does not introduce a mandatory NetCDF dependency. `git diff --check` passed.

The tests link real production modules in isolated directories, never run
DRUtES main and never touch the Rhine run outputs. Checks cover:

- independent strong-residual identity, including source, reaction, signs,
  current-time matrix and previous-time RHS;
- tau scales and temporal bound, zero flow and multiplier;
- real FEM matrix assembly in steady, lumped-Euler and consistent-Euler modes;
- constant-field preservation and null-hook no-op;
- absent/on/off configuration and rejection of negative/nonfinite multipliers;
- a synthetic steady convection-diffusion strip assembled with real
  `build_stiff_np`, capacity and SUPG dispatch, and solved as a small dense
  matrix. This is a numerical unit benchmark, not a river simulation.

Strip benchmark: 20 subdivisions in x, two rows of nodes, 40 P1 triangles,
u=(1,0), D=0.005, C(0)=0 and C(1)=1. Exact solution is
`C(x)=(exp((x-1)/D)-exp(-1/D))/(1-exp(-1/D))`.

| Method | Minimum C | Maximum C | Nodal RMS error |
| --- | ---: | ---: | ---: |
| Galerkin, SUPG off | -1.30653 | 1.0 | 0.24029 |
| SUPG, factor 1 | -0.13344 | 1.0 | 0.06049 |

SUPG strongly reduces the oscillation and error but does NOT eliminate all
undershoot on these anisotropic triangles. The regression assertion checks
that reduction, not an unsupported positivity guarantee. Sharp fronts may
require shock capturing/crosswind stabilization or a bound-preserving method.
No new Rhine simulation was performed during implementation/testing.
A separate scenario1-supg2-rI7lgx case was subsequently prepared in the
2026-10-07 simulation batch, retaining the 14-day first-scenario duration;
preparation does not authorize launching it.
