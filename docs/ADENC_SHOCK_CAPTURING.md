# ADEnc residual shock capturing

Implemented 2026-10-07 as an optional addition to SUPG, not a bound-preserving
limiter. Existing run directories, private executables and forcing data are
unchanged. The running SUPG-factor-2 case retains its original binary.

## Settings

Optional `drutes.conf/netcdf/shock.conf` has two records:

```
# enable
y
# nonnegative finite multiplier
1.0
```

Missing file, `n`, or factor zero disables shock capturing independently of
SUPG. Root defaults enable factor 1; existing configurations without the file
retain their behavior. A future isolated run needs this file and the rebuilt
binary. Never launch main in an existing output directory: it clears `out/*`.

## Formulation

Using the current Picard iterate C* (column 2), previous-time C_old (existing
`elnode_prev`), storage H and depth-integrated flux q, the normalized residual is

`r* = (C* - C_old)/dt + (q.grad(C*) - reaction*C* - source)/H`.

The temporal term is omitted for steady state. P1 gradients are formed from
nodal values and production basis gradients, avoiding large-UTM plane fitting.
At each quadrature point:

`g = grad(C*)`, `d_g = g/|g|`, `h_g = 2/sum_i |grad(N_i).d_g|`, `u=q/H`,

`nu_sc = (h_g/2) min(|u|, factor |r*|/|g|)`.

Zero gradient, flow or factor gives zero viscosity. The advective cap prevents
unbounded viscosity near flat fields, but does not guarantee a discrete maximum
principle. Factors above one cannot exceed that fixed cap. This is an isotropic
residual-viscosity variant, not a named crosswind-only algorithm. For the
general residual shock-capturing family see
[Stabilized finite element methods with shock capturing](https://www.sciencedirect.com/science/article/pii/S0045782502002220).
The particular cap/length choice here is exploratory, not calibrated.

nu_sc has units m2/s; H*nu_sc is depth-integrated artificial diffusion. Using
DRUtES's negative equation sign, the added spatial matrix is

`A_sc(i,j) = -dt integral H nu_sc grad(N_i).grad(N_j)`.

Viscosity is frozen at the current Picard iterate and diffusion acts implicitly
on the new unknown. Capacity, RHS and physical dispersion callback are unchanged.
The local correction is symmetric, negative semidefinite in this sign convention,
annihilates constants and has zero row/column sums. These properties do not prove
global conservation of the hydrological projection.

## Limitations

- Single 2D P1 ADEnc equation, existing steady/implicit-Euler assembly.
- Reuses the existing hook; no additional shared FEM/Schwarz edits.
- Broken residual omits coefficient/interface jumps and P1 diffusive Hessians;
  no conservative flux remapping is added.
- Smooth nonconstant fields may acquire viscosity. A steady thin boundary
  layer can saturate the cap; front smearing is an expected tradeoff.
- No concentration clipping, positivity guarantee or bound enforcement.
- Nonlinear convergence may require more iterations or smaller adaptive steps.
  Synthetic tests use underrelaxed Picard; production relaxation/convergence
  rules are unchanged. Rhine performance and Schwarz/subcycling are unvalidated.
- Does not fix the observation-point export discrepancy.

## Verification

Tests link production modules in temporary directories without running main.
They cover residual/time/source/reaction algebra, viscosity scaling/cap,
zero residual/gradient/flow/factor, sign, partition of unity, real-hook assembly
with SUPG off and shock on, unchanged capacity/RHS, absent/off/on settings and
rejection of negative/nonfinite factors. The existing 40-triangle steady strip
is solved by repeated real assembly and dense linear solves. The nonlinear
test checks convergence and reduced undershoot, not improved exact accuracy
over SUPG. Smearing must be assessed in a separate authorized Rhine run.

Full suite: 35 passed (including embedded Fortran checks). The steady strip
uses SUPG factor 1, shock factor 1, and 0.5 Picard underrelaxation; the shock
fixed point converged in 30 iterations to tolerance 1e-9.

| Method | Minimum C | Maximum C | Nodal RMS error |
| --- | ---: | ---: | ---: |
| Galerkin | -1.30653 | 1.0 | 0.24029 |
| SUPG | -0.13344 | 1.0 | 0.06049 |
| SUPG + shock capturing | -0.04157 | 1.0 | 0.09626 |

Thus undershoot is reduced by about 69% relative to SUPG, but RMS error grows
about 59% because of smearing. This is not proof of suitability for the Rhine.
