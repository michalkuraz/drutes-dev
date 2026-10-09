# Optional conservative ADEnc transport and inventory audit

Implemented locally 2026-10-08. Existing server runs, private run binaries,
original scenario inputs and NetCDF data were not changed. The repository's
new `drutes.conf/netcdf/conservative.conf` defaults to `n` / `n`: no change to
the default transport formulation or production output filenames in legacy mode.

## Enable in a NEW isolated run directory

Two logical records, with comments allowed above records:

```text
# Conservative transport formulation
y
# Accepted-step mass-balance output
y
```

`n` / `y` enables an audit of the original formulation without changing its
matrices, Dirichlet timing or time-step selection. A missing file or `n` / `n`
retains the original path. `riverbank.conf` must contain `y` for either audit
or conservative mode. Supported: single ADEnc equation, 2D P1 triangles,
standard Picard, transient implicit Euler with lumped or consistent capacity.
Other solvers/time methods and restart from backup are rejected, not silently
used with missing storage history. Schwarz still compiles/links; this optional
ADEnc mode does not claim support for Schwarz time stepping.

Never run `bin/drutes` inside a directory with results to preserve: main clears
that working directory's `out/*`. Keep large NetCDF inputs as symlinks to the
original read-only forcing; changing this option needs no data copy.

## Equation and local matrices

The original advective equation is

    H*C_t + q.grad(C) - div(K grad(C)) = S.

The optional conservative equation is

    d(H*C)/dt + div(q*C - K grad(C)) = S.

Here `q=Q/W*d` is depth-integrated water flux [L²/T], not velocity;
`H=Q/(W*v)` is the existing effective storage depth [L]; and `K` is the
existing depth-integrated dispersion tensor [L³/T]. This change does not
alter Q, W, v, channel directions, the fixed active-element mask or the
physical dispersivities. It does not enforce water continuity on mHM fields.

For test function Ni, integrate the COMPLETE solute flux by parts:

    integral Ni*d(H*C)/dt
      - integral grad(Ni).(q*C-K grad(C))
      + integral_boundary Ni*Jn = integral Ni*S.

DRUtES stores the negative equation. The local ADEnc hook therefore replaces
the old `-dt*integral Ni*q.grad(Nj)` with
`+dt*integral grad(Ni).q*Nj`. Diffusion retains its existing sign/operator.
The old-time load uses `-M(H_old)*C_old`, not `-M(H_new)*C_old`.
`M(H_old)` uses the SAME quadrature and capacity lumping choice as the solve.
Accepted quadrature depths are cached; rejected attempts never commit them.

This is a globally conservative continuous P1 weak formulation. It is NOT
a locally conservative finite-volume/DG method and does not guarantee a
nonnegative solution or eliminate concentration overshoots.

## Banks and external boundaries

Internal active/inactive banks have `Jn=0`, now a natural condition. The
legacy advective-form Robin bank addition must NOT also be applied in this
mode: doing both would double-count bank transport. The ADEnc dispatcher
skips that addition only when conservative mode is enabled.

Unlabelled TRUE external mesh edges remain distinct from internal banks.
At external outflow, zero normal diffusive flux leaves `Jn=(q.n)*C`. An
explicit `-dt*integral_edge Ni*(q.n)*Nj` matrix retains this solute loss.
Two-point Gauss quadrature uses the active-element-side hydrological lookup,
the stored channel direction and effective width. Tangential/dry/zero-width
edges contribute zero. Unlabelled external inflow with `q.n < -1e-12` raises
an error identifying the element/edge: it needs an explicit inlet condition.
No incoming concentration is guessed and no external edge is silently sealed.

Original physical Dirichlet port edges are excluded from those natural edge
matrices; their discrete reaction gives total solute exchange. In particular,
`C=0` at an inlet AFTER a pulse can allow dispersive loss. It is not a no-flow
condition. Only existing Dirichlet port records are supported in this version.

The first Rhine scenario still has physical inlet101 and reserved unused102,
NOT an explicit river outlet. An internal active-domain end remains sealed
unless an appropriate outlet is explicitly introduced. Conservative algebra
does not repair that geometry/physical boundary choice.

## Time and boundary history

Conservative steps are clipped to the next inlet data record and daily mHM
record boundary. They use the left trace at the end of the interval: a step
ending exactly at the six-hour pulse change still receives the interval's
`C_in=1`; the next step uses zero. Likewise, a step ending at midnight uses
the preceding hydrological record, and the next step updates storage/flux.
Inventory remains conserved across that storage change because the old H
is retained independently. The last inlet record persists beyond its time.
Series must contain exactly two finite columns and strictly increasing times.

`LSstate_time` and `LSprevious_time` preserve accepted boundary traces for
the four FEM solution columns. Assembly, Dirichlet elimination, shock/SUPG
evaluation and extraction no longer reuse the new inlet value as the old one.
No history or CSV row is committed on a rejected time step.

The inventory audit is for ACCEPTED solver states, not time-interpolated
observation records. Prefer observation times aligned with accepted steps
when comparing HC maps to its CSV, especially at discontinuities. Values at
an exact forcing transition describe the left accepted trace.

## SUPG and shock capturing

The SUPG temporal matrix still uses H_new, but its old-time load uses H_old.
Residual viscosity uses `(H_new*C_new-H_old*C_old)/(dt*H_new)` rather than
`(C_new-C_old)/dt`. Physical K and the configured stabilization factors are
unchanged. The local broken P1 residual still omits distributional coefficient
jumps on hydrological interfaces; this is not a full interface-residual SUPG.
Its test derivatives sum to zero, so its assembled temporal/spatial corrections
do not create global inventory. Artificial diffusion also has zero column sum.

## Output: out/adenc_mass_balance.csv

Only accepted steps are written. The diagnostic captures the last Picard
assembly BEFORE Dirichlet elimination, row scaling or preconditioning.
It reports:

- `time_s`, `dt_s`: accepted step endpoint and duration;
- `inventory`: signed M = integral_active H*C dA, using the solve's capacity
  weights (consistent or lumped, including represented Dirichlet-node storage);
- `initial_inventory`: represented initial M; an initially nonzero inlet gives
  a nonzero P1 initial inventory even when the interior initial concentration is zero;
- `cumulative_in`, `cumulative_out`: time-integrated natural exterior transport
  and assembled Dirichlet reactions, split by sign;
- `cumulative_source`: physical assembled source load, EXCLUDING the old-time
  capacity load (ordinary ADEnc source/reaction callbacks are zero);
- `error = M-M0-cumulative_in+cumulative_out-cumulative_source`;
- `relative_error`: error divided by `abs(M0)+cumulative_in+abs(cumulative_source)`
  (bounded below by floating-point tiny);
- `step_error`: corresponding single-step closure defect;
- `free_residual`: maximum absolute assembled unscaled free-node residual;
- `negative_nodal_inventory`: capacity-weighted negative nodal concentrations.

The final column is NOT the exact integral of the negative part of a P1 field.
Dirichlet exchange is a discrete reaction consistent with the entire final
operator, including stabilization, not a separate noisy boundary-gradient
estimate. Do not interpret it as a validated water discharge measurement.
For audit-only legacy runs, missing conservative volume/storage terms can
produce a nonzero budget error; the diagnostic does not force closure.

With C in [M/L³], M is [M]. For normalized unit concentration, report a
normalized inventory/volume-equivalent, not kilograms without an actual
concentration scale. The existing generic `do_masscheck` is not replaced;
use this ADEnc-specific CSV for this discrete inventory audit.

In conservative mode `ADEls_mass` exports `H*C` under
`solute_inventory_density` [M/L²]. The `conc_flux` filename is retained for
compatibility, but its corrected label is DEPTH-INTEGRATED WATER flux [L²/T].
It is still hydrological q, NOT the solute flux `q*C-K grad(C)`.

## Code changes

- `ncconservative.f90`: optional reader, event/trace timing, depth histories,
  conservative volume/capacity corrections and Dirichlet history.
- `ncbalance.f90`: unscaled local-equation capture, accepted-step audit/CSV.
- `ncglobvars.f90`: flags, depth caches and optional coefficient clock.
- `ncboundary.f90`: separates natural exterior edges and integrates outflow.
- `ncpointers.f90`: attaches optional hooks and selects bank/formulation path.
- `lsconstitutive.f90`, `init_netcdf.f90`: optional clocks/storage export/reader.
- `ncsupg.f90`: conservative temporal residual/load when enabled.
- `pde_objs.f90`, `fem.f90`, `femmat.f90`: null-by-default step and boundary
  history hooks; ordinary models retain the existing code paths.
- `Makefile`: new module dependencies; no blanket `-cpp` change.

No changes to generic capmat/stiffmat, linear solvers, row balancing,
preconditioning, Schwarz implementation, GUI, widths, directions or NetCDF files.

## Verification and remaining work

Production-assembly regression tests exercise both Euler capacity choices
and Galerkin/SUPG2/SUPG2+shock1 for 30 steps each with spatially discontinuous
flux and varying/jumping H. Closed-domain inventory and unscaled budget errors
must remain below 1e-11. A separate pulse test checks inlet Dirichlet history,
external outflow and cumulative budget. Rejection, daily/pulse clipping,
absent/off/audit configuration and unsupported restart/Schwarz are checked.
Audit-only intentionally detects legacy mass drift rather than hiding it.

The test driver uses real production modules and assembler with bounds checking;
a tiny dense solve isolates assembly correctness from iterative-solver tolerance.
This is not a long Rhine run or a calibrated physical validation.

The Rhine driver `check_adenc_bank_initialization --assemble` exercises full
initialization, initial export and conservative matrix assembly, then REJECTS
the diagnostic attempt without solving or time stepping. Use only a fresh
diagnostic output directory with conservative mode enabled. It succeeds on
the staged SUPG2 first-scenario inputs: 4968 active triangles, 723 banks,
46 inlet nodes, inlet series (3,2). The inherited 713 zero-width active triangles
and missing DEM coverage warnings remain. In zero-width/dry cases the existing
storage callback's fallback depth is retained; algebraic inventory there is
not evidence of a physically meaningful river volume.

Before paper-grade interpretation: identify a physically appropriate outlet,
assess the width/mask/storage reconstruction, run fresh short/full cases,
inspect the inventory CSV and compare pulse spreading/oscillations. Neither
global algebraic conservation nor this preflight establishes correct Rhine
hydraulics, local flux continuity, positivity or benchmark validation.
