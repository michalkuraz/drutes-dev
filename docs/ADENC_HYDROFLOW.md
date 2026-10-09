# Compatible hydraulic reconstruction for ADEnc

Implemented 2026-10-08, optional and **disabled in the repository configuration**.
This is a hydraulic coupling/projection layer, not a new hydrodynamic model,
positivity limiter, or validation of the Rhine model. Existing server runs and
saved inputs/results have not been changed or restarted.

## Why this layer is needed

The original ADEnc callbacks combine Q sampled at individual quadrature points,
an element's geometric width and a channel-derived direction. The algebraic
relation q = H*v*d holds at a point, but does not enforce continuity of normal
water flux across FE edges or the balance H_t + div(q) = R. In conservative
transport, that residual can change even a spatially constant concentration.

The diagnosis of the stopped Rhine conservative runs found maxima at node1153,
where one triangle sampled Q approximately0.194 and369.768m3/s at different
quadrature points. This layer removes that particular mixed-point definition
and projects a preferred hydraulic field into a compatible edge-flux space.
It does NOT establish that the hydrological input or inferred geometry is exact.

## Supported scope and assumptions

- One equation, standard Picard, 2D P1 triangles, conservative implicit Euler,
  with the existing no-flow-bank mode enabled. Schwarz still builds; this mode
  does not enable or validate its numerical use.
- Fixed active geometry, positive storage. No wetting/drying or disappearance
  of hydrological data on an active element; those conditions fail explicitly.
- One centroid hydrological cell defines each element's Q, width and P0 depth.
  Initially node-active elements with missing/below-Qmin/nonpositive centroid Q
  or zero width are removed. This is a documented homogenization/mask change,
  NOT exact clipping of triangles against hydrological cell boundaries.
- Every retained connected FE component must have explicitly listed inlet and outlet
  edges. All other exposed edges, including exterior mesh edges, are sealed
  WATER banks. No downstream outlet is inferred from a channel endpoint.
- Source policy0 means **no distributed lateral water sources/sinks**.
  Multiple resolved tributaries are supported as separate explicit inlets;
  their discharge can be read from Qrouted. This is NOT automatic inference
  of an mRM routing graph or local runoff from neighbouring raster values.
- Source policy1 (2026-10-09) reads explicitly prescribed incoming water and
  contaminant-load series from lateral.conf; see ADENC_LATERAL_INFLOWS.md.
  It does not infer sources from Qrouted differences or implement withdrawals.
- No clipping of concentrations or artificial cancellation of div(q)*C.

Qrouted alone is not a full hydrological water balance. Genuine lateral runoff,
abstraction, or routing storage needs additional input and corresponding solute
sources before this approximation can be interpreted as the complete mHM model.
Do not set R=0 silently and claim that all downstream Qrouted values are preserved.

## Mathematics

The previous empirical velocity law is shared through hydro_velocity:

    v = vref * (Q/Qref)^0.4
    H0 = Q/(W*v), q0 = (Q/W)*direction.

The implementation retains the previous algebraic evaluation with W factors;
they cancel analytically. This defines preferred hydrology, not a constraint
that the corrected field must preserve this speed exactly.

Since2026-10-09 W is the **full hydrological-cell transverse span**, not a
narrow projected neighbour contact. H is an equivalent cell-area storage
coefficient, not physical channel depth. See ADENC_CELL_WIDTH.md for the
geometry/volume convention and its limits. Edge contacts/banks remain separate
flux constraints; no lateral source or actual mHM storage is inferred.

Each unique edge has one oriented integrated flux F [m3/s], positive outwards
from its first owner; the second owner uses -F. The element incidence matrix B
therefore gives its net outward discharge. At a transport trial step:

    B F = A_e * (Rbar_e - (Hnew_e-Hold_e)/dt).

Rbar=0 for policy0. For1 it is the exact interval-average explicitly supplied
water source, with the matching solute load added to the transport equation.

Initial reconstruction uses BF=0 for daily-step forcing, or BF=-A*dH/dt for
linearly interpolated forcing. Prescribed inlet edges have negative F;
banks have F=0; interior and explicitly listed outlet edges are unknown.
Outlet edges additionally satisfy F>=0: they may close, but cannot become
unconfigured inflows. This inequality was added on 2026-10-09 after the first
long test reversed some local outlet edges despite a positive total discharge.

Among fields satisfying these constraints, minimize

    sum_edges (F-Fpreferred)^2 / weight,
    weight = edge_length^2/(Aleft+Aright),

with Aleft at a boundary. Preferred interior normal flux is the mean of the
two centroid-derived vectors. A matrix-free Jacobi-preconditioned conjugate
gradient solve of B W B^T determines the correction. This auxiliary solver is
contained entirely in nchydroflow.f90; the transport solver is unchanged.
The generic projection kernel also handles compatible closed components with
a potential gauge and rejects incompatible closed balances. Production setup
requires explicit inlet/outlet components rather than inventing open boundaries;
the optional initial filter only excludes entirely unported components.

`hydro_project_outflow` enforces the exterior-outlet inequality with a monotone
active set: solve the equality projection, fix any negative outlet to zero,
then re-solve the complete balance until all outlets are nonnegative. It is
NOT clipping a finished hydraulic field. This algorithm is restricted to
one-owner exterior edges; the graph-Laplacian M-matrix makes binding monotone.
Interior one-way constraints would require a different active-set algorithm.
An incompatible component requiring net incoming water at its outlets still
fails; no water source or incoming solute concentration is invented. The
correction-limit guard remains unchanged. Logs include bound outlet count and
number of projection passes.

The RT0 reconstruction in a triangle of area A is

    q(x) = sum_edges F_out,e * (x - opposite_vertex_e)/(2*A),
    div(q) = sum_edges F_out,e/A.

It has a shared constant normal trace across each edge. Storage is P0 per
element. Together with the existing P1 quadrature and backward-Euler storage
history this makes the water balance compatible with constant-state transport,
for source-free flow and matching initial/inlet concentrations.

Report the correction norm ||F-Fpreferred||2/||Fpreferred||2. A configurable
limit rejects excessive correction; do not increase it merely to get a run to
start. This norm includes prescribed-port and sealed-bank corrections, so it
is a diagnostic of the total reconstruction, not only its interior solve.

## Configuration and ports

Optional `drutes.conf/netcdf/hydroflow.conf` starts with y/n. Absent or n leaves
the previous ADEnc behavior unchanged. For y, the remaining records are:

1. Integer lateral source policy:0=none,1=explicit lateral.conf.
2. Integer forcing mode: 0=daily step, 1=linear interpolation (recommended).
3. Positive finite relative water-balance tolerance, at most1e-6.
4. Positive integer maximum auxiliary PCG iterations.
5. Positive finite maximum relative edge-flux correction.
6. Number of port-edge records (at least2).
7. That many records: `node1 node2 kind discharge`.
8. Optional final `y/n`: remove entire shared-edge connected components with
   **no explicit inlet AND no explicit outlet**, once at initialization.
   Absent or `n` preserves the strict historical behavior. `y` logs component
   count, removed FE count and area, then rebuilds the active assembly mask,
   node participation, banks and hydraulic topology. Components with even one
   explicit port retain the existing inlet/outlet completeness checks. An
   explicit lateral source in a component selected for removal is an error,
   not silently discarded. Vertex-only contact is not a hydraulic connection.
   This is an explicit domain restriction, NOT a fabricated source or a repair
   of disconnected routing geometry. Connected tributary arms remain active.

Mode1 linearly interpolates daily Q and the depths derived from each daily Q
using the existing velocity law. Thus H is continuous across day boundaries;
the velocity law is exact at the data timestamps, not necessarily between them.
It requires a valid subsequent daily slice to bracket the entire requested time
interval. This is an explicit temporal-reconstruction assumption, not a claim
about mHM's within-day routing. Mode0 retains daily jumps, which can demand large
storage-compensating fluxes in short time steps and can trigger the outlet or
correction guard. Do not hide such a failure by raising the correction limit.

`kind=-1` is an inlet. A positive discharge prescribes that edge's Q [m3/s].
Zero uses the adjacent centroid hydro cell's Qrouted at the current forcing
date. All auto-Q edges sharing an original concentration boundary ID form one
inlet: each receives `Qcell*edge_length/sum_port_lengths`. Thus multiple edges
do not inject a full Q each. If their centroid cells have differing Q, the
total is their length-weighted average, not an independently measured port Q.
Do not mix prescribed and auto-Q edges under the same inlet ID. Distinct
tributaries should have distinct IDs and concentration boundary series.

`kind=+1` is an outlet; its discharge record must be0. The projection determines
its discharge from the water balance, rather than simultaneously enforcing an
incompatible mHM outlet Q. Its nonnegative discharge is enforced in the
hydraulic solve; the transport reversal guard remains as an invariant check.
Actual bank solute flux is naturally zero;
outlets retain q.n*C and zero normal dispersive flux.

Node numbers are **DRUtES FE array indices**, not necessarily original Gmsh
tags. Edges must belong to exactly one active triangle. Inlet endpoint nodes
must have the same real original concentration boundary ID (>100, not addedbc).
Only listed inlet nodes retain Dirichlet DOFs. Bank/outlet participating nodes
are free, including any old outlet Dirichlet tag. Inlet and outlet edges cannot
share a node. Duplicate, interior and absent port edges are rejected.

`hydroflow.conf.example` is a FOUR-TRIANGLE TEST mesh, not valid Rhine input.
Before enabling a Rhine case, identify all upstream inflows and the real
downstream outlet on its reconstructed active domain. The existing reserved
ID102 alone is not an outlet and the example's node numbers must not be copied.

## Runtime integration and stabilization

- init_netcdf reads settings, filters the domain after width preparation,
  rebuilds bank DOFs/topology and initializes hydraulic coefficients before
  initial outputs. Thus postprocessing and assembly use the same reconstructed
  field, not a separate visualization-only field.
- lsconstitutive flux/storage callbacks use hydro_value when enabled.
  Dispersion also follows the reconstructed flux through the existing callback.
- conservative_begin prepares the trial field with the same left-limit daily
  forcing/event clipping as the conservative storage update. The old depth
  remains accepted history. conservative_end commits only accepted hydraulic
  states; a rejected step restores accepted q/H for callbacks and outputs.
- Hydrological endpoint samples are cached by forcing day. Projection is cached
  by forcing day, exact storage-balance target and (for linear forcing) time,
  not by Picard iteration. A rejected retry with changed dt/time is re-solved
  when its hydraulic target or preferred field changes.
- ncboundary uses the shared reconstructed normal flux for outlets. All
  unlisted exposed edges are water-impermeable in this opt-in mode.
- ncsupg includes div(q)*C in the conservative strong residual and div(K).grad(C)
  for the spatially varying RT0 physical dispersion tensor. With grad(q)=bI,
  div(q)=2b and K=alpha_T*|q|I+(alpha_L-alpha_T)*q*q^T/|q|:

      div(K) = b*(2*alpha_L-alpha_T)*q/|q|.

  At q=0 the routine returns0; the norm-based mechanical tensor is not
  differentiable there. SUPG/shock remain broken-element approximations; this
  does not add all interelement diffusion jump residuals or guarantee positivity.
  Legacy disabled-mode stabilization is unchanged.

## Verification and remaining work

Tests in test_adenc_hydroflow.f90 and test_ncflux_width.py exercise the real
modules in temporary directories, never main or existing results directories:

- branched flux conservation, transient storage, closed-component compatibility;
- nonnegative outlet constraints, cascading bounds, zero net discharge,
  infeasible incoming demand, disconnected components and invalid interior
  bounds; 100 weighted graph cases compared with exhaustive bound-subset search;
- RT0 edge traces/divergence, analytic div(K) against finite differences;
- conservative SUPG constant-state temporal/spatial cancellation;
- real NetCDF daily slices, daily-step/linear forcing and production callbacks;
- shared multi-edge Qrouted inlet without duplicated discharge;
- production FEM assembly with lumped/consistent capacity and Galerkin,
  SUPG2, SUPG2+shock1, including a daily storage change and accepted inventory;
- rejected trial restores accepted hydraulic fields;
- absent/off legacy behavior and invalid configurations/port geometry/corrections.

`check_adenc_hydro_window.f90` additionally checks the actual production
hydraulics over a configured forcing window in 1800s steps without running
transport/main: finite positive storage, nonnegative outlet edges, sealed banks,
RT0 divergence and local storage/edge-flux balance. It requires a fresh directory.

These checks are not a completed Rhine simulation or evidence of calibrated
velocities, width, lateral hydrology, or concentration positivity. A proper
Rhine port configuration and short real-data assessment are still required.
At the initial verification stage only a private source/test snapshot was
uploaded; no production working tree was overwritten. Subsequent explicitly
authorized short enabled-mode runs are recorded separately below.

Architecture changes: nchydroflow (new), init_netcdf, lsconstitutive,
ncboundary, ncconservative, ncsupg and Makefile/test dependencies. No new shared
FEM hooks and no edits to global FEM assembly, main transport linear solver,
or Schwarz files.

Verification completed 2026-10-08: all35 local pytest tests passed, optimized
Mac build passed, and bounds-checked Linux build plus15 standalone acceptance
checks passed on hydrocalc in /mnt/stock/ncflux-hydroflow-check-20261008-DXKlWg.
The15 checks include three test-driver compilations, legacy bank/conservative
regressions, hydraulic kernel, daily/linear real-assembly cases and configuration
guards. The test correction limit100 is deliberately large for a unit-width
synthetic strip using a4m reference width; it is NOT a recommended Rhine value.
Logs are preserved under checks/hydro-checks-g3el0fbs. Reproduce after building:

    python3 tests/run_adenc_hydroflow_checks.py --build build --output-parent /path/to/new-checks

All three original Rhine configurations (galerkin, supg2, supg2-shock1) also
passed full FEM initialization, initial export and conservative matrix assembly
on the server, without a solve or accepted time advancement. This preflight
used absent/off hydroflow configuration and therefore verifies legacy-mode
compatibility, NOT reconstructed Rhine hydraulics. Its fresh logs are under
/mnt/stock/ncflux-hydroflow-check-20261008-DXKlWg/rhine-preflight/<case>/terminal.log.
The legacy preflight reports4968 active elements,723 bank edges and713 zero-width
active elements; passing assembly does not resolve those legacy hydraulic
limitations. The enabled reconstruction's centroid filtering and explicit
ports need a separate real-data assessment.

Subsequent enabled real-data verification (2026-10-08): two deliberately short
SUPG2 Rhine runs completed six300s steps each in
/mnt/stock/ncflux-hydro-rhine-20261008-tEDEXw. Constant-state C1 is preserved
to1.021e-9; assembly budgets and independent saved-field inventories agree at
roundoff. The600s pulse still has a small negative tail. Hydraulic correction
is about0.8114, and the outlet is a restricted numerical-test opening, with no
lateral sources. These results are not physical calibration or a positivity
proof. See ADENC_HYDROFLOW_RHINE_TEST_20261008.md, including the guarded failed
outlet selection and initial/boundary mass caveats. No long benchmark run was
launched; repository defaults and earlier outputs remain untouched.

Method background: Odsæter et al., *Postprocessing of Non-Conservative Flux for
Compatibility with Transport in Heterogeneous Media*, CMAME315 (2017),799–830,
https://arxiv.org/abs/1605.04076. This implementation's weights, centroid
homogenization and source/port policy are choices for ADEnc; do not transfer
that paper's convergence or physical-validation claims to this model.
