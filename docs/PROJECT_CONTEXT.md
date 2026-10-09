# DRUtES project context

LATEST2026-10-09 mandatory mRM routing: discovered Maryam's actual Rhine L0
dem.asc/fdir.asc/facc.asc and L11 graph in restart/mRM_restart_001.nc on
hydrocalc. ncrouting.f90 requires fixed drutes.conf/netcdf/mRM_restart_001.nc
for EVERY ADEnc init (missing/invalid=>ERROR STOP), maps by actual lat/lon to
Qrouted axes rather than dimension names or IDs, validates links and cycles,
and exposes cell/FE upstream/downstream/direction queries. Network-based UTM
unit vectors replace channel directions before widths/hydro preparation;
true mRM outlets retain the explicit channel exit vector with a logged count.
No automatic lateral rates/storage/FE ports inferred. Other models and generic
FEM/solver/Schwarz numerical sources unchanged. mHM filenames are configurable;
_001/_002 are documented domain1/domain2 conventions, not a universal suffix.
VPN/DNS failed after a125KiB routing-only extract downloaded. Local fixed-name
file is that exact static-field extract with source_file provenance, NOT the
full467MiB restart and NOT usable to restart mHM. Reader also accepts fullfile.
User authorized completing tests/commit OFFLINE; no further server actions or
case2/3 input changes. Data excluded from Git; case1 symlinks installed only
after commit; existing outputs/private executables never overwritten or rerun.
Artifacts runs/routing-tests-20261009-QRW8Ce: full checked+optimized builds,
synthetic failures/mapping checks, real1438links/1outlet/476headwaters/401
confluences and independent Python graph comparison; all41pytest tests and30
hydro/bank/conservative standalone checks passed. Production case1 init/dry
assembly passed with4correct observation FE; one additional zero-width FE
removed (714), projected raw-field correction~0.9455 is LARGE. This changes
future geometry/hydraulics, not validation or retroactive correction of saved
case1 results. First old preflight guard only rejected its inherited case2
inlet amplitude0.19055; staged unit case1 inlet fixed, no numerical code fix.
See ADENC_ROUTING.md; replacing full restart/server config installation waits
for restored VPN. Prior lateral implementation below is included in the tested
working version; no source-free case2 launch or calibrated physics claim.

LATEST2026-10-09 lateral-inflow extension: case2 initialization reproduced on
Mac and freshly rebuilt hydrocalc code. Its May2 rising-wave empirical storage
growth2253.4m3/s exceeds weighted inlet1903.1m3/s at t0; first two full-day
budgets also require negative outlet volumes under source policy0. This is a
coupling/storage/input-assumption incompatibility, not a demonstrated corrupt
mHM file or ADE time-integration error. No case2 transport run launched.
The inspected NetCDF contains only Qrouted/coordinates/time; DEM contains no
routing graph. User authorized explicit tributary/lateral-inflow support:
nchydroflow now accepts policy1 with mandatory lateral.conf groups of active
FE indices and time/Q/C records. Total group Q and Q*C are distributed per area,
linearly interpolated separately and interval-averaged; source knots clip dt.
Water balance, ADEnc-local P1 load, existing mass audit and SUPG/shock residual
use the same sources. Rejected trials restore accepted sources. Existing
multiple resolved tributary inlet IDs remain supported. No automatic Q-difference
source inference, withdrawals, generic FEM/solver/Schwarz changes, benchmark
source enabling, commit or server deployment. Root/default configs unchanged.
See ADENC_LATERAL_INFLOWS.md and lateral.conf.example (toy indices only).
Bounds-checked fresh Mac full build,30 standalone checks and40 pytest tests
passed. Matching C1, clean dilution, variable contaminant loads, overlap,
integration/knots/rejection and invalid input cases tested; synthetic inventory
relative errors<=2.15e-16. These are not physical validation of case2.
Source-free Rhine hydraulic regression over the October24 window also passed
all672x1800s trials through1209600s with nonnegative outlets; no transport/main
run. Only its staged diagnostic date changed from May2 to October24. Final
coverage guard requires source series to cover the configured duration at read.
Artifacts: runs/lateral-tests-20261009-jg4qef. Missing genuine mHM routing/lateral
inputs or independently justified tributary mappings still block case2 launch.

DEPLOYMENT requested2026-10-09: user now explicitly requests a local commit,
server git pull, compilation inside ~/drutes-dev and executable copies into
fresh Galerkin/SUPG2/SUPG2+shock1 benchmarks. Parent reserved:
/mnt/stock/ncflux-inside-benchmarks-20261009-5aQSXB. Existing outputs remain
untouched and shared NetCDF/DEM stay symlinked. Server's old dirty fixes were
archived/stashed and its HEAD fast-forwarded to55e11ad; an earlier preparation
stopped before compiling/launching because AppleDouble metadata differences
were treated as source changes. No new simulation has launched at this dated
snapshot. Continue with the explicitly requested commit/pull workflow; do not
rerun the first prepare.sh or use an existing results out/. Latest deployment
evidence will be in runs/inside-benchmark-deploy-20261009-kaMjg7/REPORT.md.

LATEST observation fix2026-10-09: reproduced the export mismatch with actual
production initialization and saved SUPG2 nodal fields. init_observe/old inside
selected FE5538/10189/14710/2877 rather than true10679/6415/10417/17993; all
selected points lie outside their chosen triangles. The P1 interpolation itself
agrees with independent barycentric interpolation when the correct FE is used.
User authorized replacing containment: src/tools/geom_tools.f90 now provides
deterministic translated/scaled barycentric inside for2D triangles, finite and
degenerate guards, scale-aware tolerance and boundary flag. Original routine
retained publicly as inside_shoot (Fortran identifiers cannot contain '-');
intervals/general polygons retain legacy behavior. No FEM matrices, ADEnc
physics, Schwarz sources, configs or original results changed. New tests:
tests/test_geom_inside.py and tests/fortran/test_geom_inside.f90. Bounds-checked
full build passes;60000 random containment comparisons against independent
quad-precision half-plane oracle, vertex/edge/order/scale/UTM checks,1000 random
full-Rhine-mesh unique-location queries and4known points pass. Isolated Linux
production-init checks now select all4correct active FE, with100/100 repeated
containment and roundoff agreement on saved day2/day12 concentrations. Local
bin/drutes rebuilt; server checkout and running/production binaries NOT updated,
only isolated diagnostic module linked. No new simulation launched.15hydraulic/
bank/conservative checks plus analysis/chart/width/dispersion regressions pass.
See runs/inside-barycentric-20261009-5hxi8k/REPORT.md and logs; old mismatch
paragraphs below are historical. Existing saved nodal fields remain usable;
old obspt exports remain wrong until reconstructed or rerun with the fix.

LATEST assessment2026-10-09: all3 bounded-hydraulics variants FINISHED exit0,
14days,779 steps,17 saved concentration snapshots,397--413s wall; no drutes
processes remain. Read-only check/download, no new runs or source/config edits.
Saved active c ranges: Galerkin[-.202053,1.072579], SUPG2[-.090123,1.025660],
SUPG2+shock1[-.001441,1]. Discrete mass-budget errors<=1.50e-13; independent
all-daily-Q P0H*P1c spatial inventories match audit<=4.21e-15. Shock1 strongly
reduces oscillations but broadens pulse; not positivity/calibrated accuracy.
Reconstructed original-coordinate samples show pulse at points1/2/3 around
days2/7/12 and risingpoint4 nearend; shockpoint3 still risingday14. Exported
obspt histories disagree for ALL4 points (especially1/3), despite using the
correct column2 c. Headers/source confirm columns time,c,H*c,q-components,cumq.
User specifically challenged quantity selection; verified SUPG2point1day2
column2=-.00210647 vs spatialP1c=.14142310; point3 exports0 but spatialday12
.06972162. Exact export/location cause NOT established/fixed. Fortran ES24.16
omits E at tiny three-digit exponents; diagnostic reader accounts for this,
no actual NaN hidden. See runs/hydro-assessment-20261009-YdM7fE/REPORT.md,
assessment.json, analyze.py and scientific comparison PNGs. Old RUNNING
paragraphs below are historical. Numerical mass balance != validated hydrology.

LATEST pre-shutdown2026-10-09: bounded SUPG2 benchmark FINISHED exit0 at1209600s
(14days). Full output assessment pending; final audit relative error-1.158e-13,
negative nodal inventory8.0522e4, so no positivity/calibrated-physics claim.
Galerkin and SUPG2+shock1 remain RUNNING at357676s/349954s (>4days), detached
workers verified PPID1/ownsession/noTTY. User requested leaving computations
until morning and will turn off MacBook; do not terminate/restart existing
cases. Server processes continue without Mac or SSH; no new local automation
is claimed. THREE_VARIANTS.md now records this newer status; prior RUNNING
SUPG2 paragraphs below are dated history.

LATEST follow-up2026-10-09: user requested all three stabilization variants.
Kept the already running bounded SUPG2 (no restart); separately launched
Galerkin and SUPG2+shock1 after identical-input diff (onlysupg.conf/shock.conf
excluded), private binary hash check and production FEM assembly/flag checks.
New cases in /mnt/stock/ncflux-hydro-outflow-fix-20261008-tDsU4W/variants-W1ELRu/
galerkin and supg2-shock1. Workers/model PIDs2614690/2614705 and2614692/2614706,
detached PPID1/ownsessions; started22:22:03UTC Oct8. All three confirmed RUNNING
with accepted steps after launch; SUPG2 then1013400s (>old3dayfailure), other
two>12h. Dated startup evidence, not completion. Same14day6hpulsephysics,
bounded hydraulics/sourcepolicy0, alpha200/0.2, identical private binary and
shared NetCDF symlinks. Original results/checkout preserved. Definitions are
the previously used Galerkin / SUPG factor2 / SUPG2+residual shock factor1,
NOT AFC/FCT. See runs/hydroflow-test-20261008-GGI9Sq/THREE_VARIANTS.md for
paths/provenance, later read-only status checks and comparison caveats.

LATEST2026-10-09: user authorized fixing the long-run outlet reversal.
Only nchydroflow.f90 numerical implementation changed: hydro_project_outflow
enforces F_out>=0 using an exterior-edge monotone active set, re-solving ALL
local water balances after binding negative outlets to0, not clipping flux.
Incompatible incoming water demand, dry-data and correction-limit guards stay;
no sources, incoming concentrations, FEM, solver or Schwarz changes.38 local
pytest tests, Mac build, Linux bounds-checked build/15 checks pass;100 random
graph cases match exhaustive bound-subset reference. Full production Rhine
hydraulic-window preflight PASSED672x1800s through1209600s, all outlets>=0,
independent global-scaled max water residual1.764e-11; per-element normalization
max1.182e-8 is also reported (PCG uses global normalization). Full FEM dry
assembly passes. New detached SUPG2 transport run started22:16:41UTC Oct8 in
/mnt/stock/ncflux-hydro-outflow-fix-20261008-tDsU4W/benchmark-supg2-bounded-14d-l7PyH3,
worker2614099 PPID1/SID2614099, model2614106. RUNNING at startup inspection,
not transport completion. Inputs byte-identical to old failed14-day run;
shared NetCDF symlinks unchanged, private binary96e67cd8..., old results intact.
Read runs/hydroflow-test-20261008-GGI9Sq/OUTFLOW_FIX.md and
BOUNDED_BENCHMARK_RUN.md. Root hydroflow remains n. This fixes local outlet
admissibility, NOT concentration positivity or full mHM/calibrated hydrology.
The failed run described below is historical and remains preserved.

LATEST long-benchmark status2026-10-08: the detached benchmark-supg2-14d-RS7VLL
already FAILED exit1 at the first trial after259200s,194 accepted steps (3days,
not14days). Error: Hydroflow outlet reversed: an inflow concentration is required.
Last trial logged water residual5.693e-12, correction0.82716, so the reversal
guard, NOT PCG failure/correction-limit/memory reporting, stopped this run.
Outputs/diagnostics preserved at the directory below. Detached launch worked
and continued after SSH disconnect; prior RUNNING/PID paragraphs are historical,
not evidence of an ongoing job. No guards disabled, sources/ports changed or
replacement simulation launched. User was warned before being told the run
could continue with Mac off. Resolve hydraulic boundary/source assumptions
before any authorized fresh rerun; do not rerun main in this results directory.

Follow-up 2026-10-08 long benchmark: user explicitly authorized the long run
and will power off MacBook. One tested SUPG2/conservative/hydroflow case started
detached via nohup+setsid in NEW isolated
/mnt/stock/ncflux-hydro-rhine-20261008-tEDEXw/benchmark-supg2-14d-RS7VLL.
Original14-day end1209600s, six-hour unit inlet pulse, bulk IC0, dt300/max1800,
same tested61 ports, linear H/Q, no lateral sources, alpha200/0.2,Qmin300.
Bounds-checked private binary unchanged; full initialization/matrix preflight
passed. Worker PID2612317 PPID1/SID2612317, model PID2612324, TTYnone,
stdin/dev/null and server log redirection. Separate post-launch SSH confirmed
RUNNING through78493.70209s about23:37Rome; this is dated startup evidence,
NOT a completion claim. worker.sh has flock/one-shot launch guard and records
status.txt,exit-code.txt,model/worker PIDs and start/finish UTC; no wall timeout.
Do not rerun main in this directory: outputs would be cleared. No new other
variants, automations or source/physics changes. NetCDF symlinked, previous
results intact. Guards can still reject future forcing/outlet reversals;
hydraulic correction is large. See runs/hydroflow-test-20261008-GGI9Sq/
BENCHMARK_RUN.md for provenance and read-only inspection commands. Older
paragraph saying no model remained applies to completed SHORT tests only.

Follow-up 2026-10-08 enabled Rhine tests: user authorized a short hydrocalc test.
New isolated parent /mnt/stock/ncflux-hydro-rhine-20261008-tEDEXw contains two
completed SUPG2/conservative/hydroflow runs, constant-state and short-pulse,
each six300s steps to1800s, exit0. Linear daily H, source policy0, explicit
44 original inlet101 edges and17 west-facing test outlet edges; not a calibrated
hydrological boundary. Initial wider cap preflight rejected a reversed corner
flux; guard retained, failed artifacts preserved, successful geometry explicit.
Saved C1 deviation <=1.021e-9; relative budget errors <=4.134e-15/1.435e-15;
independent P0(H)*P1(C) spatial inventory agrees at roundoff. Short600s pulse
final Cmin=-3.154e-6; no positivity claim. Hydraulic relative correction0.8114
is LARGE and requires hydrological assessment before benchmark interpretation.
Initial discrete pulse inventory is nonzero due prescribed inlet P1 values;
post-release strong C0 inlet can export mass, which audit includes. Root
hydroflow.conf remains n; old server runs/results/working tree and shared
NetCDF targets unchanged. No drutes processes remained after tests. Read
ADENC_HYDROFLOW_RHINE_TEST_20261008.md for exact inputs, limits and artifacts.

Follow-up 2026-10-08 compatible ADEnc hydraulics: added optional nchydroflow.f90
with shared oriented FE edge fluxes, RT0 reconstruction and P0 centroid storage.
Local backward-Euler water balance is imposed by a private matrix-free PCG
projection; existing velocity law is reused. Explicit inlet/outlet FE edges,
including multiple Qrouted-driven tributary inlet IDs, are required. No hidden
lateral source inference; policy0 is source-free except listed ports. Optional
linear daily Q/H interpolation avoids instantaneous storage jumps; the original
daily mode remains selectable. Flux/storage/bank/conservative hooks and SUPG/
shock residuals use the compatible field (including div(q) and analytic div(K)).
No edits to shared FEM, primary solver or Schwarz files. Root hydroflow.conf is
n: existing model inputs retain old behavior. README/config example documents
that toy port indices are NOT Rhine configuration and node IDs are FE array
positions. All35 local pytest tests and optimized Mac build pass. Bounds-checked
Linux build and15 standalone checks pass in the NEW isolated directory
/mnt/stock/ncflux-hydroflow-check-20261008-DXKlWg; logs under
checks/hydro-checks-g3el0fbs. No server working-tree overwrite, new simulation
or existing-output reuse. Previous conservative models were explicitly STOPPED
at user request; current pgrep finds no drutes processes. Older paragraphs below
describing running PIDs are dated history, not current status. A real Rhine port
configuration and short enabled-mode assessment are still required. See
ADENC_HYDROFLOW.md; no claim of concentration positivity or full mHM replication.
All three original Rhine variants also passed Linux full initialization,
initial export and conservative matrix assembly, with no solve/time advance,
in rhine-preflight/<case>/terminal.log under that new private directory. These
real-data checks used absent/off hydroflow configuration, so they establish
legacy compatibility only, not enabled reconstruction validity. Their legacy
713 zero-width active triangles remain a documented limitation.

Follow-up 2026-10-08 server conservative launch: user explicitly authorized
stopping the three older bank-only runs and launching all three new variants.
Private source snapshot exported to /mnt/stock/ncflux-conservative-20261008;
server ~/drutes-dev working tree preserved. Linux build, bank/conservative
regressions and full initialization/export/conservative assembly passed for
ALL three cases before the old launcher2606342 was sent TERM. Old batch
20261008T151502Z-banks-69WBpi outputs preserved, statuses INTERRUPTED.
New batch /mnt/stock/ncflux-conservative-20261008/runs/
20261008T164531Z-conservative-V9Flrm runs galerkin PID2608442, supg2 PID2608444,
supg2-shock1 PID2608446; launcher2608329. Conservative and balance flags y/y;
shared forcing symlinked to prior export, original hashes verified, no recopy.
First 6/6/7 accepted steps reached 2315/2315/2846 seconds with finite CSV values
and max relative budget error about 2.3e-15. This is STARTUP evidence only,
not completed-run conservation/physical validation. See
ADENC_CONSERVATIVE_SERVER_20261008.md for exact paths/provenance.

Follow-up 2026-10-08 conservative ADEnc: implemented optional conservative
d(H*C)/dt + div(q*C-K*grad(C)) transport and accepted-step unscaled inventory
audit. New ncconservative/ncbalance modules; three null-by-default shared
step/history hooks; no generic capmat/stiffmat/solver/Schwarz implementation
changes. Optional netcdf/conservative.conf defaults n/n, so root and existing
server runs retain previous physics. Single 2D P1 Picard transient Euler only;
restart from backup/Schwarz rejected for this mode, Schwarz still builds.
Variable/jumping-H closed tests, pulse/Dirichlet/outflow budget, rejection,
event clipping and legacy audit tested through production assembly. Fresh
Rhine conservative full-init/initial-export/matrix preflight passed without a
solve. No server interaction/deployment/restart or long simulation performed.
Existing zero-width, internal sealed outlet and physical reconstruction
limitations remain; see ADENC_CONSERVATIVE.md. Prior dated paragraphs below
may describe older Git/deployment states (adjacency fix is now committed in
3e9d563); do not treat them as current run-status evidence.

Follow-up 2026-10-08 server restart: stopped the remaining old Galerkin process,
preserving all results, pulled testing to 8548d416 on hydrocalc, and rebuilt
the Linux executable including Schwarz. All three initial restart attempts
failed at initial output in active_node_element: it iterated smartarray storage
capacity instead of its valid entry count %pos. Fixed locally and deployed by
SCP (NOT yet committed); production-fill regression tests now cover invalid
spare entries and inactive-only nodes. All 35 tests pass. Extended the Rhine
preflight with --full (real callback linking, FEM initialization and initial
spatial/observation export, no stepping) in a separate output directory; it
passes on Mac and Linux. Updated run_banks.sh requires this check for all cases.
Three simulations are now stepping in new isolated server batch
/mnt/stock/ncflux-three-variants-20261008/runs/20261008T151502Z-banks-69WBpi:
galerkin PID2606471, supg2 PID2606473, supg2-shock1 PID2606475. Run state is a
dated observation, not a completion/validation claim. Failed restart logs remain
in 20261008T150735Z-banks-ByCd17. Shared NetCDF forcing is unchanged/symlinked.

Follow-up 2026-10-08: optional netcdf/riverbank.conf implements ADEnc zero-total-
solute-flux internal banks via ncboundary.f90. Root enabled; legacy mode when
absent/off. Active FE triangles only assemble (optional PDE assembly_mask);
participating bank DOFs restored after legacy missing-DEM annotations. Edge
Robin correction uses one-sided hydrological q.n; actual external mesh/ports
remain separate. Only standard 2D Picard supported; Schwarz still builds.
35 tests pass, constant-H closed-box mass checks cover 180 assembled steps,
Rhine initialization succeeds (4968 elements, 723 banks). No server deployment
or production simulation. The scenario1 mesh has only inlet101; reserved102 is
NOT an outlet, so an internal active-domain end needs explicit outlet labeling
before long-run physical interpretation. 713 zero-width active triangles and
interior continuity/storage conservation issues remain. See ADENC_NOFLOW_BANKS.md.

Follow-up 2026-10-08 restart fix: server export run_all.sh now launches into
fresh runs/<UTC timestamp>-<unique suffix>/ directories, not cases/ outputs.
Preserves old results/markers; flock plus process detection prevents duplicate
active batches. Each attempt has a private binary and regenerated shared-data
symlinks. Updated server script/README/SHA256SUMS, --check passed, previous
script backed up in script-backup-rEDEBT3E. No simulation launched by this fix.
The original transfer tar.gz contains the old launcher; installed server and
local export directory contain the fixed version.

Follow-up 2026-10-08: exported three identical-physics Rhine scenario1 variants
(Galerkin, SUPG2, SUPG2+shock1) to exports/ncflux-three-variants-20261008.tar.gz
and scp to miguel@hydrocalc.science.fzp.czu.cz:/mnt/stock/. Unpacked as
/mnt/stock/ncflux-three-variants-20261008; native Linux build and SHA256 checks
passed. Shared NetCDF inputs stored once, relative symlinks in each case.
run_all.sh builds then launches three processes concurrently with rerun guards.
No server simulation launched yet. AFC/FCT is not included/implemented.

Follow-up 2026-10-07: optional ADEnc netcdf/shock.conf adds capped isotropic
residual viscosity via the existing ncsupg element hook. Root enabled factor1;
missing settings retain old behavior. No additional shared FEM/Schwarz edits.
Current Picard nodal values drive the viscosity; physical dispersion unchanged.
The active scenario1-supg2-rI7lgx run/binary/configs remain untouched. See
ADENC_SHOCK_CAPTURING.md for equations, nonlinear cost and smearing limitations.

Follow-up 2026-10-07: user selected SUPG multiplier 2. Root supg.conf changed
to 2.0; a new isolated scenario1-supg2-rI7lgx folder in the simulation batch
was prepared with 200/0.2 m dispersion and the original first-scenario 14-day
inputs. NetCDF files are symlinks, not copies. Production initialization
checked factor=2, date, inlet series and widths. Model has NOT been launched.

Follow-up 2026-10-07 SUPG: optional netcdf/supg.conf enables ADEnc-only
residual stabilization via ncsupg.f90 and a null-by-default PDE element hook
after capacity assembly. Includes consistent temporal correction and old-time
RHS; lumped base capacity is unchanged. Root defaults enable factor 1; existing
scenario/run copies remain unchanged. Algebra/real-assembly tests pass, and a
small synthetic strip shows min C improving -1.3065 to -0.1334, RMS error
0.2403 to 0.0605. This does not establish positivity or Rhine stability.
Schwarz hooks are present but not numerically validated. See ADENC_SUPG.md.

Follow-up 2026-10-07: ADEnc now supports direction-aligned dispersion with
separate alpha_L and alpha_T on the existing dispersivity line. One scalar
preserves legacy isotropic behavior. `ncdispersion.f90` provides strict record
parsing and tensor algebra; `ADElsdisp` uses the actual depth-integrated flux
vector and returns K=|q|[alpha_T I+(alpha_L-alpha_T)dd^T]. Root defaults are
200 m / 0.2 m, exploratory and uncalibrated. Prepared scenario folders and
finished run copies retain the original 2000 m setting. No new simulation
was requested here. See docs/ADENC_DISPERSION.md for units and limitations.

Follow-up 2026-10-04: prepared independent drutes.conf3a (Rhine-only continuous
unit-concentration release) and drutes.conf3b (Moselle-only release). Both use
Qmin=100, start 2015-06-01 and the inherited one-day duration/2000 m dispersivity.
Their identical meshes add five physical inlet-102 line records on existing
western exterior triangle edges, without changing nodes or triangles. Their
original grid-centre channel2.dat is extended upstream to that inlet; the
experimental DEM-guided channel2.dat-v2 remains comparison-only. Boundaries are
101 Rhine, 102 Moselle, 103 reserved inactive. An isolated production-module
initialization check (not main/time stepping/solve_pde) passed both setups:
44/6 surviving inlet nodes, positive widths and inward Moselle direction,
correct concentration callbacks at start/middle/end. The extended line's 1414
samples lie inside the FE mesh with Q >=116.89 on the initial date. No simulation
run; original drutes.conf and drutes.conf3 unchanged. Qmin selects all qualifying
cells, not only named rivers; equal inlet concentrations are not equal masses.
The inherited northern DEM coverage warnings remain. See each new scenario's
SCENARIO_README.md and scenario.json for setup and verification details.

Follow-up 2026-10-01: ADEnc now reads `channel_count` after initial concentration,
loading `channel.dat`, `channel2.dat`, through `channelN.dat` as separate
polylines. The reader requires two finite coordinates per data line, reports
file/line errors, and rejects duplicate consecutive points within each file.
Shared coordinates between files remain valid for confluences. Segment direction
selection defensively skips zero-length segments. Boundary count is derived from
the original mesh maximum ID plus its reserved inactive boundary; configuration
records must be ordered from 101 through that ID. Separate physical inlets can
load separate concentration files. See `docs/ADENC_CHANNEL_INPUTS.md`. This is
code support; scenario 3 still needs actual tributary geometry and inlet labeling.

Follow-up 2026-10-04: scenario 3 now includes `netcdf/channel2.dat` for the
lower Moselle and sets channel_count=2. Its nine downstream-ordered UTM 32N
points follow the routed-cell corridor, ending on segment 22 of the original
Rhine line. All 808 sampled locations are inside the FE domain and valid Q cells.
At 2015-06-01 the tributary Q is approximately 130.5--136.0 m3/s; the unchanged
Qmin=300 excludes it. Scenario 3 still needs an activation threshold/date decision
and a separate release boundary. See its SCENARIO_README.md and scenario.json.
Only a channel-reader driver was executed; no transport simulation was run.

Updated 2026-09-05 at the user's request: read the numerical source first,
then the GUI, and retain working context. Repository root on this machine is
`/Users/miguel/drutes-dev`. Snapshot: branch `testing`, HEAD
`2834aad` (`improved width computation`), preceded by `25467cf`
(`ADE nc fixed`). This is an architectural orientation, not a complete
line-by-line numerical correctness audit. No model run, build, authentication
test, or automated test suite was performed during this context refresh.

Follow-up completed 2026-09-06: implemented and checked contact-limited river
widths, set the requested root `Qmin` to 35, ran all 15 tests successfully,
and rebuilt root `bin/drutes`. The width section below describes this newer
working-tree state; HEAD above is the earlier committed baseline. Existing
simulation outputs and GUI project copies were not changed.

## Working boundaries

- Preserve existing code, configurations, outputs, and untracked work. At
  inspection, tracked files were clean but many GUI modules, documentation,
  configuration extras, build products, and outputs were untracked. Untracked
  does not mean disposable. No implementation changes were requested in this
  refresh; only this context document and its `AGENTS.md` entry point were added.
- The model starts by clearing `out/*` relative to its working directory in
  `src/core/main.f90`. Inspection must not casually execute `bin/drutes` on
  existing results. A future numerical test needs an explicitly isolated run
  directory with the required inputs.
- Never copy OAuth credentials into chat, source, notes, or screenshots. Native
  authentication already exists and should be continued, not reimplemented.
- Current conversation language is Czech; GUI labels and code are English.

## Numerical source: `src/`

There are 104 `.f90` files, totaling 51,212 physical lines including comments
and blank lines at this snapshot. These totals do not imply that every source
variant is built. `src/drutes_gui.egg-info/` is Python packaging metadata, not
part of the numerical solver.

DRUtES is a modular Fortran finite-element simulator for nonlinear, potentially
coupled transport/flow PDEs, with 1D, 2D/axisymmetric, and 3D infrastructure.
Individual model/dimension combinations have their own restrictions.

| Directory | Responsibility |
| --- | --- |
| `src/core/` | Main program, numeric kinds, mesh/global types, PDE types, global state and shared callbacks. |
| `src/pointerman/` | Connect model names and solver settings to concrete procedures. |
| `src/femtools/` | FEM initialization, quadrature, local capacity/stiffness/load assembly, global assembly, nonlinear/time stepping and mass/flux calculations. |
| `src/mathtools/` | Geometry-independent numerical helpers, quadrature and linear-solver interfaces including GMRES. |
| `src/pma++/` | Matrix classes, sparse/full storage, reordering, matrix I/O, direct and iterative linear algebra. |
| `src/decompo/` | Schwarz domain decomposition, subdomains, coarse levels and subcycling variants. |
| `src/tools/` | Configuration and mesh readers, initialization, geometry/projection utilities, output, logging, timing and inverse-model objective support. |
| `src/models/` | Model-specific readers, state, constitutive relations, boundary conditions and callback linking. |

Execution route: `main` -> `parse_globals`/mesh and solver inputs ->
`set_pointers` -> observations and `feminit` -> optional decomposition ->
`solve_pde` -> outputs and final diagnostics. `fem.solve_pde` writes initial
outputs and advances time through `pde_common%treat_pde` and adaptive time-step
handling. `femmat.solve_picard` assembles the coupled system, calls
`solve_matrix`, checks convergence, and accepts/rejects iterates.

`global_objs.f90` defines `node`, `element`, observation and integration-point
types. `nodes%element(node)%data` stores adjacent FE elements.
`elements%data` is connectivity; `%gc` centroids; `%areas` measures;
`%neighbours` adjacency; `%material` material IDs. Internal array positions
must not be assumed identical to arbitrary external mesh tags.

`pde_objs.f90` defines `PDE_str`, `pde_fnc_str`, and `pde_common_str`.
Procedure pointers provide dispersion, convection, elasticity/storage,
reaction/source, flux, initial and boundary conditions, and value/gradient
evaluation. The point descriptor distinguishes `gqnd`, `obpt`, and `ndpt`.
Changing one model callback can affect assembly AND postprocessing.

Model names dispatched by `manage_pointers.f90` include `RE`, `REstd`,
`REtest`, `Re_dual`, `ADE`, `ADEnc`, `heat`, `boussi`, `kinwave`, `freeze`,
`LTNE`, `ICENE`, and `REevap`. Several are explicitly developmental. The GUI
supports only a subset; do not equate the GUI's model list with solver scope.

Richards implementations live in `models/RE/`: `RE` uses total hydraulic head,
while `REstd` uses pressure head and has different boundary support.
`re_reader.f90` reads matrix/root parameters; `re_constitutive.f90` contains
van Genuchten/Mualem, constitutive tables, root sink and related functions;
`re_total.f90` supplies total-head operations. Other model directories cover
dual porosity, evaporation/vapour, solute ADE, heat, soil freezing, Boussinesq
and kinematic-wave flow.

Heat coupling is decided from the first value in `heat.conf`:
`heat_pointers.f90` creates one heat PDE, or Richards as PDE 1 and heat as
PDE 2. `heat_reader.f90` reads convection rows ONLY when Richards coupling
is off. Initial temperature may be a per-material numeric block or the
legacy `file` input; the GUI exposes per-layer temperatures.

`tools/postpro.f90` writes observation series, spatial profiles and Gmsh
data. Observation series begin with time, then solution, printed mass
properties, flux components, and cumulative flux. Exact column counts depend
on the model and dimension. `out/solver.time` is written when enabled around
linear-system solves; multiple solves can share a simulation time.

The root `Makefile` builds `bin/drutes`, objects in `build/objs`, modules in
`build/mods`, and logs in `build/logs/compile.log`. Defaults include gfortran,
optimization, implicit-none, single-image coarrays and default-real-8.
NetCDF support is conditional on detecting `nf-config` and `nc-config`
(`HAVE_NETCDF`). Consult the Makefile before building; source presence alone
does not imply inclusion or parallel execution in the default build.

## Current work: NetCDF river flux and hypothetical width

Relevant files are in `src/models/fluxLS/`, plus geometric helpers in
`src/tools/geom_tools.f90`. Inputs include
`drutes.conf/netcdf/netcdf.conf`, `mRM_Fluxes_States.nc`, `dem.nc`, and
`channel.dat`. Keep these user data intact.

- `ncglobvars.f90`: NetCDF axes/time/cache, Qrouted variable metadata, cell
  bounds, active-element flags and flow directions; hydro mesh and FE mapping.
  Current geographic zone is 32. Do not silently change coordinate assumptions.
- `netcdfflux.f90`: public subroutines `ncflux_get_xy` (original bilinear),
  `_bilin`, `_pw` (triangular piecewise linear), `_nn` (nearest valid one of
  four bracketing sample locations), and `_cell` (piecewise constant value of
  the containing hydrological cell using bounds). `_nn` is not a global
  nearest-valid search and is not the same as `_cell`.
- Each public lookup takes `(x, y, cur_hrs, qval, ok, errmsg)`, transforms UTM
  to lat/lon, looks up an exact integer-hour time index, and caches one slice.
  NetCDF data are read `(nlon,nlat,1)` then transposed to `(nlat,nlon)`.
  Handle `ok`; do not interpret a failed lookup as physical zero automatically.
- `init_netcdf.f90` uses `_cell` at mesh nodes for initialization/activation.
  A node is assigned the added boundary if lookup fails or Q is below Qmin;
  an FE element is inactive if all of its nodes have that boundary. This mask
  is constructed at initialization, not recomputed for every later time slice.
  The DEM-based node-deactivation block in `init_netcdf` is commented out,
  but `nctools::terrain_slopes` still marks nodes adjacent to missing DEM
  elevations as `addedbc` afterwards, without rebuilding `activeel`. Do not
  describe DEM-based boundary marking as completely disabled. This pre-existing
  behavior was not changed by the width update.
- `ncmesh.f90` constructs hydrological quadrilaterals from bounds and maps
  coordinates to UTM; `ncdem.f90` loads/interpolates elevations;
  `nctools.f90` provides coordinate conversion and slopes.
- `ncmap.f90::mapel` maps FE centroids to containing hydrological quads in
  `el2ncgrid`. This is separate from discharge lookup at an arbitrary point.
- Channel direction comes from ordered consecutive point pairs in
  `channel.dat`. For each active FE centroid, initialization now chooses the
  nearest segment using clamped point-to-segment distance, then stores its
  unit direction in `ncfluxdata%fluxvct`. It no longer requires the
  perpendicular projection to fall strictly inside some segment.
- User reported the channel holes disappeared after that direction-mapping
  change. Element 8777 was a historical debugging example. Recheck data at
  the actual simulation time before any new diagnosis; an old first-slice
  inspection is not evidence about another date. Neither unusual element
  numbers nor the screenshot alone prove missing NetCDF values.
- `lsconstitutive.f90::ncflux` uses `_cell`, flow direction, and width to
  return `Q * direction / width` (or zero on certain failures). Its output is
  named `conc_flux`, but this callback does not multiply by concentration.
  Do not infer its physical formula from the filename alone.
- Current time sampling there and in the storage callback is
  `ora_di_ini + int(time / 86400) * 24`: daily steps in seconds-to-hours
  conversion, followed by an exact NetCDF-time match.
- Nodal flux evaluation uses the first adjacent element for activation,
  direction and width. Gauss/observation evaluation uses the specified
  element. These choices matter for discontinuities and displayed results.

### Width: contact-limited implementation (2026-09-06)

`src/models/fluxLS/ncfluxarea.f90::ncflux_active_width(element_number)` returns
a cached effective width in UTM coordinate length units. Initialization calls
`ncflux_prepare_widths` after FE mapping and channel-direction assignment.
The cache uses the initial NetCDF slice: valid finite Q > 0 and Q >= Qmin.
Missing/dry/below-threshold cells do not supply active contacts. Rebuild the
cache explicitly if the mesh, directions, mapping or initial mask changes;
loading another daily discharge slice does not change the fixed geometry.

`ncwidth_geometry.f90::river_contact_width` computes actual shared edge
overlaps between a hydrological quad and its active neighbouring quads.
Each opening contributes `overlap_length * abs(dot(d, outward_edge_normal))`,
with normalized flow direction d. It sums inlet openings and outlet openings
separately (parallel branches), and limits the previous full-cell transverse
span by the smaller available inlet/outlet width. Partial faces are supported
by the geometry routine; structured-grid neighbours are gathered from the
eight surrounding cells so corner-only contacts are also recognized.

An absent inlet/outlet is treated as an open channel end. A known connection
only through a corner has zero width unless there is also a positive face
opening on that side. A direction tangent to all real contacts returns zero.
An isolated cell retains its full span because no connecting width is known.
Inactive/unmapped/invalid FE entries and zero directions return zero.
Initialization logs the number of zero-width elements mapped to flowing cells.
No epsilon-width replacement is used to hide a closed contact.

This is a geometric effective-width approximation, not a measured wetted
river width, a local cross-section through the FE centroid, or a conservative
remapping of Qrouted between cells. Inlet and outlet transport conservation
for varying directions/branch discharges would require a separate numerical
assessment. Flow directions and existing FE/hydrological mapping were not
changed. Both flux and effective-depth/storage callbacks use the same width.
The storage callback now guards failed lookups and Q <= 0 before division.

User explicitly requested Qmin = 35: root
`drutes.conf/netcdf/netcdf.conf` now contains 35 (m3/s for this Qrouted data).
The initial date is 2015-06-01 (NetCDF time index 20240, zero-based), not the
first data slice. At that date, Q >= 35 selects 196 cells forming one
face-connected component; the previous Qmin=0 selected 1,439 positive cells.

Verification on 2026-09-06:

- `python -m pytest -q`: 15 passed, including new Fortran geometry and real-module
  cache tests. Cases include oblique/partial contacts, branches, corner-only
  connections, rotated UTM coordinates, reversed orientation, both descending
  axes, zero/NaN/missing discharges, invalid mapping, and explicit cache rebuild.
  Integration tests build the real solver modules in temporary directories and
  run a dedicated test driver, never DRUtES main.
- `make build_target`: successful optimized build of root `bin/drutes`.
- A separate read-only geometry driver used the real NetCDF, Gmsh mesh,
  channel polyline and production mapping/direction/width routines: 35,137 FE
  triangles, 6,873 active by the initial node criterion; 5,875 of those have
  centroids in flowing cells. None of these 5,875 has zero width; 5,141 widths
  are reduced. Positive widths range about 408.746--16,492.418 m. Another 998
  node-active FE centroids are below Qmin/unmapped and get zero width, consistent
  with the cache's centroid-cell criterion. This checks width geometry, not
  complete FEM boundary handling or a transient simulation.
- No full model simulation was run; existing `out/` results are untouched.
  No existing GUI project executable or configuration copies were replaced.

Remaining related follow-ups, not addressed by this change: `velocity(Q,w)`
still algebraically cancels width for positive w; duplicate channel points
can yield zero direction; nodal evaluation still uses the first adjacent FE
element. The initial FE activation criterion (any active node) differs from
the centroid-cell criterion used by the width cache near river banks.

## GUI: Python modules and workflow

GUI/application code is split across `drutes_gui/`, root `pages/`, and
`drutespy/config/`; reading only `drutes_gui/` misses the editors.
These directories plus the two test modules total 4,394 Python lines at this
snapshot. `pyproject.toml` requires Python >=3.12 and `streamlit[auth]`;
pytest is the test extra. The declared Streamlit minimum is not a verified
compatibility guarantee for every API currently used.

- `drutes_gui/app.py`: entry point, theme and root `logo.png`, Google native
  login, project gateway, session navigation, save-all, subprocess terminal,
  results pages, configuration/output ZIP downloads.
- `drutes_gui/projects.py`: validates account/project path components,
  creates projects by staging and renaming copied defaults, resolves and lists
  projects, and produces configuration archives.
- `drutespy/config/parameter.py`: parameter types, schema definitions, values,
  original values and line positions.
- `parser.py`: strict positional parser, skips blank/whole-line `#` comments,
  supports counted lists and y/n booleans, rejects unexpected extra records.
- `writer.py`: edits modified values bottom-up and replaces the target via a
  temporary file. Unchanged records/comments are retained; replacing a list
  replaces its full line span, so do not promise arbitrary comment retention
  inside modified list blocks.
- `configfile.py`: UI-independent load/get/set/save API combining parser/writer.
  Model schemas: `global_config.py`, `mesh_config.py`, `solver_config.py`,
  `heat_config.py`, `matrix_config.py`, `root_uptake_config.py`.
- `pages/global_configuration.py`: global controller plus UI. Supports only
  `RE` and `heat`; saving forces 1D, internal mesh, pure output, Picard and
  implicit Euler plus other `FIXED_VALUES`. Never use this as a generic
  round-trip editor for every Fortran model/configuration.
- `pages/mesh_configuration.py`: live SVG of soil, mesh spacing dx and
  observation points; range/continuity/alignment messages below the title.
- `pages/solver_configuration.py`: linear solver and solver-time toggle.
- `pages/model_configuration.py`: heat, Richards matrix and root-uptake
  editors; shared continuation controls and boundary upload validation.
- Output pages: `richards_outputs.py`, `heat_outputs.py`,
  `solver_time_output.py`, `simulation_log.py`.

Navigation: Google login -> project selection -> project home -> global ->
mesh / solver / heat or matrix -> optional coupled matrix / root uptake ->
save all -> separate run page. `visited_configuration_pages` makes labels
change from Edit to Return and enables completion when required pages were
opened; visiting is intentionally not synonymous with editing. Individual
saves persist edits. `save_all` reloads/saves files on disk; it is not a
transaction committing every outstanding widget value.

The physical path is lowercase `drutes.conf/mesh/drumesh1d.conf`; some display
labels say `drumesh1D.conf`. Preserve correct spelling on case-sensitive
servers. Mesh layer count drives read-only counts and row resizing in the
model pages; those files are updated when saved, not immediately on mesh edit.
Fortran material counts use maximum material ID, so arbitrary material IDs
need care when assuming one row per layer.

Richards UI retains hidden RCZA tokens; initial options map to `hpres`,
`H_tot`, `theta`. Bottom/top are stored IDs 101/102. Boundary choices map to
1 Dirichlet, 2 Neumann, -1 height-defined Dirichlet, 3 free drainage,
4 seepage, 5 atmospheric. Time-dependent inputs use 101.bc/102.bc beside the
model config; ordinary inputs require two columns, atmospheric three and
forced time dependence. Heat offers only Dirichlet/Neumann with two-column
files. Root uptake exposes per-layer Feddes h1, h2, h3, h4, Smax.

Heat coupling hides AND removes convection rows from the saved file;
turning it off exposes and inserts them again. This matches the Fortran
reader. Coupled results have separate heat and Richards tabs.

Richards units: time sec/min/hrs/day and length mm/cm/m. Heat: fixed sec/m.
Length is stored as `# GUI length unit:` in global.conf, not an extra
positional Fortran record. Unit selection controls metadata/labels; it does
not rescale all physical numeric inputs. Keep physical parameters consistent.

### Authentication, projects, and running

Native `st.login()`, `st.user`, `st.logout()`; required `[auth]` fields are
`redirect_uri`, `cookie_secret`, `client_id`, `client_secret`,
`server_metadata_url`. Local intended redirect is
`http://localhost:8501/oauth2callback`; Google discovery is
`https://accounts.google.com/.well-known/openid-configuration`.
Local secrets belong in `.streamlit/secrets.toml`, template in
`.streamlit/secrets.toml.example`. Credential contents were not read during
this review. Current `.gitignore` ignores secrets, `user/`, and `backups/`;
the ignore rules were checked. The ignore file itself was untracked.

Projects live under `user/<normalized-google-email>/<project-name>/` with
their own `drutes.conf/`, `bin/drutes`, `out/`, and `project.json`.
They copy current root defaults and executable at creation. Rebuilding root
`bin/drutes` does NOT update existing projects' private executable copies.
Configuration ZIP includes configs, binary and metadata; output ZIP includes
the entire project's `out/`. Archives are built in memory and skip symlinks.

The runner uses `subprocess.Popen` with project cwd and merged stdout/stderr.
A daemon reader thread accumulates lines; a 0.5-second Streamlit fragment
renders an escaped fixed-height terminal and running-only kill control.
Termination targets the process (terminate, then kill if needed), not a
separate process group. This is session-local execution, not a durable job
queue. The GUI prevents switching project/re-editing through the run-page
buttons while running and stops the process on explicit logout.

Documented local launch (paths exist, but server not started in this review):

```sh
cd /Users/miguel/drutes-dev
source /Users/miguel/.venv/bin/activate
streamlit run drutes_gui/app.py
```

Browser: `http://localhost:8501/`.

### Outputs, checks, and documentation

- Richards series: `out/obspt_RE_matrix-N.out`; column names inferred from
  comments. Multiple selected points each have independent property selectors.
  Profiles: `RE_matrix_{press_head,sat-dg,theta,flux}-N.dat` with node ID,
  coordinate, value. Index 0 is initial; configured observation times follow;
  the expected last index is observation count + 1 for end time.
- Heat series: `obspt_heat-N.out`, exactly four columns for this 1D GUI:
  time, temperature, heat flux, cumulative transfer. Profiles:
  `heat_temperature-N.dat` and `heat_heat_flux-N.dat`. Display labels use
  seconds, metres, degrees Celsius, W/m² and J/m².
- Charts allow multiple properties on one axis with comma-separated labels;
  there is no normalization into a common physical unit or independent axes.
- `solver.time`: optional separate timing display/download. `DRUtES.log`:
  parsed timestamped events, metrics and rendered report plus original-file
  download. Do not treat that rendered subset as a complete diagnostic log.
- Existing automated tests cover global parser/writer behavior and project
  creation/path validation/archives (`tests/test_global_config.py`,
  `tests/test_projects.py`), not a full numerical or browser regression suite.
- Manual source: `docs/manual/drutes_gui_manual.tex`; screenshots:
  `docs/manual/images/`; PDF copies: `docs/manual/build/drutes_gui_manual.pdf`
  and `output/pdf/drutes_gui_manual.pdf`. The TeX graphic path supports the
  repository root and manual directory. Render/layout freshness was not
  checked in this review.

Other static-review follow-ups, not changes made now: boundary validation
checks the expected column count but converts only columns 1 and 2, even
for atmospheric data; mesh alignment errors are displayed without blocking
save; coupled-heat/root-uptake navigation and completion requirements are
not handled uniformly on every route. Review these with targeted tests if
asked to change the related workflow.
