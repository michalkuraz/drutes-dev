# Explicit tributary and lateral inflows for ADEnc

Implemented 2026-10-09. Opt-in, ADEnc-local, no changes to global FEM/solver/
Schwarz files. Existing source policy0 and missing/off hydroflow configs retain
their behavior. No benchmark inputs or NetCDF targets are changed automatically.

## What the available mHM file does and does not provide

The inspected `mRM_Fluxes_States.nc` contains coordinates, bounds, daily time
and `Qrouted` [m3/s] only. `dem.nc` contains elevation, not the routing graph.
There is no lateral-runoff, routing-storage or flow-direction variable in these
files. A neighbouring discharge may belong to a tributary, downstream cell,
another branch or a different basin. Discharge differences are not automatically
local runoff, especially under transient routing/storage.

The May2 case2 empirical storage prescription exceeds its upstream supply.
This implementation does NOT fill that deficit with an invented source. Its
actual upstream/tributary network and lateral series must first be established
from mHM routing/local-runoff data or independently justified boundary inputs.
Even adding genuine inflows may require reassessing the empirical depth closure.

## Two distinct representations (avoid double counting)

1. A resolved tributary in the active FE domain: existing explicit inlet edges
   in hydroflow.conf, distinct original concentration boundary ID and BC series.
   Its auto-Qrouted discharge is split across edges sharing that ID. The channel
   geometry and active-domain inlet must actually reach this tributary.
2. An unresolved tributary or distributed lateral inflow: source policy1 with
   explicit total water/contaminant inputs in `netcdf/lateral.conf`. Injection
   is spread uniformly per receiving area over a specified group of active FE.
   This is a declared source representation, not automatic tributary discovery.

Do not inject the same water through both a port and a lateral group. Multiple
groups can overlap additively, but an element may not be repeated within a group.
No withdrawals, evaporation, incoming diffusion or wetting/drying is supported
by this new source input. Those would require separate physical terms.

## Equations and time discretization

For group g with total discharge Qg and contaminant concentration Cg, let Ag
be the summed area of its receiving triangles. Each receives water source
Rg=Qg/Ag [m/s] and solute source Sg=(Qg*Cg)/Ag. Contributions from different
groups add. Qg is never injected in full into every receiving triangle.

    H_t + div(q) = R
    (H*C)_t + div(q*C - K*grad(C)) = S
    BF = A*(Rbar - (Hnew-Hold)/dt)

There is NO additional artificial -R*C sink in the conservative solute equation:
dilution follows from its storage/flux balance. R*C would arise on expanding
the equation into a nonconservative concentration form. Clean water (Cg=0)
adds water but no solute; matching inlet/source/initial C preserves a constant
concentration when the hydraulic reconstruction is compatible.

At input timestamps construct load Lg=Qg*Cg. Qg and Lg are each interpolated
linearly; C between records is therefore Lg/Qg where Qg>0, NOT independent
linear interpolation of C. For a trial step use the exact interval means of
these linear Q/load series in BOTH hydraulic reconstruction and transport.
Steps are clipped at all source knots; no source extrapolation or last-value
holding beyond the final timestamp. Series must cover the intended duration.
Trials at a left-limit endpoint use the actual solver start for knot detection.
Only accepted steps commit source and hydraulic states; rejection restores them.

## Configuration

Set the lateral policy (the record after y) in hydroflow.conf to1. Policy0 ignores
the optional lateral file. For1, `drutes.conf/netcdf/lateral.conf` is mandatory:

    number_of_groups
    number_of_receiving_elements number_of_time_records
    element_array_index_1
    element_array_index_2
    ...
    time_seconds total_Q_m3_per_second concentration
    ...
    [repeat group block]

Times start at0 and strictly increase; at least2 records per group. Coverage of
the configured simulation duration is checked on reading, and no runtime step
may extrapolate past the final record. All fields
must be finite; time, Q and C are nonnegative. Element indices are DRUtES FE
ARRAY indices, not arbitrary Gmsh tags, and must reference active triangles.
Blank and # comment lines are allowed. Inline # comments are removed; extra
fields/trailing records and repeated elements within a group are rejected.
Concentration uses the same convention as the transported unknown; no assumed
mg/L-to-kg/m3 conversion. The toy `.example` file is not a Rhine configuration.

## Implementation and verification

- nchydroflow: reads groups, distributes water/load per area, clips knots,
  changes the equality projection target and tracks accepted/trial sources.
- ncconservative: adds the negative-sign P1 source load locally, and its positive
  physical amount to the existing mass audit. Generic zerord assembly currently
  omits Ni weighting, so the new source deliberately bypasses that shared path.
- ncsupg: includes the SAME interval-average solute source in SUPG/shock residuals.
- Boundary conditions retain sealed banks and nonnegative outlets. Existing
  tributary inlet ports still use independently configured concentrations.
- check_adenc_hydro_window: includes R in its local water-balance residual and
  reports total lateral water/load in two additional output columns.

The standalone production-assembly tests cover both capacity variants,
Galerkin/SUPG2/SUPG2+shock1, C1 matching-source preservation, clean dilution,
varying concentration/load interpolation, zero/overlapping sources, interval
integration, knot clipping, rejection and invalid input guards. Existing
source-free regressions remain. The overlapping sub-metre toy geometry uses a
local origin to avoid large-UTM-coordinate cancellation at high divergence;
other original production fixtures retain their UTM coordinates.
These checks do not establish a physically validated case2 or positivity.

Run after a bounds-checked build, always in NEW test directories:

    python3 tests/run_adenc_hydroflow_checks.py --build build --output-parent runs

Case2 remains blocked until genuine inflow series/mapping and storage consistency
are established. No old run is restarted and no missing inflow is silently invented.
