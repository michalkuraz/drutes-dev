# ADEnc changes: hydraulic reconstruction and particle retention

This note records the ADEnc update through November 2026. For the full
input format, see [ADEnc configuration](ADENC_CONFIGS.md).

## What changed

- Transport and hydraulic settings moved from `transport.conf` into the
  `transport` block of `drutes.conf/netcdf/netcdf.conf`. The standalone file
  and unused input templates were removed from the repository. Stabilization
  (SUPG and residual viscosity) still uses `stabilization.conf`.
- ADEnc now always uses conservative transport with impermeable banks and
  hydraulic reconstruction. There is no conservative-transport or
  reconstruction switch. The no-flow banks record must be `y`; the
  accepted-step mass audit is independently optional.
- A two-point outlet profile replaces the manually entered list of FE node
  pairs. The mesh supplies inlet labels; the outlet profile supplies the
  location of the downstream opening.
- Dispersivity and first-order retention can vary by mesh material. The
  retention rate has its own input line, immediately after dispersivity.

## Configure a case

The `transport` block goes after the channel count and before the
concentration boundary records. The coordinates below are the outlet for the
current Rhine mesh, not generic defaults for another catchment:

```text
transport
y                       # mandatory: no solute flux through banks
y                       # optional: write accepted-step mass audit
outlet profile
293200 5740000           # first endpoint: mesh x y [m]
293200 5751500           # second endpoint: mesh x y [m]; order does not matter
hydraulic settings
0                       # no distributed lateral source file
1                       # linear interpolation between daily mRM values
1.0                     # maximum relative L2 correction of the raw edge flows
y                       # remove disconnected components with no inlet or outlet
end transport
```

The `hydraulic settings` subsection is optional; without it the four
values are `0`, `1`, `1.0`, and `n`. An absent or unusable outlet is an
error, not a request to fall back to unreconstructed flux. A disconnected
component with only one port type remains an error even when the filter
is enabled.

The segment should run along the opening and slightly overlap the real
outlet. ADEnc considers edges with only one active neighbouring triangle.
Their midpoints must project onto the segment and lie within 2.5% of its
length on either side. The routed flow must cross the edge outward and
have the same across-profile direction as the edge's outward normal.
Where candidates overlap along the profile, the closer edge takes
priority. Neither endpoint order nor a minimum DEM elevation determines
the outflow direction. Check the selected opening whenever the mesh,
profile location, or active channel threshold changes.

## Water balance and tributary dilution

The hydraulic reconstruction produces one shared normal flux per FE edge
and balances fluxes, storage changes, and specified sources per active
triangle. It uses the shared PCG routine in `solvers`, with relative tolerance
`1e-10` and at most `20000` iterations per linear solve, fixed in
[`nchydroflow.f90`](../src/models/fluxLS/nchydroflow.f90). These are separate
from the transport solver in `solver.conf`. The remaining limit `1.0` is
a physical-input consistency guard: it allows a global relative L2
correction of at most 100%, not a 100% error in the water balance.

Concentration boundary IDs label the mesh inlets; at an inlet, the model
uses the adjacent routed discharge with the existing edge-length weighting.
For example, scenario 3a mixes contaminated Rhine water with clean Moselle
water (concentration zero at the Moselle inlet); scenario 3b reverses the
release. A tributary not resolved by the mesh can instead be represented
as a lateral source with its own discharge and concentration
(`hydraulic settings` policy `1` and `lateral.conf`).
Reconstruction does not invent missing water inputs, and an increase
in downstream `Qrouted` alone is not an explicit input of clean water.

## Material-specific dispersivity and retention

Use `uniform` to apply the same pair of lines to the entire mesh. For
different zones, use `materials` / `end materials` and repeat the pair for
each mesh material ID, in order from 1 through the largest ID:

```text
materials
200 0.2       # material 1: longitudinal and transverse dispersivity [m]
0             # material 1: no retention
500 2.0       # material 2: longitudinal and transverse dispersivity [m]
1e-5          # material 2: effective retention rate k_dep [1/s]
end materials
```

These numbers illustrate the format, not calibrated values. In a Gmsh
2.2 mesh, the material ID is the triangle's physical surface tag, not
its element number or inlet label.

The sink is `-k_dep * M`, where `M = H_eff * C` is the model's suspended
inventory per unit horizontal area. The equivalent depth `H_eff`
converts concentration to inventory; it is **not** the physical depth
particles fall through. Do not calculate `k_dep` by dividing a settling
velocity by `H_eff`. In a well-mixed physical water column an estimate
would instead use settling velocity divided by **real** depth. Without
those data, treat `k_dep` as an effective rate for calibration or
sensitivity analysis. Zero disables the sink; `1e-5 1/s` alone gives
a suspended-inventory half-life of about 19.3 hours.

Conservative assembly does not mean suspended mass remains constant in
the presence of a sink. The audit counts this physical loss separately
as `step_deposited` and `cumulative_deposited` in
`out/adenc_mass_balance.csv` when the audit is on and any material has
a positive rate. ADEnc does not evolve a sediment-bed stock or
model resuspension.

## Migration and verification

Move bank and audit records into `netcdf.conf`. Remove the old
conservative flag, explicit port list, PCG tolerance and iteration
records, and standalone `transport.conf`. Add the actual outlet profile
and confirm the mesh inlet labels. Keep the NetCDF forcing, DEM and
mRM routing inputs, and any applicable concentration series. Start each
run in a fresh working directory: DRUtES clears `out/` on startup.

The local build, 43 configuration checks, 41 hydraulic checks, and
initialization and matrix assembly for 15 example configurations passed.
Reversing the profile endpoints selected the same 17 outlet edges in
each of the six local scenario configurations. Before the final
configuration simplification, the hydraulic forcing windows also passed
over each scenario's configured duration (one or 14 days). These checks
verify input handling, edge selection and numerical water balance;
they do not validate real river hydraulics, dilution or sedimentation
against measurements.
