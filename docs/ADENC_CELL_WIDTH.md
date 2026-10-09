# Cell-scale effective width for homogenized ADEnc

Implemented2026-10-09. This replaces the contact-limited width used in ADEnc
flux/storage callbacks; no configuration option, arbitrary width floor or
physical-channel-width calibration is introduced.

## Definition

For a hydrological quadrilateral with UTM vertices x_i, normalized direction d
and transverse unit vector n=(-d_y,d_x), use the full transverse span:

    Wcell = max_i((x_i-xcentre).n) - min_i((x_i-xcentre).n).

`ncwidth_geometry::river_cell_width(vertices,direction)` implements this
definition, with finite/zero-direction/degenerate-geometry guards. It is
translation/rotation invariant and normalizes the supplied direction. For
an axis-aligned rectangle of sides dx,dy:

    Wcell = dx*abs(d_y) + dy*abs(d_x).

Thus a direction becoming tangent to a neighbour face cannot reduce the
entire cell's storage width to a few metres. This is **homogenized cell
geometry**, NOT a measured wetted river width or a cross-sectional depth.

`ncfluxarea::ncflux_prepare_widths` caches this width for each mapped active
FE, using its hydrological cell and supplied direction. Dry/missing/below-Qmin
centroid cells remain excluded. Initial activation/cache semantics and explicit
rebuild requirements are unchanged. Different FE at the same true routing
node inherit the same direction/width; true routing outlets still retain the
explicit channel exit direction. Geometry remains fixed over the simulation.

All existing consumers of `ncflux_active_width` automatically use Wcell:
raw/hydraulic flux, equivalent storage, legacy constitutive/boundary callbacks.
No alternate contact width may be used in just one of those formulas.

## Storage and flux

Preferred fields remain

    q0 = (Q/Wcell)*d,
    v0 = vref*(Q/Qref)^0.4,
    Heff = Q/(Wcell*v0).

Width cancels from the existing empirical velocity expression. This change
removes the artificial amplification of q0 and Heff from near-zero contact
widths; it does NOT directly change that empirical velocity law. Corrected
RT0 q/Heff need not equal v0: projection magnitude remains a diagnostic.

Heff is an equivalent water-volume-per-horizontal-area coefficient [m], not
physical water depth. Across a whole hydrological cell of horizontal area A:

    V = A*Heff = (Q/v0)*Leff,  Leff = A/Wcell.

This makes the volume representation consistent with the chosen effective
width/length convention, not with measured channel bathymetry or actual mHM
routing storage. Summing over centroid-classified FE approximates that area;
it is NOT exact hydrological-cell clipping/remapping.

## Contacts and conservation

`river_contact_width` remains a tested legacy/diagnostic API. It is no longer
called when preparing ADEnc's storage/flux width. Its nearly tangent-face and
corner restrictions must not inflate Heff throughout a cell.

When hydroflow is enabled, shared FE edge fluxes, explicit ports, sealed banks,
RT0 reconstruction and backward-Euler water balance remain enforced by
nchydroflow. No edge constraint, source, Qrouted data, timestep, stabilizer,
transport/FEM/Schwarz solver or conservation audit is changed. Corner-only
disconnected FE components are NOT connected by this scalar-width definition;
existing explicit-component/port/feasibility guards remain active.

Cell-scale width is not a guarantee of feasible source-free water balance,
positivity, realistic corrected velocities or calibrated transport. Genuine
lateral inflows and storage assumptions may still require separate treatment.
Do not fabricate inflows or relax guards to force benchmark completion.

## Motivation and checks

Routing-enabled case1 previously generated W~7.4m on kilometre-scale cells,
Heff~68m and ~157km3 of derived water storage. At day4+30min the surrogate
storage demand20824m3/s exceeded inlet371m3/s, making nonnegative outlets
impossible. These were coupling-derived values, not demonstrated mHM errors.
Preserved diagnostics: runs/routing-window-local-20261009-v5KQW5/.

Geometry/cache regression tests cover full versus contact spans, nearly
tangential directions, reversed flow, every degree of rectangle direction,
rotated/translated UTM coordinates, invalid cells and unchanged masking/cache
semantics. Production checks must use fresh directories, never existing out/.
Before further server transfer, require a fresh checked build, full automated
suite, standalone hydraulic/bank/conservative/lateral checks, all three variant
init/assembly checks and the full14day hydraulic forcing-window diagnostic.
These are implementation/feasibility checks, not physical validation.
