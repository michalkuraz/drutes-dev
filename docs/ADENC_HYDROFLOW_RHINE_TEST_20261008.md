# Enabled Rhine hydraulic-reconstruction tests, 2026-10-08

## Outcome

Both deliberately short SUPG2 tests completed six accepted steps to1800s,
exit0, with no solver failure. Constant concentration is preserved to about
1e-9 on this real mesh with time-dependent reconstructed storage. Independent
spatial integration agrees with the assembly inventory audit at roundoff.
This verifies the tested numerical water/solute coupling; it is NOT a
calibrated Rhine simulation, a positivity proof, or a14-day benchmark.

| Check | Constant C=1 | Ten-minute pulse |
| --- | ---: | ---: |
| Model duration |1800s |1800s |
| dt / accepted steps |300s /6 |300s /6 |
| Real runtime reported by model |4.84s |4.88s |
| Maximum relative assembly-budget error |4.134e-15 |1.435e-15 |
| Max saved active-node deviation from C=1 |1.021e-9 |not applicable |
| Final saved nodal minimum |0.9999999989785 |-3.154e-6 |
| Final saved nodal maximum |1.0000000006839 |0.231173 |
| Max saved spatial/audit inventory difference |2.384e-6 |1.164e-10 |

Inventory differences are in the model's consistent concentration-volume
units, not an independently established kg scale. The constant-case inventory
is about7.73e8; its largest independent-integration discrepancy is relative
about3.1e-15. The small negative pulse contribution at1800s is0.00873,
about3.3e-8 of total inventory262877.96. No concentration values were clipped.

## Provenance and deliberate test changes

Server: miguel@hydrocalc.science.fzp.czu.cz.
New isolated parent:

    /mnt/stock/ncflux-hydro-rhine-20261008-tEDEXw

Folders constant-state/ and short-pulse/ each have private executable,
drutes.conf/ and initially empty out/. No old run was restarted or overwritten.
Source/build from the separate bounds-checked snapshot
/mnt/stock/ncflux-hydroflow-check-20261008-DXKlWg/source, not ~/drutes-dev.
Both runs used a300s wall-clock limit; neither reached it.

Inherited settings:2015-10-24, Qmin300m3/s, alpha_L200m, alpha_T0.2m,
SUPG factor2, shock off, conservative/audit on, lumped base capacity.
Hydroflow: source policy0, linear daily Q/H interpolation, tolerance1e-10,
20000 maximum projection iterations and correction limit1.0 unchanged.
Constant case changes initial C from0 to1; original inlet remains1 throughout
the short interval. Pulse case keeps initial bulk C0 but replaces the original
six-hour benchmark release by a600s diagnostic release, followed by C0.
Both use1800s end time, maximum dt300s and saved times600/1200/1800s.

Large mRM and DEM inputs remain symlinks to
/mnt/stock/ncflux-three-variants-20261008/shared/. No data recopy or alteration.

## Geometry, port selection and guarded failed attempt

Production centroid filtering removes714 of the4968 formerly node-active
triangles, leaving4254 triangles and2578 participating nodes in ONE connected
component. All retained triangles have positive effective width. This is P0
homogenization and a changed mask, not exact cut-cell river geometry.

The original inlet101 supplies44 exposed edges (FE indices485--529), combined
using the configured length-weighted Qrouted policy, NOT a full Q on each edge.
Initial combined inlet discharge is366.7185m3/s. The controlled downstream
opening uses17 exposed west-facing edges with midpoint x<294000m,
y>5740000m and outward-normal dot final-channel-direction>0.25.
Their exact61 total port records are in hydroflow-cap.conf in the server parent
and copied in each successful test's drutes.conf/netcdf/hydroflow.conf.
All unlisted exposed edges remain impermeable water/solute banks.

An initial wider opening also included the lower end-cap corner. Its preflight
was rejected by the normal guard: edge10970--10971 gave inward flow22.314m3/s.
The failed outputs remain in ports-preflight/. The successful cap excludes
four lower-corner edges geometrically by the y threshold, rather than disabling
the guard or increasing the correction limit. cap-preflight/ then passed full
initialization, output and assembly before either simulation was started.

This opening is an explicit NUMERICAL TEST boundary, not a measured Rhine
cross-section. Additional tributary inflows/lateral sources are not inferred.

## Remaining limitations

- Auxiliary water-balance residual is about6e-12 relative, with928--929 PCG
  iterations, but relative correction ||F-Fpreferred||/||Fpreferred|| is about
  0.8114. That is a LARGE hydrological correction, not evidence of field accuracy.
- In the diagnostic domain the initial storage rate is about-235.44m3/s,
  requiring net outlet about602.16m3/s from inlet366.72m3/s. This balance follows
  the inferred P0 depths and linear forcing; it does not reproduce all mHM Q.
- SUPG is still not a positivity-preserving method. The short pulse has a small
  negative nodal tail; longer, sharper-front and convergence tests remain needed.
- Initial pulse inventory is239171.65, NOT zero: the prescribed inlet nodes are
  already C1 in the initial P1 field although initial bulk C is configured0.
  Budget comparisons must account for this discrete initial inventory.
- After the pulse, inlet C0 is still a strong Dirichlet condition. Audit outflow
  includes dispersive/algebraic transfer at that boundary, not solely the far
  downstream advective outlet. Small relative budget error alone cannot prove
  a physically correct release mass or retention.

## Reproducing analysis without rerunning the model

Read-only analyzer: tests/analyze_adenc_hydro_run.py. It reads EVERY scalar
NodeData time block in each solute file (multiple times share one Gmsh file),
uses only participating nodes, and independently integrates

    M(t) = sum_e area_e * H_e(t) * (C1+C2+C3)/3.

Geometry/next-day data are in discovery-next/{nodes,active,next}.txt, produced
by tests/fortran/export_adenc_hydro_domain.f90. The standalone analyzer assumes
the recorded Rhine velocity constants vref1.7m/s,Qref1672m3/s and first-day
linear H; do not reuse blindly for different constants or longer runs.

    python3 analyze_adenc_hydro_run.py constant-state discovery-next --constant 1
    python3 analyze_adenc_hydro_run.py short-pulse discovery-next

Run the commands from the server parent above. Complete JSON summaries are
constant-summary.json and pulse-summary.json there, and locally in
runs/hydroflow-test-20261008-GGI9Sq/. All original outputs/logs/configurations,
including the guarded failed attempt, remain on the server. No drutes process
remained after completion. No numerical implementation changed during these
tests; the repository's default hydroflow.conf remains n.
