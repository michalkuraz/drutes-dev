# ADEnc server area audit — 2026-10-08

Read-only audit of ongoing run
`/mnt/stock/ncflux-three-variants-20261008/runs/20261008T122213Z-Y6MRBb`.
No model source, configuration, or simulation was changed. Script executed on
server through stdin, without installing files there.

Script: `scripts/integrate_adenc_concentration.py`. Standard-library Python;
Gmsh 2.2 scalar NodeData and P1 triangular connectivity. Reads only complete
snapshots; uses actual node IDs, not array position assumptions.

For each triangle, the signed integral is exactly A*(C1+C2+C3)/3.
Positive and negative integrals are computed by clipping the affine field at
C=0, not by clipping nodal values and interpolating. Four analytical triangle
tests passed, and signed=positive-negative is checked on every snapshot.

## Inputs and area integral

Inlet 101 uses a piecewise-constant concentration: C=1 for 0 <= t < 21600 s,
then C=0. All three cases have the same six-hour pulse. Initial interior C=0;
Dirichlet nodes already have C=1 at t=0, so the represented initial FE field
has a nonzero area integral of 4.85642838e6.

Total mesh area: 7.64204659015e10 m2. This is the entire FE domain, not only
active-channel area. C is integrated as represented in all triangles.

| Time [days] | Galerkin integral C dA | SUPG2 integral C dA | SUPG2+shock1 integral C dA |
| --- | ---: | ---: | ---: |
| 0 | 4.85642838e6 | 4.85642838e6 | 4.85642838e6 |
| 0.25 | 2.27991090e8 | 2.20910732e8 | 2.06020548e8 |
| 0.5 | 1.04422114e8 | 9.30658798e7 | 6.61820487e7 |
| 1 | 2.70911717e7 | 2.63707244e7 | 2.53619453e7 |
| 2 | 1.94693788e7 | 2.04513300e7 | 1.77341446e7 |
| 3 | 2.36016378e7 | 2.32993987e7 | 1.55675342e7 |
| 4 | 2.30974052e7 | 2.33828643e7 | 1.19514519e7 |

At day 3, integral of negative magnitude / integral of positive concentration
is 0.176462 (Galerkin), 0.0411703 (SUPG2), and 0.000189429 (SUPG2+shock1).
All complete concentration snapshots checked are finite.

## What this does NOT establish

Integral C dA has units [mass/length] when C has [mass/volume]. It is not mass.
The physical storage proxy in current code is H=Q/(W*v), so an appropriate
inventory is M=integral H*C dA. H varies spatially and is updated daily.
Do not interpret these area-integral decreases as percentages of lost mass.

`ADEls_mass` actually returns C*|q|*element_area. Consequently the nodal
`conc_in_river` export is not the storage inventory and must not be integrated
again to obtain mass. The `conc_flux` callback exports the hydrological
depth-integrated flux, not the complete solute flux C*q-K*grad(C).

Current standard assembly uses H*dC/dt + q dot grad(C) - div(K grad(C))=0.
ADEnc does not replace the default zero `der_convect` callback. A conservative
storage equation instead reads d(H*C)/dt + div(q*C-K grad(C))=0.
These coincide if H_t+div(q)=0, with compatible interface/boundary treatment;
that continuity identity is not enforced by the current reconstruction of H
from Q, active width, and empirical velocity. This is a reason to audit, not
a measured numerical mass deficit.

For the advective formulation, the formal inventory balance contains the
additional volume contribution integral C*(H_t+div(q)) dA. Piecewise fields
also require accounting for interface jumps. SUPG/shock discrete contributions
and Dirichlet reactions must be included in a discrete balance.

To verify conservation, obtain storage at the assembly quadrature points (or
the actual lumped capacity weights) at each accepted step, the FE Dirichlet
reaction/boundary transport, source/reaction terms, and daily storage changes.
Compare M(t)-M(0) with integrated inlet minus outlet solute flux. Both upstream
and downstream boundaries are Dirichlet; after the pulse, upstream C=0 can
also allow diffusive loss. A zero downstream boundary value does not imply
zero diffusive outlet flux. Existing sparse snapshots alone do not close this
budget. No claim of verified mass conservation is warranted.

## Interpretation

A maximum concentration near 0.06 at day 3 is not by itself inconsistent:
dispersion, variable channel geometry/storage, and numerical diffusion can
lower a pulse peak. But the large early area-integral decline prevents calling
the result physically confirmed. SUPG2+shock1 has much smaller negative
contributions, yet its area integral continues decreasing between days 3 and
4 whereas SUPG2 is nearly constant. Investigate storage and boundary budgets
before attributing this difference solely to harmless pulse spreading.
