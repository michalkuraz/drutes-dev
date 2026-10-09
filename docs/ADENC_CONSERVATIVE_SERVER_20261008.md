# Conservative ADEnc server launch — 2026-10-08

User requested server launch and explicitly chose to stop the earlier three
bank-only simulations, preserving all results. No old output directory was
reused and no old result file was removed. The old verified launcher2606342
was sent TERM; its child models2606471/2606473/2606475 exited. Their statuses
are INTERRUPTED under:

    /mnt/stock/ncflux-three-variants-20261008/runs/20261008T151502Z-banks-69WBpi

## New batch

Server: miguel@hydrocalc.science.fzp.czu.cz

    /mnt/stock/ncflux-conservative-20261008/runs/20261008T164531Z-conservative-V9Flrm

Start approximately 2026-10-08 16:45 UTC / 18:45 Europe/Rome.
Launcher PID2608329; model PIDs2608442 (galerkin),2608444 (supg2),
2608446 (supg2-shock1). Each variant has its own drutes.conf, empty-at-launch
out/, terminal.log, pid.txt, status.txt, started-utc.txt. The private Linux
executable is in the batch's bin/ directory; its SHA256 is recorded alongside.
No Mac executable was transferred or used. All three have conservative and
inventory audit enabled, plus impermeable internal banks.

Same original first-scenario inputs: 2015-10-24, duration1209600 seconds,
6-hour unit-concentration pulse, Qmin300, alpha_L200m/alpha_T0.2m.
Stabilization variants: none; SUPG2; SUPG2+shock1. No outlet/width/forcing
physics was silently modified. Root Mac configuration remains n/n.

## Code and forcing provenance

Export root: /mnt/stock/ncflux-conservative-20261008
Launcher: bash run_conservative.sh
Latest batch path is recorded in latest-run.txt.

Local archive: exports/ncflux-conservative-20261008.tar.gz (approximately10MB).
Archive SHA256 verified before extraction:

    85eb7c13df74faaa71f3eaf1d53b448a9435b81a0f593d9e65fe30a29206a6b1

Private source includes base3e9d56318d71db521dcfe053487ccfe68c9af7ae plus
all current conservative working-tree changes and the two new modules.
SOURCE.md, source-working-diff.patch and SHA256SUMS identify the snapshot;
the commit alone does not identify this uncommitted implementation.
The server's existing dirty ~/drutes-dev checkout was NOT modified/pulled.

shared/ is a symlink to /mnt/stock/ncflux-three-variants-20261008/shared.
Every case links to those same existing data files. No large NetCDF file was
copied/transferred/edited. Hashes checked on the server:

- mRM_Fluxes_States.nc: 2d5635524727d16a91dc0195ee69202088cd7053f6c9c86329315ba01f5b6396
- dem.nc: 329937bf8809234e5041c5438eb57a2731af15f0d1d1142cc3bd48b8b2303dd8

Do not move/remove/modify those targets while runs are active.

## Verification before launch

Native Linux build included Schwarz. Production module bank regression and
conservative variable/jumping-storage, pulse/outflow, rejection and audit
regression passed. All three configurations passed full FEM initialization,
initial export and conservative matrix assembly without solving. Those
diagnostic outputs reside in separate preflight/ directories, not model out/.

The preparation-only batch20261008T164339Z-conservative-ZVsfM1 is retained:
its cases are READY and were NEVER launched. The running batch above is a
separate fresh launch; no duplicate simulations were started in preparation.

## Initial actual-run check (not a final result)

All three processes alive and RUNNING; first6/6/7 accepted steps reached
2314.683/2314.683/2846.1513 seconds. All recorded CSV values were finite.
Maximum absolute relative inventory errors:

- galerkin: 2.099e-15;
- supg2: 1.849e-15;
- supg2-shock1: 2.245e-15.

Each writes out/adenc_mass_balance.csv. This establishes startup and early
discrete budget closure, NOT conservation throughout14days, positivity,
correct Rhine hydraulics or calibration. Negative nodal contributions were
already nonzero; conservative assembly is not a positivity limiter.
Inherited713 zero-width elements and the missing explicit internal river
outlet remain; see ADENC_CONSERVATIVE.md. Full duration/daily transitions,
pulse transport, residuals and final budget still need assessment.

No new recurring monitoring automation was created by this launch.
