# NetCDF transport configuration — 2026-10-12

The model now uses two files under `drutes.conf/netcdf/`:
`transport.conf` and `stabilization.conf`. The former separate
`riverbank.conf`, `conservative.conf`, `hydroflow.conf`, `supg.conf` and
`shock.conf` are rejected, including when a combined file is also present.
This prevents an obsolete setting from being silently ignored.

## Transport and hydraulics

The first four data records of `transport.conf`, in order, are:

```text
# Zero total solute flux on banks
y
# Conservative transport
y
# Write accepted-step mass balance to out/adenc_mass_balance.csv
y
# Compatible hydraulic reconstruction
y
```

If the fourth flag is `n`, the file ends there. If it is `y`, append:

1. Lateral source policy: `0` = no distributed sources; `1` = explicit `lateral.conf`.
2. Hydraulic forcing: `0` = daily steps; `1` = linear daily interpolation.
3. Relative water-balance tolerance, positive and at most `1e-6`.
4. Maximum auxiliary projection PCG iterations, a positive integer.
5. Maximum relative edge-flux correction, positive.
6. Number of explicit port edges, at least two.
7. That many `node1 node2 kind discharge` records.
8. Optional final `y/n` to remove unported disconnected components at initialization.

Port nodes are FE array indices. `kind=-1` is an inlet: discharge `0` uses
Qrouted and positive values prescribe m³/s. `kind=1` is an outlet and requires
discharge `0`; its nonnegative flux is determined by the hydraulic projection.
Port lists must match the active mesh. The example file is a four-triangle
test configuration; its node indices are not Rhine boundary definitions.

Conservative/audit modes require no-flow banks and standard 2D Picard.
Hydraulic reconstruction also requires conservative mode. Conservative
time integration requires the supported transient implicit Euler formulation;
restart from backup and Schwarz execution remain unsupported in this mode.
The audit flag can be off when conservative hydraulic reconstruction is on.

Absent `transport.conf` means all four options are disabled, provided no
obsolete file is present. Root defaults preserve their previous physics:
banks on, conservative/audit/hydraulic reconstruction off. Current scenario
3a/3b input folders enable all four options.

## Stabilization

`stabilization.conf` has exactly four data records:

```text
# Enable SUPG
y
# SUPG multiplier
2.0
# Enable residual viscosity (shock capturing)
y
# Residual-viscosity multiplier
1.0
```

| Variant | SUPG flag | SUPG multiplier | Residual-viscosity flag | Multiplier |
|---|---|---|---|---|
| Galerkin | n | 2.0 | n | 1.0 |
| SUPG | y | 1.0 | n | 1.0 |
| SUPG2 | y | 2.0 | n | 1.0 |
| SUPG2 + residual viscosity | y | 2.0 | y | 1.0 |

Both multipliers must be finite and nonnegative, even when disabled.
Zero disables the corresponding contribution. The two switches are independent.
An absent file disables both stabilizations with multiplier defaults of one.
Neither option guarantees positivity; residual viscosity can increase spreading.
The physical dispersivity values remain in `netcdf.conf`.

Flags must be lowercase `y` or `n` on separate records. Blank lines and `#`
comments are allowed. Changing these files does not require recompilation.

## Migrating an existing configuration

```sh
python3 scripts/migrate_adenc_configs.py /path/to/drutes.conf --backup-parent /path/to/backups
```

The utility preserves the old parameter values and port list, writes the two
combined files, archives the five small old files before removing them, and
leaves large NetCDF files untouched. It refuses differing existing combined
files. Pass input configuration folders explicitly; do not alter the archived
inputs of completed scientific runs merely to inspect their outputs.

Every new model run needs its own empty `out/`: main clears `out/*` at startup.
Numerical calculations, boundary formulas, flux projection and stabilization
formulas were not changed by this input consolidation.
