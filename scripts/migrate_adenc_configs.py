#!/usr/bin/env python3
"""Combine NetCDF transport settings, preserving values and backing up old files.

Pass configuration directories (e.g. drutes.conf, drutes.conf3b), not run roots.
Only the five named small input files are migrated; NetCDF data are untouched.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import tarfile
import tempfile

OLD = ("riverbank.conf", "conservative.conf", "hydroflow.conf", "supg.conf", "shock.conf", "hydroflow.conf.example")


def records(path: Path, default: list[str]) -> list[str]:
    if not path.exists():
        return default
    return [line.split("#", 1)[0].strip() for line in path.read_text().splitlines()
            if line.split("#", 1)[0].strip()]


def combine(directory: Path) -> dict[str, str]:
    banks = records(directory / "riverbank.conf", ["n"])
    conservative = records(directory / "conservative.conf", ["n", "n"])
    hydro = records(directory / "hydroflow.conf", ["n"])
    supg = records(directory / "supg.conf", ["n", "1.0"])
    shock = records(directory / "shock.conf", ["n", "1.0"])
    if len(banks) != 1 or len(conservative) != 2 or len(supg) != 2 or len(shock) != 2 or not hydro:
        raise ValueError(f"Unexpected configuration record counts in {directory}")
    if any(value not in ("y", "n") for value in [*banks, *conservative, hydro[0], supg[0], shock[0]]):
        raise ValueError(f"Invalid y/n setting in {directory}")
    if hydro[0] == "n" and len(hydro) != 1:
        raise ValueError(f"Disabled hydroflow has unused records in {directory}")
    transport = ("# Zero total solute flux on banks [y/n]\n" + banks[0] + "\n"
                 "# Conservative transport [y/n]\n" + conservative[0] + "\n"
                 "# Write accepted-step mass balance [y/n]\n" + conservative[1] + "\n"
                 "# Compatible hydraulic reconstruction [y/n]\n" + hydro[0] + "\n")
    if hydro[0] == "y":
        if len(hydro) < 7:
            raise ValueError(f"Incomplete hydraulic configuration in {directory}")
        count = int(hydro[6])
        if count < 2 or len(hydro) not in (7 + count, 8 + count):
            raise ValueError(f"Invalid hydraulic port count in {directory}")
        labels = ["Lateral source policy: 0=none, 1=explicit lateral.conf",
                  "Hydraulic forcing: 0=daily steps, 1=linear interpolation",
                  "Relative water-balance tolerance", "Maximum projection PCG iterations",
                  "Maximum relative edge-flux correction", "Number of explicit port edges"]
        for label, value in zip(labels, hydro[1:7]):
            transport += f"# {label}\n{value}\n"
        transport += "# node1 node2 kind discharge; kind -1=inlet, +1=outlet\n"
        transport += "\n".join(hydro[7:7 + count]) + "\n"
        if len(hydro) == 8 + count:
            transport += "# Remove unported disconnected components at initialization [y/n]\n" + hydro[-1] + "\n"
    stabilization = ("# Enable SUPG [y/n]\n" + supg[0] + "\n"
                     "# SUPG time-scale multiplier\n" + supg[1] + "\n"
                     "# Enable residual viscosity (shock capturing) [y/n]\n" + shock[0] + "\n"
                     "# Residual-viscosity multiplier\n" + shock[1] + "\n")
    outputs = {"transport.conf": transport, "stabilization.conf": stabilization}
    if (directory / "hydroflow.conf.example").exists():
        example = records(directory / "hydroflow.conf.example", [])
        outputs["transport.conf.example"] = (
            "# TOY four-triangle mesh; replace node indices with verified ports for your mesh.\n"
            "# No-flow banks / conservative transport / mass audit / hydraulic settings\n"
            "y\ny\ny\n" + "\n".join(example) + "\n")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="+", type=Path)
    parser.add_argument("--backup-parent", required=True, type=Path)
    args = parser.parse_args()
    plans = []
    for configuration in args.directories:
        directory = configuration.resolve() / "netcdf"
        if not directory.is_dir():
            raise ValueError(f"Missing NetCDF configuration directory: {directory}")
        old = [directory / name for name in OLD if (directory / name).exists()]
        if not old and all((directory / name).exists() for name in ("transport.conf", "stabilization.conf")):
            print(f"Already migrated: {directory}")
            continue
        outputs = combine(directory)
        for name, content in outputs.items():
            if (directory / name).exists() and (directory / name).read_text() != content:
                raise ValueError(f"Refusing to overwrite differing {directory / name}")
        plans.append((directory, old, outputs))
    args.backup_parent.mkdir(parents=True, exist_ok=True)
    backup = Path(tempfile.mkdtemp(prefix="config-backup-", dir=args.backup_parent.resolve()))
    with tarfile.open(backup / "old-configs.tar.gz", "w:gz") as archive:
        for index, (directory, old, _) in enumerate(plans):
            for path in old:
                archive.add(path, arcname=f"{index}-{directory.parent.name}/netcdf/{path.name}", recursive=False)
    for directory, old, outputs in plans:
        for name, content in outputs.items():
            path = directory / name
            if not path.exists():
                path.write_text(content)
        for path in old:
            path.unlink()
        print(f"Migrated: {directory}")
    print(f"Original small inputs preserved: {backup / 'old-configs.tar.gz'}")


if __name__ == "__main__":
    main()
