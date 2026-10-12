#!/usr/bin/env python3
"""Create fresh examples from saved inputs, with shared NC data and three stabilizations."""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--case', action='append', required=True, help='label=/path/to/input/drutes.conf')
    args = parser.parse_args()
    destination = args.destination.resolve()
    if destination.exists():
        raise SystemExit(f'Refusing to overwrite existing examples: {destination}')
    specifications = []
    for specification in args.case:
        label, folder = specification.split('=',1)
        if not label or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789-' for c in label):
            raise ValueError(f'Invalid example label: {label}')
        source = Path(folder).resolve()
        assert (source/'netcdf/netcdf.conf').is_file(), source
        for nc in (source/'netcdf').glob('*.nc'):
            assert nc.resolve().is_file(), nc
        specifications.append((label,source))
    destination.mkdir(parents=True)
    directories = []
    for label, source in specifications:
        for variant, supg, shock in [('galerkin','n','n'),('supg2','y','n'),('supg2-shock1','y','y')]:
            configuration = destination/f'{label}-{variant}'/'drutes.conf'
            shutil.copytree(source,configuration,symlinks=True,ignore=shutil.ignore_patterns('*.nc'))
            for nc in (source/'netcdf').glob('*.nc'):
                (configuration/'netcdf'/nc.name).symlink_to(nc.resolve())
            old_stabilization = configuration/'netcdf/stabilization.conf'
            if old_stabilization.exists():
                old_stabilization.write_text(supg+'\n2.0\n'+shock+'\n1.0\n')
            else:
                (configuration/'netcdf/supg.conf').write_text(supg+'\n2.0\n')
                (configuration/'netcdf/shock.conf').write_text(shock+'\n1.0\n')
            directories.append(configuration)
    subprocess.run([sys.executable,str(Path(__file__).with_name('migrate_adenc_configs.py')),
                    *map(str,directories),'--backup-parent',str(destination/'migration-backups')],check=True)
    (destination/'README.md').write_text(
        '# Consolidated NetCDF transport examples\n\n'
        'Each folder contains a fresh drutes.conf and three stabilization variants are supplied per case.\n'
        'Large NetCDF inputs are symlinks to the existing data. No simulations have been launched here.\n'
        'Run from a fresh working folder with this configuration and an empty out directory.\n'
        'The archived original simulations and their inputs remain the provenance of published results.\n')
    print(f'Prepared {len(directories)} examples in {destination}')


if __name__ == '__main__':
    main()
