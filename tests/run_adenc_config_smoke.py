#!/usr/bin/env python3
"""Compare short real Rhine/Moselle runs using old and consolidated inputs.

Every invocation creates fresh run directories and symlinks large NC inputs.
The reference executable must be built from the pre-consolidation source.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile

import numpy as np

from analyze_adenc_hydro_run import nodal_blocks


def fortran_float(value):
    text = value.decode() if isinstance(value,bytes) else str(value)
    # Older Fortran exports omit E when the exponent needs three digits.
    text = re.sub(r'(?<=[0-9])([+-][0-9]{3})$',r'E\1',text.strip())
    return float(text)


def read_budget(path):
    return np.genfromtxt(path,delimiter=',',names=True,
                         converters={i:fortran_float for i in range(12)})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--baseline-binary', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--legacy-backup', type=Path, required=True)
    parser.add_argument('--output-parent', type=Path, required=True)
    args = parser.parse_args()
    args.output_parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='rhine-smoke-', dir=args.output_parent.resolve()))
    original = {}
    with tarfile.open(args.legacy_backup) as archive:
        # Migration archives prefix each configuration with its numeric index.
        for member in archive.getmembers():
            if Path(member.name).parts[0].endswith('-'+args.inputs.name):
                original[Path(member.name).name] = archive.extractfile(member).read().decode()
    assert set(original) == {'riverbank.conf','conservative.conf','hydroflow.conf','supg.conf','shock.conf'}
    results = []
    for variant, supg, shock in [('galerkin','n','n'),('supg2','y','n'),('supg2-shock1','y','y')]:
        runs = []
        for mode, binary in [('baseline',args.baseline_binary.resolve()),('combined',args.binary.resolve())]:
            case = work / f'{variant}-{mode}'
            case.mkdir()
            shutil.copytree(args.inputs, case/'drutes.conf', symlinks=True, ignore=shutil.ignore_patterns('*.nc'))
            inputs = case/'drutes.conf/netcdf'
            for path in (args.inputs/'netcdf').glob('*.nc'):
                assert path.resolve().is_file(), path
                (inputs/path.name).symlink_to(path.resolve())
            configuration = case/'drutes.conf/global.conf'
            text = configuration.read_text()
            text, count = re.subn(r'(#-----end time \[real\]\s*\n)[^\n]+',r'\g<1>3600.0',text,count=1)
            assert count == 1
            text, count = re.subn(r'(#number of observation times\s*\n).*?(?=###-----end observation time settings)',
                                 r'\g<1>3\n\n#observation time values [real]\n1200.0\n2400.0\n3600.0\n\n',text,count=1,flags=re.S)
            assert count == 1
            configuration.write_text(text)
            if mode == 'baseline':
                (inputs/'transport.conf').unlink()
                (inputs/'stabilization.conf').unlink()
                for name, content in original.items():
                    (inputs/name).write_text(content)
                (inputs/'supg.conf').write_text(supg+'\n2.0\n')
                (inputs/'shock.conf').write_text(shock+'\n1.0\n')
            else:
                (inputs/'stabilization.conf').write_text(supg+'\n2.0\n'+shock+'\n1.0\n')
            (case/'out').mkdir()
            with (case/'terminal.log').open('w') as log:
                process = subprocess.run([str(binary)],cwd=case,stdout=log,stderr=subprocess.STDOUT,timeout=240)
            (case/'exit.status').write_text(str(process.returncode)+'\n')
            budget = read_budget(case/'out/adenc_mass_balance.csv')
            assert budget.size and np.all(np.isfinite(np.array(budget.tolist())))
            assert abs(budget['time_s'][-1]-3600) < 1e-8
            assert 'F I N I S H E D' in (case/'terminal.log').read_text()
            if process.returncode:
                log = (case/'terminal.log').read_text()
                assert 'get_RAM_use' in log or ('postpro.f90' in log and 'RAM consumption' in log), case
            runs.append(case)
        baseline, combined = runs
        checked = 0
        maximum_difference = 0.0
        for path in sorted((baseline/'out').glob('obspt*.out')):
            old = np.loadtxt(path,comments='#',ndmin=2,converters=fortran_float)
            new = np.loadtxt(combined/'out'/path.name,comments='#',ndmin=2,converters=fortran_float)
            assert old.shape == new.shape and np.all(np.isfinite(new)), path
            np.testing.assert_allclose(new,old,rtol=1e-10,atol=1e-12,err_msg=str(path))
            maximum_difference = max(maximum_difference,float(np.max(np.abs(new-old))))
            checked += 1
        assert checked == 7
        count = 0
        for path in sorted((baseline/'out').glob('ADE_in_watershed_solute_concentration-*.msh')):
            if '-el_avg-' in path.name:
                continue
            old_blocks = list(nodal_blocks(path,number=fortran_float))
            new_blocks = list(nodal_blocks(combined/'out'/path.name,number=fortran_float))
            assert len(old_blocks) == len(new_blocks)
            for (old_time,_,old_values),(new_time,_,new_values) in zip(old_blocks,new_blocks):
                assert old_time == new_time and old_values.keys() == new_values.keys()
                old = np.array(list(old_values.values()))
                new = np.array([new_values[node] for node in old_values])
                assert np.all(np.isfinite(new))
                np.testing.assert_allclose(new,old,rtol=1e-10,atol=1e-12)
                maximum_difference = max(maximum_difference,float(np.max(np.abs(new-old))))
                count += 1
        assert count >= 4
        old = read_budget(baseline/'out/adenc_mass_balance.csv')
        new = read_budget(combined/'out/adenc_mass_balance.csv')
        assert old.dtype.names == new.dtype.names and old.shape == new.shape
        for field in old.dtype.names:
            np.testing.assert_allclose(new[field],old[field],rtol=1e-10,atol=1e-9,err_msg=field)
        results.append({'variant':variant,'duration_s':3600,'observation_histories':checked,
                        'spatial_snapshots':count,'maximum_field_difference':maximum_difference,
                        'max_relative_budget_error':float(np.max(np.abs(new['relative_error']))),
                        'baseline_exit':int((baseline/'exit.status').read_text()),
                        'combined_exit':int((combined/'exit.status').read_text())})
        print(f'PASS {variant}: seven histories/{count} maps/budget agree; field difference {maximum_difference:g}',flush=True)
    (work/'RESULTS.json').write_text(json.dumps(results,indent=2)+'\n')
    print(f'All three short real-data comparisons passed. Results: {work}')


if __name__ == '__main__':
    main()
