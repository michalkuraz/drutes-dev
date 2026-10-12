#!/usr/bin/env python3
"""Run two fresh style-regression simulations against preserved reference fields."""
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import tempfile

import numpy as np

from run_adenc_config_smoke import fortran_float, read_budget
from analyze_adenc_hydro_run import nodal_blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output-parent', type=Path, required=True)
    parser.add_argument('--resume', type=Path)
    args = parser.parse_args()
    work = args.resume.resolve() if args.resume else Path(tempfile.mkdtemp(prefix='two-simulations-', dir=args.output_parent.resolve()))
    results = []

    for variant in ['galerkin', 'supg2-shock1']:
        reference = args.reference.resolve() / (variant + '-combined')
        case = work / variant
        completed = (case / 'terminal.log').exists() and 'F I N I S H E D' in (case / 'terminal.log').read_text()
        case.mkdir(exist_ok=completed)
        source = reference / 'drutes.conf'
        if not completed:
            shutil.copytree(source, case / 'drutes.conf', symlinks=True,
                            ignore=shutil.ignore_patterns('*.nc'))

        for path in (source / 'netcdf').glob('*.nc'):
            assert path.resolve().is_file(), path
            if not completed:
                (case / 'drutes.conf/netcdf' / path.name).symlink_to(path.resolve())

        if not completed:
            (case / 'out').mkdir()
            with (case / 'terminal.log').open('w') as log:
                process = subprocess.run([str(args.binary.resolve())], cwd=case,
                                         stdout=log, stderr=subprocess.STDOUT, timeout=240)
            (case / 'exit.status').write_text(str(process.returncode)+'\n')
            assert process.returncode == 0, (case, process.returncode)
        assert 'F I N I S H E D' in (case / 'terminal.log').read_text()
        old_budget = read_budget(reference / 'out/adenc_mass_balance.csv')
        budget = read_budget(case / 'out/adenc_mass_balance.csv')
        assert abs(budget['time_s'][-1] - 3600) < 1.e-8

        for name in budget.dtype.names:
            np.testing.assert_allclose(budget[name], old_budget[name], rtol=1.e-10, atol=1.e-12)

        histories = 0
        maximum = 0.0
        for path in (reference / 'out').glob('obspt*.out'):
            old = np.loadtxt(path, comments='#', ndmin=2, converters=fortran_float)
            new = np.loadtxt(case / 'out' / path.name, comments='#', ndmin=2, converters=fortran_float)
            assert np.all(np.isfinite(new)), path
            np.testing.assert_allclose(new, old, rtol=1.e-10, atol=1.e-12)
            maximum = max(maximum, float(np.max(np.abs(new - old))))
            histories += 1

        maps = 0
        for path in (reference / 'out').glob('ADE_in_watershed_solute_concentration-*.msh'):
            if '-el_avg-' in path.name:
                continue
            old = list(nodal_blocks(path, number=fortran_float))
            new = list(nodal_blocks(case / 'out' / path.name, number=fortran_float))
            assert len(old) == len(new)
            for (ot, oi, ov), (nt, ni, nv) in zip(old, new):
                assert ot == nt and oi == ni
                assert ov.keys() == nv.keys()
                old_values = np.array(list(ov.values()))
                new_values = np.array([nv[node] for node in ov])
                assert np.all(np.isfinite(new_values))
                np.testing.assert_allclose(new_values, old_values, rtol=1.e-10, atol=1.e-12)
                maximum = max(maximum, float(np.max(np.abs(new_values - old_values))))
                maps += 1

        assert histories == 7 and maps == 5
        results.append(dict(variant=variant, duration_s=3600, exit=0,
                            histories=histories, maps=maps, maximum_difference=maximum))
        print('PASS', variant, results[-1], flush=True)

    (work / 'RESULTS.json').write_text(json.dumps(results, indent=2) + '\n')
    print('Results:', work, flush=True)


if __name__ == '__main__':
    main()
