#!/usr/bin/env python3
"""Test consolidated input parsing, obsolete-file rejection and migration."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output-parent', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = args.build.resolve()
    args.output_parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='config-checks-', dir=args.output_parent.resolve()))
    objects = sorted(p for p in (build / 'objs').glob('*.o') if p.name != 'main.o')
    nf_flags = shlex.split(subprocess.check_output(['nf-config', '--fflags'], text=True))
    libraries = shlex.split(subprocess.check_output(['nf-config', '--flibs'], text=True))
    libraries.insert(0, '-L' + subprocess.check_output(['nc-config', '--libdir'], text=True).strip())
    binary = work / 'test-config'
    subprocess.run(['gfortran', '-fcoarray=single', '-fdefault-real-8', '-fcheck=all',
                    '-I', str(build / 'mods'), *nf_flags,
                    str(root / 'tests/fortran/test_adenc_config.f90'),
                    *map(str, objects), *libraries, '-o', str(binary)], check=True, cwd=work)
    results = []

    def check(label, transport=None, stabilization=None, expected=None, error=None, old=None):
        case = work / label
        inputs = case / 'drutes.conf/netcdf'
        inputs.mkdir(parents=True)
        (case / 'out').mkdir()
        if transport is not None:
            (inputs / 'transport.conf').write_text(transport)
        if stabilization is not None:
            (inputs / 'stabilization.conf').write_text(stabilization)
        for name, text in (old or {}).items():
            (inputs / name).write_text(text)
        result = subprocess.run([str(binary)], cwd=case, text=True, capture_output=True, timeout=30)
        output = result.stdout + result.stderr
        (case / 'terminal.log').write_text(output)
        if error:
            assert result.returncode != 0 and error in output, (label, output)
        else:
            assert result.returncode == 0, (label, output)
            values = next(line.split()[1:] for line in output.splitlines() if line.startswith('CONFIG'))
            actual = tuple(map(int, values[:6])) + tuple(map(float, values[6:]))
            assert actual == expected, (label, actual, expected)
        results.append({'check': label, 'passed': True})
        print(f'PASS {label}', flush=True)
        return case

    check('absent', expected=(0,0,0,0,0,0,1,1))
    check('disabled', 'n\nn\nn\nn\n', 'n\n2\nn\n1\n', expected=(0,0,0,0,0,0,2,1))
    for mode, settings, flags in [('galerkin','n\n2\nn\n1\n',(0,0)),
                                   ('supg1','y\n1\nn\n1\n',(1,0)),
                                   ('supg2','y\n2\nn\n1\n',(1,0)),
                                   ('supg2-shock1','y\n2\ny\n1\n',(1,1)),
                                   ('shock-only','n\n2\ny\n1\n',(0,1))]:
        factor = 1 if mode == 'supg1' else 2
        check(mode, 'y\ny\ny\nn\n', settings, expected=(1,1,1,0,*flags,factor,1))
    check('audit-only','y\nn\ny\nn\n', expected=(1,0,1,0,0,0,1,1))
    check('comments', '# banks\ny # enabled\n\ny\ny\nn\n# end\n',
          '# SUPG\ny\n2\n# residual viscosity\ny\n1\n', expected=(1,1,1,0,1,1,2,1))
    hydro = '0\n1\n1e-12\n1000\n100\n4\n1 3 -1 0\n3 5 -1 0\n2 4 1 0\n4 6 1 0\ny\n'
    check('hydraulic-header', 'y\ny\ny\ny\n' + hydro, 'y\n2\ny\n1\n', expected=(1,1,1,1,1,1,2,1))
    for name in ('riverbank.conf','conservative.conf','hydroflow.conf','supg.conf','shock.conf'):
        check('obsolete-' + name, old={name:'n\n'}, error='no longer supported')
        check('mixed-' + name, 'n\nn\nn\nn\n', 'n\n1\nn\n1\n',
              old={name:'n\n'}, error='no longer supported')
    for label, transport, stabilization, error in [
        ('transport-short','y\ny\ny\n',None,'Missing y/n'),
        ('transport-invalid','yes\nn\nn\nn\n',None,'single y or n'),
        ('transport-extra','n\nn\nn\nn\n99\n',None,'Unexpected trailing'),
        ('stabilization-short',None,'y\n2\n','Missing y/n'),
        ('stabilization-invalid',None,'y n\n2\nn\n1\n','single y or n'),
        ('stabilization-extra',None,'n\n2\nn\n1\n99\n','Unexpected trailing'),
        ('supg-negative',None,'y\n-2\nn\n1\n','SUPG factor must be'),
        ('shock-nan',None,'n\n2\ny\nNaN\n','Shock factor must be'),
        ('hydro-without-conservative','y\nn\nn\ny\n'+hydro,None,'requires conservative'),
    ]:
        check(label,transport,stabilization,error=error)
    case = work / 'migration'
    inputs = case / 'drutes.conf/netcdf'
    inputs.mkdir(parents=True)
    (case / 'out').mkdir()
    originals = {'riverbank.conf':'y\n','conservative.conf':'y\ny\n',
                 'hydroflow.conf':'y\n'+hydro,'supg.conf':'y\n2\n','shock.conf':'y\n1\n'}
    for name, text in originals.items():
        (inputs / name).write_text(text)
    subprocess.run([sys.executable,str(root / 'scripts/migrate_adenc_configs.py'),
                    str(case / 'drutes.conf'),'--backup-parent',str(work / 'backups')],check=True)
    assert all(not (inputs / name).exists() for name in originals)
    import tarfile
    archive = next((work / 'backups').glob('*/old-configs.tar.gz'))
    with tarfile.open(archive) as backup:
        saved = {Path(member.name).name:backup.extractfile(member).read().decode() for member in backup.getmembers()}
    assert saved == originals
    result = subprocess.run([str(binary)],cwd=case,text=True,capture_output=True,check=True)
    assert 'CONFIG 1 1 1 1 1 1' in result.stdout
    results.append({'check':'migration-values-and-backup','passed':True})
    (work / 'RESULTS.json').write_text(json.dumps(results,indent=2)+'\n')
    print(f'All {len(results)} consolidated configuration checks passed. Logs: {work}')


if __name__ == '__main__':
    main()
