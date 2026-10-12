#!/usr/bin/env python3
"""Run real initialization/export/assembly checks in fresh directories for each example."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--examples',type=Path,required=True)
    parser.add_argument('--output-parent',type=Path,required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.output_parent.mkdir(parents=True,exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='example-checks-',dir=args.output_parent.resolve()))
    build = args.build.resolve()
    objects = sorted(p for p in (build/'objs').glob('*.o') if p.name != 'main.o')
    nf_flags = shlex.split(subprocess.check_output(['nf-config','--fflags'],text=True))
    libraries = shlex.split(subprocess.check_output(['nf-config','--flibs'],text=True))
    libraries.insert(0,'-L'+subprocess.check_output(['nc-config','--libdir'],text=True).strip())
    binary = work/'preflight'
    subprocess.run(['gfortran','-fcoarray=single','-fdefault-real-8','-fcheck=all',
                    '-I',str(build/'mods'),*nf_flags,
                    str(root/'tests/fortran/check_adenc_bank_initialization.f90'),
                    *map(str,objects),*libraries,'-o',str(binary)],cwd=work,check=True)
    configurations = sorted(args.examples.glob('*/drutes.conf'))
    assert configurations, args.examples
    results = []
    for source in configurations:
        case = work/source.parent.name
        case.mkdir()
        shutil.copytree(source,case/'drutes.conf',symlinks=True,ignore=shutil.ignore_patterns('*.nc'))
        for nc in (source/'netcdf').glob('*.nc'):
            (case/'drutes.conf/netcdf'/nc.name).symlink_to(nc.resolve())
        (case/'out').mkdir()
        with (case/'terminal.log').open('w') as log:
            result = subprocess.run([str(binary),'--assemble'],cwd=case,stdout=log,stderr=subprocess.STDOUT,timeout=120)
        text = (case/'terminal.log').read_text()
        assert result.returncode == 0 and 'Conservative Rhine matrix assembly passed' in text, (case,text[-4000:])
        results.append({'example':source.parent.name,'passed':True,'log':str(case/'terminal.log')})
        print(f'PASS {source.parent.name}',flush=True)
    (work/'RESULTS.json').write_text(json.dumps(results,indent=2)+'\n')
    print(f'All {len(results)} example initialization/export/assembly checks passed. Logs: {work}')


if __name__ == '__main__':
    main()
