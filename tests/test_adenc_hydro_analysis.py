import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    'hydro_analysis', Path(__file__).with_name('analyze_adenc_hydro_run.py'))
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def block(time, component=1):
    return (f'$NodeData\n1\n"C"\n1\n{time}\n3\n0\n{component}\n3\n'
            '1 0.25\n2 0.25\n3 0.25\n$EndNodeData\n')


def test_reads_multiple_times_in_one_gmsh_file(tmp_path):
    path = tmp_path/'scalar.msh'
    path.write_text(block(0)+block(300))
    records = list(analysis.nodal_blocks(path))
    assert [r[0] for r in records] == [0, 300]
    assert records[-1][2] == {1: .25, 2: .25, 3: .25}


def test_rejects_non_scalar_data(tmp_path):
    path = tmp_path/'scalar.msh'
    path.write_text(block(0, component=2))
    with pytest.raises(ValueError, match='Expected scalar'):
        list(analysis.nodal_blocks(path))


def test_independent_inventory_for_p0_depth_p1_concentration(tmp_path):
    geometry = tmp_path/'geometry'; geometry.mkdir()
    geometry.joinpath('nodes.txt').write_text('1 0 0 0\n2 1 0 0\n3 0 1 0\n')
    geometry.joinpath('active.txt').write_text('1 1 2 3 .333 .333 2 1672 1 0\n')
    geometry.joinpath('next.txt').write_text('1 1672\n')
    folder = tmp_path/'run'; folder.mkdir(); folder.joinpath('out').mkdir()
    folder.joinpath('exit.status').write_text('0\n')
    folder.joinpath('terminal.log').write_text('F I N I S H E D !!!!!\n')
    folder.joinpath('out/ADE_in_watershed_solute_concentration-0.msh').write_text(block(0)+block(300))
    expected = .5*(1672/(2*1.7))*.25
    folder.joinpath('out/adenc_mass_balance.csv').write_text(
        'time_s,inventory,initial_inventory,relative_error,free_residual\n'
        f'300,{expected},{expected},0,0\n')
    result = analysis.analyze(folder, geometry, constant=.25)
    assert result['numerical_finish']
    assert result['active_nodes'] == 3
    assert result['snapshots'][-1]['inventory'] == pytest.approx(expected)
    assert result['snapshots'][-1]['spatial_minus_audit_inventory'] == pytest.approx(0)
    assert result['snapshots'][-1]['max_constant_error'] == 0
