import glob
import logging
import os
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_parallel_shock_coarse_uniform.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    for subdir in ['bin', 'pvtk', 'rst']:
        pattern = os.path.join(_athena_exe_dir(), subdir, basename + '.*')
        for fname in glob.glob(pattern):
            os.remove(fname)


def _latest_output_file(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError('No output files found for pattern: ' + pattern)
    return matches[-1]


def _latest_restart_file(basename):
    pattern = os.path.join(_athena_exe_dir(), 'rst', basename + '.*.rst')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError('No restart files found for pattern: ' + pattern)
    return os.path.relpath(matches[-1], _athena_exe_dir())


def _read_restart_integer(restart_file, block, name):
    full_path = os.path.join(_athena_exe_dir(), restart_file)
    with open(full_path, 'rb') as handle:
        text = handle.read(65536).decode('latin1', errors='ignore')
    if '<par_end>' not in text:
        raise RuntimeError('Restart parameter header is missing <par_end>: ' +
                           full_path)
    active_block = None
    for raw_line in text.split('<par_end>', 1)[0].splitlines():
        line = raw_line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('<') and line.endswith('>'):
            active_block = line[1:-1]
            continue
        if active_block != block or '=' not in line:
            continue
        key, value = line.split('=', 1)
        if key.strip() != name:
            continue
        return int(value.split('#', 1)[0].strip())
    raise RuntimeError(f'Missing restart parameter <{block}>/{name} in ' +
                       full_path)


def _read_particle_density(basename):
    data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_d'))
    return np.asarray(data['pdens'])


def _measure_particle_count(basename):
    return float(np.sum(_read_particle_density(basename)))


def _run_command(label, nproc, arguments, restart_file=None):
    command = ['./athena']
    if restart_file is None:
        command += ['-i', _athena_input_path()]
    else:
        command += ['-r', restart_file]
    command += list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = proc.stdout + proc.stderr
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + ': ' +
                           ' '.join(command) + '\n' + output)
    return output


def _recenter_args(basename):
    return [
        'job/basename=' + basename,
        'mesh/nx1=32',
        'mesh/nx2=16',
        'mesh/nx3=1',
        'meshblock/nx1=8',
        'meshblock/nx2=8',
        'meshblock/nx3=1',
        'particles/ppc=1.0',
        'problem/ps_enable_injection=false',
        'problem/ps_enable_frame_tracking=true',
        'problem/ps_frame_mode=recenter',
        'problem/ps_recenter_x_target=0.001',
        'problem/ps_recenter_x_trigger=0.002',
        'problem/ps_recenter_shift_cells=1',
        'time/nlim=2',
        'output1/file_type=bin',
        'output1/variable=prtcl_d',
        'output1/id=prtcl_d',
        'output1/dcycle=1',
        'output1/ghost_zones=false',
        'output2/dcycle=0',
        'output3/dcycle=0',
        'output4/dcycle=0',
    ]


def _velocity_diag_args(basename):
    return [
        'job/basename=' + basename,
        'mesh/nx1=32',
        'mesh/nx2=16',
        'mesh/nx3=1',
        'meshblock/nx1=8',
        'meshblock/nx2=8',
        'meshblock/nx3=1',
        'particles/ppc=0.0',
        'particles/deposit_moments=true',
        'problem/ps_enable_injection=false',
        'problem/ps_enable_frame_tracking=true',
        'problem/ps_frame_mode=velocity',
        'problem/ps_frame_t_start=0.0',
        'problem/ps_frame_t_ramp=0.0',
        'problem/ps_frame_vfrac=0.0',
        'problem/ps_feedback_diag_dcycle=1',
        'time/nlim=1',
        'output1/dcycle=0',
        'output2/dcycle=0',
        'output3/dcycle=0',
        'output4/dcycle=0',
    ]


def _injection_restart_args(basename, nlim, rst_dcycle):
    return [
        'job/basename=' + basename,
        'mesh/nx1=32',
        'mesh/nx2=16',
        'mesh/nx3=1',
        'meshblock/nx1=8',
        'meshblock/nx2=8',
        'meshblock/nx3=1',
        'particles/ppc=0.0',
        'particles/deposit_moments=true',
        'particles/deposit_qscale=1.0e-3',
        'particles/couple_moments_to_mhd=true',
        'particles/couple_j_to_efield_representation=cell_centered',
        'particles/couple_moments_momentum_to_mhd=false',
        'particles/couple_moments_energy_to_mhd=false',
        'particles/pic_feedback_mode=coupled',
        'problem/ps_eta=0.1',
        'problem/ps_enable_injection=true',
        'problem/ps_enable_gas_subtraction=false',
        'problem/ps_enable_frame_tracking=false',
        'time/nlim=' + str(nlim),
        'time/tlim=1.0',
        'output1/file_type=bin',
        'output1/variable=prtcl_d',
        'output1/id=prtcl_d',
        'output1/dcycle=1',
        'output1/ghost_zones=false',
        'output2/dcycle=0',
        'output3/dcycle=0',
        'output4/dcycle=' + str(rst_dcycle),
    ]


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    cases = [
        {
            'name': 'serial_recenter',
            'basename': 'pic_parallel_recenter_serial',
            'nproc': 1,
        },
        {
            'name': 'mpi2_recenter',
            'basename': 'pic_parallel_recenter_mpi2',
            'nproc': 2,
        },
    ]

    for case in cases:
        _remove_outputs(case['basename'])
        output = _run_command(case['name'], case['nproc'],
                              _recenter_args(case['basename']))
        _RESULTS[case['name']] = output

    basename = 'pic_parallel_velocity_diag_zero_dv'
    _remove_outputs(basename)
    output = _run_command('velocity_diag_zero_dv', 1,
                          _velocity_diag_args(basename))
    _RESULTS['velocity_diag_zero_dv'] = {
        'kind': 'velocity_diag_zero_dv',
        'output': output,
    }

    restart_cases = [
        ('serial_injection_restart', 1),
        ('mpi2_injection_restart', 2),
    ]
    for case_name, nproc in restart_cases:
        full = 'pic_parallel_inj_' + case_name + '_full'
        seg = 'pic_parallel_inj_' + case_name + '_seg'
        rst = 'pic_parallel_inj_' + case_name + '_rst'
        for basename in [full, seg, rst]:
            _remove_outputs(basename)

        _run_command(case_name + '_full_run', nproc,
                     _injection_restart_args(full, 4, 0))
        full_density = _read_particle_density(full)
        full_count = float(np.sum(full_density))

        _run_command(case_name + '_seg_run', nproc,
                     _injection_restart_args(seg, 2, 2))
        seg_count = _measure_particle_count(seg)
        rst_path = _latest_restart_file(seg)
        full_rst_path = os.path.join(_athena_exe_dir(), rst_path)
        if not os.path.exists(full_rst_path):
            raise RuntimeError('Expected restart file not found: ' + full_rst_path)
        restart_next_tag = _read_restart_integer(rst_path, 'problem',
                                                 'ps_next_tag')

        _run_command(case_name + '_restart_run', nproc,
                     _injection_restart_args(rst, 4, 0),
                     restart_file=rst_path)
        restart_density = _read_particle_density(rst)
        restart_count = float(np.sum(restart_density))
        density_err = float(np.max(np.abs(full_density - restart_density)))
        _RESULTS[case_name] = {
            'full_count': full_count,
            'restart_count': restart_count,
            'seg_count': seg_count,
            'restart_next_tag': restart_next_tag,
            'density_err': density_err,
        }


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for case_name, output in _RESULTS.items():
        if isinstance(output, dict):
            continue
        saw_recenter = 'pic_parallel_shock frame(recenter)' in output
        saw_event = 'dev=1' in output
        clean = 'recenter particle migration failed' not in output
        logger.info('%s saw_recenter=%s saw_event=%s clean=%s',
                    case_name, saw_recenter, saw_event, clean)
        ok = saw_recenter and saw_event and clean and ok

    for case_name, result in _RESULTS.items():
        if not isinstance(result, dict):
            continue
        if result.get('kind') == 'velocity_diag_zero_dv':
            saw_diag = 'pic_parallel_shock feedback_diag:' in result['output']
            saw_cfg = 'pic_parallel_shock feedback_diag_cfg:' in result['output']
            logger.info('%s saw_diag=%s saw_cfg=%s',
                        case_name, saw_diag, saw_cfg)
            ok = saw_diag and saw_cfg and ok
            continue
        full_count = result['full_count']
        restart_count = result['restart_count']
        seg_count = result['seg_count']
        restart_next_tag = result['restart_next_tag']
        density_err = result['density_err']
        abs_err = abs(full_count - restart_count)
        logger.info('%s full_count=% .8e restart_count=% .8e abs_err=% .8e '
                    'density_err=% .8e seg_count=% .8e restart_next_tag=%d',
                    case_name, full_count, restart_count, abs_err,
                    density_err, seg_count, restart_next_tag)
        min_next_tag = max(1, int(round(seg_count)))
        ok = (abs_err <= 1.0e-12 and density_err <= 1.0e-12 and
              restart_next_tag >= min_next_tag) and ok

    return ok
