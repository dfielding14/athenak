import glob
import logging
import os
import re
import subprocess

import numpy as np
import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_relativistic_gyro_paper.athinput'
_RESULTS = {}
_C = 3.0
_B = np.array([0.0, 0.0, 1.0])
_U0 = np.array([1.0, 0.0, 0.0])
_CONTINUUM_C_VALUES = (1.5, 3.0, 10.0)
_CONTINUUM_THETA_VALUES = (0.8, 0.4, 0.2, 0.1)
_CONTINUUM_TEND = 1.2
_LINEAR_WAVE_TLIM_SCALE = 4.0


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    for pattern in [
            os.path.join(exe_dir, 'pvtk', basename + '.*.part.vtk'),
            os.path.join(exe_dir, 'rst', basename + '.*.rst*'),
            os.path.join(exe_dir, basename + '-errs.dat')]:
        for fname in glob.glob(pattern):
            os.remove(fname)


def _run_athena(label, arguments, restart_file=None, expect_fail=None):
    command = ['./athena']
    if restart_file is None:
        command += ['-i', _athena_input_path()]
    else:
        command += ['-r', restart_file]
    command += list(arguments)

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if expect_fail is None:
        if proc.returncode != 0:
            raise RuntimeError('Command failed for ' + label + '\n' + output)
        return
    if proc.returncode == 0:
        raise RuntimeError('Expected failure for ' + label + ', but command passed')
    if expect_fail not in output:
        raise RuntimeError('Unexpected failure reason for ' + label + '\n' + output)


def _latest_file(pattern, label):
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError('No ' + label + ' files found for pattern: ' + pattern)
    return matches[-1]


def _latest_pvtk_file(basename):
    return _latest_file(os.path.join(_athena_exe_dir(), 'pvtk',
                                     basename + '.prtcl_all.*.part.vtk'),
                        'particle VTK')


def _latest_restart_file(basename):
    return _latest_file(os.path.join(_athena_exe_dir(), 'rst',
                                     basename + '.*.rst'),
                        'restart')


def _read_big_endian_floats(contents, offset, count, label):
    payload = contents[offset:offset + 4*count]
    if len(payload) != 4*count:
        raise RuntimeError('Truncated ' + label + ' payload in particle VTK output')
    return np.frombuffer(payload, dtype='>f4').astype(np.float64)


def _find_marker(contents, pattern, offset, label):
    match = re.search(pattern, contents[offset:])
    if match is None:
        raise RuntimeError('Could not find ' + label + ' marker in particle VTK output')
    return offset + match.start(), offset + match.end(), match


def _read_pvtk_snapshot(basename):
    path = _latest_pvtk_file(basename)
    with open(path, 'rb') as fp:
        contents = fp.read()

    header = re.search(
        rb'# AthenaK particle data at time=\s*([^ ]+)\s+nranks=.*cycle=([0-9]+)',
        contents)
    if header is None:
        raise RuntimeError('Could not read particle VTK cycle header in ' + path)
    time = float(header.group(1))
    cycle = int(header.group(2))

    _, offset, match = _find_marker(
        contents, rb'\nPOINTS\s+([0-9]+)\s+float\n', 0, 'POINTS')
    npoint = int(match.group(1))
    offset += 4*3*npoint

    for name in ['gid', 'ptag', 'species']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' int\nLOOKUP_TABLE default\n')
        _, offset, _ = _find_marker(contents, pattern, offset, 'SCALARS ' + name)
        offset += 4*npoint

    scalars = {}
    for name in ['deltaf_f0', 'deltaf_weight']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' float\nLOOKUP_TABLE default\n')
        _, offset, _ = _find_marker(contents, pattern, offset, 'SCALARS ' + name)
        scalars[name] = _read_big_endian_floats(contents, offset, npoint, name)
        offset += 4*npoint

    _, offset, _ = _find_marker(contents, rb'\nVECTORS vel float\n',
                                offset, 'VECTORS vel')
    velocity = _read_big_endian_floats(contents, offset, 3*npoint, 'vel')
    return {
        'path': path,
        'time': time,
        'cycle': cycle,
        'npoint': npoint,
        'velocity': velocity.reshape(npoint, 3),
        'deltaf_f0': scalars['deltaf_f0'],
        'deltaf_weight': scalars['deltaf_weight'],
    }


def _boris_expected_velocity(cycle, time, light_speed=_C):
    if cycle <= 0:
        return _U0/np.sqrt(1.0 + np.dot(_U0, _U0)/(light_speed*light_speed))

    dt = time/cycle
    state = _U0.copy()
    for _ in range(cycle):
        gamma = np.sqrt(1.0 + np.dot(state, state)/(light_speed*light_speed))
        t = 0.5*dt*_B/gamma
        s = 2.0*t/(1.0 + np.dot(t, t))
        state_prime = state + np.cross(state, t)
        state += np.cross(state_prime, s)
    gamma = np.sqrt(1.0 + np.dot(state, state)/(light_speed*light_speed))
    return state/gamma


def _continuum_expected_velocity(time, light_speed):
    gamma = np.sqrt(1.0 + np.dot(_U0, _U0)/(light_speed*light_speed))
    phase = np.linalg.norm(_B)*time/gamma
    state = np.array([
        _U0[0]*np.cos(phase) + _U0[1]*np.sin(phase),
        _U0[1]*np.cos(phase) - _U0[0]*np.sin(phase),
        _U0[2],
    ])
    return state/gamma


def _momentum_from_velocity(velocity, light_speed):
    velocity2 = np.sum(velocity*velocity, axis=1)
    gamma = 1.0/np.sqrt(1.0 - velocity2/(light_speed*light_speed))
    return velocity*gamma[:, None]


def _kinetic_energy_from_momentum(momentum, light_speed):
    momentum2 = np.sum(momentum*momentum, axis=1)
    gamma = np.sqrt(1.0 + momentum2/(light_speed*light_speed))
    return (gamma - 1.0)*light_speed*light_speed


def _phase_error(measured, expected):
    measured_phase = np.arctan2(-measured[1], measured[0])
    expected_phase = np.arctan2(-expected[1], expected[0])
    delta = measured_phase - expected_phase
    return float(abs(np.arctan2(np.sin(delta), np.cos(delta))))


def _continuum_tag(light_speed, theta_max):
    return ('c' + str(light_speed).replace('.', 'p')
            + '_theta' + str(theta_max).replace('.', 'p'))


def _measure_continuum_case(basename, light_speed):
    snapshot = _read_pvtk_snapshot(basename)
    velocity = snapshot['velocity']
    momentum = _momentum_from_velocity(velocity, light_speed)
    pnorm = np.linalg.norm(momentum, axis=1)
    energy = _kinetic_energy_from_momentum(momentum, light_speed)
    expected_velocity = _continuum_expected_velocity(snapshot['time'], light_speed)
    pnorm0 = float(np.linalg.norm(_U0))
    energy0 = float((np.sqrt(1.0 + pnorm0*pnorm0/(light_speed*light_speed)) - 1.0)
                    * light_speed*light_speed)
    snapshot.update({
        'dt': snapshot['time']/snapshot['cycle'],
        'phase_error': _phase_error(np.mean(velocity, axis=0), expected_velocity),
        'pnorm_rel_error': float(np.max(np.abs(pnorm - pnorm0))/pnorm0),
        'energy_rel_error': float(np.max(np.abs(energy - energy0))/energy0),
        'velocity_spread': float(np.max(np.abs(velocity - velocity[0]))),
    })
    return snapshot


def _measure_case(basename):
    snapshot = _read_pvtk_snapshot(basename)
    snapshot['expected_velocity'] = _boris_expected_velocity(
        snapshot['cycle'], snapshot['time'])
    return snapshot


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    base_full = 'pic_rel_gyro_full'
    base_seg = 'pic_rel_gyro_seg'
    base_rst = 'pic_rel_gyro_rst'
    base_deltaf_box = 'pic_rel_gyro_deltaf_box'
    continuum_basenames = []
    for light_speed in _CONTINUUM_C_VALUES:
        for theta_max in _CONTINUUM_THETA_VALUES:
            tag = _continuum_tag(light_speed, theta_max)
            continuum_basenames.append('pic_rel_gyro_continuum_' + tag)
    basenames = [base_full, base_seg, base_rst, base_deltaf_box] + continuum_basenames
    for basename in basenames:
        _remove_outputs(basename)

    _run_athena('full', ['job/basename=' + base_full, 'time/nlim=2'])
    _RESULTS['full'] = _measure_case(base_full)

    _run_athena('segment', ['job/basename=' + base_seg, 'time/nlim=1'])
    _RESULTS['segment'] = _measure_case(base_seg)
    restart_file = os.path.relpath(_latest_restart_file(base_seg),
                                   _athena_exe_dir())

    _run_athena('restart',
                ['job/basename=' + base_rst, 'time/nlim=2',
                 'output1/file_number=0', 'output2/dcycle=0'],
                restart_file=restart_file)
    _RESULTS['restart'] = _measure_case(base_rst)

    _run_athena(
        'guard_restart_physical_model',
        ['particles/pic_physical_mode=engineering',
         'particles/pic_cr_light_speed=1.0',
         'particles/pic_cr_initial_state=velocity',
         'time/nlim=2', 'output1/dcycle=0', 'output2/dcycle=0'],
        restart_file=restart_file,
        expect_fail='Particle restart physical-model metadata mismatch',
    )
    _run_athena(
        'guard_velocity_initial_state_below_c',
        ['particles/pic_cr_initial_state=velocity',
         'particles/cr_vx0=3.1', 'time/nlim=0'],
        expect_fail='must define |v| < <particles>/pic_cr_light_speed',
    )
    _run_athena(
        'paper_deltaf_expanding_box',
        ['job/basename=' + base_deltaf_box,
         'particles/pic_deltaf_mode=physical',
         'particles/pic_deltaf_f0=kappa_aniso',
         'particles/pic_deltaf_p0=1.0',
         'particles/pic_deltaf_kappa=1.25',
         'particles/pic_deltaf_aniso_x1=0.75',
         'particles/pic_deltaf_aniso_x2=1.25',
         'particles/pic_expanding_box_mode=on',
         'particles/pic_expansion_law=linear',
         'particles/pic_expansion_rate_x1=0.1',
         'time/nlim=1', 'output2/dcycle=0'],
    )
    _RESULTS['deltaf_box'] = _measure_case(base_deltaf_box)
    _run_athena(
        'guard_restart_extension_model',
        ['particles/pic_deltaf_mode=physical',
         'particles/pic_deltaf_f0=kappa_iso',
         'time/nlim=2', 'output1/dcycle=0', 'output2/dcycle=0'],
        restart_file=restart_file,
        expect_fail='Particle restart physical-model metadata mismatch',
    )
    _run_athena(
        'guard_nonpositive_box_scale',
        ['particles/pic_expanding_box_mode=on',
         'particles/pic_expansion_rate_x1=-2.0', 'time/nlim=0'],
        expect_fail='must produce finite, positive scale factors',
    )
    _run_athena(
        'guard_paper_hall_requires_extension_mode',
        ['particles/pic_cr_hall_mode=current_to_ct_experimental',
         'time/nlim=0'],
        expect_fail='requires <particles>/pic_physical_mode=extended_mhd_pic',
    )
    _run_athena(
        'guard_negative_particle_load_cost',
        ['particles/pic_load_balance_cost_per_particle=-1.0',
         'time/nlim=0'],
        expect_fail='<particles>/pic_load_balance_cost_per_particle must be >= 0',
    )
    _RESULTS['continuum'] = {}
    for light_speed in _CONTINUUM_C_VALUES:
        cases = []
        for theta_max in _CONTINUUM_THETA_VALUES:
            tag = _continuum_tag(light_speed, theta_max)
            basename = 'pic_rel_gyro_continuum_' + tag
            _run_athena(
                'continuum_' + tag,
                ['job/basename=' + basename,
                 'particles/pic_cr_light_speed=' + str(light_speed),
                 'particles/pic_theta_max=' + str(theta_max),
                 # The inherited linear-wave pgen scales tlim by lambda/|v_wave|.
                 'time/tlim=' + str(_CONTINUUM_TEND/_LINEAR_WAVE_TLIM_SCALE),
                 'time/nlim=1000',
                 'output1/dcycle=1000',
                 'output2/dcycle=0'],
            )
            cases.append(_measure_continuum_case(basename, light_speed))
        _RESULTS['continuum'][light_speed] = cases


def _check_velocity(label, measured, expected, abs_tol=1.0e-6):
    err = float(np.max(np.abs(measured - expected[None, :])))
    logger.info('%s measured=%s expected=%s max_abs_err=% .8e',
                label, measured[0], expected, err)
    return err <= abs_tol


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for label in ['segment', 'full', 'restart']:
        result = _RESULTS[label]
        logger.info('%s path=%s time=% .8e cycle=%d npoint=%d',
                    label, result['path'], result['time'],
                    result['cycle'], result['npoint'])
        ok = result['npoint'] == 64 and ok
        ok = bool(np.all(np.isfinite(result['velocity']))) and ok
        ok = _check_velocity(label + ':boris',
                             result['velocity'],
                             result['expected_velocity']) and ok

    ok = _check_velocity('restart_vs_full',
                         _RESULTS['restart']['velocity'],
                         _RESULTS['full']['velocity']) and ok
    deltaf_box = _RESULTS['deltaf_box']
    logger.info('deltaf_box f0=[% .8e,% .8e] weight=[% .8e,% .8e]',
                np.min(deltaf_box['deltaf_f0']), np.max(deltaf_box['deltaf_f0']),
                np.min(deltaf_box['deltaf_weight']),
                np.max(deltaf_box['deltaf_weight']))
    ok = bool(np.all(np.isfinite(deltaf_box['deltaf_f0']))) and ok
    ok = bool(np.all(deltaf_box['deltaf_f0'] > 0.0)) and ok
    ok = bool(np.all(np.isfinite(deltaf_box['deltaf_weight']))) and ok
    ok = bool(np.max(np.abs(deltaf_box['deltaf_weight'])) > 1.0e-6) and ok
    for light_speed, cases in _RESULTS['continuum'].items():
        phase_errors = []
        timesteps = []
        for case in cases:
            logger.info(
                'continuum C=% .8e dt=% .8e cycles=%d phase_error=% .8e '
                'pnorm_rel_error=% .8e energy_rel_error=% .8e velocity_spread=% .8e',
                light_speed, case['dt'], case['cycle'], case['phase_error'],
                case['pnorm_rel_error'], case['energy_rel_error'],
                case['velocity_spread'])
            ok = abs(case['time'] - _CONTINUUM_TEND) <= 1.0e-12 and ok
            ok = case['npoint'] == 64 and ok
            ok = case['pnorm_rel_error'] <= 2.0e-6 and ok
            ok = case['energy_rel_error'] <= 2.0e-6 and ok
            ok = case['velocity_spread'] <= 1.0e-7 and ok
            phase_errors.append(case['phase_error'])
            timesteps.append(case['dt'])
        orders = [
            np.log(phase_errors[n]/phase_errors[n + 1])
            / np.log(timesteps[n]/timesteps[n + 1])
            for n in range(len(cases) - 1)
        ]
        logger.info('continuum C=% .8e phase_orders=%s', light_speed, orders)
        ok = bool(np.all(np.diff(timesteps) < 0.0)) and ok
        ok = bool(np.all(np.diff(phase_errors) < 0.0)) and ok
        ok = min(orders) >= 1.8 and ok
    return ok
