import hashlib
import json
import logging
import os
import subprocess

import numpy as np

from scripts.particles import pic_entity_deposit_mink
from scripts.particles import pic_entity_deposit_reflect

logger = logging.getLogger('athena' + __name__[7:])

_ATHENAK_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..'))
_SNAPSHOT_PATH = os.path.join(
    _ATHENAK_ROOT, 'tst', 'publication', 'readiness', 'entity_snapshot.json')
_ENTITY_SHAPES_PATH = 'src/kernels/particle_shapes.hpp'
_ENTITY_SHAPES_SHA256 = (
    'ca6e1acbd182d44c057ac8286777106e9fc7192e6f5fffe8d42effb64192641c')
_RESULTS = {}


def _run_output(command, cwd):
    proc = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed: ' + ' '.join(command) + '\n' + output)
    return (proc.stdout or '').strip()


def _sha256(contents):
    return hashlib.sha256(contents).hexdigest()


def _read_bytes(path):
    with open(path, 'rb') as fp:
        return fp.read()


def _entity_blob(entity_path, commit, path):
    return subprocess.check_output(
        ['git', 'show', commit + ':' + path], cwd=entity_path)


def _require_markers(contents, markers, label):
    text = contents.decode('utf-8')
    missing = [marker for marker in markers if marker not in text]
    if missing:
        raise RuntimeError(label + ' is missing source markers: ' + repr(missing))


def _verify_entity_snapshot():
    with open(_SNAPSHOT_PATH, encoding='utf-8') as fp:
        snapshot = json.load(fp)

    entity_path = snapshot['path']
    commit = snapshot['git_commit']
    tree = snapshot['git_tree']
    observed_commit = _run_output(['git', 'rev-parse', 'HEAD'], entity_path)
    observed_tree = _run_output(['git', 'show', '-s', '--format=%T', 'HEAD'],
                                entity_path)
    worktree_status = _run_output(['git', 'status', '--porcelain'], entity_path)

    file_hashes_match = True
    for entry in snapshot['files']:
        path = os.path.join(entity_path, entry['path'])
        file_hashes_match = (
            _sha256(_read_bytes(path)) == entry['sha256'] and file_hashes_match)

    shapes_blob = _entity_blob(entity_path, commit, _ENTITY_SHAPES_PATH)
    _require_markers(
        _entity_blob(entity_path, commit, 'src/kernels/pushers/sr.hpp'),
        ['Inline void velocityEMPush_Boris', 'COEFF *= ONE / math::sqrt',
         'u0[0] += CROSS_x1'],
        'Entity SR Boris pusher')
    _require_markers(
        shapes_blob,
        ['Inline void order', 'Inline void for_deposit',
         'S[0]  = ONE - di', 'S[1]  = di'],
        'Entity particle shapes')
    _require_markers(
        _read_bytes(os.path.join(
            _ATHENAK_ROOT, 'src', 'particles', 'particles_pushers.cpp')),
        ['// Magnetic rotation', 'Real tx = qdt_2m*Bx*inv_gamma_minus;',
         'state_x += state_py*rot_z - state_pz*rot_y;'],
        'AthenaK SR Boris pusher')
    _require_markers(
        _read_bytes(os.path.join(
            _ATHENAK_ROOT, 'src', 'particles', 'particles_moments.cpp')),
        ['void ShapeOrder', 'void ShapeForDeposit',
         'S[0] = static_cast<Real>(1.0) - di;', 'S[1] = di;'],
        'AthenaK particle shapes')

    return {
        'commit': commit,
        'tree': tree,
        'commit_match': observed_commit == commit,
        'tree_match': observed_tree == tree,
        'worktree_clean': worktree_status == '',
        'declared_file_hashes_match': file_hashes_match,
        'particle_shapes_blob_hash_match': _sha256(shapes_blob) ==
                                           _ENTITY_SHAPES_SHA256,
    }


def _entity_boris_step(state, electric, magnetic, q_over_m, dt):
    coeff = 0.5*q_over_m*dt
    e0 = coeff*np.array(electric, dtype=np.float64)
    u0 = np.array(state, dtype=np.float64) + e0
    coeff /= np.sqrt(1.0 + np.dot(u0, u0))
    b0 = coeff*np.array(magnetic, dtype=np.float64)
    coeff = 2.0/(1.0 + np.dot(b0, b0))
    u1 = (u0 + np.cross(u0, b0))*coeff
    return u0 + np.cross(u1, b0) + e0


def _athenak_boris_step(state, electric, magnetic, q_over_m, dt):
    qdt_2m = 0.5*q_over_m*dt
    state_minus = np.array(state, dtype=np.float64)
    state_minus += qdt_2m*np.array(electric, dtype=np.float64)
    inv_gamma_minus = 1.0/np.sqrt(1.0 + np.dot(state_minus, state_minus))
    tvec = qdt_2m*np.array(magnetic, dtype=np.float64)*inv_gamma_minus
    rot = 2.0*tvec/(1.0 + np.dot(tvec, tvec))
    state_prime = state_minus + np.cross(state_minus, tvec)
    state_plus = state_minus + np.cross(state_prime, rot)
    return state_plus + qdt_2m*np.array(electric, dtype=np.float64)


def _entity_test_expected_momentum(state0, magnetic, time):
    bmag = np.linalg.norm(magnetic)
    bhat = magnetic/bmag
    udotb = np.dot(state0, bhat)
    phase = bmag*time/np.sqrt(1.0 + np.dot(state0, state0))
    return (bhat*udotb + (state0 - bhat*udotb)*np.cos(phase) +
            np.cross(state0, bhat)*np.sin(phase))


def _measure_boris_overlap():
    initial = np.array([1.0, -2.0, 0.1])
    magnetic = 0.2*np.array([0.66, 0.55, 0.44])
    electric = np.zeros(3)
    dt = 0.01
    entity_state = initial.copy()
    athenak_state = initial.copy()
    max_diff = 0.0
    max_continuum_error = 0.0
    max_pnorm_drift = 0.0

    # Reuse the frozen Entity pusher regression's 2,000-step B-only trajectory.
    for step in range(2000):
        entity_state = _entity_boris_step(
            entity_state, electric, magnetic, 1.0, dt)
        athenak_state = _athenak_boris_step(
            athenak_state, electric, magnetic, 1.0, dt)
        expected = _entity_test_expected_momentum(initial, magnetic,
                                                  (step + 1)*dt)
        max_diff = max(max_diff, float(np.max(np.abs(
            entity_state - athenak_state))))
        max_continuum_error = max(
            max_continuum_error,
            float(np.max(np.abs(entity_state - expected))))
        max_pnorm_drift = max(
            max_pnorm_drift,
            float(abs(np.linalg.norm(entity_state) - np.linalg.norm(initial))))

    generic_entity = np.array([0.45, -0.25, 0.15])
    generic_athenak = generic_entity.copy()
    generic_electric = np.array([0.07, -0.03, 0.02])
    generic_magnetic = np.array([0.30, 0.10, -0.20])
    for _ in range(64):
        generic_entity = _entity_boris_step(
            generic_entity, generic_electric, generic_magnetic, -0.6, 0.025)
        generic_athenak = _athenak_boris_step(
            generic_athenak, generic_electric, generic_magnetic, -0.6, 0.025)

    return {
        'entity_regression_steps': 2000,
        'entity_vs_athenak_max_abs_diff': max_diff,
        'entity_vs_continuum_max_abs_error': max_continuum_error,
        'entity_b_only_max_pnorm_drift': max_pnorm_drift,
        'generic_eb_64_step_max_abs_diff': float(np.max(np.abs(
            generic_entity - generic_athenak))),
    }


def _entity_shape_order(staggered, order, index, frac):
    if order == 1:
        if not staggered:
            return index, np.array([1.0 - frac, frac])
        if frac < 0.5:
            return index - 1, np.array([0.5 - frac, 0.5 + frac])
        return index, np.array([1.5 - frac, frac - 0.5])
    if not staggered:
        if frac < 0.5:
            shape0 = 0.5*(0.5 - frac)*(0.5 - frac)
            shape1 = 0.75 - frac*frac
            return index - 1, np.array([shape0, shape1, 1.0 - shape0 - shape1])
        shape0 = 0.5*(1.5 - frac)*(1.5 - frac)
        shape1 = 0.75 - (1.0 - frac)*(1.0 - frac)
        return index, np.array([shape0, shape1, 1.0 - shape0 - shape1])
    shape0 = 0.5*(1.0 - frac)*(1.0 - frac)
    shape2 = 0.5*frac*frac
    return index - 1, np.array([shape0, 1.0 - shape0 - shape2, shape2])


def _athenak_shape_order(staggered, order, index, frac):
    if order == 1:
        if not staggered:
            index_min = index
            shape = np.array([1.0 - frac, frac])
        elif frac < 0.5:
            index_min = index - 1
            shape0 = 0.5 - frac
            shape = np.array([shape0, 1.0 - shape0])
        else:
            index_min = index
            shape0 = 1.5 - frac
            shape = np.array([shape0, 1.0 - shape0])
        return index_min, shape
    if not staggered:
        if frac < 0.5:
            index_min = index - 1
            shape0 = 0.5*(0.5 - frac)*(0.5 - frac)
            shape1 = 0.75 - frac*frac
        else:
            index_min = index
            shape0 = 0.5*(1.5 - frac)*(1.5 - frac)
            d1 = 1.0 - frac
            shape1 = 0.75 - d1*d1
        return index_min, np.array([shape0, shape1, 1.0 - shape0 - shape1])
    index_min = index - 1
    shape0 = 0.5*(1.0 - frac)*(1.0 - frac)
    shape2 = 0.5*frac*frac
    return index_min, np.array([shape0, 1.0 - shape0 - shape2, shape2])


def _shape_for_deposit(shape_order, order, initial, final):
    initial_min, initial_shape = shape_order(False, order, *initial)
    final_min, final_shape = shape_order(False, order, *final)
    initial_out = np.zeros(order + 2)
    final_out = np.zeros(order + 2)
    if initial_min < final_min:
        index_min = initial_min
        index_max = index_min + order + 1
        initial_out[:order + 1] = initial_shape
        final_out[1:order + 2] = final_shape
    elif initial_min > final_min:
        index_min = final_min
        index_max = index_min + order + 1
        initial_out[1:order + 2] = initial_shape
        final_out[:order + 1] = final_shape
    else:
        index_min = initial_min
        index_max = index_min + order
        initial_out[:order + 1] = initial_shape
        final_out[:order + 1] = final_shape
    return index_min, index_max, initial_out, final_out


def _measure_shape_overlap():
    max_shape_diff = 0.0
    max_partition_error = 0.0
    shape_evaluations = 0
    for order in [1, 2]:
        for staggered in [False, True]:
            for frac in [0.0, 0.1, 0.49, 0.5, 0.65, 0.99, 1.0]:
                entity_min, entity_shape = _entity_shape_order(
                    staggered, order, 7, frac)
                athenak_min, athenak_shape = _athenak_shape_order(
                    staggered, order, 7, frac)
                if entity_min != athenak_min:
                    raise RuntimeError('ShapeOrder index differential mismatch')
                max_shape_diff = max(
                    max_shape_diff,
                    float(np.max(np.abs(entity_shape - athenak_shape))))
                max_partition_error = max(
                    max_partition_error,
                    float(abs(np.sum(entity_shape) - 1.0)))
                shape_evaluations += 1

    trajectories = [
        ((25, 0.65), (24, 0.99)),
        ((21, 0.65), (20, 0.80)),
        ((2, 0.20), (2, 0.80)),
        ((2, 0.80), (3, 0.10)),
        ((3, 0.10), (2, 0.80)),
    ]
    trajectory_evaluations = 0
    for order in [1, 2]:
        for initial, final in trajectories:
            entity = _shape_for_deposit(
                _entity_shape_order, order, initial, final)
            athenak = _shape_for_deposit(
                _athenak_shape_order, order, initial, final)
            if entity[:2] != athenak[:2]:
                raise RuntimeError('ShapeForDeposit index differential mismatch')
            max_shape_diff = max(
                max_shape_diff,
                float(np.max(np.abs(entity[2] - athenak[2]))),
                float(np.max(np.abs(entity[3] - athenak[3]))))
            trajectory_evaluations += 1

    return {
        'shape_order_evaluations': shape_evaluations,
        'trajectory_support_evaluations': trajectory_evaluations,
        'entity_vs_athenak_max_abs_diff': max_shape_diff,
        'max_partition_of_unity_error': max_partition_error,
    }


def _run_existing_deposit_regression(label, module):
    logger.info('Running existing AthenaK regression: %s', label)
    module._RESULTS.clear()
    module.run()
    return bool(module.analyze())


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    _RESULTS['snapshot'] = _verify_entity_snapshot()
    _RESULTS['boris'] = _measure_boris_overlap()
    _RESULTS['shapes'] = _measure_shape_overlap()
    _RESULTS['athenak_regressions'] = {
        'pic_entity_deposit_mink': _run_existing_deposit_regression(
            'pic_entity_deposit_mink', pic_entity_deposit_mink),
        'pic_entity_deposit_reflect': _run_existing_deposit_regression(
            'pic_entity_deposit_reflect', pic_entity_deposit_reflect),
    }


def analyze():
    logger.debug('Analyzing test ' + __name__)
    snapshot = _RESULTS['snapshot']
    boris = _RESULTS['boris']
    shapes = _RESULTS['shapes']
    regressions = _RESULTS['athenak_regressions']

    logger.info('entity snapshot commit=%s tree=%s', snapshot['commit'],
                snapshot['tree'])
    logger.info('boris exact-overlap metrics=%s', boris)
    logger.info('shape-support exact-overlap metrics=%s', shapes)
    logger.info('existing AthenaK Entity-style regressions=%s', regressions)

    ok = all([
        snapshot['commit_match'],
        snapshot['tree_match'],
        snapshot['worktree_clean'],
        snapshot['declared_file_hashes_match'],
        snapshot['particle_shapes_blob_hash_match'],
        boris['entity_vs_athenak_max_abs_diff'] <= 1.0e-14,
        boris['entity_vs_continuum_max_abs_error'] <= 2.0e-7,
        boris['entity_b_only_max_pnorm_drift'] <= 1.0e-13,
        boris['generic_eb_64_step_max_abs_diff'] <= 1.0e-14,
        shapes['entity_vs_athenak_max_abs_diff'] == 0.0,
        shapes['max_partition_of_unity_error'] <= 1.0e-15,
        all(regressions.values()),
    ])
    return bool(ok)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    run()
    if not analyze():
        raise SystemExit('pic_entity_differential: FAIL')
    print('pic_entity_differential: PASS')
