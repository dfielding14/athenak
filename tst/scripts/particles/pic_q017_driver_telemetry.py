"""Focused regression for bounded Q-017 driver telemetry."""

import glob
import logging
import math
import os
import re
import subprocess

import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_q017_driver_telemetry.athinput'
_BASENAME = 'pic_q017_driver_telemetry'
_RESULTS = {}

_REQUIRED_NAMES = {
    'schema_version',
    'mpi.ranks',
    'timer.driver.seconds_rank_max',
    'timer.task_lists.seconds_rank_max',
    'timer.task_lists.seconds_rank_mean',
    'timer.task_lists.calls_rank_max',
    'timer.task_list.before_timeintegrator.seconds_rank_max',
    'timer.task_list.before_stagen.seconds_rank_max',
    'timer.task_list.stagen.seconds_rank_max',
    'timer.task_list.after_stagen.seconds_rank_max',
    'timer.task_list.after_timeintegrator.seconds_rank_max',
    'timer.output_publication.seconds_rank_max',
    'timer.output_publication.calls_rank_max',
    'timer.amr_load_balance.seconds_rank_max',
    'timer.amr_load_balance.calls_rank_max',
    'timer.particle.adaptive_deltaf.seconds_rank_max',
    'timer.particle.adaptive_deltaf.calls_rank_max',
    'timer.particle.push.seconds_rank_max',
    'timer.particle.push.calls_rank_max',
    'timer.particle.deposition.seconds_rank_max',
    'timer.particle.deposition.calls_rank_max',
    'timer.particle.migration.seconds_rank_max',
    'timer.particle.migration.calls_rank_max',
    'cycles',
    'meshblocks.total',
    'meshblocks.rank_min',
    'meshblocks.rank_max',
    'mesh.cells_per_meshblock',
    'mesh.active_cells',
    'particles.total',
    'particles.rank_min',
    'particles.rank_max',
    'updates.meshblock_cycles',
    'updates.particle_updates',
    'throughput.zone_cycles_per_second',
    'throughput.particle_updates_per_second',
    'load.meshblock_efficiency',
    'load.particle_efficiency',
    'load.cost.total',
    'load.cost.rank_min',
    'load.cost.rank_max',
    'load.cost.efficiency',
    'load.cost.invalid_meshblocks',
    'amr.enabled',
    'amr.meshblocks_created',
    'amr.meshblocks_deleted',
    'amr.meshblocks_communicated',
    'particle_memory.sync_kernel_timers',
    'particle_memory.record_bytes',
    'particle_memory.root_level',
    'particle_memory.max_level',
    'particle_memory.resident_records.bytes_total',
    'particle_memory.direct_views.allocated_snapshot_bytes_total',
    'particle_memory.direct_views.allocated_snapshot_bytes_rank_max',
    'particle_memory.invalid_records',
    'particle_memory.species.0.count',
    'particle_memory.species.0.resident_bytes',
    'particle_memory.species.1.count',
    'particle_memory.species.1.resident_bytes',
    'particle_memory.level.1.count',
    'particle_memory.level.1.resident_bytes',
    'particle_memory.level.2.count',
    'particle_memory.level.2.resident_bytes',
    'particle_memory.species.0.level.1.count',
    'particle_memory.species.0.level.2.count',
    'particle_memory.species.1.level.1.count',
    'particle_memory.species.1.level.2.count',
}


def _athena_exe_dir():
    return os.environ.get('ATHENA_Q017_EXE_DIR',
                          os.path.join(os.getcwd(), 'build', 'src'))


def _athena_input_path():
    source_root = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               '..', '..', '..'))
    return os.path.join(source_root, 'inputs', _INPUT_DECK)


def _remove_outputs():
    pattern = os.path.join(_athena_exe_dir(), 'bin', _BASENAME + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _parse_telemetry(output):
    telemetry = {}
    pattern = re.compile(r'^q017\.telemetry\.([A-Za-z0-9_.]+)=([^\s]+)$',
                         re.MULTILINE)
    for match in pattern.finditer(output):
        name = match.group(1)
        if name in telemetry:
            raise RuntimeError('Duplicate Q-017 telemetry name: ' + name)
        telemetry[name] = float(match.group(2))
    return telemetry


def _check_equal(name, expected):
    measured = _RESULTS.get(name, float('nan'))
    logger.info('%s measured=% .8e expected=% .8e', name, measured, expected)
    return measured == expected


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    _remove_outputs()

    command = ['./athena', '-i', _athena_input_path()]
    logger.info('Executing Q-017 telemetry regression: %s', ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Q-017 telemetry command failed\n' + output)
    _RESULTS.update(_parse_telemetry(output))


def analyze():
    missing = sorted(_REQUIRED_NAMES - set(_RESULTS))
    if missing:
        logger.warning('Missing Q-017 telemetry names: %s', ', '.join(missing))
        return False

    ok = True
    for name, value in sorted(_RESULTS.items()):
        logger.info('%s=% .8e', name, value)
        ok = math.isfinite(value) and value >= 0.0 and ok

    ok = _check_equal('schema_version', 2.0) and ok
    ok = _check_equal('mpi.ranks', 1.0) and ok
    ok = _check_equal('cycles', 1.0) and ok
    ok = _check_equal('meshblocks.total', 15.0) and ok
    ok = _check_equal('mesh.cells_per_meshblock', 64.0) and ok
    ok = _check_equal('mesh.active_cells', 960.0) and ok
    ok = _check_equal('particles.total', 240.0) and ok
    ok = _check_equal('updates.meshblock_cycles', 15.0) and ok
    ok = _check_equal('updates.particle_updates', 240.0) and ok
    ok = _check_equal('timer.task_lists.calls_rank_max', 8.0) and ok
    ok = _check_equal('timer.output_publication.calls_rank_max', 3.0) and ok
    ok = _check_equal('timer.amr_load_balance.calls_rank_max', 1.0) and ok
    ok = _check_equal('load.cost.invalid_meshblocks', 0.0) and ok
    ok = _check_equal('amr.enabled', 1.0) and ok
    ok = _check_equal('timer.particle.adaptive_deltaf.calls_rank_max', 0.0) and ok
    ok = _check_equal('timer.particle.push.calls_rank_max', 1.0) and ok
    ok = _check_equal('timer.particle.deposition.calls_rank_max', 1.0) and ok
    ok = _check_equal('timer.particle.migration.calls_rank_max', 7.0) and ok
    ok = _check_equal('particle_memory.sync_kernel_timers', 1.0) and ok
    ok = _check_equal('particle_memory.record_bytes', 224.0) and ok
    ok = _check_equal('particle_memory.root_level', 1.0) and ok
    ok = _check_equal('particle_memory.max_level', 2.0) and ok
    ok = _check_equal('particle_memory.resident_records.bytes_total', 53760.0) and ok
    ok = _check_equal('particle_memory.invalid_records', 0.0) and ok
    ok = _check_equal('particle_memory.species.0.count', 120.0) and ok
    ok = _check_equal('particle_memory.species.1.count', 120.0) and ok
    ok = _check_equal('particle_memory.level.1.count', 112.0) and ok
    ok = _check_equal('particle_memory.level.2.count', 128.0) and ok
    ok = _check_equal('particle_memory.species.0.level.1.count', 56.0) and ok
    ok = _check_equal('particle_memory.species.0.level.2.count', 64.0) and ok
    ok = _check_equal('particle_memory.species.1.level.1.count', 56.0) and ok
    ok = _check_equal('particle_memory.species.1.level.2.count', 64.0) and ok
    return ok
