import logging
import os
from pathlib import Path
import subprocess

from scripts.particles import pic_paper_coupling_conservation as conservation

logger = logging.getLogger('athena' + __name__[7:])

_REPO_ROOT = Path(__file__).resolve().parents[3]
_EXE_DIR = _REPO_ROOT / 'tst' / 'build' / 'src'
_INPUT_DECK = _REPO_ROOT / 'inputs' / 'tests' / (
    'pic_paper_coupling_conservation.athinput')
_RESULTS = {}

_PARTICLE_CYCLE_CHAIN = [
    'id.adapt_deltaf = tl["before_timeintegrator"]->AddTask('
    '&Particles::AdaptDeltaF, this, none);',
    'id.save_old = tl["before_timeintegrator"]->AddTask('
    '&Particles::SaveOldPositions, this, id.adapt_deltaf);',
    'id.push = tl["before_timeintegrator"]->AddTask('
    '&Particles::Push, this, id.save_old);',
    'id.zero_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::ZeroMoments, this, id.push);',
    'id.irecv_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::InitRecvMoments, this, id.zero_mom);',
    'id.dep_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::DepositMoments, this, id.irecv_mom);',
    'id.rest_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::RestrictMoments, this, id.dep_mom);',
    'id.send_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::SendMoments, this, id.rest_mom);',
    'id.recv_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::RecvMoments, this, id.send_mom);',
    'id.crecv_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::ClearRecvMoments, this, id.recv_mom);',
    'id.csend_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::ClearSendMoments, this, id.crecv_mom);',
    'id.bcs_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::ApplyMomentPhysicalBCs, this, id.csend_mom);',
    'id.prol_mom = tl["before_timeintegrator"]->AddTask('
    '&Particles::ProlongateMoments, this, id.bcs_mom);',
]

_PAPER_STAGE_INSERT_CHAIN = [
    'TaskID insert_dep = (feedback_in_mhd_src ? pmhd->id.rkupdt : '
    'pmhd->id.efld);',
    'TaskID insert_loc = (feedback_in_mhd_src ? pmhd->id.srctrms : '
    'pmhd->id.efldsrc);',
    'stagen_tl->InsertTask(&Particles::SaveOldPositions, this, '
    'insert_dep, insert_loc);',
    'stagen_tl->InsertTask(&Particles::ZeroMoments, this, sid, insert_loc);',
    'stagen_tl->InsertTask(&Particles::InitRecvMoments, this, sid, insert_loc);',
    'stagen_tl->InsertTask(&Particles::DepositMoments, this, sid, insert_loc);',
    'stagen_tl->InsertTask(&Particles::RestrictMoments, this, sid, insert_loc);',
    'stagen_tl->InsertTask(&Particles::SendMoments, this, sid, insert_loc);',
    'stagen_tl->InsertTask(&Particles::RecvMoments, this, sid, insert_loc);',
    'stagen_tl->InsertTask(&Particles::ClearRecvMoments, this, sid, insert_loc);',
    'stagen_tl->InsertTask(&Particles::ClearSendMoments, this, sid, insert_loc);',
    'stagen_tl->InsertTask(&Particles::ApplyMomentPhysicalBCs, this, sid, '
    'insert_loc);',
    'stagen_tl->InsertTask(&Particles::ProlongateMoments, this, sid, '
    'insert_loc);',
]

_MHD_STAGE_CHAIN = [
    'id.copyu = tl["stagen"]->AddTask(&MHD::CopyCons, this, none);',
    'id.flux = tl["stagen"]->AddTask(&MHD::Fluxes, this, id.copyu);',
    'id.sendf = tl["stagen"]->AddTask(&MHD::SendFlux, this, id.flux);',
    'id.recvf = tl["stagen"]->AddTask(&MHD::RecvFlux, this, id.sendf);',
    'id.rkupdt = tl["stagen"]->AddTask(&MHD::RKUpdate, this, id.recvf);',
    'id.srctrms = tl["stagen"]->AddTask(&MHD::MHDSrcTerms, this, id.rkupdt);',
    'id.picdampu = tl["stagen"]->AddTask(&MHD::ApplyPICWaveDamping, this, '
    'id.srctrms);',
    'id.efld = tl["stagen"]->AddTask(&MHD::CornerE, this, id.picdampu);',
    'id.efldsrc = tl["stagen"]->AddTask(&MHD::EFieldSrc, this, id.efld);',
    'id.sende = tl["stagen"]->AddTask(&MHD::SendE, this, id.efldsrc);',
    'id.recve = tl["stagen"]->AddTask(&MHD::RecvE, this, id.sende);',
    'id.ct = tl["stagen"]->AddTask(&MHD::CT, this, id.recve);',
    'id.expboxb = tl["stagen"]->AddTask(&MHD::ApplyPICExpandingBoxB, this, '
    'id.ct);',
    'id.expboxu = tl["stagen"]->AddTask(&MHD::ApplyPICExpandingBoxU, this, '
    'id.expboxb);',
    'id.expboxfb = tl["stagen"]->AddTask(&MHD::ApplyPICExpandingBoxFeedback, '
    'this, id.expboxu);',
    'id.expboxdampu = tl["stagen"]->AddTask(&MHD::ApplyPICExpandingBoxWaveDamping, this, '
    'id.expboxfb);',
    'id.sendu_oa = tl["stagen"]->AddTask(&MHD::SendU_OA, this, id.expboxdampu);',
    'id.recvu_oa = tl["stagen"]->AddTask(&MHD::RecvU_OA, this, id.sendu_oa);',
    'id.restu = tl["stagen"]->AddTask(&MHD::RestrictU, this, id.recvu_oa);',
    'id.sendu = tl["stagen"]->AddTask(&MHD::SendU, this, id.restu);',
    'id.recvu = tl["stagen"]->AddTask(&MHD::RecvU, this, id.sendu);',
    'id.sendu_shr = tl["stagen"]->AddTask(&MHD::SendU_Shr, this, id.recvu);',
    'id.recvu_shr = tl["stagen"]->AddTask(&MHD::RecvU_Shr, this, id.sendu_shr);',
]

_RUNTIME_TOKENS = [
    'physical_mode=paper_mhd_pic',
    'state=momentum_p_over_m',
    'C=3',
    'background=coupled',
    'feedback=coupled',
    'induction=ideal_mhd_only',
    'deposition=tsc',
    'wave_damping=off',
    'restart_schema=8',
]


def _normalized(path):
    return ' '.join(path.read_text(encoding='ascii').split())


def _require_ordered(source, snippets, label):
    position = -1
    for snippet in snippets:
        normalized = ' '.join(snippet.split())
        next_position = source.find(normalized, position + 1)
        if next_position < 0:
            raise RuntimeError(label + ' missing ordered source fragment:\n'
                               + normalized)
        position = next_position


def _check_static_task_graph():
    particles_tasks = _normalized(
        _REPO_ROOT / 'src' / 'particles' / 'particles_tasks.cpp')
    particles_hpp = _normalized(
        _REPO_ROOT / 'src' / 'particles' / 'particles.hpp')
    mhd_tasks = _normalized(_REPO_ROOT / 'src' / 'mhd' / 'mhd_tasks.cpp')
    deck = _normalized(_INPUT_DECK)

    _require_ordered(particles_tasks, _PARTICLE_CYCLE_CHAIN,
                     'particle cycle task chain')
    _require_ordered(particles_tasks, _PAPER_STAGE_INSERT_CHAIN,
                     'paper stage insertion chain')
    _require_ordered(mhd_tasks, _MHD_STAGE_CHAIN, 'MHD stage task chain')

    required_fragments = [
        ('particle task graph',
         particles_tasks,
         'couple_fluid_feedback_order == '
         'CoupledFluidFeedbackOrder::mhd_src_terms'),
        ('particle task graph',
         particles_tasks,
         'auto comm_tl = (couple_moments_to_mhd ? '
         'tl["after_timeintegrator"] : tl["before_timeintegrator"]);'),
        ('CT source gate',
         mhd_tasks,
         'if ((ppart != nullptr) && ppart->AddsCRCurrentToCT()) {'),
        ('CT source identity',
         particles_hpp,
         '(pic_physical_mode == PICPhysicalMode::engineering) || '
         '((pic_physical_mode == PICPhysicalMode::extended_mhd_pic) && '
         '(pic_cr_hall_mode == '
         'PICCRHallMode::current_to_ct_experimental))'),
        ('paper fixture',
         deck,
         'pic_physical_mode = paper_mhd_pic'),
        ('paper fixture',
         deck,
         'couple_fluid_feedback_order = mhd_src_terms'),
        ('paper fixture',
         deck,
         'pic_cr_hall_mode = off'),
        ('paper fixture',
         deck,
         'pic_wave_damping_mode = off'),
        ('paper fixture',
         deck,
         'pic_expanding_box_mode = off'),
    ]
    for label, source, fragment in required_fragments:
        if fragment not in source:
            raise RuntimeError(label + ' missing source fragment:\n' + fragment)

    _RESULTS['particle_cycle_tasks'] = len(_PARTICLE_CYCLE_CHAIN)
    _RESULTS['paper_stage_insert_tasks'] = len(_PAPER_STAGE_INSERT_CHAIN) - 2
    _RESULTS['mhd_stage_tasks_through_expanding_box_u'] = len(_MHD_STAGE_CHAIN)
    _RESULTS['static_graph'] = 'pass'


def _run_runtime_identity():
    command = [
        './athena', '-i', os.path.relpath(_INPUT_DECK, _EXE_DIR),
        'time/nlim=0',
        'particles/couple_j_to_efield_coeff=7.0',
    ]
    logger.info('Executing runtime identity trace: %s', ' '.join(command))
    proc = subprocess.run(command, cwd=_EXE_DIR, capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Runtime identity trace failed\n' + output)
    trace_lines = [
        line for line in output.splitlines() if 'PIC runtime model:' in line
    ]
    if len(trace_lines) != 1:
        raise RuntimeError('Expected one PIC runtime identity trace\n' + output)
    trace = trace_lines[0]
    for token in _RUNTIME_TOKENS:
        if token not in trace:
            raise RuntimeError('Runtime identity trace missing token: ' + token
                               + '\n' + trace)
    _RESULTS['runtime_identity'] = trace


def _run_dynamic_induction_isolation():
    original_cwd = os.getcwd()
    try:
        os.chdir(_REPO_ROOT / 'tst')
        conservation.run()
        if not conservation.analyze():
            raise RuntimeError('pic_paper_coupling_conservation: FAIL')
    finally:
        os.chdir(original_cwd)
    _RESULTS['dynamic_induction_isolation'] = 'pass'


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _check_static_task_graph()
    _run_runtime_identity()
    _run_dynamic_induction_isolation()


def analyze():
    logger.info('PIC paper task-stage trace: %s', _RESULTS)
    return (_RESULTS.get('static_graph') == 'pass'
            and _RESULTS.get('dynamic_induction_isolation') == 'pass'
            and 'induction=ideal_mhd_only' in _RESULTS.get('runtime_identity', ''))
