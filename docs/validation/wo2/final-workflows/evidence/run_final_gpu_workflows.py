#!/usr/bin/env python3
"""Prepare or run the exact final binary through full and paper-smoke workflows.

Default is preparation only. --execute requires the parent's assigned allocation.
No builds, submissions, scientific overrides or source edits are performed.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
WO2 = HERE.parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def contained(path):
    path = Path(path).resolve()
    if not path.is_relative_to(WO2):
        raise ValueError(f'Final validation path must be below {WO2}: {path}')
    return path

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--build-dir', type=Path, required=True)
    ap.add_argument('--run-name', default='acceptance-hip')
    ap.add_argument('--expected-binary-sha256')
    ap.add_argument('--git-metadata-source', type=Path)
    ap.add_argument('--execute', action='store_true')
    args = ap.parse_args()
    source, binary, build = map(contained, (args.source, args.binary, args.build_dir))
    script = source/'scripts/cgl_lf_workflow.py'
    if not script.is_file() or not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError('An existing final source tree and executable are required')
    if args.expected_binary_sha256 and sha(binary) != args.expected_binary_sha256:
        raise ValueError('Supplied final binary does not match the expected SHA256')
    output = (HERE/args.run_name).resolve()
    if not output.is_relative_to(HERE) or output == HERE:
        raise ValueError('Run output must be a child of final-research')
    if output.exists():
        raise ValueError(f'Refusing to overwrite retained run: {output}')
    output.mkdir(parents=True)
    os.chdir(HERE)
    module_spec = importlib.util.spec_from_file_location('wo2_final_workflow', script)
    workflow = importlib.util.module_from_spec(module_spec)
    sys.modules[module_spec.name] = workflow
    module_spec.loader.exec_module(workflow)
    if workflow.ROOT_DIR.resolve() != source:
        raise ValueError('Workflow ROOT_DIR must be the supplied WO2 source tree; stage a real script, not an external symlink')
    plans = {}
    for name in ('full', 'paper-smoke'):
        cases = workflow.workflow_cases(name)
        plans[name] = {'case_names':[c.name for c in cases],
            'inputs':{c.input_path:sha(source/c.input_path) for c in cases},
            'command':[sys.executable, str(script), name, '--no-build',
                '--build-dir',str(build), '--athena-bin',str(HERE/'launch_final.py'),
                '--output-dir',str(output/name)]}
    if len(plans['full']['case_names']) != 29:
        raise ValueError('Final full workflow no longer contains the expected 29 cases')
    source_hashes = {}
    for directory, _, files in os.walk(source/'src', followlinks=True):
        for name in files:
            path=Path(directory)/name
            if path.suffix in ('.cpp','.hpp','.h') or path.name=='CMakeLists.txt':
                source_hashes[str(path.relative_to(source))]=sha(path)
    environment = dict(os.environ)
    for key in list(environment):
        if key.startswith('ATHENAK_CGL_LF_') or key.startswith('KOKKOS_PROFILE'):
            del environment[key]
    environment.update(WO2_FINAL_BINARY=str(binary), WO2_FINAL_SHA256=sha(binary),
        WO2_FINAL_LAUNCH_LOG=str(output/'launches.jsonl'),
        PATH=str(HERE)+os.pathsep+environment['PATH'], PYTHONDONTWRITEBYTECODE='1',
        TMPDIR=str(HERE/'tmp'), MPLCONFIGDIR=str(HERE/'runtime/mpl-cache'),
        XDG_CACHE_HOME=str(HERE/'runtime/cache'))
    for key in ('TMPDIR','MPLCONFIGDIR','XDG_CACHE_HOME'):
        Path(environment[key]).mkdir(parents=True,exist_ok=True)
    git_metadata = None
    if args.git_metadata_source:
        metadata_source=args.git_metadata_source.resolve()
        gitdir=subprocess.check_output(['git','-C',str(metadata_source),
            'rev-parse','--absolute-git-dir'],text=True).strip()
        environment.update(GIT_DIR=gitdir,GIT_WORK_TREE=str(source),GIT_OPTIONAL_LOCKS='0')
        config_index=int(environment.get('GIT_CONFIG_COUNT','0'))
        environment['GIT_CONFIG_COUNT']=str(config_index+1)
        environment[f'GIT_CONFIG_KEY_{config_index}']='submodule.kokkos.ignore'
        environment[f'GIT_CONFIG_VALUE_{config_index}']='all'
        revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE,env=environment,text=True).strip()
        status=subprocess.check_output(['git','status','--short','--untracked-files=all'],
            cwd=HERE,env=environment,text=True)
        (output/'source-git-status.txt').write_text(status)
        git_metadata={'metadata_source':str(metadata_source),'git_dir':gitdir,
            'work_tree':str(source),'optional_locks':False,'revision':revision,
            'submodule_status':'kokkos ignored only for Git status because snapshot .git pointer is relocated; binary/source provenance retained separately',
            'status_sha256':sha(output/'source-git-status.txt')}
    report = {'status':'prepared', 'created_utc':datetime.now(timezone.utc).isoformat(),
        'git_metadata':git_metadata, 'source':str(source), 'binary':str(binary), 'binary_sha256':sha(binary),
        'workflow_sha256':sha(script), 'launcher_sha256':sha(HERE/'launch_final.py'),
        'runner_sha256':sha(__file__), 'build_dir':str(build), 'plans':plans,
        'source_hashes':source_hashes,
        'cmake_cache_sha256':sha(build/'CMakeCache.txt') if (build/'CMakeCache.txt').exists() else None,
        'allocation':environment.get('SLURM_JOB_ID'), 'results':[]}
    def save():
        (output/'runner-manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    if not args.execute:
        print(json.dumps({'prepared':str(output),'cases':{k:len(v['case_names']) for k,v in plans.items()}},indent=2))
        return 0
    if not environment.get('SLURM_JOB_ID'):
        raise ValueError('--execute requires a coordinated Slurm allocation')
    report['status'] = 'running'; save()
    for name, plan in plans.items():
        if sha(binary) != report['binary_sha256'] or sha(script) != report['workflow_sha256']:
            raise ValueError('Final binary or workflow changed during acceptance')
        with (output/f'{name}.log').open('w') as stream:
            result = subprocess.run(plan['command'], cwd=HERE, env=environment,
                                    stdout=stream, stderr=subprocess.STDOUT)
        manifest_file=output/name/'manifest.json'
        manifest=json.loads(manifest_file.read_text()) if manifest_file.exists() else {}
        actual=[c['name'] for c in manifest.get('cases',[])]
        passed=(result.returncode==0 and manifest.get('status')=='passed'
                and actual==plan['case_names'] and not manifest.get('disabled_cases',[]))
        report['results'].append({'workflow':name,'returncode':result.returncode,
            'passed':passed,'actual_cases':actual,'expected_cases':plan['case_names'],
            'disabled_cases':manifest.get('disabled_cases',[]),
            'manifest_sha256':sha(manifest_file) if manifest_file.exists() else None})
        save(); print(name, 'passed' if passed else 'FAILED', flush=True)
    report['binary_unchanged']=sha(binary)==report['binary_sha256']
    report['workflow_unchanged']=sha(script)==report['workflow_sha256']
    report['source_unchanged']=all(sha(source/path)==value for path,value in source_hashes.items())
    report['inputs_unchanged']=all(sha(source/path)==value for plan in plans.values() for path,value in plan['inputs'].items())
    report['status']='passed' if (all(x['passed'] for x in report['results']) and
        report['binary_unchanged'] and report['workflow_unchanged'] and
        report['source_unchanged'] and report['inputs_unchanged']) else 'failed'
    save()
    return 0 if report['status']=='passed' else 1

if __name__=='__main__':
    raise SystemExit(main())
