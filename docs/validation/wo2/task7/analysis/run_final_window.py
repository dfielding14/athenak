#!/usr/bin/env python3
"""Exclusive actual-GPU correctness and unprofiled paired timing for Tasks7/3."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import struct
import subprocess
import sys
import time

import numpy as np

ROOT=Path('/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2')
BASE=Path('/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2-baseline')
sys.path.insert(0,str(ROOT/'p1-research'))
sys.path.insert(0,str(BASE/'vis/python'))
from restart_compare import restart_hashes, normalized_restart
from athena_read import hst

p=argparse.ArgumentParser()
p.add_argument('--phase',choices=['task7-check','task7-time'],required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--manifest',type=Path,required=True)
p.add_argument('--extra-input',action='append',default=[],help='Additional check fixture NAME=PATH')
p.add_argument('--repeats',type=int,default=3)
p.add_argument('--cases',nargs='+')
a=p.parse_args()
final_manifest=json.loads(a.manifest.read_text())
assert final_manifest['backend']=='hip', 'This runner reserves one actual GPU per rank'
for group in ('binaries','trace_binaries'):
    for item in final_manifest[group].values():
        assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['sha256']
FINAL_SOURCE=Path(final_manifest['source'])
assert os.environ.get('SLURM_JOB_ID'),'requires parent-managed exclusive allocation window'
assert a.repeats>=3
out=a.output.resolve();assert str(out).startswith(str(ROOT)+'/')
out.mkdir(parents=True,exist_ok=False);(out/'inputs').mkdir()

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def input_file(name, source, parameters, timing=False):
    # Give both binaries exactly the same saved input, including all overrides.
    text=source.read_text();blocks={};current=None
    for line in text.splitlines():
        match=re.match(r'\s*<([^>]+)>',line)
        if match:
            current=match[1];blocks.setdefault(current,[])
        elif current:blocks[current].append(line)
    for key,value in parameters.items():
        block,param=key.split('/',1);lines=blocks.setdefault(block,[])
        pattern=re.compile(r'^\s*'+re.escape(param)+r'\s*=')
        indices=[i for i,line in enumerate(lines) if pattern.match(line)]
        assert len(indices)<=1,key
        line=f'{param} = {value}'
        if indices:lines[indices[0]]=line
        else:lines.append(line)
    for block in list(blocks):
        if re.fullmatch('output[0-9]+',block):del blocks[block]
    blocks['output1']=['file_type = hst','data_format = %24.16e','dcycle = 1']
    blocks['output2']=['file_type = bin','variable = mhd_w_bcc','dcycle = 1000000']
    blocks['output3']=['file_type = rst','dcycle = 1000000']
    path=out/'inputs'/(name+'.athinput')
    path.write_text('\n'.join('<'+b+'>\n'+'\n'.join(lines)+'\n' for b,lines in blocks.items()))
    return path

def compare(left,right,allow_work_roundoff):
    result={'failures':[],'files':{}}
    inventory=lambda d:{str(p.relative_to(d)):p for p in d.rglob('*')
                        if p.is_file() and p.suffix in ('.bin','.rst','.hst')}
    l,r=inventory(left),inventory(right)
    if l.keys()!=r.keys():result['failures'].append('file inventory differs')
    for name in sorted(l.keys()&r.keys()):
        row={'left_raw_sha256':sha(l[name]),'right_raw_sha256':sha(r[name])}
        if name.endswith('.rst'):
            x,y=restart_hashes(l[name]),restart_hashes(r[name])
            row.update(left_restart=x,right_restart=y)
            row['equal']=x['normalized_sha256']==y['normalized_sha256']
            if not row['equal'] and allow_work_roundoff:
                lx,_=normalized_restart(l[name]);ry,_=normalized_restart(r[name])
                header=x['binary_header_offset']
                assert header==y['binary_header_offset'] and len(lx)==len(ry)
                params=lx[:header].decode()
                assert '<turb_driving>' not in params and '<z4c>' not in params
                nmb,=struct.unpack_from('=i',lx,header)
                offset=header+252+nmb*20+13*8
                vx=struct.unpack_from('=2d',lx,offset);vy=struct.unpack_from('=2d',ry,offset)
                physical=(lx[:offset]==ry[:offset] and lx[offset+16:]==ry[offset+16:])
                close=bool(np.allclose(vx,vy,rtol=256*np.finfo(float).eps,atol=0.0))
                different=np.flatnonzero(np.frombuffer(lx,dtype='u1')!=np.frombuffer(ry,dtype='u1'))
                row['named_diagnostic_comparison']={
                    'members':['lf_diag[13]:qpar_work','lf_diag[14]:qperp_work'],
                    'range':[offset,offset+16],'left':vx,'right':vy,
                    'different_offsets':different.tolist(),
                    'all_other_normalized_bytes_identical':physical,
                    'within_256epsilon_relative':close,
                    'source':'src/outputs/restart.cpp:271-292,342',
                    'normalization_applied':False}
                row['equal']=physical and close
        elif name.endswith('.hst'):
            x,y=hst(str(l[name])),hst(str(r[name]))
            row['equal']=x.keys()==y.keys();row['columns']={}
            for key in x.keys()&y.keys():
                exact=np.array_equal(x[key],y[key])
                same_shape=x[key].shape==y[key].shape
                allowed=(same_shape and allow_work_roundoff and key in ('lf_qprwrk','lf_qpewrk')
                         and np.allclose(x[key],y[key],rtol=256*np.finfo(float).eps,atol=0.0))
                row['columns'][key]={'exact':exact,'allowed_reduction_roundoff':bool(allowed),
                    'max_abs':float(np.max(np.abs(x[key]-y[key]))) if x[key].shape==y[key].shape else None}
                row['equal'] &= exact or allowed
        else:row['equal']=row['left_raw_sha256']==row['right_raw_sha256']
        if not row['equal']:result['failures'].append(name)
        result['files'][name]=row
    return result

def compare_faces(left,right):
    result={'snapshots':0,'values':0,'failures':[]}
    l={p.name:p for p in (left/'trace').glob('*.json')}
    r={p.name:p for p in (right/'trace').glob('*.json')}
    if not l and not r:
        for directory in (left,right):
            history=hst(str(next(directory.glob('*.mhd.hst'))))
            assert np.all(history['lf_nstage']==0),'missing trace despite active LF stages'
        result['reason']='LF disabled throughout below-field-floor case; zero LF stages'
        return result
    assert l and l.keys()==r.keys(),(l.keys(),r.keys())
    for name in sorted(l):
        x,y=json.loads(l[name].read_text()),json.loads(r[name].read_text())
        assert x['gids']==y['gids'] and x['active_kji']==y['active_kji']
        for axis in (1,2,3):
            bounds=[list(z) for z in x['active_kji']]
            if axis>1 and bounds[3-axis][0]==bounds[3-axis][1]:continue
            bounds[3-axis][1]+=1
            sl=(slice(None),slice(4,6))+tuple(slice(lo,hi+1) for lo,hi in bounds)
            arrays=[]
            for file,d in ((l[name],x),(r[name],y)):
                meta=d['arrays']['flux'+str(axis)]
                arrays.append(np.fromfile(file.parent/meta['file'],dtype='=f8').reshape(meta['shape'])[sl])
            result['values']+=arrays[0].size
            bits=[np.ascontiguousarray(value).view('=u8') for value in arrays]
            if not np.array_equal(*bits):
                result['failures'].append({'snapshot':name,'axis':axis,
                    'count':int(np.count_nonzero(bits[0]!=bits[1])),
                    'max_abs':float(np.max(np.abs(arrays[0]-arrays[1])))})
        result['snapshots']+=1
    return result

def launch(binary,inp,ranks,work,trace=False):
    work.mkdir(parents=True);(work/'tmp').mkdir()
    env=os.environ.copy();env['TMPDIR']=str(work/'tmp')
    for key in ('ATHENAK_CGL_LF_ARITHMETIC','ATHENAK_CGL_LF_DIAGNOSTICS',
                'ATHENAK_CGL_LF_STS_FLUX','ATHENAK_CGL_LF_P1_TRACE_DIR',
                'ATHENAK_CGL_LF_PROFILE','ATHENAK_CGL_LF_PROFILE_DETAIL',
                'ATHENAK_CGL_LF_TASK_TRACE','KOKKOS_TOOLS_LIBS','KOKKOS_PROFILE_LIBRARY'):
        env.pop(key,None)
    if trace:
        (work/'trace').mkdir()
        env.update(ATHENAK_CGL_LF_P1_TRACE_DIR=str(work/'trace'),
                   ATHENAK_CGL_LF_P1_MAX_CYCLE='0',ATHENAK_CGL_LF_P1_MAX_STAGE='1',
                   ATHENAK_CGL_LF_P1_POINTS='STSFluxes_end')
    cmd=['srun','--exact','-N1','-n',str(ranks),'--ntasks-per-node',str(ranks),
         '--threads-per-core=1','--cpu-bind=threads','-c7','--gpus-per-task=1',
         '--gpu-bind=closest',str(binary),'-i',str(inp)]
    start=time.perf_counter()
    proc=subprocess.run(cmd,cwd=work,env=env,capture_output=True,text=True,timeout=900)
    elapsed=time.perf_counter()-start
    stdout=proc.stdout+proc.stderr;(work/'stdout.log').write_text(stdout)
    row=dict(command=cmd,returncode=proc.returncode,launch_inclusive_seconds=elapsed,
             binary_sha256=sha(binary),input_sha256=sha(inp),cwd=str(work))
    if proc.returncode==0:
        row['cycles']=int(re.findall(r'^time=\S+ cycle=(\d+)',stdout,re.M)[-1])
        row['solver_seconds']=float(re.findall(r'^cpu time used\s*=\s*(\S+)',stdout,re.M)[-1])
        row['solver_seconds_per_cycle']=row['solver_seconds']/row['cycles']
        history=hst(str(next(work.glob('*.mhd.hst'))))
        row['lf_cell_stages']=float(history['lf_nstage'][-1])
        row['repair_counters']={k:float(history[k][-1]) for k in
            ('lf_dfloor','lf_pfloor','lf_nonfin','lf_nonpos','lf_hardbd','lf_hwproj')}
    (work/'command.json').write_text(json.dumps(row,indent=2)+'\n')
    return row

record={'manifest':str(a.manifest.resolve()),'manifest_sha256':sha(a.manifest),'source':str(FINAL_SOURCE),'phase':a.phase,'slurm_job':os.environ['SLURM_JOB_ID'],'cases':{},'failures':[],
        'timing':'Unprofiled driver timer starts after setup/initial outputs and includes evolution, later outputs and final diagnostics; launch time separately retained. Alternating paired order, one warmup per binary and at least three timed repeats.',
        'histogram_roundoff':'Only full-diagnostic q-work sums permit256epsilon relative reduction-order error. All field/restart state and other history columns exact.'}
def save(): (out/'results.json').write_text(json.dumps(record,indent=2)+'\n')

if a.phase=='task7-check':
    binaries={key:Path(final_manifest['trace_binaries'][key]['path']) for key in ('unfused','fused')}
    source=FINAL_SOURCE
    specs=[]
    for case in ('lf1d','lf2d','lf3d','low_b','smr_varied'):
        for arithmetic,diagnostics in (('safe','full'),('fast','none')):
            specs.append((case,BASE/'docs/validation/wo1/inputs'/(case+'.athinput'),
                          arithmetic,diagnostics,1,{}))
    for case in ('lf2d','lf3d'):
        specs.append((case,BASE/'docs/validation/wo1/inputs'/(case+'.athinput'),
                      'safe','none',1,{'mhd/lf_coefficient_mode':'local'}))
        specs.append((case,BASE/'docs/validation/wo1/inputs'/(case+'.athinput'),
                      'fast','full',4,{'mesh/x2max':'1.5','mesh/x3max':'0.75'}))
    for case in ('cgl_lf_staggered_checkerboard','cgl_lf_field_reversal_1d',
                 'cgl_lf_field_reversal_2d'):
        specs.append((case,source/'inputs/unit_tests'/(case+'.athinput'),'safe','full',1,{}))
    for option in a.extra_input:
        name,path=option.split('=',1)
        for arithmetic,diagnostics,ranks in (('safe','full',1),('fast','none',1),('safe','none',4)):
            specs.append((name,Path(path),arithmetic,diagnostics,ranks,{}))
    for case,source,arithmetic,diagnostics,ranks,extra in specs:
        if a.cases and case not in a.cases:continue
        name=f'{case}-{arithmetic}-{diagnostics}-{ranks}rank'
        parameters={'mhd/cgl_lf_profile':'false','mhd/cgl_lf_profile_detail':'false',
                    'mhd/cgl_lf_arithmetic':arithmetic,'mhd/cgl_lf_diagnostics':diagnostics,
                    'time/nlim':'20' if 'field_reversal' in case else '4',**extra}
        inp=input_file(name,source,parameters)
        runs={key:launch(binary,inp,ranks,out/name/key,True) for key,binary in binaries.items()}
        entry={'runs':runs}
        if all(r['returncode']==0 for r in runs.values()):
            entry['comparison']=compare(out/name/'unfused',out/name/'fused',True)
            entry['face_comparison']=compare_faces(out/name/'unfused',out/name/'fused')
            if entry['comparison']['failures'] or entry['face_comparison']['failures']:
                record['failures'].append(name)
        else:record['failures'].append(name)
        record['cases'][name]=entry;save();print(name,'failures',len(record['failures']),flush=True)
else:
    task7=a.phase=='task7-time'
    binaries={key:Path(final_manifest['binaries'][key]['path']) for key in ('unfused','fused')}
    cases=('lf2d','lf3d') if task7 else ('smr','smr_varied','smr_outflow')
    modes=(('safe','full'),('fast','none')) if task7 else (('safe','full'),)
    keys=list(binaries)
    for case in cases:
        if a.cases and case not in a.cases:continue
        for arithmetic,diagnostics in modes:
            for ranks in (1,4):
                name=f'{case}-{arithmetic}-{diagnostics}-{ranks}rank'
                parameters={'mhd/cgl_lf_profile':'false','mhd/cgl_lf_profile_detail':'false',
                    'mhd/cgl_lf_arithmetic':arithmetic,'mhd/cgl_lf_diagnostics':diagnostics,
                    'time/nlim':'40' if task7 else '12','time/tlim':'1000'}
                inp=input_file(name,BASE/'docs/validation/wo1/inputs'/(case+'.athinput'),parameters,True)
                entry={'runs':[]};record['cases'][name]=entry
                for repeat in range(-1,a.repeats):
                    order=keys if repeat%2 else keys[::-1]
                    for key in order:
                        work=out/name/key/('warmup' if repeat<0 else f'r{repeat}')
                        row=launch(binaries[key],inp,ranks,work)
                        row.update(binary=key,repeat=repeat);entry['runs'].append(row)
                        if row['returncode']:record['failures'].append(f'{name}/{key}/{repeat}')
                        save();print(name,key,repeat,row.get('solver_seconds'),flush=True)
                    if all(r['returncode']==0 for r in entry['runs'][-2:]):
                        tag='warmup' if repeat<0 else f'r{repeat}'
                        comparison=compare(out/name/keys[0]/tag,out/name/keys[1]/tag,task7)
                        entry.setdefault('comparisons',{})[tag]=comparison
                        if comparison['failures']:record['failures'].append(f'{name}/comparison/{tag}')
                    save()
                successful={key:[r for r in entry['runs'] if r['binary']==key and
                                r['repeat']>=0 and r['returncode']==0] for key in keys}
                valid=all(len(rows)==a.repeats for rows in successful.values())
                valid &= all(not c['failures'] for c in entry.get('comparisons',{}).values())
                if valid:
                    entry['medians']={key:statistics.median(r['solver_seconds_per_cycle']
                        for r in successful[key]) for key in keys}
                    entry['candidate_speedup']=entry['medians'][keys[0]]/entry['medians'][keys[1]]
                else:
                    entry['timing_accepted']=False
                    entry['timing_rejection']='Incomplete repeats or failed physical-state comparison'
                save()
print('Failures:',len(record['failures']),flush=True)
raise SystemExit(bool(record['failures']))
