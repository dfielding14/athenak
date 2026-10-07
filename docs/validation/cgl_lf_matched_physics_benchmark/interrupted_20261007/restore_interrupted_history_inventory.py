#!/usr/bin/env python3
"""Add existing histories to canceled-launch analysis; never edit raw provenance.
Run only after both original analysis commands have finished. Retains pre-amendment
metrics and uses repository readers/diagnostic functions without recomputing FFTs.
"""
import copy,datetime,hashlib,json,pathlib,shutil,subprocess,sys,time
S=pathlib.Path('/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2');sys.path.insert(0,str(S/'scripts'))
import analyze_cgl_lf_physics_benchmark as single
import compare_cgl_lf_physics_benchmark as pair
M=pathlib.Path('/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched');O=M/'interrupted-analysis';run=M/'active-to14'
for mode in ['paired_startup','active_interrupted']:
 rec=json.loads((O/(mode+'.execution.json')).read_text())
 if rec.get('returncode') != 0:raise RuntimeError('Original analysis unfinished/failed: '+mode)
source=run/'benchmark_metadata.json';metadata=json.loads(source.read_text());source_record=single.retained_file(source)
assert metadata['outputs']['user_history']==[] and metadata['outputs']['mhd_history']==[]
metadata['outputs']['user_history']=[p.name for p in sorted(run.glob('*.user.hst'))]
metadata['outputs']['mhd_history']=[p.name for p in sorted(run.glob('*.mhd.hst'))]
assert len(metadata['outputs']['user_history'])==len(metadata['outputs']['mhd_history'])==1
metadata['analysis_inventory_amendment']={'purpose':'Read existing histories omitted from the canceled launcher initial metadata; raw launch metadata is unchanged.','original_metadata':source_record,'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'fields_changed':['outputs.user_history','outputs.mhd_history'],'termination_note':'User requested stopping the allocation. The original launcher did not finalize and its missing ended_utc/returncode remain missing; no successful completion is claimed.'}
metadata_path=O/'active-analysis-metadata.json';metadata_path.write_text(json.dumps(metadata,indent=2)+'\n')
segments=single.load_segments(run,metadata,metadata_path)
user_paths,ub,_=single.segment_paths(segments,'user_history','**/*.user.hst');mhd_paths,mb,_=single.segment_paths(segments,'mhd_history','**/*.mhd.hst')
user,ui=single.merge_histories(user_paths,ub);mhd,mi=single.merge_histories(mhd_paths,mb)
record={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'method':'Existing snapshot/PDF/FFT/force products retained exactly; only history discovery, history products, numerical-health evidence, dependent adequacy reasons, provenance and rendered figures/reports refreshed with existing analyzer functions.','script':single.retained_file(pathlib.Path(__file__)),'source_metadata':source_record,'analysis_metadata':single.retained_file(metadata_path),'single_run_source':single.retained_file(pathlib.Path(single.__file__)),'comparison_source':single.retained_file(pathlib.Path(pair.__file__)),'analysis_revision':subprocess.check_output(['git','-C',str(S),'rev-parse','HEAD'],text=True).strip(),'outputs':[]}
startclock=time.monotonic()
for relative in ['paired-startup/active','active-interrupted']:
 dest=O/relative;path=dest/'metrics.json';before=single.retained_file(path);data=json.loads(path.read_text())
 backup=dest/'metrics-before-history-inventory-amendment.json';shutil.copyfile(path,backup)
 old_reasons=single.energy_coverage_reasons(data['history'],data['model']);start,end=data['requested_window'];block=data['sampling']['block_duration']
 data['history']=single.json_value(single.history_products(user,mhd,data['model'],start,end,block))
 data['simulation_integrity']=single.solver_health(segments,user,mhd,start,end)
 data['adequacy']['reasons']=[r for r in data['adequacy']['reasons'] if r not in old_reasons]+single.energy_coverage_reasons(data['history'],data['model'])
 data['adequacy']['sampling_prerequisites_satisfied']=not data['adequacy']['reasons']
 data['provenance']['history_inventory_amendment']={'original_metrics':single.retained_file(backup),'amendment':metadata['analysis_inventory_amendment'],'amendment_script':record['script']}
 data['provenance']['metadata']=single.retained_file(metadata_path);data['provenance']['retained_metadata']=metadata
 data['provenance']['segments']=[{k:str(v) if isinstance(v,pathlib.Path) else v for k,v in seg.items()} for seg in segments]
 data['provenance']['user_history']=ui;data['provenance']['mhd_history']=mi
 old=json.loads(backup.read_text())
 unchanged_hashes={}
 for key in ['snapshots','scalars','spectra','PDFs','pressure_balance_by_scale','forcing_decomposition']:
  assert data[key]==old[key],key
  old_hash=hashlib.sha256(json.dumps(old[key],sort_keys=True,separators=(',',':')).encode()).hexdigest()
  new_hash=hashlib.sha256(json.dumps(data[key],sort_keys=True,separators=(',',':')).encode()).hexdigest()
  assert old_hash==new_hash,key
  unchanged_hashes[key]={'old_sha256':old_hash,'new_sha256':new_hash}
 path.write_text(json.dumps(single.json_value(data),indent=2,allow_nan=False)+'\n')
 single.make_figures(data,dest);single.write_report(data,dest)
 record['outputs'].append({'directory':str(dest),'old_metrics':before,'new_metrics':single.retained_file(path),'unchanged_products':unchanged_hashes})
 print('Updated histories and figures:',dest,flush=True)
# Main paired figures/reports consume the corrected cached single-run products.
dest=O/'paired-startup';path=dest/'metrics.json';backup=dest/'metrics-before-history-inventory-amendment.json';shutil.copyfile(path,backup)
result=json.loads(path.read_text());data={mode:json.loads((dest/mode/'metrics.json').read_text()) for mode in pair.MODES}
failed=[mode for mode in pair.MODES if any(r.get('returncode') not in (0,None) for r in data[mode]['simulation_integrity']['segments'])]
health_classes=[data[mode]['simulation_integrity']['classification'] for mode in pair.MODES]
integrity={'classification':'concerning' if failed or 'concerning' in health_classes else 'consistent' if all(c=='consistent' for c in health_classes) else 'inconclusive','failed_runs':failed,'scope':'simulation completion and retained numerical health; separate from parameter matching and descriptive contrast'}
result['simulation_integrity']=integrity;result['physical_evidence']=pair.physical_evidence(data,integrity);result['figure_status']=pair.figure_status(data)
for mode in pair.MODES:
 result['runs'][mode]['metrics']=single.retained_file(dest/mode/'metrics.json');result['runs'][mode]['simulation_integrity']=data[mode]['simulation_integrity']
result['provenance']['history_inventory_amendment']={'original_metrics':single.retained_file(backup),'analysis_metadata':single.retained_file(metadata_path),'script':record['script'],'scope':record['method']}
path.write_text(json.dumps(single.json_value(result),indent=2,allow_nan=False)+'\n');audit=pair.make_figures(data,dest);(dest/'figure-audit.json').write_text(json.dumps(audit,indent=2)+'\n');pair.write_report(result,data,dest)
assert single.retained_file(source)==source_record,'Raw metadata changed unexpectedly'
record.update(ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),wall_seconds=time.monotonic()-startclock,raw_metadata_unchanged=True)
(O/'history-inventory-amendment.json').write_text(json.dumps(record,indent=2)+'\n')
print('Completed metadata-preserving history augmentation',flush=True)
