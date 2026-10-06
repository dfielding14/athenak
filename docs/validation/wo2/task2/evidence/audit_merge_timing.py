"""Read-only repeatability and health audit of the completed merge timings."""
from pathlib import Path
import hashlib,json,statistics,sys
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here.parent/'p1-research'))
from restart_compare import normalized_restart
root=here/'timing/task1-refill-hip'
original=(root/'results.json').read_bytes();data=json.loads(original)
assert len(data['results'])==16 and not data['failures'] and 'medians' in data
report={'source_results_sha256':hashlib.sha256(original).hexdigest(),
        'normalization':'Only the documented36 unused root RegionIndcs coarse-index bytes; no forcing or other normalization.',
        'history_health':[],'repeatability':[],'ranges':{},'failures':[]}
for sample in data['results']:
 label=f"ratio{sample['ratio']}-repeat{sample['repeat']}-merge{sample['enabled']}"
 counters={key:sample['history_final'][key] for key in
           ('lf_dfloor','lf_pfloor','lf_nonfin','lf_nonpos','lf_hardbd','lf_hwproj')}
 report['history_health'].append({'run':label,'counters':counters})
 if any(value!=0 for value in counters.values()):report['failures'].append({'run':label,'health':counters})
for ratio in (2,20):
 for enabled in (False,True):
  samples=[x for x in data['results'] if x['ratio']==ratio and x['enabled']==enabled and x['timed']]
  assert len(samples)==3
  timing=[x['solver_seconds_per_cycle'] for x in samples]
  report['ranges'][f'ratio{ratio}-merge{enabled}']={'median':statistics.median(timing),'min':min(timing),'max':max(timing),'samples':timing}
  reference=next(x for x in samples if x['repeat']==0)
  refdir=root/f'ratio{ratio}-repeat0-merge{enabled}'
  for sample in samples:
   if sample['repeat']==0:continue
   directory=root/f"ratio{ratio}-repeat{sample['repeat']}-merge{enabled}"
   comparison={'reference':str(refdir),'run':str(directory),'files':{},'physical_equal':True}
   assert reference['output_hashes'].keys()==sample['output_hashes'].keys()
   for name in reference['output_hashes']:
    a,b=refdir/name,directory/name
    raw_a,raw_b=a.read_bytes(),b.read_bytes()
    assert hashlib.sha256(raw_a).hexdigest()==reference['output_hashes'][name]
    assert hashlib.sha256(raw_b).hexdigest()==sample['output_hashes'][name]
    record={'raw_equal':raw_a==raw_b,'reference_raw_sha256':reference['output_hashes'][name],
            'raw_sha256':sample['output_hashes'][name]}
    if name.endswith('.rst'):
     norm_a,meta_a=normalized_restart(a);norm_b,meta_b=normalized_restart(b)
     assert meta_a['normalized_byte_count']==meta_b['normalized_byte_count']==36
     equal=norm_a==norm_b
     record.update(reference_metadata=meta_a,metadata=meta_b,
                   raw_difference_offsets=[i for i,(x,y) in enumerate(zip(raw_a,raw_b)) if x!=y],
                   raw_lengths=[len(raw_a),len(raw_b)])
    else:equal=raw_a==raw_b
    record['physical_equal']=equal;comparison['files'][name]=record
    if not equal:comparison['physical_equal']=False;report['failures'].append({'run':str(directory),'file':name,'reason':'physical repeatability'})
   report['repeatability'].append(comparison)
assert (root/'results.json').read_bytes()==original
(here/'timing-repeatability.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'physical_comparisons':len(report['repeatability']),'all_health_zero':all(not any(x['counters'].values()) for x in report['history_health']),
                  'ranges':report['ranges'],'failures':report['failures']},indent=2))
assert not report['failures']
