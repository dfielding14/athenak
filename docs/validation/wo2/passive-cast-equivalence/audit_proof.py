from pathlib import Path
import hashlib,json,re,subprocess
root=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
report={'scope':'Real=double at final CPU/HIP flags; only two static_cast<Real> wrappers; no application execution','checks':{},'normalization':{}}
a=(root/'cpu/old.o').read_bytes();b=(root/'cpu/new.o').read_bytes()
report['checks']['cpu_complete_object']={'equal':a==b,'bytes':len(a),'sha256':sha(root/'cpu/old.o')}
for backend in ['cpu','hip']:
 texts=[(root/backend/(label+'.ll')).read_text() for label in ['old','new']]
 if backend=='cpu':
  def norm(text):
   text=re.sub(r'(; include/eos/cgl_passive.hpp:(?:19|20):)\d+',r'\g<1>COLUMN',text)
   return re.sub(r'(!DILocation\(line: (?:19|20), column: )\d+',r'\g<1>COLUMN',text)
  report['normalization']['cpu_ir']='Only source-column numbers on changed header lines19/20, in comments and DILocation metadata; full machine object already raw-identical.'
 else:
  norm=lambda text:re.sub(r'__hip_cuid_[0-9a-f]+','__hip_cuid_COMPILER_ID',text)
  report['normalization']['hip_ir']='Only compiler-generated __hip_cuid_* identifier names (two occurrences per IR); no instructions or arithmetic attributes removed.'
  assert all(len(re.findall(r'__hip_cuid_[0-9a-f]+',t))==2 for t in texts)
 normalized=[norm(t) for t in texts]
 for label,text in zip(['old','new'],normalized):(root/backend/(label+'-metadata-normalized.ll')).write_text(text)
 report['checks'][backend+'_ir']={'equal':normalized[0]==normalized[1],'raw_equal':texts[0]==texts[1]}
 for function in ['wo2_passive_specific','wo2_passive_encode']+(['wo2_passive_kernel'] if backend=='hip' else []):
  bodies=[]
  for text in normalized:
   match=re.search(r'^define[^\n]*@'+function+r'\([^\n]*\{\n.*?^\}',text,re.M|re.S);assert match,function
   bodies.append(match.group())
  report['checks'][backend+'_'+function]={'equal':bodies[0]==bodies[1],'sha256':hashlib.sha256(bodies[0].encode()).hexdigest()}
for section in ['text','rodata']:
 a=root/f'hip/old-device.{section}';b=root/f'hip/new-device.{section}'
 report['checks']['hip_device_'+section]={'equal':a.read_bytes()==b.read_bytes(),'bytes':a.stat().st_size,'sha256':sha(a)}
report['passed']=all(x['equal'] for x in report['checks'].values());assert report['passed']
(root/'equivalence-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
