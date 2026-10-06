from pathlib import Path
import argparse,difflib,hashlib,json,os,shlex,shutil,subprocess
here=Path(__file__).resolve().parent;final=here.parent
ap=argparse.ArgumentParser();ap.add_argument('backend',choices=['cpu','hip']);args=ap.parse_args()
backend=args.backend;out=here/backend
old=(final/'source-fused/src/eos/cgl_passive.hpp').read_text()
new=Path('/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2/src/eos/cgl_passive.hpp').read_text()
expected=old.replace('return {log(ppar) + 2.0*logb - 3.0*logr,\n          log(pperp) - log(ppar) + 2.0*logr - 3.0*logb};',
    'return {static_cast<Real>(log(ppar) + 2.0*logb - 3.0*logr),\n          static_cast<Real>(log(pperp) - log(ppar) + 2.0*logr - 3.0*logb)};')
assert new==expected and new!=old, 'Proof is limited to exactly the two cast wrappers'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
flags_file=final/f'build-{backend}/src/CMakeFiles/athena.dir/flags.make'
variables={}
for line in flags_file.read_text().splitlines():
 if line.startswith('CXX_') and ' = ' in line:
  key,value=line.split(' = ',1);variables[key]=shlex.split(value)
flags=[*variables['CXX_DEFINES'],'-I'+str(out/'include'),*variables['CXX_INCLUDES'],*variables['CXX_FLAGS']]
compiler=shutil.which('CC'); assert compiler
report={'backend':backend,'scope':'Real=double, exactly two complete-expression static_cast<Real> wrappers',
        'compiler':compiler,'compiler_version':subprocess.check_output([compiler,'--version'],text=True),
        'final_flags_file':str(flags_file),'final_flags_sha256':sha(flags_file),
        'probe_sha256':sha(here/'probe.cpp'),'commands':[],'objects':{},'comparisons':{}}
for label,text in [('old',old),('new',new)]:
 (out/'include/eos/cgl_passive.hpp').write_text(text)
 (out/(label+'-cgl_passive.hpp')).write_text(text)
 for product,options in [('object',['-c']),('ir',['-S','-emit-llvm']+(['--cuda-device-only'] if backend=='hip' else []))]:
  target=out/f'{label}.{ "o" if product=="object" else "ll"}'
  command=[compiler,*flags,*options,str(here/'probe.cpp'),'-o',str(target)]
  report['commands'].append({'header':label,'product':product,'argv':command})
  proc=subprocess.run(command,cwd=out,env=os.environ.copy(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  (out/f'{label}-{product}.log').write_text(proc.stdout)
  report['commands'][-1]['returncode']=proc.returncode
  (out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
  if proc.returncode:raise RuntimeError(proc.stdout)
  report['objects'][target.name]={'sha256':sha(target),'bytes':target.stat().st_size}
 for attr in ['old','new']:
  if (out/(attr+'-cgl_passive.hpp')).exists():report[attr+'_header_sha256']=sha(out/(attr+'-cgl_passive.hpp'))
for suffix in ['o','ll']:
 a=(out/f'old.{suffix}').read_bytes();b=(out/f'new.{suffix}').read_bytes()
 report['comparisons'][suffix]={'raw_byte_equal':a==b,'old_sha256':sha(out/f'old.{suffix}'),'new_sha256':sha(out/f'new.{suffix}')}
 if suffix=='ll' and a!=b:
  (out/'ir.diff').write_text(''.join(difflib.unified_diff(a.decode().splitlines(True),b.decode().splitlines(True))))
report['passed']=all(x['raw_byte_equal'] for x in report['comparisons'].values())
(out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'backend':backend,'comparisons':report['comparisons'],'passed':report['passed']},indent=2))
