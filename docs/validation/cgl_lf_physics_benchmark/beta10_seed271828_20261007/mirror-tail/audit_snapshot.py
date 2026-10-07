import sys,json,hashlib,re
from pathlib import Path
import numpy as np
R=Path('/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2')
sys.path.insert(0,str(R/'scripts'))
import analyze_cgl_lf_physics_benchmark as a
B=Path('/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark');D=B/'analysis-checks/mirror-tail-t14p25'
p=B/'canonical-fixed-to18/bin/cgl_lf_physics_benchmark_beta10.mhd_w_bcc.00057.bin'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(p)=='41232a8f3d2000935522d9e1cdfa6f0a0dddc729823ba88ac480428268819e17'
f,lengths,t=a.paper.read_snapshot(p)
rho,pp,pq=[f[k] for k in ('dens','eint','p_perp')];b=[f[k] for k in ('bcc1','bcc2','bcc3')]
b2=sum(v*v for v in b);x=2*(pq-pp)/b2;en=a.float32_x_envelope(pp,pq,b);inds=np.argwhere(x-en>1)
rows=[]
for ix in inds:
 ix=tuple(int(i) for i in ix);X=float(x[ix]);e=float(en[ix]);rows.append({'global_kji':ix,'rho':float(rho[ix]),'p_parallel':float(pp[ix]),'p_perp':float(pq[ix]),'B':[float(v[ix]) for v in b],'B2':float(b2[ix]),'X':X,'X_minus_1':X-1,'float32_X_rounding_envelope':e,'excess_over_envelope':X-1-e,'excess_to_envelope_ratio':(X-1)/e,'delta_p_minus_B2_over_2':float(pq[ix]-pp[ix]-.5*b2[ix]),'beta_scalar':float(2*(pp[ix]+2*pq[ix])/3/b2[ix])})
m={'snapshot':str(p),'sha256':sha(p),'time':t,'shape_zyx':list(x.shape),'lengths_xyz':list(lengths),'cell_count':int(x.size),'strict_mirror_count':int(np.count_nonzero(x>1)),'rounding_robust_mirror_count':len(rows),'rounding_robust_firehose_count':int(np.count_nonzero(x+en < -2)),'robust_volume_fraction':len(rows)/x.size,'max_X':float(x.max()),'cells':rows,'provenance':{str(q):sha(q) for q in [Path(__file__),R/'scripts/analyze_cgl_lf_physics_benchmark.py',R/'scripts/analyze_cgl_lf_paper.py',R/'src/eos/ideal_c2p_mhd.hpp',R/'src/mhd/mhd_sts.cpp',R/'src/driver/driver.cpp',B/'canonical-fixed-to18/benchmark_metadata.json']},'runtime':{'job':5629677,'solver_run':False,'numpy':np.__version__}}
(D/'cells.json').write_text(json.dumps(m,indent=2)+'\n');print(json.dumps(m,indent=2))
