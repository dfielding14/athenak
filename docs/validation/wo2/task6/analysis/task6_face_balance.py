#!/usr/bin/env python3
"""Match periodic block boundary fluxes and isolate same/coarse-fine residuals."""
import argparse
import json
from pathlib import Path
import numpy as np
from read_cgl_restart import read_restart

p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args()
rst=read_restart(sorted((a.run/'rst').glob('*.rst'))[0]);ng,nx,ny,nz=rst['indcs'][:4]
root_level=rst['root_level'];fine_level=int(rst['loc'][:,3].max())
n=32*2**(fine_level-root_level)
file=next(f for f in sorted((a.run/'trace').glob('*RecvFlux_end.json'))
          if json.loads(f.read_text())['sweep']=='pre')
d=json.loads(file.read_text());segments={}
for axis in (1,2):
    meta=d['arrays']['flux'+str(axis)]
    flux=np.fromfile(file.parent/meta['file'],dtype='=f8').reshape(meta['shape'])
    for m,gid in enumerate(d['gids']):
        loc=rst['loc'][gid];scale=2**(fine_level-int(loc[3]))
        for side in (0,1):
            fixed=(int(loc[axis-1])*nx+side*nx)*scale%n
            transverse_axis=2-axis
            for j in range(nx):
                transverse=(int(loc[transverse_axis])*nx+j)*scale
                value=(flux[m,4:6,0,ng+j,ng+side*nx] if axis==1 else
                       flux[m,4:6,0,ng+side*nx,ng+j])
                for sub in range(scale):
                    key=(axis,fixed,(transverse+sub)%n)
                    segments.setdefault(key,[]).append(dict(gid=gid,side=side,
                        level=int(loc[3]),scale=scale,flux=value.tolist()))
groups={}
for (axis,fixed,transverse),pair in segments.items():
    assert len(pair)==2 and pair[0]['side']!=pair[1]['side'],((axis,fixed,transverse),pair)
    kind='same_level' if pair[0]['level']==pair[1]['level'] else 'coarse_fine'
    scale=max(x['scale'] for x in pair)
    key=(kind,axis,fixed,transverse//scale,scale)
    group=groups.setdefault(key,dict(kind=kind,axis=axis,plane=fixed/n,
        transverse_cell=transverse//scale,scale=scale,gids=sorted(x['gid'] for x in pair),
        energy=0.,mu=0.,pairs=[]))
    for label,index in [('energy',0),('mu',1)]:
        group[label]+=sum((1 if x['side'] else -1)*x['flux'][index]/n for x in pair)
    group['pairs'].append(pair)
summary={kind:{label:sum(g[label] for g in groups.values() if g['kind']==kind)
               for label in ('energy','mu')} for kind in ('same_level','coarse_fine')}
result=dict(file=file.name,summary=summary,faces=sorted(groups.values(),key=lambda g:abs(g['energy']),reverse=True))
(a.run/'face-balance.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(summary,indent=2))
for g in result['faces'][:12]:print({k:v for k,v in g.items() if k!='pairs'})
