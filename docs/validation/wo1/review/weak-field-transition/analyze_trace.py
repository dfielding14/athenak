from pathlib import Path
import csv,math,json
base=Path('/tmp/cgl-b4-transition-20261001')
allcases={}
for name in ['plm','dc','reverse','floor6','all_magnetized']:
 rows=[]
 for r in csv.DictReader((base/name/'b4-state.csv').open()):
  rows.append({k:v if k=='phase' else float(v) for k,v in r.items()})
 index={(int(r['cycle']),int(r['stage']),r['phase'],int(r['cell'])):r for r in rows}
 summaries=[]
 print('\nCASE',name)
 for cycle in range(5):
  cp=[r for r in rows if r['cycle']==cycle and r['phase']=='coll_post']
  pre=[r for r in rows if r['cycle']==cycle and r['phase']=='coll_pre']
  best=max(cp,key=lambda r:r['pperp']/r['ppar'])
  cell=64 if name!='reverse' else 63
  t=next(r for r in cp if r['cell']==cell)
  walls=[];c2p=[]
  for r in cp:
   old=index[(cycle,int(r['stage']),'coll_pre',int(r['cell']))]
   if r['A']!=old['A']:
    walls.append({'cell':int(r['cell']),'B':r['By'],'dA':r['A']-old['A'],'ratio_before':old['pperp']/old['ppar'],'ratio_after':r['pperp']/r['ppar']})
  for stage in (1,2):
   for i in range(128):
    old=index[(cycle,stage,'c2p_pre',i)];new=index[(cycle,stage,'c2p_post',i)]
    if old['A']!=new['A']:
     c2p.append({'cell':i,'stage':stage,'B':new['By'],'dA':new['A']-old['A']})
  result={'cycle':cycle+1,'max_ratio':best['pperp']/best['ppar'],'max_cell':int(best['cell']),'B_at_max':best['By'],'target_ratio':t['pperp']/t['ppar'],'target_A':t['A'],'target_B':t['By'],'target_rho':t['rho'],'walls':walls,'c2p':c2p}
  summaries.append(result)
  print('cycle',cycle+1,'maxr',result['max_ratio'],'cell',result['max_cell'],'B',result['B_at_max'],'walls',walls,'c2p',c2p)
 allcases[name]=summaries
(base/'trace-summary.json').write_text(json.dumps(allcases,indent=2)+'\n')
print('\nTARGET CELL 64: PLM')
with (base/'plm/b4-state.csv').open() as f:
 for r in csv.DictReader(f):
  if int(r['cell'])!=64:continue
  rho=float(r['rho']);a=float(r['A']);b=abs(float(r['By']))
  q=a/rho-2*math.log(rho)+3*math.log(b)
  print(r['cycle'],r['stage'],r['phase'],'rho',rho,'A',a,'B',b,'q',q,'r_calc',math.exp(q),'FA',r['FAL'],r['FAR'],'FB',r['FByL'],r['FByR'])
