"""Independent passive invariants, flow identity, heating, linear and restart tests.

The restart decoder checks the Linux/double ABI and keeps the complete raw state.
Every invocation writes only to its explicit output directory.
"""
from pathlib import Path
import argparse, cmath, hashlib, json, math, os, re, shlex, struct, subprocess, time
import numpy as np
HERE=Path(__file__).resolve().parent

def input_text(*,passive=True,nx=128,mode=0,tlim=.05,nlim=-1,lf=None,forcing=False,
               recon='plm',amp=.1,thermal=.2,cfl=.3,nu=0.0,dfloor=1e-12,
               pfloor=1e-12,weak_field=False,acceleration=False,forcing_policy=None,strict=True):
    dims=(16,16,16) if forcing else (nx,1,1)
    blocks=(8,8,8) if forcing else (min(nx//4,64),1,1)
    d={'job':{'basename':'passive_validation'},'mesh':{'nghost':4},'meshblock':{},
       'time':{'evolution':'dynamic','integrator':'rk2','sts_integrator':'rkl2' if passive and lf=='sts' else 'none',
               'sts_max_dt_ratio':-1,'cfl_number':cfl,'nlim':nlim,'tlim':tlim,'ndiag':1},
       'mhd':{'eos':'cgl' if passive else 'isothermal','passive':str(passive).lower(),
              'iso_sound_speed':1.0,'reconstruct':recon,'rsolver':'hlle','fofc':'true',
              'dfloor':dfloor,'pfloor':pfloor,'bfloor':1e-10,'nu_coll':nu,
              'mirror_limiter':'false','firehose_limiter':'false','backup_limiters':'false'},
       'problem':{'pgen_name':'cgl_passive_validation','case':mode,'amplitude':amp,'thermal_amplitude':thermal,'p0':2.0},
       'output1':{'file_type':'hst','data_format':'%25.17e','dcycle':1},
       'output2':{'file_type':'rst','dcycle':100000}}
    for n,(total,block) in enumerate(zip(dims,blocks),1):
        d['mesh'].update({f'nx{n}':total,f'x{n}min':0,f'x{n}max':1,
                          f'ix{n}_bc':'periodic',f'ox{n}_bc':'periodic'})
        d['meshblock'][f'nx{n}']=block
    if passive and lf:
        d['mhd'].update(cgl_heat_flux='landau_fluid',cgl_heat_flux_integrator=lf,
                        lf_k_parallel=2*math.pi,lf_coefficient_mode='background',
                        lf_c_parallel0=math.sqrt(2),cgl_lf_strict_admissibility=str(strict).lower())
    if weak_field:
        d['problem'].update(bx=1e-12,by=0,bz=0)
    if acceleration:
        d['mhd_srcterms']={'const_accel':'true','const_accel_dir':2,'const_accel_val':.5}
    if forcing:
        d['turb_driving']={'driving_type':1,'tcorr':1,'dedt':.02,'nlow':1,'nhigh':2,
                           'exp_prp':5/3,'exp_prl':0,'rseed':519,'record_injected_work':'true'}
        if forcing_policy:d['turb_driving']['projection_policy']=forcing_policy
    return ''.join('\n<'+block+'>\n'+''.join(f'{k} = {v}\n' for k,v in values.items())
                   for block,values in d.items())

def decode_restart(path,passive):
    raw=Path(path).read_bytes();marker=b'<par_end>\n';h=raw.index(marker)+len(marker)
    nmb,root=struct.unpack_from('=ii',raw,h)
    ng,n1,n2,n3=struct.unpack_from('=4i',raw,h+8+72+76)
    assert nmb>0 and root>=0 and ng==4 and min(n1,n2,n3)>0,'unexpected restart ABI'
    nx,ny,nz=n1+2*ng,n2+2*ng if n2>1 else 1,n3+2*ng if n3>1 else 1
    nv=6 if passive else 4;cc=nx*ny*nz
    sizes=[nv*cc,(nx+1)*ny*nz,nx*(ny+1)*nz,nx*ny*(nz+1)]
    count=sum(sizes);off=len(raw)-nmb*count*8
    assert off>h+252 and struct.unpack_from('=Q',raw,off-8)[0]==count*8
    data=np.frombuffer(raw,dtype='=f8',offset=off).reshape(nmb,count)
    assert np.isfinite(data).all(),'nonfinite full-precision restart state'
    u=data[:,:sizes[0]].reshape(nmb,nv,nz,ny,nx)
    f=[];start=sizes[0]
    for count,shape in zip(sizes[1:],[(nz,ny,nx+1),(nz,ny+1,nx),(nz+1,ny,nx)]):
        f.append(data[:,start:start+count].reshape(nmb,*shape));start+=count
    slices=(slice(None),slice(ng,ng+n3) if n3>1 else slice(None),
            slice(ng,ng+n2) if n2>1 else slice(None),slice(ng,ng+n1))
    uactive=u[(slice(None),slice(None),*slices[1:])]
    bcc=np.stack([.5*(f[0][...,:-1]+f[0][...,1:]),.5*(f[1][:,:,:-1,:]+f[1][:,:,1:,:]),
                  .5*(f[2][:,:-1,:,:]+f[2][:,1:,:,:])],axis=1)
    bactive=bcc[(slice(None),slice(None),*slices[1:])]
    t,dt,cycle=struct.unpack_from('=ddi',raw,h+8+72+2*76)
    ret={'u':u,'active':uactive,'faces':f,'bcc':bactive,'time':t,'dt':dt,'cycle':cycle,
         'sha256':hashlib.sha256(raw).hexdigest(),'payload':data}
    if passive:
        rho=uactive[:,0];bm=np.maximum(np.sqrt(np.sum(bactive*bactive,axis=1)),1e-10)
        j=uactive[:,4]/rho;a=uactive[:,5]/rho
        pp=np.exp(j+3*np.log(rho)-2*np.log(bm));pt=np.exp(a+j+np.log(rho)+np.log(bm))
        ret.update(pp=pp,pt=pt,U=.5*pp+pt)
        assert np.all(pp>0) and np.all(pt>0)
    return ret

def history(path):
    lines=Path(path).read_text().splitlines()
    labels=re.findall(r'\[\d+\]=([^\s]+)',next(x for x in lines if '[1]=' in x))
    data=np.loadtxt(path,ndmin=2)
    assert len(labels)==data.shape[1],(labels,data.shape)
    return {key:data[:,n] for n,key in enumerate(labels)}

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--launcher',default='')
    p.add_argument('--suite',choices=['all','identity','heating','linear','advection','restart','fences'],default='all')
    args=p.parse_args(argv)
    args.execute=True
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=True)
    if list(root.glob('*/stdout.log')):
        raise RuntimeError('Output directory already contains retained solver results')
    records=[];states={};checks=[]
    def save():
        temporary=root/'results.json.tmp'
        temporary.write_text(json.dumps({'execute':args.execute,'binary':str(args.binary),
           'binary_sha256':hashlib.sha256(args.binary.read_bytes()).hexdigest() if args.binary.exists() else None,
           'records':records,'checks':checks},indent=2)+'\n')
        temporary.replace(root/'results.json')
    def run(name,passive=True,restart=None,overrides=(),expect_failure=None,binary=None,**kw):
        out=root/name;out.mkdir(exist_ok=True);(out/'tmp').mkdir(exist_ok=True)
        inp=out/'input.athinput';input_data=input_text(passive=passive,**kw)
        if not restart:
            for option in overrides:
                block,keyvalue=option.split('/',1);key=keyvalue.split('=',1)[0]
                marker='<'+block+'>\n'
                if marker not in input_data:
                    input_data+='\n'+marker+keyvalue+'\n'
                else:
                    segment=input_data.split(marker,1)[1].split('\n<',1)[0]
                    if not re.search(r'^'+re.escape(key)+r'\s*=',segment,re.M):
                        input_data=input_data.replace(marker,marker+keyvalue+'\n',1)
        inp.write_text(input_data)
        command=shlex.split(args.launcher)
        selected_binary=(binary or args.binary).resolve()
        command += [str(selected_binary),'-r' if restart else '-i',str(restart or inp),*overrides]
        record={'name':name,'passive':passive,'binary':str(selected_binary),
                'binary_sha256':hashlib.sha256(selected_binary.read_bytes()).hexdigest() if selected_binary.exists() else None,'configuration':kw,'restart_input':str(restart) if restart else None,'expected_failure':expect_failure,'command':command,'cwd':str(out),
                'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest()};records.append(record);save()
        if not args.execute:return None
        print('START '+name,flush=True)
        env=os.environ.copy();env['TMPDIR']=str(out/'tmp');start=time.monotonic()
        with (out/'stdout.log').open('w') as log:
            proc=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        record.update(exit_code=proc.returncode,seconds=time.monotonic()-start);save()
        if expect_failure:
            text=(out/'stdout.log').read_text()
            passed=proc.returncode!=0 and expect_failure in text
            checks.append({'name':name,'passed':passed,'expected_error':expect_failure});save()
            print(('PASS ' if passed else 'FAIL ')+name+' expected rejection',flush=True)
            assert passed,(name,proc.returncode,text[-1000:])
            return None
        assert proc.returncode==0,(name,proc.returncode)
        rst=sorted((out/'rst').glob('*.rst'))[-1];state=decode_restart(rst,passive)
        hst=next(out.glob('*.hst'));state['hist']=history(hst);states[name]=state
        record.update(restart=str(rst),restart_sha256=state['sha256'],time=state['time'],cycle=state['cycle'])
        state['restart']=rst
        save();print(f'DONE {name} cycle={state["cycle"]} time={state["time"]:.17g}',flush=True)
        return state
    def check(name,passed,**evidence):
        checks.append({'name':name,'passed':bool(passed),**evidence});save()
        print(('PASS ' if passed else 'FAIL ')+name+' '+json.dumps(evidence),flush=True)
        assert passed,(name,evidence)
    if args.suite in ['all','identity']:
        for label,kw in [('plm',{}),('dc',{'recon':'dc'}),('ppm4',{'recon':'ppm4'}),('ppmx',{'recon':'ppmx'}),('wenoz',{'recon':'wenoz'}),
                         ('density-floor',{'amp':.95,'dfloor':.2}),
                         ('thermal-floor',{'pfloor':3.0,'strict':False}),('weak-field',{'weak_field':True}),
                         ('accelerated',{'acceleration':True}),('forced',{'forcing':True,'nlim':12,'tlim':1}),
                         ('forced-3d',{'forcing':True,'nlim':12,'tlim':1,
                                       'forcing_policy':'mks24_alfvenic_perpendicular'})]:
            ref=run('iso-'+label,False,**kw)
            variants=['none','sts']
            for method in variants:
                trial=run('passive-'+label+'-'+method,True,lf=None if method=='none' else method,
                          nu=.3 if method!='none' else 0,**kw)
                if args.execute:
                    pairs=[(ref['u'][:,:4],trial['u'][:,:4]),*zip(ref['faces'],trial['faces'])]
                    mismatch=sum(int(np.count_nonzero(x.view('=u8')!=y.view('=u8'))) for x,y in pairs)
                    time_equal=np.array_equal(ref['hist']['time'],trial['hist']['time'])
                    check(label+'-'+method,mismatch==0 and time_equal and ref['dt']==trial['dt'],
                          flow_bit_mismatches=mismatch,time_sequence_equal=time_equal,
                          next_dt_equal=ref['dt']==trial['dt'])
        # Both budgets scale with cfl_number. Use weak but nonzero LF so the
        # explicit stability bound is above the unchanged advective bound.
        ref=run('iso-explicit-matched',False)
        trial=run('passive-explicit-matched',True,lf='explicit',nu=.3,
                  overrides=['mhd/lf_k_parallel=100000'])
        if args.execute:
            pairs=[(ref['u'][:,:4],trial['u'][:,:4]),*zip(ref['faces'],trial['faces'])]
            mismatch=sum(int(np.count_nonzero(x.view('=u8')!=y.view('=u8'))) for x,y in pairs)
            time_equal=np.array_equal(ref['hist']['time'],trial['hist']['time'])
            check('explicit-matched-dt',mismatch==0 and time_equal and ref['dt']==trial['dt'],
                  flow_bit_mismatches=mismatch,time_sequence_equal=time_equal)
    if args.suite in ['all','heating']:
        oracle=.2*.1*2*math.pi/4;errors=[]
        for nx in [128,256,512]:
            h=.004*128/nx
            a=run(f'heat-{nx}-h',nx=nx,mode=1,amp=.1,thermal=.2,tlim=h)
            b=run(f'heat-{nx}-2h',nx=nx,mode=1,amp=.1,thermal=.2,tlim=2*h)
            if args.execute:
                initial=float(a['hist']['thermal-U'][0]);ua=float(a['U'].mean());ub=float(b['U'].mean())
                measured=(-3*initial+4*ua-ub)/(2*h);error=abs(measured/oracle-1);errors.append(error)
                check(f'periodic-heating-{nx}',error<.01,oracle=oracle,measured=measured,relative_error=error,
                      initial_U=initial,final_U=ua)
        if args.execute:
            order=math.log(errors[0]/errors[-1])/math.log(4)
            check('periodic-heating-refinement',order>1.7,order=order,errors=errors)
    if args.suite in ['all','linear']:
        amp=1e-6;t=.04;k=2*math.pi;omega=k;p0=2;c=math.sqrt(p0)
        for method in [None,'sts','explicit']:
            errors=[]
            for nx in ([128,256,512] if method!='explicit' else [256]):
                out=run(f'linear-{method or "none"}-{nx}',nx=nx,mode=2,amp=amp,thermal=amp,tlim=t,lf=method)
                if args.execute:
                    xx=(np.arange(nx)+.5)/nx;phase=np.exp(-1j*k*xx)
                    measured=[complex(2*np.mean((out[key].reshape(-1)-p0)*phase)) for key in ['pp','pt']]
                    expected=[]
                    for multiplier,diff in [(3,math.sqrt(8/math.pi)*c/k),(1,math.sqrt(2/math.pi)*c/k)]:
                        damping=diff*k*k if method else 0
                        particular=amp*p0*(damping-1j*multiplier*omega)/(damping-1j*omega)
                        expected.append(particular*cmath.exp(-1j*omega*t)+
                                        (multiplier*p0*amp-particular)*math.exp(-damping*t))
                    error=max(abs(m-e)/(p0*amp) for m,e in zip(measured,expected));errors.append(error)
                    check(f'linear-response-{method}-{nx}',error<.02,error_normalized_by_p0_amp=error,
                          measured=[[z.real,z.imag] for z in measured],expected=[[z.real,z.imag] for z in expected])
            if args.execute and len(errors)>1:
                order=math.log(errors[0]/errors[-1])/math.log(4)
                check('linear-refinement-'+str(method),order>1.7,order=order,errors=errors)
    if args.suite in ['all','advection']:
        errors=[]
        for nx in [64,128,256]:
            out=run(f'advection-{nx}',nx=nx,mode=3,thermal=.1,tlim=.2)
            if args.execute:
                x=(np.arange(nx)+.5)/nx-.4*.2
                targets=[2*(1+.1*np.cos(2*math.pi*x)),2*(1+.05*np.sin(2*math.pi*x))]
                error=max(float(np.mean(np.abs(out[key].reshape(-1)-target))) for key,target in zip(['pp','pt'],targets))
                errors.append(error);check(f'advection-{nx}',error<.003,l1_error=error)
        if args.execute:
            order=math.log(errors[0]/errors[-1])/math.log(4)
            check('advection-refinement',order>1.7,order=order,errors=errors)
    if args.suite in ['all','restart']:
        for label,kw in [('zero-rate',{'lf':'sts','nu':0.0}),
                         ('forced',{'forcing':True,'lf':'sts','nu':.3})]:
            ref=run('restart-direct-'+label,tlim=1,nlim=16,**kw)
            first=run('restart-initial-'+label,tlim=1,nlim=8,**kw)
            path=first['restart'] if args.execute else root/('restart-initial-'+label)/'rst/last.rst'
            continued=run('restart-continued-'+label,restart=path,
                          overrides=['time/nlim=16','output1/file_number=0','output1/last_time=-1',
                                     'output2/file_number=0','output2/last_time=-1'])
            if args.execute:
                mismatch=int(np.count_nonzero(ref['payload'].view('=u8')!=continued['payload'].view('=u8')))
                check('restart-full-state-'+label,mismatch==0 and ref['time']==continued['time'] and ref['dt']==continued['dt'],
                      full_state_bit_mismatches=mismatch,time_equal=ref['time']==continued['time'],
                      next_dt_equal=ref['dt']==continued['dt'])
        unencoded=root/'passive-without-encoding.rst'
        if args.execute:
            raw=path.read_bytes();boundary=raw.index(b'<par_end>\n')+len(b'<par_end>\n')
            header,count=re.subn(rb'^passive_restart_encoding[ \t]*=.*\n',b'',raw[:boundary],flags=re.M)
            assert count==1,'Expected exactly one passive encoding marker'
            unencoded.write_bytes(header+raw[boundary:])
        run('restart-reject-missing-encoding',restart=unencoded,
            expect_failure='restart EOS must agree with passive J/A encoding')
        run('restart-reject-active-mode',restart=path,overrides=['mhd/passive=false'],
            expect_failure='restart EOS must agree with passive J/A encoding')
        run('restart-reject-version',restart=path,overrides=['mhd/passive_restart_encoding=2'],
            expect_failure='restart EOS must agree with passive J/A encoding')
        run('restart-reject-isothermal-mode',restart=path,
            overrides=['mhd/eos=isothermal','mhd/passive=false'],
            expect_failure='<mhd>/cgl_heat_flux requires <mhd>/eos = cgl')
    if args.suite in ['all','fences']:
        for name,flags,error in [
            ('sgs-output',['output3/file_type=bin','output3/variable=mhd_sgs','output3/dt=.1'],
             'passive CGL mhd_sgs output requires a physical-U definition'),
            ('amr',['mesh_refinement/refinement=static'],'passive J/A has not yet validated'),
            ('viscosity',['mhd/viscosity=.001'],'passive J/A has not yet validated'),
            ('resistivity',['mhd/ohmic_resistivity=.001'],'passive J/A has not yet validated'),
            ('hyperviscosity',['mhd/hyperviscosity=.001'],'passive J/A has not yet validated'),
            ('cooling',['mhd_srcterms/ism_cooling=true','mhd_srcterms/hrate=0'],
             'passive J/A has not yet validated'),
            ('rel-cooling',['mhd_srcterms/rel_cooling=true','mhd_srcterms/crate_rel=.001'],
             'passive J/A has not yet validated'),
            ('shearing',['shearing_box/qshear=1.5','shearing_box/omega0=1'],
             'passive J/A has not yet validated'),
            ('coupled-fluid',['ion-neutral/drag_coefficient=1'],
             'passive J/A has not yet validated'),
            ('radiation',['radiation/nangles=4'],'passive J/A has not yet validated'),
            ('outflow',['mesh/ix1_bc=outflow','mesh/ox1_bc=outflow'],'currently requires periodic physical boundaries'),
            ('reflect',['mesh/ix1_bc=reflect','mesh/ox1_bc=reflect'],'currently requires periodic physical boundaries'),
            ('inflow',['mesh/ix1_bc=inflow','mesh/ox1_bc=inflow'],'currently requires periodic physical boundaries'),
            ('user',['mesh/ix1_bc=user','mesh/ox1_bc=user'],'currently requires periodic physical boundaries'),
            ('kinematic',['time/evolution=kinematic'],'requires dynamic evolution'),
            ('ordinary-conduction',['mhd/conductivity=.001'],'Ordinary <mhd>/conductivity is disabled'),
            ('nonpositive-isothermal-speed',['mhd/iso_sound_speed=0'],'passive CGL requires finite positive'),
            ('negative-isothermal-speed',['mhd/iso_sound_speed=-1'],'passive CGL requires finite positive')]:
            run('fence-'+name,overrides=flags,expect_failure=error)
        builtin=args.binary
        run('fence-unconverted-pgen',binary=builtin,overrides=['problem/pgen_name=advection'],
            expect_failure='this built-in pgen has no passive J/A initializer')
    save();print(json.dumps({'root':str(root),'planned_runs':len(records),'checks':checks},indent=2))
def run_validation(directory, suite, launcher=''):
    binary = Path(os.environ.get('ATHENAK_CGL_PASSIVE_BINARY', './athena')).resolve()
    main(['--binary', str(binary), '--output', str(directory),
          '--suite', suite, '--launcher', launcher])


if __name__ == '__main__':
    main()
