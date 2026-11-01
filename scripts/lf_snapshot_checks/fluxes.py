"""Sampled periodic CGL face fluxes; raw fields are z,y,x,(Bx,By,Bz,rho,ppar,pperp)."""
import numpy as np
from numba import njit, prange
from scipy.integrate import quad_vec
from scipy.optimize import linprog

C=np.sqrt(np.array([8.,2.])/np.pi)
ALPHA=np.array([(3*np.pi-8)/np.sqrt(8*np.pi),np.sqrt(2/np.pi)])

def integrated_kernel():
    grid=np.linspace(-16,16,4097)
    value,_=quad_vec(lambda t:2/np.pi*t*np.exp(-t)*np.arctan2(np.exp(grid),t),
                     0,np.inf,epsabs=1e-13,epsrel=1e-12)
    return grid,value

@njit
def primitive(s,a,grid,values):
    if a==0:return 2*s/np.pi
    x=s*a
    if x<np.exp(grid[0]):return 2*s/np.pi
    if x>np.exp(grid[-1]):return (1-4/(np.pi*x))/a
    pos=(np.log(x)-grid[0])/(grid[1]-grid[0])
    j=min(int(pos),len(grid)-2)
    return (values[j]+(pos-j)*(values[j+1]-values[j]))/a

@njit
def cell(raw,k,j,i):
    nz,ny,nx,_=raw.shape
    v=raw[k%nz,j%ny,i%nx]
    bx,by,bz=np.float64(v[0]),np.float64(v[1]),np.float64(v[2])
    rho=np.float64(v[3])
    return np.array([bx,by,bz,v[4]/rho,v[5]/rho,np.sqrt(bx*bx+by*by+bz*bz)])

@njit
def sample(raw,point):
    x,y,z=point
    i,j,k=int(np.floor(x)),int(np.floor(y)),int(np.floor(z))
    fx,fy,fz=x-i,y-j,z-k
    out=np.zeros(6)
    for dk in range(2):
        for dj in range(2):
            for di in range(2):
                weight=(fx if di else 1-fx)*(fy if dj else 1-fy)*(fz if dk else 1-fz)
                if weight==0:continue
                v=raw[(k+dk)%raw.shape[0],(j+dj)%raw.shape[1],(i+di)%raw.shape[2]]
                for n in range(3):out[n]+=weight*v[n]
                bx,by,bz=np.float64(v[0]),np.float64(v[1]),np.float64(v[2])
                rho=np.float64(v[3])
                out[3]+=weight*(v[4]/rho)
                out[4]+=weight*(v[5]/rho)
                out[5]+=weight*np.log(np.sqrt(bx*bx+by*by+bz*bz))
    return out

@njit
def vl(a,b):
    return 2*a*b/(a+b) if a*b>0 else 0.

@njit
def receivers(raw,ids,dx,cfg):
    # cfg: floors, local flag/background c, background nu, mirror/firehose flags,
    # soft rate, backup flag/rate, mirror/firehose thresholds and backup factors.
    out=np.empty((len(ids),8)); fixed=np.empty((len(ids),2))
    for f in range(len(ids)):
        axis,k,j,i=ids[f]; pos=np.array([i,j,k])
        left=pos.copy();left[axis]-=1
        vr=raw[k,j,i].astype(np.float64);vlr=raw[left[2]%raw.shape[0],left[1]%raw.shape[1],left[0]%raw.shape[2]].astype(np.float64)
        rho=max(.5*(vr[3]+vlr[3]),cfg[0]);pp=max(.5*(vr[4]+vlr[4]),cfg[1])
        pt=max(.5*(vr[5]+vlr[5]),cfg[1])
        right=cell(raw,k,j,i);l=cell(raw,left[2],left[1],left[0])
        bf=.5*(right[5]+l[5]);bv=.5*(right[:3]+l[:3]);norm=np.sqrt(np.sum(bv*bv))
        bh=bv/max(bf,norm);att=np.sqrt(np.sum(bh*bh))
        c=np.sqrt(max(pp/rho,cfg[3])) if cfg[4] else cfg[5]
        delta=pt-pp;b2=bf*bf
        soft=(cfg[7] and delta>=.5*cfg[12]*b2) or (cfg[8] and delta<=-.5*cfg[13]*b2)
        hard=cfg[10] and ((cfg[7] and delta>=.5*cfg[14]*cfg[12]*b2)
                or (cfg[8] and delta<=max(-.5*cfg[15]*cfg[13]*b2,-b2)))
        nu=cfg[6]+(cfg[11] if hard else cfg[9] if soft else 0.)
        gradients=np.zeros(3)
        for d in range(3):
            unit=np.zeros(3,dtype=np.int64);unit[d]=1
            if d==axis:g=(right[3:6]-l[3:6])/dx[d]
            else:
                rp=cell(raw,k+unit[2],j+unit[1],i+unit[0]);rm=cell(raw,k-unit[2],j-unit[1],i-unit[0])
                lp=cell(raw,left[2]+unit[2],left[1]+unit[1],left[0]+unit[0]);lm=cell(raw,left[2]-unit[2],left[1]-unit[1],left[0]-unit[0])
                g=np.empty(3)
                for ch in range(2):
                    g[ch]=vl(vl(rp[3+ch]-right[3+ch],right[3+ch]-rm[3+ch]),
                             vl(lp[3+ch]-l[3+ch],l[3+ch]-lm[3+ch]))/dx[d]
                g[2]=.25*(rp[5]-rm[5]+lp[5]-lm[5])/dx[d]
            gradients+=bh[d]*g
        for ch in range(2):
            pressure=pp if ch==0 else pt
            drive=rho*gradients[ch]/pressure
            if ch==1:drive-=(1-pt/pp)*gradients[2]/bf
            ratio=-drive/(cfg[16]+ALPHA[ch]*nu/c)
            fixed[f,ch]=C[ch]*c*pressure*ratio/(1+abs(ratio))
        out[f]=np.array([rho,c,pp,pt,att,bh[axis],bf,nu])
    return out,fixed

@njit(parallel=True)
def traces(raw,ids,recv,dx,ds,stops,a,rational_bins,grid,values,bfloor):
    out=np.empty((len(stops),3,len(ids),2))
    for f in prange(len(ids)):
        axis,k,j,i=ids[f]
        arms=np.empty((2,3));arms[0]=np.array([float(i),float(j),float(k)])
        arms[0,axis]-=.5;arms[1]=arms[0]
        rho,c,pp,pt,att=recv[f,:5]
        g=pt/rho*(1-pt/pp)
        acc=np.zeros((3,2));prev=np.zeros((3,2));save=0;valid=True
        for step in range(1,stops[-1]+1):
            v=np.empty((2,6))
            for side in range(2):
                sign=1 if side==0 else -1
                first=sample(raw,arms[side]);norm=np.sqrt(np.sum(first[:3]**2))
                if not norm>bfloor:valid=False;break
                mid=arms[side]+.5*sign*ds*first[:3]/(dx*norm)
                v[side]=sample(raw,mid);norm=np.sqrt(np.sum(v[side,:3]**2))
                if not norm>bfloor:valid=False;break
                arms[side]+=sign*ds*v[side,:3]/(dx*norm)
            if not valid:
                out[:,:,f,:]=np.nan
                break
            s=step*ds
            for ch in range(2):
                jump=v[0,3+ch]-v[1,3+ch]
                if ch==1:jump-=g*(v[0,5]-v[1,5])
                for method in range(3):
                    if method<2:
                        rate=ALPHA[ch]*recv[f,7]/c if method==0 else a[ch]
                        p=primitive(s,rate,grid,values)
                    else:p=rational_bins[step-1,ch]
                    acc[method,ch]+=(p-prev[method,ch])*jump/((2*step-1)*ds)
                    prev[method,ch]=p
            if step==stops[save]:
                for method in range(3):
                    for ch in range(2):
                        pressure=pp if ch==0 else pt
                        ratio=-rho*att*acc[method,ch]/pressure
                        out[save,method,f,ch]=C[ch]*c*pressure*ratio/(1+abs(ratio))
                save+=1
    return out

def config_array(config):
    """Normalize CGL input parameters for the shared sampled recovery kernels."""
    config=dict(config)
    if 'lf_k_parallel' not in config:raise ValueError('Snapshot must specify lf_k_parallel')
    if 'cgl_firehose_threshold' in config:
        alias=config['cgl_firehose_threshold']
        if alias not in ('parallel','oblique'):raise ValueError('Invalid legacy firehose threshold')
        value=1.4 if alias=='oblique' else 2.
        if 'firehose_threshold' in config and float(config['firehose_threshold'])!=value:
            raise ValueError('Conflicting firehose threshold and legacy alias')
        config['firehose_threshold']=value
    cfg=np.array([config.get(k,v) for k,v in [
        ('dfloor',1e-12),('pfloor',1e-12),('bfloor',1e-12),('tfloor',0.),
        ('_local',config.get('lf_coefficient_mode','local')=='local'),('lf_c_parallel0',1.),
        ('nu_coll',0.),('mirror_limiter',False),('firehose_limiter',False),
        ('limiter_nu_coll',0.),('backup_limiters',False),('limiter_backup_nu',1e10),
        ('mirror_threshold',1.),('firehose_threshold',2.),('mirror_backup_factor',2.),
        ('firehose_backup_factor',1.),('lf_k_parallel',2*np.pi)]],dtype=float)
    if not np.isfinite(cfg).all() or cfg[16]<=0 or cfg[6]<0 or cfg[9]<0 or cfg[11]<0:
        raise ValueError('Invalid LF coefficients')
    return cfg

def evaluate_faces(fields,face_ids,spacing,config,radii,ds,k_band,poles=8):
    """Return method-major q, shared face geometry and explicitly controlled references."""
    raw=fields; ids=np.asarray(face_ids,dtype=np.int64);dx=np.asarray(spacing,dtype=float)
    if ids.ndim!=2 or ids.shape[1]!=4 or not np.all((ids[:,0]>=0)&(ids[:,0]<3)):
        raise ValueError('face_ids must be (axis,z,y,x)')
    if raw.ndim!=4 or raw.shape[-1]!=6 or np.any(dx<=0) or ds<=0:
        raise ValueError('Expected raw zyx6, positive xyz spacing and step')
    if np.any(ids[:,1:]<0) or np.any(ids[:,1:]>=np.array(raw.shape[:3])):
        raise ValueError('Face cell indices out of range')
    radii=np.asarray(radii,dtype=float)
    if not np.isfinite(radii).all() or not np.isfinite(ds):raise ValueError('Nonfinite reaches/step')
    stops=np.rint(radii/ds).astype(np.int64)
    if not len(stops) or stops[0]<1 or np.any(np.diff(stops)<=0):
        raise ValueError('Reaches must yield distinct increasing positive step counts')
    kmin,kmax=map(float,k_band)
    if not np.isfinite([kmin,kmax]).all() or not 0<kmin<kmax:raise ValueError('k_band requires finite 0<kmin<kmax')
    cfg=config_array(config)
    recv,fixed=receivers(raw,ids,dx,cfg)
    if np.any(recv[:,6]<=cfg[2]) or not np.all(np.isfinite(recv)):
        raise ValueError('Weak/nonfinite sampled face state needs separate treatment')
    a=ALPHA*cfg[6]/np.median(recv[:,1])
    shifts=np.geomspace(.2*kmin,4*kmax,poles);weights=[];errors=[]
    k=np.geomspace(kmin,kmax,601)
    for rate in a:
        A=(k[:,None]+rate)*shifts/(k[:,None]**2+shifts**2);one=np.ones(len(k))
        fit=linprog(np.r_[np.zeros(poles),1.],A_ub=np.vstack(
            (np.column_stack((A,-one)),np.column_stack((-A,-one)))),
            b_ub=np.r_[one,-one],bounds=(0,None),method='highs')
        if not fit.success:raise ValueError(fit.message)
        weights.append(fit.x[:-1]*shifts);errors.append(float(fit.x[-1]))
    grid,values=integrated_kernel()
    weights=np.array(weights)
    distances=ds*np.arange(1,stops[-1]+1)
    t=distances[:,None]*shifts
    basis=(-np.expm1(-t)-t*np.exp(-t))/shifts**2
    rational_bins=basis@weights.T
    flux=traces(raw,ids,recv,dx,ds,stops,a,rational_bins,grid,values,cfg[2])
    q=np.concatenate((fixed[None],flux[:,0],flux[:,1],flux[:,2]),axis=0)
    if not np.all(np.isfinite(q)):raise ValueError('Ray crossed a magnetic null/nonfinite state')
    n=len(stops)
    labels=['fixed_k']+[f'SHD_receiver_R={r:g}' for r in stops*ds]+[
        f'SHD_uniform_R={r:g}' for r in stops*ds]+[f'rational_uniform_R={r:g}' for r in stops*ds]
    return dict(q=q,bn=recv[:,5],bf=recv[:,6],nu=recv[:,7],labels=labels,
        radii=stops*ds,receiver_reference=n,uniform_reference=2*n,
        uniform_inverse_scattering_length=a,poles=shifts,
        fit_relative_error=errors,fit_band=[kmin,kmax],
        assumptions=['bcc substitutes for CT normal field','receiver thermal/scattering prefactors frozen',
        'rational and uniform SHD exclude intermittent scattering and use identical finite ray lengths',
        'finite traces do not establish infinite-reach Helmholtz solve accuracy'])
