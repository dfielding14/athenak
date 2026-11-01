"""Positive-background, receiving-prefactor optical-depth sensitivity.

The field-line coordinate is cell-center xyz. Rates remain finite, including
backup rates. Compare step and terminal reach before interpreting any result.
"""
import numpy as np
from numba import njit
from scipy.special import sici
from fluxes import sample, receivers, config_array, C, ALPHA


@njit
def samples(raw, positions):
    out=np.empty((len(positions),7))
    for n in range(len(positions)):
        p=positions[n];out[n,:6]=sample(raw,p)
        x,y,z=p;i,j,k=int(np.floor(x)),int(np.floor(y)),int(np.floor(z))
        fx,fy,fz=x-i,y-j,z-k;rho=0.
        for dz in range(2):
            for dy in range(2):
                for dx in range(2):
                    w=(fx if dx else 1-fx)*(fy if dy else 1-fy)*(fz if dz else 1-fz)
                    rho+=w*raw[(k+dz)%raw.shape[0],(j+dy)%raw.shape[1],(i+dx)%raw.shape[2],3]
        out[n,6]=rho
    return out


def primitives(x):
    G=np.zeros_like(x);M=np.zeros_like(x)
    use=(x>0)&(x<30);v=x[use]
    si,ci=sici(v)
    F=np.sin(v)*ci-np.cos(v)*(si-np.pi/2)
    G[use]=np.cos(v)*ci+np.sin(v)*(si-np.pi/2)
    M[use]=v*G[use]-F+np.pi/2
    large=x>=30;v=x[large]
    G[large]=-1/v**2+6/v**4-120/v**6+5040/v**8
    M[large]=np.pi/2-2/v+8/v**3-144/v**5+5760/v**7
    return G,M


def segment(lo,hi,flo,fhi):
    gl,ml=primitives(lo);gh,mh=primitives(hi)
    slope=np.divide(fhi-flo,hi-lo,out=np.zeros_like(lo),where=hi>lo)
    out=flo*(gh-gl)+slope*((mh-ml)-lo*(gh-gl))
    return np.where(lo==0,slope*mh,out)


def optical_faces(fields,face_ids,spacing,config,reach,ds):
    """Return q, q_frozen [faces,2], bn,bf and path-event fractions."""
    if config.get('nu_coll',0)<=0:
        raise ValueError('Optical coordinate requires strictly positive background scattering')
    ids=np.asarray(face_ids,dtype=np.int64);dx=np.asarray(spacing,dtype=float)
    if reach<=0 or ds<=0 or ids.ndim!=2 or ids.shape[1]!=4:
        raise ValueError('Positive reach/step and (axis,z,y,x) faces required')
    cfg=config_array(config)
    rec,_=receivers(fields,ids,dx,cfg)
    rho,c,pp,pt,att,bn,bf,nu=rec.T
    if np.any(bf<=cfg[2]) or not np.isfinite(rec).all():
        raise ValueError('Invalid/weak receiving field')
    start=ids[:,1:][:,::-1].astype(float)
    start[np.arange(len(ids)),ids[:,0]]-=.5
    g=pt/rho*(1-pt/pp)
    def source(v):return np.column_stack((v[:,3],v[:,4]-g*v[:,5]))
    def direction(p):
        b=samples(fields,p)[:,:3];norm=np.linalg.norm(b,axis=1)
        if np.any(norm<=0):raise ValueError('Magnetic null along optical ray')
        return b/norm[:,None]
    def rates(v):
        delta=v[:,6]*(v[:,4]-v[:,3]);b2=np.exp(2*v[:,5])
        soft=(bool(cfg[7])&(delta>=.5*cfg[12]*b2))|(bool(cfg[8])&(delta<=-.5*cfg[13]*b2))
        hard=bool(cfg[10])&((bool(cfg[7])&(delta>=.5*cfg[14]*cfg[12]*b2))|
              (bool(cfg[8])&(delta<=np.maximum(-.5*cfg[15]*cfg[13]*b2,-b2))))
        return cfg[6]+np.where(hard,cfg[11],np.where(soft,cfg[9],0.)),soft,hard
    initial=source(samples(fields,start))
    pos=[start.copy(),start.copy()]
    eta=np.zeros((2,len(ids),2));prev=np.zeros_like(eta)
    optical=np.zeros((len(ids),2));frozen=optical.copy()
    steps=int(np.ceil(reach/ds));h=reach/steps
    cross=np.zeros(len(ids),bool);hardseen=cross.copy();events=0
    for it in range(steps):
        for side,sign in enumerate((1,-1)):
            old=pos[side].copy();mid=old+.5*sign*h*direction(old)/dx
            end=old+sign*h*direction(mid)/dx
            values=[samples(fields,p) for p in (old,mid,end)]
            rate=[rates(v) for v in values]
            f=[source(v)-initial for v in values]
            for a,b in ((0,1),(1,2)):
                left=[old,mid][a];right=[mid,end][a]
                fraction=np.zeros(len(ids));current=rate[a][0].copy()
                oldf=f[a].copy();pending=np.ones(len(ids),bool)
                cross|=rate[a][1]|rate[b][1];hardseen|=rate[a][2]|rate[b][2]
                for event in range(6):
                    if not pending.any():break
                    changed=pending&(current!=rate[b][0])
                    finish=np.ones(len(ids))
                    ix=np.flatnonzero(changed)
                    if len(ix):
                        low=fraction[ix].copy();high=np.ones(len(ix))
                        for _ in range(30):
                            t=.5*(low+high)
                            rr=rates(samples(fields,left[ix]+t[:,None]*(right[ix]-left[ix])))[0]
                            same=rr==current[ix]
                            low=np.where(same,t,low);high=np.where(same,high,t)
                        finish[ix]=high;events+=len(ix)
                    indices=np.flatnonzero(pending);target=f[b][indices].copy()
                    if len(ix):
                        vv=samples(fields,left[ix]+finish[ix,None]*(right[ix]-left[ix]))
                        target[np.searchsorted(indices,ix)]=np.column_stack(
                            (vv[:,3],vv[:,4]-g[ix]*vv[:,5]))-initial[ix]
                    next_eta=eta[side,indices]+(.5*h*(finish[indices]-fraction[indices]))[:,None]*current[indices,None]*ALPHA/c[indices,None]
                    optical[indices]-=sign*segment(eta[side,indices],next_eta,oldf[indices],target)/np.pi
                    eta[side,indices]=next_eta;oldf[indices]=target
                    pending=changed;fraction=finish
                    if len(ix):current[ix]=rates(vv)[0]
                if pending.any():raise ValueError('Unresolved multiple limiter transitions; halve ds')
            a=nu[:,None]*ALPHA/c[:,None]
            frozen-=sign*segment(a*it*h,a*(it+1)*h,prev[side],f[-1])/np.pi
            prev[side]=f[-1];pos[side]=end%np.array(fields.shape[:3][::-1])
    pref=rho[:,None]*c[:,None]*C
    caps=c[:,None]*np.column_stack((pp,pt))*C
    def cap(value):
        q=value*pref*att[:,None]
        return q/(1+np.abs(q)/caps)
    result=dict(q=cap(optical),q_frozen=cap(frozen),bn=bn,bf=bf,
                path_limiter_fraction=float(cross.mean()),path_backup_fraction=float(hardseen.mean()),
                transitions=events,actual_ds=h,
                assumptions=['receiving thermal speed, density and magnetic driver coefficient frozen',
                'positive background optical coordinate; finite backup rates retained',
                'source piecewise linear in optical depth; midpoint rays and bisected switches',
                'attenuation applied before saturation, matching main flux calculator',
                'older pilot optical artifact applied attenuation after saturation',
                'bcc face geometry; compare step and terminal reach'])
    if not np.isfinite(result['q']).all():raise ValueError('Nonfinite optical result')
    return result
