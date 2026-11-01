"""Focused common-cap and uniform-scattering convergence check."""
import numpy as np
from fluxes import evaluate_faces,config_array,receivers
from optical import optical_faces


def test_optical_cap_matches_main():
    n=12;z,y,x=np.meshgrid(*[(np.arange(n)+.5)/n]*3,indexing='ij')
    angle=1.45*np.sin(4*np.pi*x)
    raw=np.ones((n,n,n,6))
    raw[...,0]=np.sin(angle);raw[...,1]=0.;raw[...,2]=np.cos(angle)
    raw[...,4]=np.exp(2*np.sin(2*np.pi*z))
    raw[...,5]=np.exp(1.5*np.sin(2*np.pi*z))
    ids=np.array([[0,k,3,i] for k in (1,3,6,9) for i in (1,3,5,8)])
    cfg=dict(nu_coll=.25,lf_k_parallel=2*np.pi)
    dx=np.ones(3)/n
    receiver,_=receivers(raw,ids,dx,config_array(cfg))
    assert receiver[:,4].min()<.85
    err=[];uniform=[]
    for h in (1/96,1/192):
        q=optical_faces(raw,ids,dx,cfg,.5,h)
        main=evaluate_faces(raw,ids,dx,cfg,[.5],h,[.1,100])
        uniform.append(np.sqrt(np.mean((q['q']-q['q_frozen'])**2))/
                       np.sqrt(np.mean(q['q_frozen']**2)))
        err.append(np.sqrt(np.mean((q['q']-main['q'][1])**2))/
                   np.sqrt(np.mean(main['q'][1]**2)))
    assert uniform[-1]<.004 and uniform[-1]<.7*uniform[0],uniform
    assert err[-1]<.003 and err[-1]<.7*err[0],err
    print('Attenuated saturated main/optical convergence',err)

if __name__=='__main__':
    test_optical_cap_matches_main()
