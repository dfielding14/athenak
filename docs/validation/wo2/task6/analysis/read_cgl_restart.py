"""Read pure CGL+LF restart states for these pinned double-precision probes.

This intentionally handles only one MHD fluid, one passive scalar, and no other
state modules. It validates the per-block byte count before reading the tail.
"""
from pathlib import Path
import struct
import numpy as np


def read_restart(path):
    raw = Path(path).read_bytes()
    header = raw.index(b'<par_end>\n')+len(b'<par_end>\n')
    text = raw[:header].decode()
    assert '<mhd>' in text and '<hydro>' not in text and '<radiation>' not in text
    nmb,root_level = struct.unpack_from('=ii',raw,header)
    indcs = struct.unpack_from('=19i',raw,header+8+72+76)
    ng,nx,ny,nz = indcs[:4]
    time,dt,cycle = struct.unpack_from('=ddi',raw,header+8+72+152)
    loc = np.frombuffer(raw,dtype='=i4',count=nmb*4,offset=header+252).reshape(nmb,4)
    n1=nx+2*ng; n2=ny+2*ng if ny>1 else 1; n3=nz+2*ng if nz>1 else 1
    shapes=[(7,n3,n2,n1),(n3,n2,n1+1),(n3,n2+1,n1),(n3+1,n2,n1)]
    sizes=[int(np.prod(s)) for s in shapes]
    nreal=sum(sizes); offset=len(raw)-nmb*nreal*8
    assert struct.unpack_from('=Q',raw,offset-8)[0]==nreal*8,'unsupported restart state layout'
    blocks=np.frombuffer(raw,dtype='=f8',count=nmb*nreal,offset=offset).reshape(nmb,nreal)
    arrays=[];start=0
    for shape,size in zip(shapes,sizes):
        arrays.append(blocks[:,start:start+size].reshape((nmb,)+shape));start+=size
    u,b1,b2,b3=arrays
    b=np.stack((.5*(b1[:,:,:,:-1]+b1[:,:,:,1:]),
                .5*(b2[:,:,:-1,:]+b2[:,:,1:,:]),
                .5*(b3[:,:-1,:,:]+b3[:,1:,:,:])),axis=1)
    rho=u[:,0];bsqr=(b*b).sum(axis=1);bmag=np.sqrt(bsqr)
    internal=u[:,4]-.5*(u[:,1:4]**2).sum(axis=1)/rho-.5*bsqr
    ratio=np.exp(u[:,5]/rho-2*np.log(rho)+3*np.log(bmag))
    ppar=internal/(.5+ratio);pperp=ratio*ppar
    fields={'rho':rho,'energy':u[:,4],'A':u[:,5],'mu':pperp/bmag,
            'ppar':ppar,'pperp':pperp,'scalar':u[:,6]/rho}
    fields.update({f'B{n+1}':b[:,n] for n in range(3)})
    return dict(nmb=nmb,root_level=root_level,indcs=indcs,time=time,dt=dt,
                cycle=cycle,loc=loc,fields=fields,u=u,bcc=b)
