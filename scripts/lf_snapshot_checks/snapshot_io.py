"""Read complete periodic uniform-leaf AthenaK v1.1 primitive snapshots lazily."""
from pathlib import Path
import glob
import numpy as np


class NativeSnapshot:
    """Index shared output or a glob containing every rank of one output."""
    def __init__(self, paths):
        paths = sorted(glob.glob(str(paths))) if isinstance(paths, (str, Path)) else list(paths)
        if not paths:
            raise ValueError('No snapshot files matched')
        self.blocks = []
        reference = None
        for name in paths:
            path = Path(name).resolve()
            with path.open('rb') as f:
                if f.readline().split()[-1] != b'version=1.1':
                    raise ValueError('Expected AthenaK binary version 1.1')
                count = int(f.readline().split(b'=')[-1])
                pre = dict(line.decode().strip().split('=', 1) for line in [f.readline() for _ in range(count-1)])
                pre = {k.strip(): v.strip() for k, v in pre.items()}
                nvar = int(f.readline().split(b'=')[-1])
                names = f.readline().decode().split()[1:]
                size = int(f.readline().split(b'=')[-1])
                header = f.read(size).decode()
                pars, block = {}, ''
                for line in header.splitlines():
                    line = line.split('#', 1)[0].strip()
                    if line.startswith('<'): block = line.strip('<>')
                    elif '=' in line:
                        key, val = line.split('=', 1)
                        pars[block+'/'+key.strip()] = val.strip()
                signature = (pre['time'], pre['cycle'], tuple(names), header)
                if reference is not None and signature != reference:
                    raise ValueError('Rank files differ in time, cycle, fields or parameter header')
                reference = signature
                locbytes, varbytes = int(pre['size of location']), int(pre['size of variable'])
                if locbytes not in (4,8) or varbytes not in (4,8) or nvar != len(names):
                    raise ValueError('Invalid binary sizes or variable count')
                locdtype, dtype = np.dtype('f'+str(locbytes)), np.dtype('f'+str(varbytes))
                while f.tell() < path.stat().st_size:
                    ind = np.frombuffer(f.read(24), dtype=np.int32)
                    logical = np.frombuffer(f.read(16), dtype=np.int32)
                    geom = np.frombuffer(f.read(6*locbytes), dtype=locdtype)
                    if len(ind)!=6 or len(logical)!=4 or len(geom)!=6:
                        raise ValueError('Truncated block header')
                    shape = tuple((ind[1::2]-ind[::2]+1)[::-1])
                    count = int(np.prod(shape))*nvar
                    offset = f.tell()
                    if offset+count*varbytes > path.stat().st_size:
                        raise ValueError('Truncated block data')
                    mbshape = np.array([int(pars[f'meshblock/nx{d}']) for d in (1,2,3)])
                    ng = int(pars['mesh/nghost'])
                    if not np.array_equal(ind[::2], np.full(3,ng)) or not np.array_equal(ind[1::2], ng+mbshape-1):
                        raise ValueError('Reject ghost-included or sliced output; require complete active blocks')
                    descriptor = (path, dtype, offset, (nvar,)+shape)
                    self.blocks.append((logical.copy(), geom.copy(), descriptor))
                    f.seek(count*varbytes, 1)
        passive = pars.get('mhd/passive','false').lower()
        if passive not in ('false','0','true','1'):
            raise ValueError('Invalid passive boolean')
        if pars.get('mhd/eos') != 'cgl' or passive in ('true','1'):
            raise ValueError('Require active CGL primitive output')
        for d in (1,2,3):
            for side in ('i','o'):
                if pars.get(f'mesh/{side}x{d}_bc') != 'periodic':
                    raise ValueError('Require periodic boundaries in all directions')
        self.origin = np.array([float(pars[f'mesh/x{d}min']) for d in (1,2,3)])
        self.lengths = np.array([float(pars[f'mesh/x{d}max']) for d in (1,2,3)])-self.origin
        levels = {int(b[0][3]) for b in self.blocks}
        shapes = {b[2][3][1:] for b in self.blocks}
        if len(levels)!=1 or len(shapes)!=1:
            raise ValueError('Require one uniform leaf level and block shape')
        self.block_shape = np.array(next(iter(shapes))[::-1])
        rootshape = np.array([int(pars[f'mesh/nx{d}']) for d in (1,2,3)])
        level = next(iter(levels))
        # The exporter records loc.level - root_level, already a relative level.
        if level < 0 or level > 30:
            raise ValueError('Invalid physical refinement level')
        self.shape = rootshape*2**level
        self.spacing = self.lengths/self.shape
        if np.any(self.shape % self.block_shape):
            raise ValueError('Block shape does not divide domain')
        spacing = self.spacing*self.block_shape
        coord_tol = 8*np.finfo(locdtype).eps*np.maximum.reduce(
            [np.abs(self.origin),np.abs(self.origin+self.lengths),self.lengths])
        if np.any(coord_tol > .05*self.spacing):
            raise ValueError('Coordinate precision cannot resolve the cell lattice')
        self.lookup = {}
        for logical, geom, data in self.blocks:
            start = np.rint((geom[::2]-self.origin)/self.spacing).astype(int)
            key = tuple(start//self.block_shape)
            if not np.allclose(geom[::2],self.origin+start*self.spacing,rtol=0,atol=coord_tol):
                raise ValueError('Block start lies off the cell lattice')
            if tuple(logical[:3]) != key:
                raise ValueError('Block logical location disagrees with geometry')
            if np.any(start % self.block_shape) or key in self.lookup:
                raise ValueError('Overlapping or misaligned blocks')
            if not np.allclose(geom[1::2]-geom[::2], spacing, rtol=0, atol=2*coord_tol):
                raise ValueError('Block geometry has inconsistent spacing')
            self.lookup[key] = data
        expected = int(np.prod(self.shape//self.block_shape))
        if len(self.lookup)!=expected or any(np.any(np.array(k)<0) or np.any(np.array(k)>=self.shape//self.block_shape) for k in self.lookup):
            raise ValueError('Snapshot files do not cover the full domain')
        parallel = 'p_parallel' if 'p_parallel' in names else 'eint'
        required = ['dens', parallel, 'p_perp', 'bcc1','bcc2','bcc3']
        if not all(k in names for k in required):
            raise ValueError('Require mhd_w_bcc primitive output with both pressures')
        self.indices = [names.index(k) for k in required]
        self.time = float(pre['time'])
        self.metadata = dict(source=[str(Path(p).resolve()) for p in paths], header=header.splitlines(), cycle=int(pre['cycle']), parameters=pars)


def load_fields(path, cache_dir=None):
    """Return [z,y,x,6] Bx,By,Bz,rho,ppar,pperp plus parsed metadata."""
    src=path if isinstance(path,NativeSnapshot) else NativeSnapshot(path)
    shape=tuple(src.shape[::-1])+(6,)
    if cache_dir is None:
        fields=np.empty(shape,dtype=np.float32)
    else:
        dest=Path(cache_dir); dest.mkdir(parents=True,exist_ok=True)
        fields=np.lib.format.open_memmap(dest/'fields.npy',mode='w+',dtype=np.float32,shape=shape)
    for key,descriptor in src.lookup.items():
        path,dtype,offset,blockshape=descriptor
        data=np.memmap(path,mode='r',dtype=dtype,offset=offset,shape=blockshape)
        x,y,z=np.array(key)*src.block_shape
        nx,ny,nz=src.block_shape
        for j,index in enumerate(src.indices[3:]+src.indices[:3]):
            values=data[index]
            if not np.isfinite(values).all() or (j>=3 and np.any(values<=0)):
                raise ValueError('Nonfinite primitive or nonpositive density/pressure in block')
            fields[z:z+nz,y:y+ny,x:x+nx,j]=values
        del values
        data._mmap.close()
    config={}
    for key,value in src.metadata['parameters'].items():
        if not key.startswith('mhd/'): continue
        name=key.split('/',1)[1]
        if value.lower() in ('true','false'): config[name]=value.lower()=='true'
        else:
            try: config[name]=float(value)
            except ValueError: config[name]=value
    if 'firehose_threshold' not in config and 'cgl_firehose_threshold' in config:
        legacy=config['cgl_firehose_threshold']
        if legacy not in ('parallel','oblique'):
            raise ValueError('Unknown legacy firehose threshold')
        config['firehose_threshold']=2.0 if legacy=='parallel' else 1.4
    meta=dict(src.metadata,spacing=src.spacing.tolist(),lengths=src.lengths.tolist(),
              origin=src.origin.tolist(),time=src.time,shape=[int(n) for n in shape[:3]],config=config)
    return fields,meta
