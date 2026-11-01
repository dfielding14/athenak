"""Portable-reader fixtures for shared, per-rank and uniformly refined output."""
from pathlib import Path
from tempfile import TemporaryDirectory
import itertools
import unittest
import numpy as np
from snapshot_io import load_fields


def fixture(path, locations, level=1, passive='false', shift=0., ghosts=False, root_n=2, origin=(0.,0.,0.), lengths=(1.,1.,1.), locbytes=8):
    header='<mesh>\nnghost=2\n'
    for d in (1,2,3):
        header+=f'nx{d}={root_n}\nx{d}min={origin[d-1]}\nx{d}max={origin[d-1]+lengths[d-1]}\nix{d}_bc=periodic\nox{d}_bc=periodic\n'
    header+='<meshblock>\nnx1=2\nnx2=2\nnx3=2\n<mhd>\neos=cgl\n'
    header+=f'passive={passive}\nnu_coll=20\ncgl_firehose_threshold=oblique\n'
    with path.open('wb') as f:
        f.write(b'Athena binary output version=1.1\n')
        f.write(f'pheader_count=5\ntime=3\ncycle=10\nsize of location={locbytes}\nsize of variable=4\n'.encode())
        f.write(b'nvars=6\nvariables: dens eint p_perp bcc1 bcc2 bcc3\n')
        f.write(f'header offset={len(header.encode())}\n'.encode());f.write(header.encode())
        for loc in locations:
            ind=np.array([2,3,2,3,2,3],dtype=np.int32)
            if ghosts: ind+=1
            f.write(ind.tobytes()); f.write(np.array([*loc,level],dtype=np.int32).tobytes())
            width=np.array(lengths)*2/(root_n*2**level)
            lo=np.array(origin)+np.array(loc)*width;lo[0]+=shift
            geom=np.column_stack((lo,lo+width)).ravel().astype('f'+str(locbytes))
            f.write(geom.tobytes())
            vals=np.empty((6,2,2,2),dtype=np.float32)
            vals[0]=1+loc[0]+2*loc[1]+4*loc[2]
            for j in range(1,6):vals[j]=j+1
            f.write(vals.tobytes())


class ReaderTests(unittest.TestCase):
    def test_uniform_refined_shared_and_perrank(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);locs=list(itertools.product(range(2),repeat=3))
            fixture(root/'shared.bin',locs)
            fixture(root/'rank0.bin',locs[:4]);fixture(root/'rank1.bin',locs[4:])
            a,m=load_fields(root/'shared.bin')
            b,n=load_fields(str(root/'rank*.bin'),root/'cache')
            self.assertEqual(a.shape,(4,4,4,6));self.assertIsInstance(b,np.memmap)
            np.testing.assert_array_equal(a,b)
            self.assertEqual(float(a[0,0,0,3]),1.)
            self.assertEqual(float(a[-1,-1,-1,3]),8.)
            self.assertEqual(m['config']['firehose_threshold'],1.4)
            self.assertEqual(m['config']['nu_coll'],20.)
            self.assertEqual(n['time'],3.)
    def test_multiblock_root_and_float_coordinates(self):
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/'multi.bin'
            for level in (0,1):
                nb=2*2**level
                locs=list(itertools.product(range(nb),repeat=3))
                fixture(p,locs,level=level,root_n=4,origin=(-2.,.3,8.),
                        lengths=(2.,3.,.7),locbytes=4)
                fields,meta=load_fields(p)
                self.assertEqual(fields.shape,(nb*2,nb*2,nb*2,6))
                np.testing.assert_allclose(meta['spacing'],np.array([2,3,.7])/(nb*2))
                self.assertEqual(float(fields[-1,-1,-1,3]),1+7*(nb-1))
                fixture(p,locs,level=level,root_n=4,shift=.001,
                        origin=(-2.,.3,8.),lengths=(2.,3.,.7),locbytes=4)
                with self.assertRaises(ValueError):load_fields(p)

    def test_invalid_files(self):
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/'bad.bin';locs=list(itertools.product(range(2),repeat=3))
            for opts in ({'passive':'1'},{'shift':.001},{'ghosts':True},{'level':-1}):
                fixture(p,locs,**opts)
                with self.assertRaises(ValueError):load_fields(p)
            for positions in (locs[:-1],locs+[locs[0]]):
                fixture(p,positions)
                with self.assertRaises(ValueError):load_fields(p)
            fixture(p,locs)
            data=p.read_bytes();p.write_bytes(data[:-1])
            with self.assertRaises(ValueError):load_fields(p)


if __name__=='__main__':unittest.main()
