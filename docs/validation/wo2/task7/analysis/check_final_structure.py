"""Host-only source and integer face-ownership checks; no AthenaK launch."""
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path

import argparse
p=argparse.ArgumentParser();p.add_argument('--artifact-root',type=Path,required=True);args=p.parse_args()
here=args.artifact_root.resolve()
manifest = json.loads((here / 'extraction-manifest.json').read_text())
original = Path(manifest['baseline']).read_text()
candidate = (here / 'source/src/diffusion/cgl_landau_fluid.cpp').read_text()
digest = lambda s: hashlib.sha256(s.encode()).hexdigest()
assert digest(original) == manifest['baseline_sha256']
assert digest(candidate) == manifest['candidate_sha256']
marker = '  // Ordinary execution owns each cell\'s lower x/y/z face'
start = candidate.index(marker)
end = candidate.index('\n  auto f1 = f.x1f;\n', start)
assert candidate[:start] + candidate[end+1:] == original

copied = candidate[start:end]
for direction, diagnostic in itertools.product((1, 2, 3), ('full', 'none')):
    kind = 'parallel_reduce' if diagnostic == 'full' else 'parallel_for'
    pos = original.index(f'Kokkos::{kind}("cgl_lf_flux{direction}"')
    begin = original.index('\n', original.index('KOKKOS_LAMBDA', pos)) + 1
    ending = f'    }}, Kokkos::Sum<array_sum::GlobalSum>(qstats{direction}));' if diagnostic == 'full' else '    });'
    finish = original.index(ending, begin)
    body = ''.join(original[begin:finish].splitlines(keepends=True)[4:])
    assert digest(body) == manifest['bodies'][f'{direction}-{diagnostic}']['sha256']
    assert ''.join('    ' + line for line in body.splitlines(keepends=True)) in copied

cases = []
for nx, ny, nz in itertools.product((1, 2, 8, 32), (1, 2, 7), (1, 2, 5)):
    if ny == 1 and nz != 1:
        continue
    active = (nx, ny, nz)
    extents = (nx+1, ny+1 if ny>1 else 1, nz+1 if nz>1 else 1)
    original_faces = set()
    for direction in range(3):
        if direction == 1 and ny == 1 or direction == 2 and nz == 1:
            continue
        dims = list(active)
        dims[direction] += 1
        original_faces.update((direction, *index) for index in itertools.product(*(range(n) for n in dims)))
    fused_faces = Counter()
    ni, nj, nk = extents
    for idx in range(ni*nj*nk):
        # Match production integer division and i-contiguous storage.
        k = idx//(ni*nj)
        j = (idx-k*ni*nj)//ni
        i = idx-k*ni*nj-j*ni
        if j < ny and k < nz:
            fused_faces[0, i, j, k] += 1
        if ny > 1 and i < nx and k < nz:
            fused_faces[1, i, j, k] += 1
        if nz > 1 and i < nx and j < ny:
            fused_faces[2, i, j, k] += 1
    assert set(fused_faces) == original_faces, active
    assert all(count == 1 for count in fused_faces.values()), active
    cases.append({'cells': active, 'padded_threads': ni*nj*nk,
                  'faces': len(original_faces)})
(here / 'structure-check.json').write_text(json.dumps({'body_identity': True,
    'original_source_unchanged': True, 'cases': cases}, indent=2))
print(f'All six extracted bodies unchanged; {len(cases)} exact face-ownership cases pass.')
