"""Fuse existing face bodies without rewriting a floating-point expression."""
from pathlib import Path
import argparse
import difflib
import hashlib
import json

here = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path,
                    default=here.parent / 'task1-research/source-task1-refill',
                    help='Source root containing src/; default is the frozen Task1+refill overlay.')
parser.add_argument('--output', type=Path, default=here,
                    help='Artifact directory for source/, patch, and extraction manifest.')
args = parser.parse_args()
baseline = args.source.resolve()
output = args.output.resolve()
if baseline == output / 'source':
    parser.error('The input source root must differ from the generated source/ overlay.')
output.mkdir(parents=True, exist_ok=True)
relative = Path('src/diffusion/cgl_landau_fluid.cpp')
original = (baseline / relative).read_text()
if 'cgl_lf_fluxes_fused' in original:
    parser.error('Input already contains the fused kernel; start from the unfused integration.')
sha = lambda text: hashlib.sha256(text.encode()).hexdigest()

def body(direction, reduction):
    kind = 'parallel_reduce' if reduction else 'parallel_for'
    marker = f'Kokkos::{kind}("cgl_lf_flux{direction}"'
    start = original.index(marker)
    begin = original.index('    ', original.index('\n', original.index('KOKKOS_LAMBDA', start)) + 1)
    end_marker = f'    }}, Kokkos::Sum<array_sum::GlobalSum>(qstats{direction}));' if reduction else '    });'
    end = original.index(end_marker, begin)
    lines = original[begin:end].splitlines(keepends=True)
    assert len(lines) > 10 and all('const int' in line for line in lines[:4])
    extracted = ''.join(lines[4:])
    assert 'idx' not in extracted
    return extracted

bodies = {(direction, reduction): body(direction, reduction)
          for reduction in (True, False) for direction in (1, 2, 3)}
guards = {
    1: 'j <= je && k <= ke',
    2: 'multi_d && i <= ie && k <= ke',
    3: 'three_d && i <= ie && j <= je',
}

def scopes(reduction):
    result = ''
    for direction in (1, 2, 3):
        result += f'      if ({guards[direction]}) {{\n'
        # Only indentation changes; the exact body bytes are recorded below.
        result += ''.join('    ' + line for line in bodies[direction, reduction].splitlines(keepends=True))
        result += '      }\n'
    return result

decode = '''      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji + ks;
      const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
      const int i = idx - m*nkji - (k - ks)*nji - (j - js)*ni + is;
'''
fused = '''  // Ordinary execution owns each cell's lower x/y/z face, including padded
  // high caps. Directional guards precede every stencil read. Detailed profiling
  // retains the original directional kernels and their separate replay buckets.
  if (!profile_detail_enabled_) {
    Kokkos::Profiling::pushRegion("cgl_lf_fluxes_fused");
    auto f1 = f.x1f;
    auto f2 = f.x2f;
    auto f3 = f.x3f;
    const int ni = ie - is + 2;
    const int nj = multi_d ? je - js + 2 : 1;
    const int nk = three_d ? ke - ks + 2 : 1;
    const int nji = nj*ni;
    const int nkji = nk*nji;
    const int nmkji = (nmb1 + 1)*nkji;
    if (collect_heat_flux_diagnostics) {
      array_sum::GlobalSum qstats_fused;
      Kokkos::parallel_reduce("cgl_lf_fluxes_fused",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &qstats) {
''' + decode + scopes(True) + '''      }, Kokkos::Sum<array_sum::GlobalSum>(qstats_fused));
      AccumulateHeatFluxDiagnostics(qstats_fused);
    } else {
      Kokkos::parallel_for("cgl_lf_fluxes_fused",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx) {
''' + decode + scopes(False) + '''      });
    }
    Kokkos::Profiling::popRegion();
    return;
  }
'''

marker = '  auto f1 = f.x1f;\n'
start = original.index('void CGLLandauFluid::AddHeatFluxes(')
offset = original.index(marker, start)
candidate = original[:offset] + fused + original[offset:]

overlay = output / 'source'
(overlay / 'src').mkdir(parents=True, exist_ok=True)
for item in (baseline / 'src').rglob('*'):
    dest = overlay / item.relative_to(baseline)
    if item.is_dir():
        dest.mkdir(parents=True, exist_ok=True)
    elif not dest.exists() and not dest.is_symlink():
        dest.symlink_to(item)
changed = overlay / relative
if changed.is_symlink():
    changed.unlink()
changed.write_text(candidate)
(output / 'task7-fusion-candidate.patch').write_text(''.join(difflib.unified_diff(
    original.splitlines(keepends=True), candidate.splitlines(keepends=True),
    fromfile='a/' + str(relative), tofile='b/' + str(relative))))
(output / 'extraction-manifest.json').write_text(json.dumps({
    'baseline': str(baseline / relative),
    'baseline_sha256': sha(original),
    'candidate_sha256': sha(candidate),
    'bodies': {f'{direction}-' + ('full' if reduction else 'none'):
               {'sha256': sha(value), 'lines': len(value.splitlines())}
               for (direction, reduction), value in bodies.items()},
    'notes': 'Each copied face body changes indentation only. No header/layout changes. Profile-detail retains original path.'
}, indent=2))
print(json.dumps({'patch_lines': len((output / 'task7-fusion-candidate.patch').read_text().splitlines()),
                  'candidate_sha256': sha(candidate)}))
