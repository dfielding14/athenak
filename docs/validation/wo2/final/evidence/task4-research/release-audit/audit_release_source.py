#!/usr/bin/env python3
"""Read-only audit of accepted fused sources versus final clean release sources."""
from pathlib import Path
import hashlib
import json
import re

W = Path('/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2')
F = W / 'final-research'
O = Path(__file__).resolve().parent
old, new = F / 'source-fused', F / 'source-release'
sha = lambda data: hashlib.sha256(data).hexdigest()
# Preserve literal spellings and operator tokens; discard only comments/whitespace.
lexer = re.compile(r'//[^\n]*|/\*.*?\*/|(?:u8|u|U|L)?"(?:\\.|[^"\\])*"|(?:u|U|L)?\'(?:\\.|[^\'\\])*\'|[A-Za-z_][A-Za-z_0-9]*|(?:\d[\w.]*|\.\d[\w.]*)(?:[+-]\d[\w.]*)?|>>=|<<=|->\*|\.\.\.|::|\+\+|--|->|&&|\|\||<=|>=|==|!=|\+=|-=|\*=|/=|%=|<<|>>|&=|\|=|\^=|\S', re.S)
def tokens(text):
    return [m.group() for m in lexer.finditer(text)
            if not m.group().startswith(('//', '/*'))]
def func(text, name):
    start = text.index('KOKKOS_INLINE_FUNCTION\nvoid '+name)
    opening = text.index('{', start)
    depth = 1
    end = opening+1
    while depth:
        if text[end] == '{': depth += 1
        if text[end] == '}': depth -= 1
        end += 1
    return text[start:end]

a = {str(p.relative_to(old)): p for p in (old/'src').rglob('*') if p.is_file()}
b = {str(p.relative_to(new)): p for p in (new/'src').rglob('*') if p.is_file()}
new_header = 'src/mhd/rsolvers/cgl_pressure_traction.hpp'
assert set(b)-set(a) == {new_header}
assert not set(a)-set(b)
report = {'scope': 'Final numerical C++ source and build flags; no applications launched',
          'accepted_source': str(old), 'release_source': str(new),
          'accepted_source_files': len(a), 'release_source_files': len(b),
          'removed_files': [], 'added_files': [new_header], 'changed_files': {},
          'source_sha256': {}, 'flags': {}}
old_helper = func(a['src/mhd/rsolvers/llf_mhd_singlestate.hpp'].read_text(),
                  'SingleStateLLF_CGLPressureTraction')
new_helper = func(b[new_header].read_text(), 'SingleStateLLF_CGLPressureTraction')
assert old_helper == new_helper
report['pressure_traction_literal_function_identical'] = True
report['pressure_traction_function_sha256'] = sha(old_helper.encode())
allowed = {
 'src/diffusion/cgl_landau_fluid.cpp','src/eos/cgl_mhd.cpp',
 'src/eos/cgl_passive.hpp','src/eos/isothermal_c2p_mhd.hpp','src/mhd/mhd.cpp',
 'src/mhd/rsolvers/hlle_cgl.hpp','src/mhd/rsolvers/llf_mhd_singlestate.hpp',
 'src/outputs/history.cpp','src/pgen/pgen.cpp',
 'src/pgen/tests/cgl_lf_boundary.cpp','src/pgen/tests/cgl_passive_validation.cpp',
 'src/srcterms/srcterms.cpp'}
for name in sorted(a):
    x, y = a[name].read_bytes(), b[name].read_bytes()
    report['source_sha256'][name] = {'accepted': sha(x), 'release': sha(y)}
    if x == y: continue
    assert name in allowed, name
    x, y = x.decode(), y.decode()
    kind = 'Only comments/whitespace; exact remaining token spellings'
    if name == 'src/eos/cgl_passive.hpp':
        before = 'return {log(ppar) + 2.0*logb - 3.0*logr,\n          log(pperp) - log(ppar) + 2.0*logr - 3.0*logb};'
        after = 'return {static_cast<Real>(log(ppar) + 2.0*logb - 3.0*logr),\n          static_cast<Real>(log(pperp) - log(ppar) + 2.0*logr - 3.0*logb)};'
        assert x.count(before) == 1
        x = x.replace(before, after)
        kind = 'Exactly two complete-expression static_cast<Real> wrappers plus comments'
    elif name == 'src/mhd/rsolvers/hlle_cgl.hpp':
        source = '#include "mhd/rsolvers/llf_mhd_singlestate.hpp"'
        target = '#include "mhd/rsolvers/cgl_pressure_traction.hpp"'
        assert x.count(source) == 1
        x = x.replace(source, target)
        kind = 'Only include of shared pressure-traction header replaces full LLF include'
    elif name == 'src/mhd/rsolvers/llf_mhd_singlestate.hpp':
        assert x.count(old_helper) == 1
        x = x.replace(old_helper, '')
        include = '#include "eos/cgl_passive.hpp"'
        assert x.count(include) == 1
        x = x.replace(include, include+'\n#include "mhd/rsolvers/cgl_pressure_traction.hpp"')
        kind = 'Only literal helper removal and shared header include'
    assert tokens(x) == tokens(y), name
    report['changed_files'][name] = {'classification': kind, 'audited': True}
assert set(report['changed_files']) == allowed
report['source_sha256'][new_header] = {'release': sha(b[new_header].read_bytes())}
for name in ['CMakeLists.txt']:
    assert (old/name).read_bytes() == (new/name).read_bytes()
    report['source_sha256'][name] = {'accepted': sha((old/name).read_bytes()),
                                     'release': sha((new/name).read_bytes())}
for backend in ['cpu', 'hip']:
    paths = [F/f'{kind}-{backend}/src/CMakeFiles/athena.dir/flags.make'
             for kind in ['build', 'release-build']]
    flags = [{line.split(' = ', 1)[0]:line.split(' = ', 1)[1]
              for line in p.read_text().splitlines() if line.startswith(('CXX_FLAGS = ', 'CXX_DEFINES = '))}
             for p in paths]
    assert flags[0] == flags[1]
    report['flags'][backend] = {'definitions_and_compile_flags_identical': True,
                                'values': flags[0], 'paths': [str(p) for p in paths]}
proof_path = F/'passive-cast-proof/equivalence-audit.json'
proof = json.loads(proof_path.read_text())
assert proof['passed']
proof_dir = proof_path.parent
for backend in ['cpu', 'hip']:
    assert tokens((proof_dir/backend/'old-cgl_passive.hpp').read_text()) == tokens(a['src/eos/cgl_passive.hpp'].read_text())
    assert tokens((proof_dir/backend/'new-cgl_passive.hpp').read_text()) == tokens(b['src/eos/cgl_passive.hpp'].read_text())
assert (proof_dir/'cpu/old.o').read_bytes() == (proof_dir/'cpu/new.o').read_bytes()
for section in ['text', 'rodata']:
    assert (proof_dir/f'hip/old-device.{section}').read_bytes() == (proof_dir/f'hip/new-device.{section}').read_bytes()
report['prior_cast_proof_chain_reverified'] = True
report['prior_cast_machine_code_proof'] = {'path': str(proof_path),
 'sha256': sha(proof_path.read_bytes()), 'passed': True,
 'scope': proof['scope'], 'normalization': proof['normalization']}
report['limits'] = [
 'This static audit does not replace final old/new application-state comparisons.',
 'Source line wrapping changes diagnostic __LINE__ values/debug locations.',
 'The 25-case release matrix includes eight passive cases, all PLM, with two forced 3D twelve-cycle cases; it compares all output fields, histories, and restart states, including live RNG.',
 'The release matrix is not a fresh same-binary isothermal-reference test, reconstruction matrix, or resumed-run test. Those acceptance gates already passed on the accepted fused binaries; no reconstruction/restart-reader/forcing numerical source changed.',
 'Full float passive applications are not newly claimed; casts repair focused float header compatibility, whose prior compile tests are documented separately.'
]
report['passed'] = True
(O/'source-audit.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({'passed': True, 'source_files': [len(a),len(b)],
                  'changed_existing_files':len(report['changed_files']),
                  'added_headers':1, 'output':str(O/'source-audit.json')},indent=2))
