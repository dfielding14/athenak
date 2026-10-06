"""Append completed final HIP followup evidence without changing original runs."""
from pathlib import Path
import csv, hashlib, json, shutil, xml.etree.ElementTree as ET
here = Path(__file__).resolve().parent
wo2 = here.parent
out = here / 'delivery/docs/validation/wo2/task1-final-followup'
ev = out / 'evidence'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
report = {
    'classification': 'Functional final GPU regression validation; no performance claim',
    'allocation': '5629018',
    'binary': str(here / 'bin/athena-hip'),
    'binary_sha256': sha(here / 'bin/athena-hip'),
    'runtime_script': str(here / 'run_pytest.sh'),
    'runtime_script_sha256': sha(here / 'run_pytest.sh'),
    'source': str(here / 'regression-source'),
    'runs': [],
}
assert report['binary_sha256'] == 'ce444cd1c34992fbb9cb8e1b1ebb1d9ec569b4d0a93f1f969f7470ca395ca027'
files = [(here / 'run_pytest.sh', 'final-run_pytest.sh'),
         (wo2 / 'scripts/pytest_baseline.py', 'final-pytest_baseline.py'),
         (wo2 / 'scripts/launch.py', 'final-launch.py')]
for label, suite, expected, ranks in [
    ('final-task1-16gpu', 'mpi-gpu', 1, [16]),
    ('final-task1-followup-hip', 'cpu', 11, [1, 4]),
]:
    layout = wo2 / 'baseline' / f'{label}-{suite}'
    xml_path = layout / 'results.xml'
    parsed = ET.parse(xml_path)
    testsuite = parsed.getroot().find('testsuite')
    assert testsuite is not None
    assert int(testsuite.attrib['tests']) == expected
    assert all(int(testsuite.attrib[k]) == 0 for k in ('errors', 'failures', 'skipped'))
    case_names = [p.attrib['name'] for p in testsuite.findall('testcase')]
    report['runs'].append({'label': label, 'suite_name': suite,
                          'backend': 'hip', 'ranks': ranks,
                          'passed': expected, 'case_names': case_names,
                          'xml': str(xml_path), 'xml_sha256': sha(xml_path)})
    files.extend([(xml_path, f'{label}-results.xml'),
                  (here / f'{label}-{suite}.log', f'{label}-{suite}.log'),
                  (here / f'{label}-launch.log', f'{label}-launch.log')])
    for p in sorted((layout / 'pytest-tmp').glob('test_decay_resolves_closure_co[0-9]/*.csv')):
        files.append((p, 'hip-' + p.name))
report_path = here / 'task1-final-gpu-followup.json'
report_path.write_text(json.dumps(report, indent=2) + '\n')
files.extend([(report_path, report_path.name),
              (Path(__file__), Path(__file__).name)])
manifest_path = out / 'evidence-manifest.json'
manifest = json.loads(manifest_path.read_text())
for source, name in files:
    shutil.copy2(source, ev / name)
    manifest['files'][name] = {'original_path': str(source),
                               'sha256': sha(source), 'bytes': source.stat().st_size}
manifest['classification'] = ('Task1 final test adaptation; exact final CPU and HIP '
    'validation passed, including two-node 16-GPU AMR churn')
manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
readme_path = out / 'README.md'
text = readme_path.read_text()
old = '''GPU counterparts and the two-node MPI churn remain pending the coordinator's
release after unrelated GPU-failure isolation. This report does not count CPU
execution of GPU-named tests as GPU validation.'''
new = '''The exact final HIP binary also passed all 11 focused checks: both capped
2D/3D churn tests, the wall test, low field, two single-rank heating cases,
one/four-rank heating agreement, and four quantitative decay cases. The separate
two-node 16-GPU churn test passed its original topology, conservation, div-B and
repair assertions. These runs used the final `run_pytest.sh` runtime in allocation
5629018, with binary SHA `ce444cd1c34992fbb9cb8e1b1ebb1d9ec569b4d0a93f1f969f7470ca395ca027`.
The [GPU audit](evidence/task1-final-gpu-followup.json) lists each passed case and
its JUnit evidence. GPU execution is confirmed independently of the earlier CPU
runs of GPU-named tests. All these checks are functional; concurrent tests make
their elapsed times unsuitable for performance claims.'''
assert old in text or new in text
readme_path.write_text(text.replace(old, new))
for name, record in manifest['files'].items():
    assert sha(ev / name) == record['sha256'], name
print(json.dumps(report, indent=2))
print('evidence manifest sha256', sha(manifest_path))
