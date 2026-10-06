#!/usr/bin/env python3
"""Retain every test observation; aggregate explicit ordered final reruns by node id."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--note', required=True)
parser.add_argument('xml', type=Path, nargs='+')
args = parser.parse_args()
results, evidence, original_ids = {}, [], None
for source in args.xml:
    root = ET.parse(source)
    records = []
    for case in root.iter('testcase'):
        key = case.get('classname').split('.')[-1] + '::' + case.get('name')
        state = 'passed'
        detail = ''
        if case.find('failure') is not None or case.find('error') is not None:
            state = 'failed'
            item = case.find('failure')
            if item is None:
                item = case.find('error')
            detail = item.get('message', '')
        elif case.find('skipped') is not None:
            item = case.find('skipped')
            state = 'xfailed' if item.get('type') == 'pytest.xfail' else 'skipped'
            detail = item.get('message', '')
        observation = {'state': state, 'evidence': str(source), 'detail': detail}
        results.setdefault(key, {'observations': []})['observations'].append(observation)
        results[key]['final'] = state
        records.append(key)
    if original_ids is None:
        original_ids = set(records)
    evidence.append({'path': str(source), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'tests': len(records)})
counts = {state: sum(results[key]['final'] == state for key in original_ids)
          for state in ('passed', 'failed', 'xfailed', 'skipped')}
record = {'classification': args.note, 'original_suite_counts': counts,
          'original_test_count': len(original_ids), 'evidence': evidence,
          'results': [{'test': k, 'in_original_suite': k in original_ids, **v}
                      for k, v in sorted(results.items())]}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(counts))
print('Unresolved:', [k for k in sorted(original_ids) if results[k]['final'] == 'failed'])
