#!/usr/bin/env python3
"""
Scan AthenaK source code for all input parameters
"""

import re
import os
from pathlib import Path
from collections import defaultdict

def scan_for_parameters(src_dir):
    """Scan source code for all parameter reads"""
    params = defaultdict(lambda: defaultdict(dict))
    
    # Patterns to match parameter reading
    patterns = [
        (r'pin->GetString\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)', 'string', None),
        (r'pin->GetInteger\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)', 'int', None),
        (r'pin->GetReal\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)', 'Real', None),
        (r'pin->GetBoolean\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)', 'bool', None),
        (r'pin->GetOrAddString\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*([^)]+)\)', 'string', 2),
        (r'pin->GetOrAddInteger\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*([^)]+)\)', 'int', 2),
        (r'pin->GetOrAddReal\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*([^)]+)\)', 'Real', 2),
        (r'pin->GetOrAddBoolean\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*([^)]+)\)', 'bool', 2),
    ]
    
    for cpp_file in Path(src_dir).rglob('*.cpp'):
        if 'test' in str(cpp_file):
            continue
            
        try:
            with open(cpp_file, 'r') as f:
                lines = f.readlines()
                for line_num, line in enumerate(lines, 1):
                    for pattern, param_type, default_idx in patterns:
                        matches = re.findall(pattern, line)
                        for match in matches:
                            block = match[0]
                            param = match[1]
                            default = match[default_idx] if default_idx and len(match) > default_idx else None
                            
                            # Clean up default value
                            if default:
                                default = default.strip().strip('"')
                            
                            # Store parameter info
                            params[block][param] = {
                                'type': param_type,
                                'default': default if default else 'required',
                                'file': str(cpp_file).replace(str(src_dir) + '/', ''),
                                'line': line_num
                            }
        except Exception as e:
            print(f"Error reading {cpp_file}: {e}")
    
    return params

def generate_markdown_docs(params):
    """Generate markdown documentation for all parameters"""
    output = []
    output.append("# AthenaK Input Parameters Reference\n")
    output.append("Complete list of all input parameters by block, extracted from source code.\n\n")
    
    # Sort blocks alphabetically
    for block in sorted(params.keys()):
        output.append(f"## Input Block: `<{block}>`\n")
        
        # Find which files use this block
        files = set()
        for param_info in params[block].values():
            files.add(param_info['file'].split('/')[-1])
        
        output.append(f"**Used by**: {', '.join(sorted(files))}\n")
        output.append("\n| Parameter | Type | Default | Source |\n")
        output.append("|-----------|------|---------|--------|\n")
        
        # Sort parameters alphabetically
        for param in sorted(params[block].keys()):
            info = params[block][param]
            file_ref = f"{info['file'].split('/')[-1]}:L{info['line']}"
            output.append(f"| `{param}` | {info['type']} | {info['default']} | {file_ref} |\n")
        
        output.append("\n")
    
    # Summary statistics
    total_blocks = len(params)
    total_params = sum(len(block_params) for block_params in params.values())
    output.append(f"\n## Summary\n")
    output.append(f"- Total input blocks: {total_blocks}\n")
    output.append(f"- Total parameters: {total_params}\n")
    
    return ''.join(output)

if __name__ == "__main__":
    src_dir = Path("../src")
    
    print("Scanning source code for parameters...")
    params = scan_for_parameters(src_dir)
    
    print(f"Found {len(params)} input blocks")
    
    # Generate documentation
    docs = generate_markdown_docs(params)
    
    # Write to file
    output_file = "source/reference/input_parameters.md"
    with open(output_file, 'w') as f:
        f.write(docs)
    
    print(f"Documentation written to {output_file}")
    
    # Print summary
    print("\nTop 10 blocks by parameter count:")
    block_counts = [(block, len(params[block])) for block in params]
    for block, count in sorted(block_counts, key=lambda x: x[1], reverse=True)[:10]:
        print(f"  <{block}>: {count} parameters")