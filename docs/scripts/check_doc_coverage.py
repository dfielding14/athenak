#!/usr/bin/env python3
"""
Check documentation coverage for AthenaK source files.

This script scans the source tree and verifies that documentation
exists for each source file.
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

# Colors for terminal output
class Colors:
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def find_source_files(src_dir: Path) -> Set[Path]:
    """Find all C++ source and header files."""
    source_files = set()
    extensions = {'.cpp', '.hpp', '.h', '.c'}
    
    for ext in extensions:
        source_files.update(src_dir.rglob(f'*{ext}'))
    
    # Filter out test files and generated files
    source_files = {
        f for f in source_files 
        if 'test' not in f.parts and 
           'build' not in f.parts and
           '.dox' not in str(f)
    }
    
    return source_files

def find_documented_files(docs_dir: Path) -> Set[str]:
    """Find all files referenced in documentation."""
    documented = set()
    
    # Look for file references in markdown and rst files
    for doc_file in docs_dir.rglob('*.md'):
        content = doc_file.read_text()
        # Look for patterns like `src/path/file.cpp`
        import re
        pattern = r'`(src/[^`]+\.(cpp|hpp|h|c))`'
        matches = re.findall(pattern, content)
        documented.update(m[0] for m in matches)
    
    for doc_file in docs_dir.rglob('*.rst'):
        content = doc_file.read_text()
        # Look for patterns in rst files
        pattern = r'``(src/[^`]+\.(cpp|hpp|h|c))``'
        matches = re.findall(pattern, content)
        documented.update(matches)
    
    return documented

def check_file_documentation(src_file: Path) -> Tuple[bool, str]:
    """Check if a source file has adequate documentation."""
    has_header = False
    has_description = False
    
    try:
        with open(src_file, 'r') as f:
            lines = f.readlines()[:50]  # Check first 50 lines
            
        # Check for file header comment
        for line in lines:
            if '\\file' in line or '@file' in line:
                has_header = True
            if '\\brief' in line or '@brief' in line:
                has_description = True
            # Also accept general header comments
            if line.startswith('//!') or line.startswith('//='):
                has_header = True
                
    except Exception as e:
        return False, f"Error reading file: {e}"
    
    if not has_header:
        return False, "Missing file header documentation"
    if not has_description and src_file.suffix == '.hpp':
        return False, "Missing brief description"
    
    return True, "OK"

def generate_coverage_report(src_dir: Path, docs_dir: Path) -> Dict:
    """Generate documentation coverage report."""
    source_files = find_source_files(src_dir)
    documented_files = find_documented_files(docs_dir)
    
    report = {
        'total_files': len(source_files),
        'documented': 0,
        'missing_docs': [],
        'missing_headers': [],
        'coverage_by_dir': {}
    }
    
    # Check each source file
    for src_file in sorted(source_files):
        rel_path = src_file.relative_to(src_dir.parent)
        str_path = str(rel_path)
        
        # Check if file is referenced in docs
        is_referenced = str_path in documented_files
        
        # Check if file has internal documentation
        has_docs, message = check_file_documentation(src_file)
        
        # Update statistics
        dir_name = str(rel_path.parent)
        if dir_name not in report['coverage_by_dir']:
            report['coverage_by_dir'][dir_name] = {
                'total': 0, 'documented': 0, 'files': []
            }
        
        report['coverage_by_dir'][dir_name]['total'] += 1
        
        if has_docs or is_referenced:
            report['documented'] += 1
            report['coverage_by_dir'][dir_name]['documented'] += 1
        else:
            if not has_docs:
                report['missing_headers'].append((str_path, message))
            if not is_referenced:
                report['missing_docs'].append(str_path)
            report['coverage_by_dir'][dir_name]['files'].append(str_path)
    
    report['coverage_percent'] = (
        report['documented'] / report['total_files'] * 100 
        if report['total_files'] > 0 else 0
    )
    
    return report

def print_report(report: Dict):
    """Print coverage report to console."""
    print(f"\n{Colors.BOLD}AthenaK Documentation Coverage Report{Colors.ENDC}")
    print("=" * 60)
    
    # Overall statistics
    coverage = report['coverage_percent']
    color = Colors.GREEN if coverage >= 80 else Colors.YELLOW if coverage >= 60 else Colors.RED
    
    print(f"\nOverall Coverage: {color}{coverage:.1f}%{Colors.ENDC}")
    print(f"Total Files: {report['total_files']}")
    print(f"Documented: {report['documented']}")
    print(f"Missing: {report['total_files'] - report['documented']}")
    
    # Coverage by directory
    print(f"\n{Colors.BOLD}Coverage by Directory:{Colors.ENDC}")
    for dir_name, stats in sorted(report['coverage_by_dir'].items()):
        dir_coverage = (
            stats['documented'] / stats['total'] * 100 
            if stats['total'] > 0 else 0
        )
        color = Colors.GREEN if dir_coverage >= 80 else Colors.YELLOW if dir_coverage >= 60 else Colors.RED
        print(f"  {dir_name:30} {color}{dir_coverage:5.1f}%{Colors.ENDC} "
              f"({stats['documented']}/{stats['total']})")
    
    # Files missing documentation
    if report['missing_headers']:
        print(f"\n{Colors.YELLOW}Files Missing Header Documentation:{Colors.ENDC}")
        for file_path, message in report['missing_headers'][:10]:
            print(f"  - {file_path}: {message}")
        if len(report['missing_headers']) > 10:
            print(f"  ... and {len(report['missing_headers']) - 10} more")
    
    # Files not referenced in docs
    if report['missing_docs']:
        print(f"\n{Colors.YELLOW}Files Not Referenced in Documentation:{Colors.ENDC}")
        for file_path in report['missing_docs'][:10]:
            print(f"  - {file_path}")
        if len(report['missing_docs']) > 10:
            print(f"  ... and {len(report['missing_docs']) - 10} more")
    
    print("\n" + "=" * 60)
    
    # Return status
    if coverage >= 80:
        print(f"{Colors.GREEN}✓ Documentation coverage is good!{Colors.ENDC}")
        return 0
    elif coverage >= 60:
        print(f"{Colors.YELLOW}⚠ Documentation coverage needs improvement{Colors.ENDC}")
        return 1
    else:
        print(f"{Colors.RED}✗ Documentation coverage is poor{Colors.ENDC}")
        return 2

def main():
    """Main entry point."""
    # Find project directories
    script_dir = Path(__file__).parent
    docs_dir = script_dir.parent / 'source'
    src_dir = script_dir.parent.parent / 'src'
    
    if not src_dir.exists():
        print(f"Error: Source directory not found at {src_dir}")
        sys.exit(1)
    
    if not docs_dir.exists():
        print(f"Error: Documentation directory not found at {docs_dir}")
        sys.exit(1)
    
    # Generate and print report
    report = generate_coverage_report(src_dir, docs_dir)
    status = print_report(report)
    
    # Optional: Write detailed report to file
    report_file = script_dir.parent / 'build' / 'coverage_report.txt'
    report_file.parent.mkdir(exist_ok=True)
    
    with open(report_file, 'w') as f:
        f.write("AthenaK Documentation Coverage Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Overall Coverage: {report['coverage_percent']:.1f}%\n")
        f.write(f"Total Files: {report['total_files']}\n")
        f.write(f"Documented: {report['documented']}\n\n")
        
        f.write("Files Missing Documentation:\n")
        for file_path in report['missing_docs']:
            f.write(f"  - {file_path}\n")
    
    print(f"\nDetailed report written to: {report_file}")
    
    sys.exit(status)

if __name__ == "__main__":
    main()