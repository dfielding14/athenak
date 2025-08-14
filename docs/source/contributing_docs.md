# Contributing to AthenaK Documentation

This guide explains how to add to and maintain the AthenaK documentation.

## Documentation Structure

```
docs/
├── source/               # Sphinx source files (Markdown/RST)
│   ├── modules/         # Module documentation
│   ├── examples/        # Example problems
│   ├── flowcharts/      # System diagrams
│   ├── reference/       # API and parameter reference
│   └── _static/         # CSS and static files
├── build/               # Generated HTML (not in git)
├── scan_parameters.py   # Auto-extract parameters
└── validate_all_docs.sh # Validation script
```

## Adding New Documentation

### 1. Document a New Module

Create `docs/source/modules/your_module.md`:

```markdown
# Module: Your Module Name

## Overview
Brief description of what the module does.

## Source Location
`src/your_module/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `main.cpp` | Core implementation | `Initialize()`, `Execute()` |

## Configuration Parameters

From `<your_block>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `param1` | Real | 1.0 | What it does |

## Usage Example

\```ini
<your_block>
param1 = 2.0
\```

## See Also
- [Related Module](other_module.md)
- Source: `src/your_module/main.cpp`
```

### 2. Add an Example

Create `docs/source/examples/your_example.md`:

```markdown
# Example: Your Example

Brief description of the physics problem.

## Problem Generator
**Source**: `src/pgen/your_problem.cpp`

## Available Input Files
- **Basic**: `inputs/hydro/your_problem.athinput`
- **With AMR**: `inputs/hydro/your_problem_amr.athinput`

## Physics
Explain the physics being demonstrated.

## Running the Simulation

\```bash
# Build
cmake -B build -DPROBLEM=your_problem
make -C build -j8

# Run
./build/src/athena -i inputs/your_problem.athinput
\```

## Complete Input File

\```ini
<time>
tlim = 1.0

<mesh>
nx1 = 256
# ... etc
\```

## Analysis
How to analyze the results.

## See Also
- [Module Documentation](../modules/your_module.md)
```

### 3. Update Navigation

Add your new page to `docs/source/index.rst`:

```rst
.. toctree::
   :maxdepth: 1
   :caption: 📦 All Modules
   
   modules/index
   modules/your_module  # Add this line
```

## Writing Style Guidelines

### Use Clear Headers
```markdown
# Main Title (Module Name)
## Major Sections
### Subsections
#### Details
```

### Include Code References
Always include file and line references:
```markdown
**Source**: `src/hydro/hydro.cpp:L125-L250`
```

### Add Practical Examples
Show actual usage, not just theory:
```markdown
## Usage Example
\```bash
./athena -i inputs/test.athinput
\```
```

### Cross-Reference
Link to related documentation:
```markdown
See [Hydro Module](../modules/hydro.md) for details.
```

## Documenting Code Changes

When you modify code, update documentation:

1. **New Feature**: Add to relevant module doc
2. **New Parameter**: Update parameter tables
3. **Bug Fix**: Note in "Recent Improvements" section
4. **New Problem Generator**: Create example page

### Example: Adding a Parameter

If you add a parameter to the code:
```cpp
// In hydro.cpp
Real new_param = pin->GetReal("hydro", "new_param");
```

Update documentation:
```markdown
# In modules/hydro.md

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `new_param` | Real | 1.0 | Controls new feature |
```

## Using Documentation Tools

### Extract Parameters Automatically

```bash
cd docs
python scan_parameters.py
# Updates source/reference/input_parameters.md
```

### Validate Documentation

```bash
cd docs
bash validate_all_docs.sh
# Checks for:
# - Broken links
# - Missing modules
# - Parameter coverage
```

### Build HTML Documentation

```bash
cd docs
make clean
make html
# View at build/html/index.html
```

## Mermaid Flowcharts

Use Mermaid for diagrams:

```markdown
\```{mermaid}
flowchart TD
    A[Start] --> B{Decision}
    B -->|Yes| C[Action 1]
    B -->|No| D[Action 2]
\```
```

Ensure proper configuration in `conf.py`:
```python
mermaid_init_js = """
mermaid.initialize({
    startOnLoad: true,
    theme: 'neutral'
});
"""
```

## Common Documentation Patterns

### Module Documentation Template

```markdown
# Module: [Name]

## Overview
[What it does]

## Source Location
`src/[directory]/`

## Key Components
[Table of files]

## Configuration Parameters
[Table of parameters]

## Execution Flow
[Mermaid diagram]

## Key Methods
[Important functions with descriptions]

## Usage Examples
[Practical examples]

## Common Issues
[Troubleshooting]

## See Also
[Related modules]
```

### Parameter Documentation

Always document:
- **Name**: Exact parameter name
- **Type**: Real, int, string, bool
- **Default**: Default value or "required"
- **Source**: File:Line reference
- **Description**: What it controls

## Testing Your Documentation

1. **Build and View**
   ```bash
   make html
   open build/html/index.html
   ```

2. **Check Links**
   ```bash
   # In browser, click all links
   # Or use link checker tool
   ```

3. **Verify Code Examples**
   ```bash
   # Copy examples and test they work
   ./athena -i your_example.athinput
   ```

## Documentation Checklist

Before committing documentation:

- [ ] Spell check completed
- [ ] All links work
- [ ] Code examples tested
- [ ] Parameters match source code
- [ ] Cross-references added
- [ ] Added to index/navigation
- [ ] Mermaid diagrams render
- [ ] Validation script passes

## Getting Help

- **Sphinx Issues**: Check [Sphinx docs](https://www.sphinx-doc.org/)
- **Mermaid Diagrams**: Use [Mermaid live editor](https://mermaid.live/)
- **Markdown**: Follow [CommonMark spec](https://commonmark.org/)

## Advanced Topics

### Custom CSS
Edit `docs/source/_static/custom.css` for styling.

### API Documentation
If Doxygen is installed, API docs are auto-generated.

### Version-Specific Docs
Use git branches for version-specific documentation.

## Quick Commands

```bash
# Full rebuild
cd docs
make clean && make html

# Quick parameter scan
python scan_parameters.py

# Validate everything
bash validate_all_docs.sh

# Open docs
open build/html/index.html
```

## Remember

Good documentation:
- **Helps users** understand how to use the code
- **Helps developers** understand how code works
- **Stays current** with code changes
- **Includes examples** that actually work
- **Cross-references** related topics

Thank you for contributing to AthenaK documentation!