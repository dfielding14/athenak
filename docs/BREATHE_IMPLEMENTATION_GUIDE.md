# Breathe Implementation Guide for AthenaK

## Overview
This guide outlines a pragmatic, incremental approach to adding Breathe/Doxygen integration to AthenaK documentation. The goal is minimal initial setup for cross-referencing benefits, with the option to deepen API documentation over time.

## Philosophy: Index-First, Comments Later
- **Phase 1**: Basic structural index without comments (immediate navigation benefits)
- **Phase 2**: Selective `@brief` comments in hot areas as code is touched
- **Phase 3**: Deeper API documentation only where it adds value

## Prerequisites
```bash
pip install breathe
apt-get install doxygen graphviz  # or brew install on macOS
```

## Step 1: Configure Doxygen

### Create/Update Doxyfile
```ini
# Essential settings
PROJECT_NAME           = "AthenaK"
PROJECT_BRIEF          = "Performance-portable astrophysical simulation framework"
OUTPUT_DIRECTORY       = docs/_doxygen
GENERATE_HTML          = NO        # We only need XML for Breathe
GENERATE_XML           = YES
XML_OUTPUT             = xml
EXTRACT_ALL            = YES        # Index everything, even without comments
EXTRACT_PRIVATE        = NO         # Skip private members
EXTRACT_STATIC         = YES
BUILTIN_STL_SUPPORT    = YES

# Input configuration
INPUT                  = ../src
FILE_PATTERNS          = *.cpp *.hpp *.h
RECURSIVE              = YES
EXCLUDE_PATTERNS       = */tests/* */old/* *_old_*

# Kokkos preprocessing (critical for parsing)
ENABLE_PREPROCESSING   = YES
MACRO_EXPANSION        = YES
EXPAND_ONLY_PREDEF     = YES
PREDEFINED             = KOKKOS_INLINE_FUNCTION= \
                         KOKKOS_FUNCTION= \
                         KOKKOS_FORCEINLINE_FUNCTION= \
                         KOKKOS_LAMBDA(...)=[](auto&&... args) \
                         "par_for(name,space,...)=par_for" \
                         "KOKKOS_CLASS_LAMBDA(...)=[](auto&&... args)"

# Skip detailed parsing of problem areas
EXCLUDE_SYMBOLS        = Kokkos::* \
                         impl::* \
                         detail::*

# Graphviz for diagrams (optional but valuable)
HAVE_DOT               = YES
DOT_NUM_THREADS        = 4
CALL_GRAPH             = YES
CALLER_GRAPH           = YES
DOT_IMAGE_FORMAT       = svg
INTERACTIVE_SVG        = YES
```

## Step 2: Configure Sphinx

### Update conf.py
```python
# Add to extensions list
extensions = [
    'myst_parser',
    'breathe',
    'sphinx.ext.graphviz',
    # ... existing extensions
]

# Breathe configuration
breathe_projects = {
    "athenak": "_doxygen/xml"
}
breathe_default_project = "athenak"

# Optional: Set default flags for all directives
breathe_default_members = ('members', 'undoc-members')
breathe_show_include = False  # Hide #include statements

# Suppress warnings for missing references (initially)
nitpicky = False
```

## Step 3: Build Workflow

### Create Makefile target
```makefile
# In docs/Makefile
doxygen:
	doxygen Doxyfile

html: doxygen
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf $(BUILDDIR)/* _doxygen/
```

## Step 4: Minimal Initial Integration

### Start with index pages only
Create `docs/source/reference/api_index.md`:
```markdown
# API Index

This is an automatically generated index of AthenaK's C++ API.

## Core Classes

```{eval-rst}
.. doxygenclass:: Mesh
   :members:
   :undoc-members:
   :project: athenak

.. doxygenclass:: MeshBlock
   :members:
   :undoc-members:
   :project: athenak

.. doxygenclass:: Hydro
   :members:
   :undoc-members:
   :project: athenak
```

## Namespaces

```{eval-rst}
.. doxygennamespace:: hydro
   :project: athenak
   :members:
   :undoc-members:
```
```

### Add selective cross-references in existing docs
In module documentation, add links to API:
```markdown
The {cpp:class}`Mesh` class manages domain decomposition...

See {cpp:func}`par_for` for parallel execution patterns...
```

## Step 5: Handle Kokkos-Specific Issues

### Create wrapper macros for documentation
In a header file visible to Doxygen:
```cpp
#ifdef DOXYGEN_PROCESSING
  // Simplified signatures for Doxygen
  #define KOKKOS_LAMBDA(...) [=]
  #define par_for(name, policy, ...) \
    template<typename F> void par_for(const std::string& name, F kernel)
#endif
```

### Document complex templates minimally
```cpp
/// @brief 5D array for device data
/// First dimension is MeshBlock, second is variable, rest are spatial
using DvceArray5D = Kokkos::View<Real*****, LayoutWrapper, DevMemSpace>;
```

## Step 6: Incremental Enhancement

### Priority order for adding comments:
1. **Public API entry points** - Problem generators, main classes
2. **Frequently referenced types** - Core data structures
3. **Complex algorithms** - Where the code isn't self-documenting
4. **Internal helpers** - Only if confusion arises

### Minimal comment style
```cpp
/// @brief Evolve system by one timestep
void Execute();

/// @brief Calculate fluxes using HLLC solver
/// @param prim Primitive variables [in]
/// @param flux Computed fluxes [out]
void CalculateFluxes(const DvceArray5D& prim, DvceFaceFld5D& flux);
```

## Step 7: Managing Warnings

### Suppress unhelpful warnings
In conf.py:
```python
# Suppress specific warnings
suppress_warnings = ['breathe.function_undocumented']

# Or selectively for specific directives
breathe_debug_trace_directives = False
breathe_debug_trace_doxygen_ids = False
```

### Use `:no-link:` for problematic symbols
```rst
.. doxygenclass:: ComplexTemplate
   :no-link:
   :project: athenak
```

## Benefits Without Full Documentation

Even with minimal/no comments, you get:
1. **Clickable cross-references** between prose and code
2. **File dependency graphs** 
3. **Class hierarchy diagrams**
4. **Include graphs**
5. **Namespace organization**
6. **Symbol search**

## Maintenance Strategy

1. **Don't force documentation** - Let it grow organically
2. **Document on modification** - Add comments when touching code
3. **Focus on interfaces** - Document public API, not implementation
4. **Use existing comments** - Convert existing `//!` comments gradually
5. **Link over duplicate** - Reference API from narrative rather than re-explaining

## Common Pitfalls to Avoid

1. **Don't document getters/setters** - Self-evident from names
2. **Don't repeat parameter names** - `@param x The x coordinate` is noise
3. **Don't parse all of Kokkos** - Use EXCLUDE_SYMBOLS
4. **Don't generate HTML from Doxygen** - Let Sphinx handle presentation
5. **Don't aim for 100% coverage** - Focus on value-added documentation

## Testing the Setup

```bash
# Generate Doxygen XML
doxygen Doxyfile

# Check XML was created
ls docs/_doxygen/xml/

# Build Sphinx with Breathe
cd docs && make clean html

# Look for successful API pages
open build/html/reference/api_index.html
```

## Future Enhancements

Once basic setup works:
- Add search functionality with `breathe_separate_member_pages`
- Generate UML diagrams with PlantUML
- Add example code snippets with `@code` blocks
- Create module-specific API pages
- Add performance annotations with custom tags

## Conclusion

Start minimal: get cross-references working with zero comments. Add documentation only where it provides value beyond what the code already expresses. The infrastructure enables future enhancement without mandatory up-front investment.