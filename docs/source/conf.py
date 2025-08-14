# Configuration file for the Sphinx documentation builder.

import os
import sys
import subprocess
from datetime import datetime

# -- Project information -----------------------------------------------------
project = 'AthenaK'
copyright = f'{datetime.now().year}, AthenaK Development Team'
author = 'AthenaK Development Team'
release = '1.0'
version = '1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'myst_parser',
    'sphinx_copybutton',
    'sphinxcontrib.mermaid',
]

# Optional: Breathe configuration for Doxygen integration (if installed)
try:
    import breathe
    extensions.append('breathe')
    breathe_projects = {
        "AthenaK": "../build/xml"
    }
    breathe_default_project = "AthenaK"
    breathe_default_members = ('members', 'undoc-members')
except ImportError:
    print("Note: breathe not installed, API documentation will be limited")

# Optional: Graphviz support (if installed)
try:
    import sphinx.ext.graphviz
    extensions.append('sphinx.ext.graphviz')
except ImportError:
    pass

# MyST configuration for Markdown support
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# Add any paths that contain templates here
templates_path = ['_templates']

# List of patterns to exclude
exclude_patterns = []

# The suffix(es) of source filenames
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------

# Use Furo theme for modern, responsive design with dark mode
html_theme = 'furo'

# Theme options
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#0066CC",
        "color-brand-content": "#0066CC",
        "font-stack": "Inter, -apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif",
        "font-stack--monospace": "SFMono-Regular, Menlo, Monaco, Consolas, Liberation Mono, Courier New, monospace",
    },
    "dark_css_variables": {
        "color-brand-primary": "#4A9EFF",
        "color-brand-content": "#4A9EFF",
    },
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "top_of_page_button": "edit",
}

# Logo
html_logo = None  # Add logo path if available
html_favicon = None  # Add favicon path if available

# Add any paths that contain custom static files
html_static_path = ['_static']

# Force cache refresh for custom.js
html_js_files = [
    ('custom.js', {'defer': 'defer', 'data-version': str(int(datetime.now().timestamp()))}),
]

# Custom CSS and JS
html_css_files = [
    'custom.css',
]

html_js_files = [
    'custom.js',
]

# HTML title
html_title = "AthenaK Documentation"

# Show source links
html_show_sourcelink = True
html_copy_source = True

# Output file base name for HTML help builder
htmlhelp_basename = 'AthenaKdoc'

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': r'''
\usepackage{amsmath,amssymb}
\setcounter{secnumdepth}{3}
''',
}

latex_documents = [
    ('index', 'AthenaK.tex', 'AthenaK Documentation',
     'AthenaK Development Team', 'manual'),
]

# -- Options for manual page output ------------------------------------------

man_pages = [
    ('index', 'athenak', 'AthenaK Documentation',
     ['AthenaK Development Team'], 1)
]

# -- Options for Texinfo output ----------------------------------------------

texinfo_documents = [
    ('index', 'AthenaK', 'AthenaK Documentation',
     'AthenaK Development Team', 'AthenaK', 
     'High-performance AMR astrophysical simulation framework',
     'Miscellaneous'),
]

# -- Extension configuration -------------------------------------------------

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

# Todo extension
todo_include_todos = True

# Graphviz configuration (if available)
try:
    graphviz_output_format = 'svg'
except:
    pass

# Mermaid configuration for proper rendering
mermaid_version = "10.6.0"
mermaid_init_js = """
mermaid.initialize({
    startOnLoad: true,
    theme: 'neutral',
    themeVariables: {
        primaryColor: '#667eea',
        primaryTextColor: '#fff',
        primaryBorderColor: '#7C0000',
        lineColor: '#764ba2',
        secondaryColor: '#006100',
        tertiaryColor: '#fff'
    },
    flowchart: {
        useMaxWidth: true,
        htmlLabels: true
    }
});
"""
mermaid_output_format = 'raw'  # Render in browser, not pre-rendered
mermaid_cli = 'mmdc'
mermaid_cmd_shell = False

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Copy button for code blocks
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

# Optional: Function to run Doxygen if available
def run_doxygen(app):
    """Run Doxygen to generate XML files for Breathe (if available)."""
    try:
        import shutil
        if shutil.which('doxygen'):
            doxygen_dir = os.path.abspath(os.path.join(app.srcdir, '..'))
            os.chdir(doxygen_dir)
            try:
                subprocess.run(['doxygen', 'Doxyfile'], check=True)
                print("Doxygen XML generated successfully")
            except subprocess.CalledProcessError as e:
                print(f"Doxygen failed (non-critical): {e}")
            finally:
                os.chdir(app.srcdir)
        else:
            print("Doxygen not found - skipping API documentation generation")
    except Exception as e:
        print(f"Doxygen check failed (non-critical): {e}")

def setup(app):
    # Only run doxygen if breathe is available
    try:
        import breathe
        app.connect('builder-inited', run_doxygen)
    except ImportError:
        pass