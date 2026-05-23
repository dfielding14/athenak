"""Sphinx configuration for AthenaK documentation."""

from __future__ import annotations

import os
import sys
from datetime import datetime

# -- Paths --------------------------------------------------------------------

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, PROJECT_ROOT)

# -- Project information ------------------------------------------------------

project = "AthenaK"
author = "AthenaK Developers"
copyright = f"{datetime.now():%Y}, {author}"

# Use git metadata or fall back to a default version tag
release = os.environ.get("ATHENAK_DOC_VERSION", "latest")
version = release

# -- General configuration ----------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.githubpages",
    "sphinxcontrib.mermaid",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "linkify",
    "substitution",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build", "Thumbs.db", ".DS_Store"]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

todo_include_todos = True

nitpicky = True
nitpick_ignore: list[tuple[str, str]] = []

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

# -- HTML output --------------------------------------------------------------

# Theme selection (default: Press). Set ATHENAK_DOC_THEME=readthedocs to switch.
_THEME_ALIASES = {
    "press": "press",
    "press-theme": "press",
    "press_theme": "press",
    "readthedocs": "sphinx_rtd_theme",
    "read_the_docs": "sphinx_rtd_theme",
    "rtd": "sphinx_rtd_theme",
    "sphinx_rtd_theme": "sphinx_rtd_theme",
}
_theme_key = os.environ.get("ATHENAK_DOC_THEME", "readthedocs").strip().lower()
_resolved_theme = _THEME_ALIASES.get(_theme_key, "press")

if _resolved_theme == "sphinx_rtd_theme":
  html_theme = "sphinx_rtd_theme"
  html_theme_options: dict[str, object] = {
      "collapse_navigation": False,
      "navigation_depth": 3,
  }
else:
  html_theme = "press"
  html_theme_options: dict[str, object] = {}

html_static_path = ["_static"]

# -- MyST substitutions -------------------------------------------------------

myst_substitutions = {
    "project": project,
}
