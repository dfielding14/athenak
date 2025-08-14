# Makefile for Sphinx documentation

# You can set these variables from the command line
# Remove -W to allow building with warnings for now
SPHINXOPTS    = --keep-going
SPHINXBUILD   = sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help"
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Run Doxygen first to generate XML (optional)
doxygen:
	@if command -v doxygen >/dev/null 2>&1; then \
		echo "Running Doxygen to generate XML..."; \
		doxygen Doxyfile; \
	else \
		echo "Doxygen not found - skipping API documentation generation"; \
		echo "Install with: brew install doxygen"; \
	fi

# Clean build directory
clean:
	@echo "Cleaning build directory..."
	rm -rf $(BUILDDIR)/*

# Build HTML documentation
html:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

# Build with optional Doxygen
html-full: doxygen html

# Build with auto-reload for development
livehtml:
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" \
		--watch ../src \
		--ignore "*.swp" \
		--ignore "*.pyc" \
		--ignore ".git/*" \
		--port 8000

# Check for broken links
linkcheck:
	@$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck" $(SPHINXOPTS) $(O)
	@echo
	@echo "Link check complete; look for any errors in $(BUILDDIR)/linkcheck/output.txt."

# Check documentation coverage
coverage:
	@$(SPHINXBUILD) -b coverage "$(SOURCEDIR)" "$(BUILDDIR)/coverage" $(SPHINXOPTS) $(O)
	@echo
	@echo "Coverage check complete; look for any errors in $(BUILDDIR)/coverage/python.txt."

# Spell checking
spelling:
	@$(SPHINXBUILD) -b spelling "$(SOURCEDIR)" "$(BUILDDIR)/spelling" $(SPHINXOPTS) $(O)
	@echo
	@echo "Spell check complete; look for any errors in $(BUILDDIR)/spelling/output.txt."

# LaTeX/PDF output
latexpdf: doxygen
	@$(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex" $(SPHINXOPTS) $(O)
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

# JSON output for search
json:
	@$(SPHINXBUILD) -b json "$(SOURCEDIR)" "$(BUILDDIR)/json" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished; the JSON files are in $(BUILDDIR)/json."

# Install documentation dependencies
install-deps:
	pip install -r requirements.txt

# Quick build (no Doxygen regeneration)
quick:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
	@echo
	@echo "Quick build finished. The HTML pages are in $(BUILDDIR)/html."

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)