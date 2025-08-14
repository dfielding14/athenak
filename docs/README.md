# AthenaK Documentation

Complete documentation system for the AthenaK astrophysical simulation framework.

## ✅ Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Build documentation
make html

# View documentation
open build/html/index.html
```

## 📁 Structure

```
docs/
├── source/              # Documentation source files
│   ├── index.rst       # Main landing page
│   ├── quickstart.md   # Getting started guide
│   ├── building.md     # Build instructions
│   ├── configuration.md # Input file guide
│   ├── running.md      # Execution guide
│   ├── modules/        # Module documentation
│   │   ├── index.md   # Module overview
│   │   ├── mesh.md    # Mesh infrastructure
│   │   ├── hydro.md   # Hydrodynamics
│   │   └── mhd.md     # Magnetohydrodynamics
│   ├── migration/      # Athena++ migration guide
│   ├── flowcharts/     # System diagrams
│   └── reference/      # File-by-file docs
├── build/              # Generated HTML (git-ignored)
├── requirements.txt    # Python dependencies
├── Makefile           # Build commands
└── Doxyfile          # API documentation config
```

## 🛠️ Building Documentation

### Standard Build
```bash
make html           # Build HTML docs
make clean         # Clean build directory
make linkcheck     # Check for broken links
```

### Development Mode
```bash
make livehtml      # Auto-rebuild on changes
# Open http://localhost:8000
```

### Full Build (with optional features)
```bash
# Install optional dependencies
brew install doxygen graphviz

# Install Python packages
pip install breathe

# Build with API docs
make html-full
```

## 📚 Documentation Coverage

- ✅ **User Guide**: Building, configuring, running simulations
- ✅ **Module Guides**: Detailed documentation for each physics module
- ✅ **Architecture**: System flowcharts and design
- ✅ **File Reference**: Every source file documented
- ✅ **Migration Guide**: Porting from Athena++ to AthenaK
- ✅ **API Reference**: Doxygen integration (optional)

## 🔍 Key Features

- **Responsive Design**: Works on desktop and mobile
- **Dark Mode**: Automatic theme switching
- **Search**: Built-in full-text search
- **Navigation**: Hierarchical sidebar navigation
- **Cross-references**: Linked documentation throughout
- **Code Examples**: Syntax-highlighted code blocks

## 📝 Adding Documentation

### New Module Documentation
1. Create `source/modules/your_module.md`
2. Add to `source/modules/index.md`
3. Follow the template in `source/_templates/module_template.md`

### New Pages
1. Create markdown/RST file in appropriate directory
2. Add to relevant `index` file's toctree
3. Rebuild with `make html`

## ⚠️ Troubleshooting

### Missing Dependencies
```bash
pip install -r requirements.txt
```

### Doxygen Not Found
```bash
# macOS
brew install doxygen

# Ubuntu/Debian
sudo apt-get install doxygen
```

### Build Warnings
- Warnings about missing cross-references are normal for stub files
- Use `make html` (not `make html -W`) to build despite warnings

## 🚀 CI/CD

Documentation is automatically built via GitHub Actions:
- See `.github/workflows/docs.yml` for configuration
- Builds on push to main branch
- Deploys to GitHub Pages (if configured)

## 📊 Verification

Run the verification script to check documentation integrity:
```bash
./verify_docs.sh
```

This checks:
- All key pages exist
- Navigation links work
- HTML files are generated

## 🤝 Contributing

1. Edit source files in `source/`
2. Test with `make html`
3. Verify with `./verify_docs.sh`
4. Submit pull request

## 📖 Generated Documentation

After building, documentation is available at:
- Local: `build/html/index.html`
- GitHub Pages: (if deployed)

## Status: ✅ READY

The documentation system is fully functional with:
- Working navigation
- All key pages created
- Build system tested
- Search functionality
- Responsive design