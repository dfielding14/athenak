#!/bin/bash
# Test script for AthenaK documentation build

echo "=== AthenaK Documentation Build Test ==="
echo

# Check Python
echo "1. Checking Python..."
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 not found"
    exit 1
fi
echo "✅ Python3 found: $(python3 --version)"
echo

# Install requirements
echo "2. Installing Python dependencies..."
pip install -q -r requirements.txt
if [ $? -eq 0 ]; then
    echo "✅ Dependencies installed"
else
    echo "❌ Failed to install dependencies"
    exit 1
fi
echo

# Check optional dependencies
echo "3. Checking optional dependencies..."
if command -v doxygen &> /dev/null; then
    echo "✅ Doxygen found: $(doxygen --version | head -1)"
else
    echo "⚠️  Doxygen not found (optional - install with: brew install doxygen)"
fi

if command -v dot &> /dev/null; then
    echo "✅ Graphviz found"
else
    echo "⚠️  Graphviz not found (optional - install with: brew install graphviz)"
fi
echo

# Clean build directory
echo "4. Cleaning build directory..."
make clean > /dev/null 2>&1
echo "✅ Build directory cleaned"
echo

# Build documentation
echo "5. Building documentation..."
make html > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Documentation built successfully"
else
    echo "❌ Build failed"
    echo "Run 'make html' to see detailed error"
    exit 1
fi
echo

# Check output
echo "6. Checking output..."
if [ -f "build/html/index.html" ]; then
    echo "✅ HTML files generated"
    file_count=$(find build/html -name "*.html" | wc -l | tr -d ' ')
    echo "   Generated $file_count HTML files"
else
    echo "❌ No HTML output found"
    exit 1
fi
echo

echo "=== Test Complete ==="
echo
echo "Documentation successfully built!"
echo "View it with: open build/html/index.html"
echo
echo "For development with live reload:"
echo "  make livehtml"
echo "  Then open http://localhost:8000"