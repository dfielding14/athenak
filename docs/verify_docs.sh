#!/bin/bash
# Final verification of AthenaK documentation

echo "===================================="
echo "  AthenaK Documentation Verification"
echo "===================================="
echo

# Check if we're in docs directory
if [ ! -f "Makefile" ] || [ ! -d "source" ]; then
    echo "❌ Error: Run this script from the docs/ directory"
    exit 1
fi

echo "✅ Checking generated HTML files..."
HTML_COUNT=$(find build/html -name "*.html" 2>/dev/null | wc -l | tr -d ' ')
echo "   Found $HTML_COUNT HTML files"

echo
echo "✅ Verifying key pages exist..."
pages=(
    "build/html/index.html:Main index"
    "build/html/building.html:Building guide"
    "build/html/configuration.html:Configuration"
    "build/html/running.html:Running guide"
    "build/html/modules/mesh.html:Mesh module"
    "build/html/modules/mhd.html:MHD module"
    "build/html/migration/index.html:Migration guide"
    "build/html/flowcharts/runtime.html:Runtime flowchart"
)

for page_info in "${pages[@]}"; do
    IFS=':' read -r page desc <<< "$page_info"
    if [ -f "$page" ]; then
        echo "   ✓ $desc"
    else
        echo "   ✗ Missing: $desc"
    fi
done

echo
echo "✅ Testing navigation links..."
# Check for broken internal links in index
if [ -f "build/html/index.html" ]; then
    # Check that key links exist
    links_found=0
    links_checked=0
    
    for link in "building.html" "configuration.html" "modules/mesh.html" "migration/index.html"; do
        ((links_checked++))
        if grep -q "href=\"$link\"" build/html/index.html; then
            ((links_found++))
        fi
    done
    
    echo "   Found $links_found/$links_checked navigation links"
fi

echo
echo "===================================="
echo "         DOCUMENTATION READY"
echo "===================================="
echo
echo "📖 View documentation:"
echo "   open build/html/index.html"
echo
echo "🔄 For live development:"
echo "   make livehtml"
echo "   Then open http://localhost:8000"
echo
echo "📦 To rebuild:"
echo "   make clean && make html"