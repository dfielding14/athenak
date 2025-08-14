#!/bin/bash
# Comprehensive documentation test for AthenaK

echo "========================================="
echo "     AthenaK Documentation Test Suite    "
echo "========================================="
echo

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counters
PASSED=0
FAILED=0
WARNINGS=0

# Function to check test result
check_result() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $2"
        ((PASSED++))
    else
        echo -e "${RED}✗${NC} $2"
        ((FAILED++))
    fi
}

# Function for warnings
warning() {
    echo -e "${YELLOW}⚠${NC} $1"
    ((WARNINGS++))
}

echo "1. Checking directory structure..."
[ -f "requirements.txt" ] && check_result 0 "requirements.txt found" || check_result 1 "requirements.txt missing"
[ -f "Makefile" ] && check_result 0 "Makefile found" || check_result 1 "Makefile missing"
[ -d "source" ] && check_result 0 "source directory found" || check_result 1 "source directory missing"
[ -f "Doxyfile" ] && check_result 0 "Doxyfile found" || check_result 1 "Doxyfile missing"
echo

echo "2. Checking Python dependencies..."
pip install -q -r requirements.txt 2>/dev/null
check_result $? "Python dependencies installed"
echo

echo "3. Checking documentation source files..."
[ -f "source/index.rst" ] && check_result 0 "Main index found" || check_result 1 "Main index missing"
[ -f "source/conf.py" ] && check_result 0 "Sphinx config found" || check_result 1 "Sphinx config missing"
[ -d "source/modules" ] && check_result 0 "Modules directory found" || check_result 1 "Modules directory missing"
[ -d "source/migration" ] && check_result 0 "Migration guide found" || check_result 1 "Migration guide missing"
[ -d "source/flowcharts" ] && check_result 0 "Flowcharts found" || check_result 1 "Flowcharts missing"
echo

echo "4. Checking key documentation files..."
[ -f "source/quickstart.md" ] && check_result 0 "Quickstart guide exists" || check_result 1 "Quickstart missing"
[ -f "source/modules/mesh.md" ] && check_result 0 "Mesh module docs exist" || check_result 1 "Mesh docs missing"
[ -f "source/modules/hydro.md" ] && check_result 0 "Hydro module docs exist" || check_result 1 "Hydro docs missing"
[ -f "source/migration/index.md" ] && check_result 0 "Migration guide exists" || check_result 1 "Migration guide missing"
[ -f "source/reference/index.md" ] && check_result 0 "File reference exists" || check_result 1 "File reference missing"
echo

echo "5. Checking optional tools..."
if command -v doxygen &> /dev/null; then
    check_result 0 "Doxygen installed"
else
    warning "Doxygen not installed (API docs limited)"
    echo "   Install with: brew install doxygen"
fi

if command -v dot &> /dev/null; then
    check_result 0 "Graphviz installed"
else
    warning "Graphviz not installed (diagrams limited)"
    echo "   Install with: brew install graphviz"
fi
echo

echo "6. Building documentation..."
make clean > /dev/null 2>&1
make html > build.log 2>&1
BUILD_RESULT=$?

if [ $BUILD_RESULT -eq 0 ]; then
    check_result 0 "Documentation built successfully"
    
    # Count generated files
    HTML_COUNT=$(find build/html -name "*.html" 2>/dev/null | wc -l | tr -d ' ')
    echo "   Generated $HTML_COUNT HTML files"
    
    # Check for key pages
    [ -f "build/html/index.html" ] && echo "   ✓ Main index page"
    [ -f "build/html/quickstart.html" ] && echo "   ✓ Quickstart guide"
    [ -f "build/html/modules/index.html" ] && echo "   ✓ Module documentation"
    [ -f "build/html/migration/index.html" ] && echo "   ✓ Migration guide"
else
    check_result 1 "Documentation build failed"
    echo "   Check build.log for details"
fi
echo

echo "7. Checking for broken internal links..."
if [ -f "build/html/index.html" ]; then
    # Simple check for common broken link patterns
    BROKEN_LINKS=$(grep -r 'href="#"' build/html/*.html 2>/dev/null | wc -l | tr -d ' ')
    if [ "$BROKEN_LINKS" -gt 0 ]; then
        warning "Found $BROKEN_LINKS potential broken links"
    else
        check_result 0 "No obvious broken links"
    fi
else
    echo "   Skipping (no HTML output)"
fi
echo

echo "========================================="
echo "              TEST SUMMARY               "
echo "========================================="
echo -e "Passed:   ${GREEN}$PASSED${NC}"
echo -e "Failed:   ${RED}$FAILED${NC}"
echo -e "Warnings: ${YELLOW}$WARNINGS${NC}"
echo

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}SUCCESS!${NC} Documentation is ready."
    echo
    echo "View documentation:"
    echo "  open build/html/index.html"
    echo
    echo "For live development:"
    echo "  make livehtml"
    echo "  # Then open http://localhost:8000"
else
    echo -e "${RED}FAILURE!${NC} Please fix the issues above."
    exit 1
fi

# Coverage report (optional)
echo
echo "Documentation coverage:"
echo "----------------------"
SRC_FILES=$(find ../src -name "*.cpp" -o -name "*.hpp" 2>/dev/null | wc -l | tr -d ' ')
echo "Source files in src/: $SRC_FILES"

if [ -d "source/reference" ]; then
    DOC_REFS=$(grep -r "\.cpp\|\.hpp" source/reference/ 2>/dev/null | wc -l | tr -d ' ')
    echo "File references in docs: ~$DOC_REFS"
fi

exit 0