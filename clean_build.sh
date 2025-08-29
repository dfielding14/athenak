#!/bin/bash
# Script to safely clean build directories and handle NFS lock files

BUILD_DIR="${1:-build_turb}"

echo "Cleaning build directory: $BUILD_DIR"

# First, try normal removal
if rm -rf "$BUILD_DIR" 2>/dev/null; then
    echo "Successfully removed $BUILD_DIR"
    exit 0
fi

# If that fails, check for NFS lock files
echo "Normal removal failed, checking for NFS lock files..."

# Find and list any NFS lock files
NFS_FILES=$(find "$BUILD_DIR" -name ".nfs*" 2>/dev/null)

if [ -n "$NFS_FILES" ]; then
    echo "Found NFS lock files:"
    echo "$NFS_FILES"
    
    # Try to find processes holding these files
    echo "Checking for processes holding these files..."
    
    # Look for common build-related processes
    PROCS=$(ps aux | grep -E "(cpptools|cmake|make|ninja|ccache)" | grep -v grep)
    
    if [ -n "$PROCS" ]; then
        echo "Found potentially problematic processes:"
        echo "$PROCS"
        
        # Extract PIDs of cpptools processes (most common culprit)
        CPPTOOLS_PIDS=$(echo "$PROCS" | grep cpptools | awk '{print $2}')
        
        if [ -n "$CPPTOOLS_PIDS" ]; then
            echo "Killing cpptools processes: $CPPTOOLS_PIDS"
            kill -9 $CPPTOOLS_PIDS 2>/dev/null
            sleep 2
        fi
    fi
    
    # Try removal again
    echo "Attempting removal again..."
    if rm -rf "$BUILD_DIR" 2>/dev/null; then
        echo "Successfully removed $BUILD_DIR after killing processes"
        exit 0
    fi
    
    # If still failing, try moving the directory aside
    echo "Still unable to remove, trying to move directory aside..."
    BACKUP_DIR="${BUILD_DIR}_old_$(date +%s)"
    if mv "$BUILD_DIR" "$BACKUP_DIR" 2>/dev/null; then
        echo "Moved $BUILD_DIR to $BACKUP_DIR"
        echo "You can manually remove $BACKUP_DIR later"
        exit 0
    fi
    
    echo "ERROR: Unable to clean build directory"
    echo "Try closing any IDEs or editors that might be accessing these files"
    exit 1
else
    echo "No NFS lock files found, but directory still can't be removed"
    echo "Check permissions or if directory is mounted"
    exit 1
fi