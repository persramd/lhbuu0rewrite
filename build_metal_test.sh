#!/bin/bash

# Build script for Metal grid density test
# M4 Mac - Metal compute shaders

set -e

echo "=== Building Metal Grid Density Test ==="
echo ""

# Paths
SRC_DIR="../src"
INCLUDE_DIR="../include"
TEST_DIR="../tests"
SHADER_DIR="../shaders"

# Check if we're running on macOS
if [[ "$(uname)" != "Darwin" ]]; then
    echo "Error: This script requires macOS for Metal support"
    exit 1
fi

# Compile Metal shaders
echo "Compiling Metal shaders..."
xcrun -sdk macosx metal -c ${SHADER_DIR}/grid_density.metal -o grid_density.air
xcrun -sdk macosx metallib grid_density.air -o grid_density.metallib
echo "  -> grid_density.metallib"
echo ""

# Build test executable
echo "Compiling source files..."

# Grid implementation (C)
clang -c ${SRC_DIR}/grid.c -I${INCLUDE_DIR} -O3 -fopenmp -march=native -o grid.o
echo "  -> grid.o"

# Utils (C)
clang -c ${SRC_DIR}/utils.c -I${INCLUDE_DIR} -O3 -o utils.o
echo "  -> utils.o"

# GPU compute (Objective-C with Metal)
clang -c ${SRC_DIR}/gpu_compute_metal.m -I${INCLUDE_DIR} -O3 -fobjc-arc \
    -framework Metal -framework Foundation -o gpu_compute_metal.o
echo "  -> gpu_compute_metal.o"

# Test program (C)
clang -c ${TEST_DIR}/test_metal_density.c -I${INCLUDE_DIR} -O3 -o test_metal_density.o
echo "  -> test_metal_density.o"

echo ""
echo "Linking executable..."

# Link everything
clang -o test_metal_density \
    test_metal_density.o \
    grid.o \
    utils.o \
    gpu_compute_metal.o \
    -framework Metal -framework Foundation -fopenmp -lm

echo "  -> test_metal_density"
echo ""
echo "=== Build Complete ==="
echo ""
echo "To run: ./test_metal_density"
echo ""
