cmake_minimum_required(VERSION 3.10)
project(LHBUU C)

# Enable Objective-C only on macOS (for Metal)
if(APPLE)
    enable_language(OBJC)
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find OpenMP
find_package(OpenMP)

# GPU compute backend detection
option(USE_GPU "Enable GPU acceleration" ON)

if(USE_GPU)
    if(APPLE)
        # Metal GPU compute on macOS
        option(USE_METAL "Enable Metal GPU acceleration" ON)
        if(USE_METAL)
            find_library(METAL_LIBRARY Metal)
            find_library(FOUNDATION_LIBRARY Foundation)
            if(METAL_LIBRARY AND FOUNDATION_LIBRARY)
                set(GPU_ENABLED TRUE)
                set(GPU_BACKEND "Metal")
                add_compile_definitions(USE_GPU_COMPUTE=1)
                add_compile_definitions(USE_METAL=1)
                message(STATUS "GPU backend: Metal (Apple)")
            else()
                message(WARNING "Metal requested but not found - CPU fallback only")
                set(GPU_ENABLED FALSE)
            endif()
        endif()
    elseif(UNIX AND NOT APPLE)
        # HIP/ROCm GPU compute on Linux (placeholder for future)
        set(GPU_ENABLED FALSE)
        set(GPU_BACKEND "HIP")
        message(STATUS "GPU backend: HIP/ROCm (not yet implemented - CPU fallback)")
    else()
        set(GPU_ENABLED FALSE)
        message(STATUS "GPU backend: None (CPU fallback only)")
    endif()
else()
    set(GPU_ENABLED FALSE)
    message(STATUS "GPU acceleration disabled (CPU only)")
endif()

# Compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O3 -march=native")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 -DDEBUG")

if(OpenMP_C_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/main.c
    src/config.c
    src/grid.c
    src/grid_tag.c
    src/potential.c
    src/potential_mdyi.c
    src/potential_gradient.c
    src/collision.c
    src/collision_full.c
    src/lorentz.c
    src/integrator.c
    src/initialize.c
    src/initialize_nuclei.c
    src/diagnostics.c
    src/utils.c
    src/performance.c
    src/neighbors.c
    src/cross_sections.c
    src/medium_corrections.c
    src/pauli_blocking.c
    src/form_factor.c
    src/drho_dr.c
    src/mdyi_cache.c
    src/half_impulse.c
    src/arch_optimize.c
    src/coulomb_fft.c
    src/inside_sigma.c
    src/p_avg.c
)

# Add GPU compute sources if enabled
if(GPU_ENABLED)
    list(APPEND SOURCES src/gpu_compute_metal.m)
else()
    # CPU-only fallback implementation
    list(APPEND SOURCES src/gpu_compute_cpu.c)
endif()

# Main executable
add_executable(lhbuu ${SOURCES})
target_link_libraries(lhbuu m)  # Link math library

if(OpenMP_C_FOUND)
    target_link_libraries(lhbuu OpenMP::OpenMP_C)
endif()

if(GPU_ENABLED)
    target_link_libraries(lhbuu ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
endif()

# Install
install(TARGETS lhbuu DESTINATION bin)

# GPU compute test executable
if(GPU_ENABLED)
    add_executable(test_gpu_compute
        tests/test_gpu_compute.c
        src/grid.c
        src/utils.c
        src/gpu_compute_metal.m
    )
    target_link_libraries(test_gpu_compute m ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
    if(OpenMP_C_FOUND)
        target_link_libraries(test_gpu_compute OpenMP::OpenMP_C)
    endif()
endif()

# Tests (optional)
option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
