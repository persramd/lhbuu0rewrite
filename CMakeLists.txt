cmake_minimum_required(VERSION 3.10)
project(LHBUU C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find OpenMP
find_package(OpenMP)

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

# Main executable
add_executable(lhbuu ${SOURCES})
target_link_libraries(lhbuu m)  # Link math library

if(OpenMP_C_FOUND)
    target_link_libraries(lhbuu OpenMP::OpenMP_C)
endif()

# Install
install(TARGETS lhbuu DESTINATION bin)

# Tests (optional)
option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
