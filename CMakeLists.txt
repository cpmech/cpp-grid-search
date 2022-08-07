cmake_minimum_required(VERSION 3.16)
project(Project)
enable_testing()

# options ################################################################################

option(A_MAKE_VERBOSE "Show additional messages during compilation/link?" OFF)
option(A_MAKE_ALL_WARNINGS "Make with all warnings (-Wall)" ON)
option(A_MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g)" ON)
option(A_MAKE_OPTIMIZED "Make optimized (-O2)" OFF)
option(A_MAKE_PROFILING "Make with profiling information (-pg)" OFF)

if(A_MAKE_VERBOSE)
    set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif(A_MAKE_VERBOSE)

if(A_MAKE_ALL_WARNINGS)
    add_definitions(-Wall)
endif(A_MAKE_ALL_WARNINGS)

if(A_MAKE_DEBUG_SYMBOLS)
    add_definitions(-g)
endif(A_MAKE_DEBUG_SYMBOLS)

if(A_MAKE_PROFILING)
    set(CMAKE_CXX_FLAGS -pg)
endif(A_MAKE_PROFILING)

if(A_MAKE_OPTIMIZED)
    add_definitions(-O2)
endif(A_MAKE_OPTIMIZED)

# cpp flags ##############################################################################

set(CMAKE_CXX_STANDARD 17)

# library ################################################################################

set(TRILIB_NAME "trilib")

add_library(${TRILIB_NAME} SHARED
    src/triangle/interface_triangle.c
    src/triangle/triangle.c
    src/triangle/tricall_report.c
)

# main ###################################################################################

# subdirectories #############################################################################

add_subdirectory(src)
add_subdirectory(examples)
