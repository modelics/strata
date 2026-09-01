# Author: Shashwat Sharma

# Copyright 2021 Shashwat Sharma and Piero Triverio

# This file is part of Strata.

# Strata is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Strata is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Strata.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.10.2)
include(ExternalProject)

project(Strata)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optionally build documentation
option(BUILD_DOC "Build the documentation" ON)

# ====== Compiler settings ======

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-march=native)
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

# ------------------------------------------------------------------------------
# Fast complex arithmetic (opt-in numerical relaxation)
#
# When ON, the strata library is compiled with -fcx-limited-range. That flag makes
# the compiler evaluate complex multiplication/division with the simple textbook
# formulas inline, instead of calling libgcc's __muldc3 / __divdc3 helpers, which
# implement the full C99 Annex G rules (extra inf/NaN/overflow scaling).
#
# Why it matters: in low-frequency multilayer Green's function (MGF) evaluation,
# complex division via __divdc3 was measured (perf) at ~14% of total runtime.
# Enabling this flag cut the "Convolution matrix initialization" wall time by ~17%
# on a 1 kHz microstrip test (1935 s -> 1609 s), with the solver taking the same
# number of iterations and S-parameters matching to the numerical floor.
#
# Trade-off: -fcx-limited-range does NOT robustly handle complex values near
# overflow/underflow or involving inf/NaN. It is safe when magnitudes are well
# scaled (as in the validated MGF cases). Turn it OFF if you ever see NaN/inf in
# results at extreme scales.
#
# Scope: this is applied ONLY to the `strata` target (see target_compile_options
# below), so it never changes how a parent project such as Rebel compiles its own
# sources -- only strata's objects get the flag.
#
# NOTE for Rebel builds: Rebel compiles strata as an isolated CMake ExternalProject,
# so a -D flag given to Rebel's top-level configure does NOT propagate in here; a
# Rebel build uses this option's DEFAULT value. To enable it for a Rebel build,
# either flip the default below to ON, or add -DSTRATA_FAST_COMPLEX_ARITHMETIC=ON
# to the external_strata CMAKE_ARGS in Rebel's CMakeLists.
# ------------------------------------------------------------------------------
option(STRATA_FAST_COMPLEX_ARITHMETIC
       "Compile strata with -fcx-limited-range (faster complex mul/div; relaxes inf/NaN/overflow handling)"
       OFF)


# ====== Users' settings ======

set(WITH_OPENBLAS_INC "" CACHE PATH "Path to an existing OpenBLAS header location.")
set(WITH_OPENBLAS_LIB "" CACHE PATH "Path to an existing OpenBLAS static library location.")


# ====== External libraries ======

set(LIB_EXT ${CMAKE_BINARY_DIR}/external)
set(BUILD_EXT ${LIB_EXT}/builds)

# ====== Find MKL ======

set(MKL_INTERFACE lp64)
set(MKL_MPI openmpi)
set(MKL_THREADING gnu_thread)

find_package(MKL CONFIG QUIET PATHS $ENV{MKLROOT})

set(MKL_INC ${MKL_INCLUDE})
set(MKL_LIB MKL::MKL)

if(NOT MKL_FOUND)
	# ------ OpenBLAS 0.3.29 ------
	set(BLA_VENDOR OpenBLAS)
	if ((NOT WITH_OPENBLAS_INC STREQUAL "") AND
	    (NOT WITH_OPENBLAS_LIB STREQUAL ""))
	    set(OPENBLAS_LIB
		"${WITH_OPENBLAS_LIB}")
	    set(OPENBLAS_INC "${WITH_OPENBLAS_INC}")
	    add_custom_target(external_openblas)
	else()
	    message(STATUS "Could NOT find OpenBLAS")
	    set(OPENBLAS_ROOT ${BUILD_EXT}/openblas)
	    ExternalProject_Add(external_openblas
		URL https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29.tar.gz
		URL_HASH MD5=853a0c5c0747c5943e7ef4bbb793162d
		PREFIX ${OPENBLAS_ROOT}
		BUILD_IN_SOURCE ON
		CONFIGURE_COMMAND ""
		BUILD_COMMAND make
		INSTALL_COMMAND make PREFIX=${OPENBLAS_ROOT} install
	    )
	    set(OPENBLAS_LIB
		"${OPENBLAS_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX}")
	    set(OPENBLAS_INC "${OPENBLAS_ROOT}/include")
	endif()
endif()

# ------ Boost ------

# ./bcp math/quadrature math/special_functions/bessel.hpp ../../../../../../../boost_1_76_0_selected --boost=../../../../../../../boost_1_76_0

find_package(Boost QUIET)

if (NOT Boost_FOUND)
    set(BOOST_ROOT ${BUILD_EXT}/boost_1_76_0_selected)
    ExternalProject_Add(external_boost
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_76_0_selected.tar.gz
        URL_HASH SHA256=a3b5990412456c7a528c255132677644e8add4f3d6b50e864667813a9bd43afb
        PREFIX ${BOOST_ROOT}
        SOURCE_DIR ${LIB_EXT}/boost_1_76_0_selected
        CONFIGURE_COMMAND tar "xzf" ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_76_0_selected.tar.gz --strip-components=1 -C ${BOOST_ROOT}
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )
    set(BOOST_INC ${BOOST_ROOT})
else()
    set(BOOST_INC ${Boost_INCLUDE_DIRS})
endif()

# ------ AmosBessel ------

set(AMOSBESSEL_ROOT ${BUILD_EXT}/AmosBessel)
ExternalProject_Add(external_amosbessel
    URL ${CMAKE_CURRENT_SOURCE_DIR}/external/AmosBessel.tar.gz
    # URL_HASH SHA256=a3b7e2d7a3bc2cab5e5a69e1790c2ca72b898d6a561255f56e154927751ffae5
    PREFIX ${AMOSBESSEL_ROOT}
    SOURCE_DIR ${LIB_EXT}/AmosBessel
    # DOWNLOAD_DIR ${LIB_EXT}
    # BUILD_IN_SOURCE ON
    # CONFIGURE_COMMAND cmake
    # BUILD_COMMAND make
    # INSTALL_COMMAND ""
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${AMOSBESSEL_ROOT}
)

set(AmosBessel_LIB
    "${AMOSBESSEL_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}AmosBessel${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(AmosBessel_INC "${AMOSBESSEL_ROOT}/include/AmosBessel")

# ------ yaml-cpp 0.6.3 ------

find_package(yaml-cpp QUIET)

if (NOT yaml_cpp_FOUND)
	set(YAML_CPP_ROOT ${BUILD_EXT}/yaml-cpp)
	ExternalProject_Add(external_yaml_cpp
	    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
	    GIT_TAG 9a3624205e8774953ef18f57067b3426c1c5ada6
	    PREFIX ${YAML_CPP_ROOT}
	    CMAKE_ARGS -DYAML_BUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=${YAML_CPP_ROOT}
	)
	
	set(YAML_LIB "${YAML_CPP_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}yaml-cpp${CMAKE_SHARED_LIBRARY_SUFFIX}")
	set(YAML_INC "${YAML_CPP_ROOT}/include")
else()
	set(YAML_LIB ${yaml-cpp_INCLUDE_DIRS})
	set(YAML_INC ${yaml-cpp_LIBRARIES})	
endif()

# ------ Documentation ------

if (BUILD_DOC)
    add_custom_target(doc COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/doc && ./generate_docs.sh)
endif()

# ====== Strata ======

if(NOT MKL_FOUND)
	set(BLAS_LIB ${OPENBLAS_LIB})
	set(BLAS_INC ${OPENBLAS_INC})
else()
	add_compile_definitions(USE_MKL)
	set(BLAS_LIB ${MKL_LIB})
	set(BLAS_INC ${MKL_INC})	
endif()

get_filename_component(BLAS_LIB_DIR ${BLAS_LIB} DIRECTORY)
get_filename_component(AmosBessel_LIB_DIR ${AmosBessel_LIB} DIRECTORY)
get_filename_component(YAML_LIB_DIR ${YAML_LIB} DIRECTORY)

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH ${BLAS_LIB_DIR} ${AmosBessel_LIB_DIR} ${YAML_LIB_DIR} ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

SET(INCLUDE_DIR
    ${BLAS_INC}
    ${BOOST_INC}
    ${AmosBessel_INC}
    ${YAML_INC}
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

SET(MGF_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/DCIM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/layers.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MGF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/quasistatic_MGF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/spectral_MGF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_interface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/coordinate_system.cpp
)

add_library(strata SHARED ${MGF_SRC})
target_link_libraries(strata ${BLAS_LIB} ${AmosBessel_LIB} ${YAML_LIB})
target_include_directories(strata PRIVATE ${INCLUDE_DIR})

# Fast complex arithmetic: apply -fcx-limited-range to strata objects only.
# See the STRATA_FAST_COMPLEX_ARITHMETIC option (with rationale) near the top.
if (STRATA_FAST_COMPLEX_ARITHMETIC)
    target_compile_options(strata PRIVATE -fcx-limited-range)
    message(STATUS "strata: -fcx-limited-range ENABLED (fast complex arithmetic; relaxed inf/NaN handling)")
else ()
    message(STATUS "strata: -fcx-limited-range disabled (set -DSTRATA_FAST_COMPLEX_ARITHMETIC=ON to enable)")
endif ()

# Add OpenBLAS dependency if externally downloaded
if (((WITH_OPENBLAS_INC STREQUAL "") OR
    (WITH_OPENBLAS_LIB STREQUAL "")) AND NOT MKL_FOUND)
    add_dependencies(strata external_openblas external_amosbessel)
else ()
    add_dependencies(strata external_amosbessel)
endif()

# Add YAML-CPP dependency if externally downloaded
if (NOT yaml_cpp_FOUND)
    add_dependencies(strata external_yaml_cpp)
endif()

# Add Boost dependency if externally downloaded
if (NOT Boost_FOUND)
    add_dependencies(strata external_boost)
endif()

install(TARGETS strata DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/inc/"
        DESTINATION "include"
        FILES_MATCHING
        PATTERN "*.hpp"
)

# ------ Test files ------

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)

