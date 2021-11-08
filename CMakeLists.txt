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


# ====== Users' settings ======

set(WITH_OPENBLAS_INC "" CACHE PATH "Path to an existing OpenBLAS header location.")
set(WITH_OPENBLAS_LIB "" CACHE PATH "Path to an existing OpenBLAS static library location.")


# ====== External libraries ======

set(LIB_EXT ${CMAKE_BINARY_DIR}/external)
set(BUILD_EXT ${LIB_EXT}/builds)

# ------ OpenBLAS 0.3.15 ------

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
        URL https://github.com/xianyi/OpenBLAS/releases/download/v0.3.15/OpenBLAS-0.3.15.tar.gz
        URL_HASH SHA256=30a99dec977594b387a17f49904523e6bc8dd88bd247266e83485803759e4bbe
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

# ------ Boost ------

# ./bcp math/quadrature math/special_functions/bessel.hpp ../../../../../../../boost_1_76_0_selected --boost=../../../../../../../boost_1_76_0

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

set(YAML_CPP_ROOT ${BUILD_EXT}/yaml-cpp)
ExternalProject_Add(external_yaml_cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG 9a3624205e8774953ef18f57067b3426c1c5ada6
    PREFIX ${YAML_CPP_ROOT}
    CMAKE_ARGS -DYAML_BUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=${YAML_CPP_ROOT}
)

set(YAML_LIB
"${YAML_CPP_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}yaml-cpp${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(YAML_INC "${YAML_CPP_ROOT}/include")

# ------ Documentation ------

add_custom_target(doc COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/doc && ./generate_docs.sh)


# ====== Strata ======

get_filename_component(OPENBLAS_LIB_DIR ${OPENBLAS_LIB} DIRECTORY)
get_filename_component(AmosBessel_LIB_DIR ${AmosBessel_LIB} DIRECTORY)
get_filename_component(YAML_LIB_DIR ${YAML_LIB} DIRECTORY)

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH ${OPENBLAS_LIB_DIR} ${AmosBessel_LIB_DIR} ${YAML_LIB_DIR} ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

SET(INCLUDE_DIR
    ${OPENBLAS_INC}
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
target_link_libraries(strata ${OPENBLAS_LIB} ${AmosBessel_LIB} ${YAML_LIB})
target_include_directories(strata PRIVATE ${INCLUDE_DIR})
if ((WITH_OPENBLAS_INC STREQUAL "") OR
    (WITH_OPENBLAS_LIB STREQUAL ""))
    add_dependencies(strata external_openblas external_amosbessel external_yaml_cpp external_boost doc)
else ()
    add_dependencies(strata external_amosbessel external_yaml_cpp external_boost doc)
endif()

install(TARGETS strata DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/inc/"
        DESTINATION "include"
        FILES_MATCHING
        PATTERN "*.hpp"
)

# ------ Test files ------

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)



