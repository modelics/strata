cmake_minimum_required(VERSION 3.10.2)

include(ExternalProject)
project(StrataTest)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-march=native)
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

message(STATUS ${CMAKE_BINARY_DIR})
message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})

# Add the Strata library to the project
# add_library(strata STATIC IMPORTED)
# set_property(TARGET strata PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
# set_property(TARGET strata PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/../build/libstrata.a)

set(STRATA_LIB
        "${CMAKE_CURRENT_SOURCE_DIR}/../build/${CMAKE_STATIC_LIBRARY_PREFIX}strata${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(STRATA_INC "${CMAKE_CURRENT_SOURCE_DIR}/../inc")

SET(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# ------ Executables ------

add_executable(testMGF ${SRC_DIR}/testMGF.cpp)
target_include_directories(testMGF PRIVATE ${INC_DIR} ${STRATA_INC})
target_link_libraries(testMGF ${STRATA_LIB})



