# Author: Shashwat Sharma

set(STRATA_LIB
        "${CMAKE_CURRENT_SOURCE_DIR}/../build/${CMAKE_SHARED_LIBRARY_PREFIX}strata${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(STRATA_INC "${CMAKE_CURRENT_SOURCE_DIR}/../inc")

SET(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# ------ Executables ------

add_executable(dev_sandbox ${SRC_DIR}/dev_sandbox.cpp)
target_include_directories(dev_sandbox PRIVATE ${INC_DIR} ${STRATA_INC})
target_link_libraries(dev_sandbox ${STRATA_LIB})

add_executable(testMGF ${SRC_DIR}/testMGF.cpp)
target_include_directories(testMGF PRIVATE ${INC_DIR} ${STRATA_INC})
target_link_libraries(testMGF ${STRATA_LIB})

add_executable(testInterp ${SRC_DIR}/testInterp.cpp)
target_include_directories(testInterp PRIVATE ${INC_DIR} ${STRATA_INC})
target_link_libraries(testInterp ${STRATA_LIB})

add_executable(testDCIM ${SRC_DIR}/testDCIM.cpp)
target_include_directories(testDCIM PRIVATE ${INC_DIR} ${STRATA_INC})
target_link_libraries(testDCIM ${STRATA_LIB})

add_executable(testSingularity ${SRC_DIR}/testSingularity.cpp)
target_include_directories(testSingularity PRIVATE ${INC_DIR} ${STRATA_INC})
target_link_libraries(testSingularity ${STRATA_LIB})


