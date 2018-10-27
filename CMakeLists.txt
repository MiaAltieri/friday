CMAKE_MINIMUM_REQUIRED(VERSION 3.8)

# name of the project is VISION
PROJECT(friday)
set (CMAKE_CXX_FLAGS "-fPIC")
set (CMAKE_C_FLAGS "-fPIC")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

FIND_PACKAGE(PythonInterp 3 REQUIRED)
FIND_PACKAGE(PythonLibs 3 REQUIRED)

if (PYTHONINTERP_FOUND)
    message("Python found")
else()
    message("Python not found")
endif()

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "CURRENT DIRECTORY = ${CMAKE_CURRENT_SOURCE_DIR}")

# enable installing dependencies
option(INSTALL_DEPENDENCIES
        "Install project dependencies"
        ON)

INCLUDE(modules/htslib.cmake)

# pybind11 to interface
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/modules/dep/pybind11)
pybind11_add_module(FRIDAY ${CMAKE_CURRENT_SOURCE_DIR}/modules/src/pybind_api.cpp)
# pybind11_add_module(FRIDAY modules/src/local_reassembly/bam_handler.cpp)

add_dependencies(FRIDAY htslib)

# add all the external libraries
target_link_libraries(FRIDAY PRIVATE z)
target_link_libraries(FRIDAY PRIVATE bz2)
target_link_libraries(FRIDAY PRIVATE curl)
target_link_libraries(FRIDAY PRIVATE lzma)
target_link_libraries(FRIDAY PRIVATE ${HTSLIB_SRC_DIR}/libhts.a)

#-------------------
# Test
#-------------------
#ENABLE_TESTING()
#ADD_SUBDIRECTORY(modules/unit_tests)

