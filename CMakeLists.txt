CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

# name of the project is VISION
PROJECT(friday)
set (CMAKE_CXX_FLAGS "-fPIC")
set (CMAKE_C_FLAGS "-fPIC")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable installing dependencies
option(INSTALL_DEPENDENCIES
        "Install project dependencies"
        ON)

INCLUDE(modules/htslib.cmake)
INCLUDE(modules/pybind11.cmake)

# pybind11 to interface
pybind11_add_module(FRIDAY ${CMAKE_CURRENT_SOURCE_DIR}/modules/src/pybind_api.cpp)

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

