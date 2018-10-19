CMAKE_MINIMUM_REQUIRED(VERSION 3.8)

# name of the project is VISION
PROJECT(vision)


SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

FIND_PACKAGE(PythonInterp 3 REQUIRED)
FIND_PACKAGE(PythonLibs 3 REQUIRED)

if (PYTHONINTERP_FOUND)
    message("Python found")
else()
    message("Python not found")
endif()

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
#target_link_libraries(bam_handler ${HTSLIB_SRC_DIR}/libhts.a)

#pybind11_add_module(vision ${CMAKE_CURRENT_SOURCE_DIR}/modules/src/pybind_src.cpp)

## FIND AND SET THE LIBRARY FILES
#SET(LIBRARY_FILES
#    ${CMAKE_CURRENT_SOURCE_DIR}/modules/src/bam_handler.cpp
#)
#
## SET UP THE LIBRARY
#MESSAGE("LIBRARY_FILES FILES: ${LIBRARY_FILES}")
#ADD_LIBRARY(vision_lib ${LIBRARY_FILES})
#
#
#
#ADD_EXECUTABLE(vision main.cpp)
#
#
### ADD THE MODULE DIRECTORY
##INCLUDE_DIRECTORIES(modules)
#
## LINK TO ZLIB LIBRARY
#add_dependencies(vision htslib)
#add_dependencies(vision pybind11)
#TARGET_LINK_LIBRARIES(vision vision_lib)
#TARGET_LINK_LIBRARIES(vision z)
#TARGET_LINK_LIBRARIES(vision ${HTSLIB_SRC_DIR}/libhts.a)


#-------------------
# Test
#-------------------
#ENABLE_TESTING()
#ADD_SUBDIRECTORY(modules/unit_tests)

