include(FetchContent)
FetchContent_Declare(
        boost_cmake
        GIT_REPOSITORY https://github.com/Orphis/boost-cmake.git
        GIT_TAG        master
)

FetchContent_GetProperties(boost_cmake)
if(NOT boost_cmake_POPULATED)
    FetchContent_Populate(boost_cmake)
    add_subdirectory(${boost_cmake_SOURCE_DIR})
endif()