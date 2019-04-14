include(FetchContent)
FetchContent_Declare(
        absl
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG        master
        GIT_PROGRESS   TRUE
)

FetchContent_GetProperties(absl)
if(NOT absl_POPULATED)
    FetchContent_Populate(absl)
    add_subdirectory(${absl_SOURCE_DIR})
endif()