cmake_minimum_required(VERSION 3.11)
project(Matrix)

find_package(GTest REQUIRED)
enable_testing()

add_executable(matrix)

target_include_directories(matrix PUBLIC ${GTEST_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)

target_sources(matrix PRIVATE
    src/double_compare.cpp
    src/main.cpp
)

target_compile_options(matrix PRIVATE 
    $<$<CONFIG:DEBUG>:
        -Wall
        -Wextra
        -Wpedantic
        -g
        -O0
        -DDEBUG
        -Wmissing-declarations
        -Wcast-align
        -Wunused
        -Wshadow
        -Wsign-conversion
        -fsanitize=address
        -fsanitize=undefined>
)

add_subdirectory(tests)

# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#   target_compile_options(unit_tests PRIVATE -g -O0 -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined)
#   target_link_options(unit_tests    PRIVATE -fsanitize=address -fsanitize=undefined)
# endif()

target_compile_features(matrix PUBLIC cxx_std_20)

target_link_libraries(matrix PRIVATE
    $<$<CONFIG:Debug>:
        -fsanitize=address
        -fsanitize=undefined>
    ${GTEST_BOTH_LIBRARIES}
)