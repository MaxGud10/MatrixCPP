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

target_compile_features(matrix PUBLIC cxx_std_20)

target_link_libraries(matrix PRIVATE
    $<$<CONFIG:Debug>:
        -fsanitize=address
        -fsanitize=undefined>
    ${GTEST_BOTH_LIBRARIES}
)