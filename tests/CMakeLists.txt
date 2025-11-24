cmake_minimum_required(VERSION 3.11)

project(gtest)

find_package(GTest REQUIRED)

add_executable(unit_tests)

target_sources(unit_tests PRIVATE
    unit/unit_test.cpp
)

target_include_directories(unit_tests PRIVATE 
    ${GTEST_INCLUDE_DIRS} 
    ../include 
    unit
)

enable_testing()

target_compile_features(unit_tests PUBLIC cxx_std_20)
target_link_libraries(unit_tests ${GTEST_BOTH_LIBRARIES})

add_test(NAME unit_tests COMMAND unit_tests)