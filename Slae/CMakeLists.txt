cmake_minimum_required(VERSION 3.24)
project(slae)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(slae_lib STATIC 
    Source/Slae.cxx
)
target_include_directories(slae_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Source)

add_executable(${PROJECT_NAME} Charts/Cxx/main.cxx)

target_link_libraries(${PROJECT_NAME} PRIVATE slae_lib)

enable_testing()
find_package(GTest REQUIRED)
add_executable(slae_test 
    Test/operations.cxx
    Test/sweep_method.cxx
)
target_link_libraries(slae_test PRIVATE slae_lib GTest::gtest_main)
target_include_directories(slae_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Source ${GTEST_INCLUDE_DIRS})
include(GoogleTest)
gtest_discover_tests(slae_test)
