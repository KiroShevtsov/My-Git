cmake_minimum_required(VERSION 3.24)
project(slae)
set(CMAKE_CXX_STANDARD 20)
add_executable(${PROJECT_NAME} Cxx/s.cxx)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Source")
enable_testing()

add_executable(
  slae_test
  Test/test.cxx
)
target_link_libraries(
  slae_test
  GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(slae_test)