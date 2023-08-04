cmake_minimum_required(VERSION 3.26)
project(TestTriangulation)

set(CMAKE_CXX_STANDARD 17)
include_directories(../src)
add_subdirectory(googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(TestTriangulation TestTriangulation.cpp
        ../src/Triangulation.cpp)
target_link_libraries(TestTriangulation gtest gtest_main)