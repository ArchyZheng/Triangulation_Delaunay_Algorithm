cmake_minimum_required(VERSION 3.26)
project(LinearInterpolation)

set(CMAKE_CXX_STANDARD 17)
add_executable(LinearInterpolation main.cpp
        src/Triangulation.cpp
        src/Triangulation.h
        test/TestTriangulation.cpp)
target_link_libraries(LinearInterpolation)
