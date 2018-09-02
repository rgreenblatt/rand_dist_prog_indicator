# Project's name
PROJECT ( g_test C CXX )

# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

add_definitions(-std=c++14 -Wall -Wextra -Wconversion)

find_package(GTest REQUIRED)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}/include" ${GTEST_INCLUDE_DIRS})

add_executable(tests 
    "${PROJECT_SOURCE_DIR}/src/tests.cpp"
    "${PROJECT_SOURCE_DIR}/src/funcs.cpp"
)

add_executable(app 
    "${PROJECT_SOURCE_DIR}/src/main.cpp"
    "${PROJECT_SOURCE_DIR}/src/funcs.cpp"
)

target_link_libraries(tests ${GTEST_LIBRARIES} libtanhsinh.a pthread) 
target_link_libraries(app tanhsinh pthread) 
