cmake_minimum_required(VERSION 2.8.12)
project(cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -flto -fPIC")


find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
#find_package(pybind11 REQUIRED)
add_subdirectory(pybind11)

include_directories(src)
file(GLOB SOURCES src/main.cpp src/other.cpp)
file(GLOB INCLUDE src/main.hpp src/other.hpp)

#pybind11_add_module(cpp src/main_py.cpp src/main.cpp src/other.cpp)

add_library(cpp MODULE src/main_py.cpp src/main.cpp src/other.cpp)
target_link_libraries(cpp PRIVATE pybind11::module GSL::gsl OpenMP::OpenMP_CXX)
set_target_properties(cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")
