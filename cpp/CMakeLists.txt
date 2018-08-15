cmake_minimum_required(VERSION 2.8.12)
project(cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -flto -fPIC")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g")


find_package(GSL REQUIRED)
#find_package(OpenMP REQUIRED)
#find_package(pybind11 REQUIRED)
add_subdirectory(pybind11)

include_directories(src)
file(GLOB SOURCES src/main.cpp src/vsh_translation.cpp)
file(GLOB INCLUDE src/main.hpp src/vsh_translation.hpp)

#pybind11_add_module(cpp src/main_py.cpp src/main.cpp)

add_library(cpp SHARED src/main_py.cpp src/main.cpp src/vsh_translation.cpp src/interactions.cpp)
#target_link_libraries(cpp PRIVATE pybind11::module GSL::gsl OpenMP::OpenMP_CXX)
target_link_libraries(cpp PRIVATE pybind11::module GSL::gsl)
set_target_properties(cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_executable(debug test/debug.cpp src/vsh_translation.cpp src/main.cpp src/interactions.cpp)
#target_link_libraries(debug PRIVATE cpp pybind11::embed GSL::gsl OpenMP::OpenMP_CXX)
target_link_libraries(debug PRIVATE cpp pybind11::embed GSL::gsl)
