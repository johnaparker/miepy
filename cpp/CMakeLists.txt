cmake_minimum_required(VERSION 2.8.12)
project(cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mavx -std=c++17 -flto -fPIC")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
endif()

find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
#find_package(pybind11 REQUIRED)
add_subdirectory(pybind11)

include_directories(src)

set(SOURCES ${SOURCES} 
    src/special.cpp
    src/vsh_functions.cpp
    src/vsh_translation.cpp
    src/interactions.cpp
    src/forces.cpp
    src/flux.cpp
    src/indices.cpp
)

set(BINDINGS
    src/main_pyb.cpp
    src/special_pyb.cpp
    src/vsh_translation_pyb.cpp
    src/interactions_pyb.cpp
    src/forces_pyb.cpp
    src/flux_pyb.cpp
    src/vsh_functions_pyb.cpp
)

add_library(cpp SHARED "${SOURCES}" "${BINDINGS}")
target_link_libraries(cpp PUBLIC pybind11::module GSL::gsl OpenMP::OpenMP_CXX)
set_target_properties(cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    add_executable(debug test/debug.cpp "${SOURCES}")
    target_link_libraries(debug cpp pybind11::embed)
endif()
