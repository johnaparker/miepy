cmake_minimum_required(VERSION 3.15)

# Set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -pedantic -Werror")

# Architecture-specific optimization flags
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx -flto -fPIC")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -flto -fPIC")
endif()

# Find dependencies
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(pybind11 CONFIG REQUIRED)

include_directories(src)

set(SOURCES
    src/special.cpp
    src/vsh_functions.cpp
    src/vsh_translation.cpp
    src/interactions.cpp
    src/forces.cpp
    src/flux.cpp
    src/indices.cpp
    src/decomposition.cpp
)

set(BINDINGS
    src/main_pyb.cpp
    src/special_pyb.cpp
    src/vsh_translation_pyb.cpp
    src/interactions_pyb.cpp
    src/forces_pyb.cpp
    src/flux_pyb.cpp
    src/vsh_functions_pyb.cpp
    src/decomposition_pyb.cpp
)

# Create the Python module using pybind11
pybind11_add_module(cpp MODULE ${SOURCES} ${BINDINGS})
target_link_libraries(cpp PRIVATE GSL::gsl OpenMP::OpenMP_CXX Eigen3::Eigen)

# Install the module to miepy/
install(TARGETS cpp DESTINATION miepy COMPONENT python)

#if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    #add_executable(debug test/debug.cpp "${SOURCES}")
    #set_target_properties(debug PROPERTIES COMPILE_FLAGS "-w")
    #target_link_libraries(debug cpp pybind11::embed)
#endif()
