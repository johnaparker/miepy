cmake_minimum_required(VERSION 3.15)
project(cpp)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -pedantic -Werror")

# Set architecture-appropriate optimization flags
if(WIN32)
  # Windows: No LTO (breaks Python symbol imports with MinGW), no -fPIC (not applicable)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx")
  else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx -flto -fPIC")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -flto -fPIC")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -flto -fPIC")
endif()

find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(pybind11 REQUIRED)

include_directories(src)

set(SOURCES ${SOURCES}
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

add_library(cpp MODULE "${SOURCES}" "${BINDINGS}")
target_link_libraries(cpp PRIVATE pybind11::module GSL::gsl OpenMP::OpenMP_CXX Eigen3::Eigen)
# Windows with MinGW: use generated import library
if(WIN32)
  # Depend on the custom target that generates the import library
  add_dependencies(cpp python_import_lib)
  # Link against the generated MinGW-compatible import library
  target_link_libraries(cpp PRIVATE ${PYTHON_MINGW_LIBRARY})
endif()
# Force vcpkg Eigen3 to be used by adding it BEFORE any other includes, and demote system eigen3
get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
target_compile_options(cpp PRIVATE
    "-isystem" "${EIGEN3_INCLUDE_DIR}"
    "-idirafter" "/usr/include/eigen3"
)
# macOS: Statically link GCC runtime libraries
if(APPLE)
    target_link_options(cpp PRIVATE "-static-libgcc")
endif()
set_target_properties(cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS cpp DESTINATION miepy)

#if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    #add_executable(debug test/debug.cpp "${SOURCES}")
    #set_target_properties(debug PROPERTIES COMPILE_FLAGS "-w")
    #target_link_libraries(debug cpp pybind11::embed)
#endif()
