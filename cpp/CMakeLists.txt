cmake_minimum_required(VERSION 3.15)
project(cpp)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
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
    src/tmatrix/tmatrix_special.cpp
    src/tmatrix/tmatrix_svwf.cpp
    src/tmatrix/tmatrix_ebcm.cpp
    src/tmatrix/tmatrix_ebcm_nonaxial.cpp
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
    src/tmatrix/tmatrix_pyb.cpp
)

option(MIEPY_VALIDATE_GAUNT "Validate Gaunt recursion against direct Wigner 3j" OFF)

# Detect __float128 / libquadmath support for extended-precision T-matrix
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
    #include <quadmath.h>
    int main() { __float128 x = 1.0Q; x = sqrtq(x); return 0; }
" HAVE_QUADMATH)
if(HAVE_QUADMATH)
    find_library(QUADMATH_LIBRARY quadmath)
endif()

add_library(cpp MODULE "${SOURCES}" "${BINDINGS}")
if(MIEPY_VALIDATE_GAUNT)
    target_compile_definitions(cpp PRIVATE MIEPY_VALIDATE_GAUNT)
endif()
target_link_libraries(cpp PRIVATE pybind11::module GSL::gsl OpenMP::OpenMP_CXX Eigen3::Eigen)
# Link libquadmath for extended-precision T-matrix support (Linux/GCC only)
# MIEPY_HAS_QUAD compile definition controls whether C++ code uses __float128;
# only set when the library is actually found, to avoid undefined symbols at runtime.
if(HAVE_QUADMATH AND QUADMATH_LIBRARY)
    target_compile_definitions(cpp PRIVATE MIEPY_HAS_QUAD=1)
    target_link_libraries(cpp PRIVATE ${QUADMATH_LIBRARY})
endif()
# Link BLAS on macOS only: Accelerate benefits the LU direct solver.
# On Linux, BLAS dispatch overhead hurts the iterative solver (many small matvecs)
# and Eigen's internal kernels are competitive, so BLAS is not linked.
if(APPLE)
    target_compile_definitions(cpp PRIVATE EIGEN_USE_BLAS)
    target_link_libraries(cpp PRIVATE "-framework Accelerate")
endif()
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
