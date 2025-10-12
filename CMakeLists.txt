cmake_minimum_required(VERSION 3.15)
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake"
  CACHE STRING "Vcpkg toolchain file")
set(VCPKG_OVERLAY_TRIPLETS "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/triplets"
  CACHE STRING "Vcpkg overlay triplets directory")
project(miepy LANGUAGES CXX)

# Find Python
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build the C++ extension
add_subdirectory(cpp)

# Build the Fortran tmatrix executable
# This uses the existing Makefile in src/miepy/tmatrix/nfmds
find_program(MAKE_EXECUTABLE NAMES gmake make REQUIRED)
find_program(GFORTRAN_EXECUTABLE NAMES gfortran gfortran-14 gfortran-13 gfortran-12 gfortran-11 REQUIRED)

set(NFMDS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/miepy/tmatrix/nfmds)
set(NFMDS_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/nfmds)
set(NFMDS_EXECUTABLE ${NFMDS_BUILD_DIR}/tmatrix)

# Create build directory for nfmds
file(MAKE_DIRECTORY ${NFMDS_BUILD_DIR})

# Custom command to build the Fortran executable
add_custom_command(
    OUTPUT ${NFMDS_EXECUTABLE}
    COMMAND ${MAKE_EXECUTABLE} -j1
        FC=${GFORTRAN_EXECUTABLE}
        objdir=${NFMDS_BUILD_DIR}
        exedir=${NFMDS_BUILD_DIR}
        compilemode=release
    WORKING_DIRECTORY ${NFMDS_SOURCE_DIR}
    COMMENT "Building Fortran tmatrix executable"
    VERBATIM
)

# Custom target to ensure the Fortran build happens
add_custom_target(nfmds ALL DEPENDS ${NFMDS_EXECUTABLE})

# Install the Fortran executable to miepy/bin
install(
    PROGRAMS ${NFMDS_EXECUTABLE}
    DESTINATION miepy/bin
    COMPONENT fortran
)

# Ensure miepy/bin directory exists in the installation
install(
    DIRECTORY
    DESTINATION miepy/bin
    COMPONENT fortran
)

# Unzip material database during build
set(MATERIALS_ZIP ${CMAKE_CURRENT_SOURCE_DIR}/src/miepy/materials/database.zip)
set(MATERIALS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/miepy/materials/database)

# Check if database directory exists, if not unzip it
if(NOT EXISTS ${MATERIALS_DIR})
    find_program(UNZIP_EXECUTABLE NAMES unzip REQUIRED)
    execute_process(
        COMMAND ${UNZIP_EXECUTABLE} -q ${MATERIALS_ZIP} -d ${CMAKE_CURRENT_SOURCE_DIR}/src/miepy/materials
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE UNZIP_RESULT
    )
    if(NOT UNZIP_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to unzip material database")
    endif()
endif()

# Install the unzipped material database
install(
    DIRECTORY ${MATERIALS_DIR}
    DESTINATION miepy/materials
    COMPONENT data
)
