cmake_minimum_required(VERSION 3.15)
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake"
  CACHE STRING "Vcpkg toolchain file")
set(VCPKG_OVERLAY_TRIPLETS "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/triplets"
  CACHE STRING "Vcpkg overlay triplets directory")
project(miepy LANGUAGES CXX)

# Find Python
# Windows needs full Development to link against Python DLL with MinGW
# Unix uses Development.Module which doesn't link libpython
if(WIN32)
  find_package(Python REQUIRED COMPONENTS Interpreter Development)

  # Generate MinGW-compatible import library from Python DLL
  # Python 3.8+ doesn't ship with libpythonXX.a for MinGW
  find_program(GENDEF_EXECUTABLE gendef REQUIRED HINTS "C:/msys64/mingw64/bin")
  find_program(DLLTOOL_EXECUTABLE dlltool REQUIRED HINTS "C:/msys64/mingw64/bin")

  # Find Python DLL (handles virtual environments)
  set(PYTHON_DLL_NAME "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}.dll")

  # Method 1: Try Python_RUNTIME_LIBRARY_DIRS (CMake 3.15+)
  if(DEFINED Python_RUNTIME_LIBRARY_DIRS)
    foreach(DIR ${Python_RUNTIME_LIBRARY_DIRS})
      if(EXISTS "${DIR}/${PYTHON_DLL_NAME}")
        set(PYTHON_DLL "${DIR}/${PYTHON_DLL_NAME}")
        break()
      endif()
    endforeach()
  endif()

  # Method 2: Get base installation from Python (for venvs)
  if(NOT PYTHON_DLL)
    execute_process(
      COMMAND ${Python_EXECUTABLE} -c "import sys; print(sys.base_prefix)"
      OUTPUT_VARIABLE PYTHON_BASE_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _result
    )
    if(_result EQUAL 0 AND EXISTS "${PYTHON_BASE_PREFIX}/${PYTHON_DLL_NAME}")
      set(PYTHON_DLL "${PYTHON_BASE_PREFIX}/${PYTHON_DLL_NAME}")
    endif()
  endif()

  # Method 3: Check executable directory and parent (standard installs)
  if(NOT PYTHON_DLL)
    get_filename_component(PYTHON_DIR ${Python_EXECUTABLE} DIRECTORY)
    if(EXISTS "${PYTHON_DIR}/${PYTHON_DLL_NAME}")
      set(PYTHON_DLL "${PYTHON_DIR}/${PYTHON_DLL_NAME}")
    else()
      get_filename_component(PYTHON_PARENT_DIR ${PYTHON_DIR} DIRECTORY)
      if(EXISTS "${PYTHON_PARENT_DIR}/${PYTHON_DLL_NAME}")
        set(PYTHON_DLL "${PYTHON_PARENT_DIR}/${PYTHON_DLL_NAME}")
      endif()
    endif()
  endif()

  # Method 4: Search in Windows system directories
  if(NOT PYTHON_DLL)
    find_file(PYTHON_DLL_FOUND ${PYTHON_DLL_NAME}
      PATHS ENV PATH
      NO_DEFAULT_PATH
    )
    if(PYTHON_DLL_FOUND)
      set(PYTHON_DLL ${PYTHON_DLL_FOUND})
    endif()
  endif()

  if(NOT PYTHON_DLL OR NOT EXISTS ${PYTHON_DLL})
    message(FATAL_ERROR "Could not find Python DLL (${PYTHON_DLL_NAME}). Searched in:\n"
      "  - Python_RUNTIME_LIBRARY_DIRS\n"
      "  - sys.base_prefix: ${PYTHON_BASE_PREFIX}\n"
      "  - Executable directory: ${PYTHON_DIR}\n"
      "  - System PATH")
  endif()

  message(STATUS "Found Python DLL: ${PYTHON_DLL}")

  # Generate import library in build directory
  set(PYTHON_DEF_FILE "${CMAKE_CURRENT_BINARY_DIR}/python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}.def")
  set(PYTHON_MINGW_LIB "${CMAKE_CURRENT_BINARY_DIR}/libpython${Python_VERSION_MAJOR}${Python_VERSION_MINOR}.a")

  add_custom_command(
    OUTPUT ${PYTHON_DEF_FILE}
    COMMAND ${GENDEF_EXECUTABLE} ${PYTHON_DLL} -o ${PYTHON_DEF_FILE}
    COMMENT "Generating Python .def file for MinGW"
    VERBATIM
  )

  add_custom_command(
    OUTPUT ${PYTHON_MINGW_LIB}
    COMMAND ${DLLTOOL_EXECUTABLE} --dllname python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}.dll --def ${PYTHON_DEF_FILE} --output-lib ${PYTHON_MINGW_LIB}
    DEPENDS ${PYTHON_DEF_FILE}
    COMMENT "Generating MinGW import library for Python"
    VERBATIM
  )

  add_custom_target(python_import_lib ALL DEPENDS ${PYTHON_MINGW_LIB})

  # Make this library available to cpp subdirectory
  set(PYTHON_MINGW_LIBRARY ${PYTHON_MINGW_LIB} CACHE INTERNAL "MinGW Python import library")
else()
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
endif()

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
set(NFMDS_EXECUTABLE ${NFMDS_BUILD_DIR}/tmatrix${CMAKE_EXECUTABLE_SUFFIX})

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
