##### Aurora/CMakeLists.txt ####################################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

PROJECT(Aurora)

# Add our cmake scripts to the module path
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

# Google Test
ADD_SUBDIRECTORY(Libs/gtest)

IF(CMAKE_BUILD_TYPE STREQUAL "")
  # Force the user to set a build type
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "None Release Debug")
ENDIF()

OPTION(AURORA_SINGLE_PRECISION "Single Precision" FALSE)
OPTION(AURORA_NO_RELATIVISTICS "Use only non-relativistic formulae" FALSE)
OPTION(AURORA_DYNAMIC_LIBRARY "Dynamic library" TRUE)
OPTION(AURORA_USE_LAPACK "Use LAPACK libraries." TRUE)
OPTION(AURORA_USE_FFTW "Use FFTW library" TRUE)
OPTION(AURORA_UNITY_BUILD "UnityBuild" FALSE)

INCLUDE(AddCompilerFlag)

# Configure the compiler
IF(MSVC)
  # Microsoft Visual C++

  # Build with multiple processes
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
  # Enable AVX instructions
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
  # Maximize speed
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
  # Whole program optimization
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /GL")
  # Disable security checks
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /GS-")

  # Link-time code generation
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE
    "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /LTCG /INCREMENTAL:NO")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE
    "${CMAKE_EXE_FLAGS_RELEASE} /LTCG /INCREMENTAL:NO")

ELSE()
  # GCC compatible compilers

  # C++14 support must be explicitly activated
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  # Deactivate runtime type information (RTTI)
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "fno-rtti")
  # Optimize the generated code for the compiling machine
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "march=native")
  # Add support for multithreading with pthread
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "pthread")
  # Warn on language extensions
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "pedantic")
  # Enable the standard set of warnings
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wall")
  # Enable some additional warnings
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wextra")
  # All warnings as errors
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Werror")
  # Look for errors in documentation comments
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wdocumentation")
  # Give a warning if a float is implicitly promoted to a double
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wdouble-promotion")
  # Warn for unknown identifiers in #if statements
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wundef")
  # Warn whenever a pointer is cast such that the required alignment of the
  # target is increased.
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wcast-align")
  # implicit conversion from enumeration type to different enumeration type
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wenum-conversion")

  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wdeprecated")
#  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wshorten-64-to-32")
#  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wcovered-switch-default")
#  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wsign-conversion")
#  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wconversion")
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wfloat-conversion")

  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "fcolor-diagnostics")

  # Release Mode: Omit frame pointers
  AURORA_ADD_FLAG(CMAKE_CXX_FLAGS_RELEASE "fomit-frame-pointer")

  # Use the GOLD linker if available
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -fuse-ld=gold -Wl,--version ERROR_QUIET
                  OUTPUT_VARIABLE LD_VERSION)
  IF("${LD_VERSION}" MATCHES "GNU gold")
    AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "fuse-ld=gold")
  ELSE()
    MESSAGE(STATUS "GOLD linker deactivated.")
  ENDIF()

  IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # Intel C++ Compiler
    AURORA_ADD_FLAG(CMAKE_CXX_FLAGS_RELEASE "parallel")
    AURORA_ADD_FLAG(CMAKE_CXX_FLAGS_RELEASE "ipo")
    AURORA_ADD_FLAG(CMAKE_CXX_FLAGS_RELEASE "O2")

  ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Clang
#    AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wshadow")
    AURORA_ADD_FLAG(CMAKE_CXX_FLAGS_RELEASE "O3")

  ELSEIF(CMAKE_COMPILER_IS_GNUCXX)
    # GNU C++ Compiler
    # Disable warning when ending a comment line with '\'
    AURORA_ADD_FLAG(CMAKE_CXX_FLAGS "Wno-comment")
    AURORA_ADD_FLAG(CMAKE_CXX_FLAGS_RELEASE "O3")

  ENDIF()
ENDIF()

# Deactivate OpenMP for MSVC since it does not support OpenMP3
IF(NOT MSVC)
  FIND_PACKAGE(OpenMP QUIET)
ENDIF()

# Activate OpenMP support of the compiler, if available
IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

IF(AURORA_UNITY_BUILD)
  INCLUDE(EnableUnityBuild)
ENDIF()

# Use the git revision as version number
FIND_PACKAGE(Git REQUIRED)

FIND_PACKAGE(TIFF REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${TIFF_INCLUDE_DIR})

IF(AURORA_USE_FFTW)
  FIND_PACKAGE(FFTW)
  IF(NOT FFTW_FOUND)
    MESSAGE(STATUS "FFTW support deactivated.")
    SET(AURORA_USE_FFTW FALSE)
  ENDIF()
ENDIF()

IF(AURORA_SINGLE_PRECISION)
  SET(AURORA_USE_LAPACK FALSE)
ENDIF()

# BLAS and LAPACK support
IF(AURORA_USE_LAPACK)
  FIND_PACKAGE(IntelBLAS QUIET)

  IF(BLAS_FOUND)
    MESSAGE(STATUS "BLAS Libraries: " ${BLAS_LIBRARIES})
    INCLUDE_DIRECTORIES(SYSTEM ${BLAS_INCLUDE_DIR})
    FIND_PACKAGE(IntelLAPACK QUIET)
  ENDIF()

  IF(LAPACK_FOUND)
    MESSAGE(STATUS "LAPACK Libraries: " ${LAPACK_LIBRARIES})
  ELSE()
    MESSAGE(STATUS "LAPACK support deactivated.")
    SET(AURORA_USE_LAPACK FALSE)
  ENDIF()
ENDIF()

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
INCLUDE_DIRECTORIES(SYSTEM ${CMAKE_SOURCE_DIR}/Libs/Eigen)

# Write the configuration file
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/AuroraLib/Config.hpp.in
               ${CMAKE_BINARY_DIR}/AuroraLib/Config.hpp)

ADD_SUBDIRECTORY(AuroraLib)
ADD_SUBDIRECTORY(Docs)
ADD_SUBDIRECTORY(UnitTest)

# QT
FIND_PACKAGE(Qt5Widgets)

IF(Qt5Widgets_FOUND)
  FIND_PACKAGE(Qt5PrintSupport)
ENDIF()

IF(Qt5PrintSupport_FOUND)
  FIND_PACKAGE(Qwt)
ENDIF()

IF(QWT_FOUND)
  SET(AURORA_BUILD_QT_CLIENTS TRUE)
  ADD_SUBDIRECTORY(QtSupport)
ELSE()
  MESSAGE(STATUS "Qt Support deactivated.")
  SET(AURORA_BUILD_QT_CLIENTS FALSE)
ENDIF()

ADD_SUBDIRECTORY(Clients)
