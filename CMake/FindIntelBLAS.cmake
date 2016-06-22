#.rst:
# FindBLAS
# --------
#
# Find BLAS library
#
# This module finds an installed fortran library that implements the
# BLAS linear-algebra interface (see http://www.netlib.org/blas/).  The
# list of libraries searched for is taken from the autoconf macro file,
# acx_blas.m4 (distributed at
# http://ac-archive.sourceforge.net/ac-archive/acx_blas.html).
#
# This module sets the following variables:
#
# ::
#
#   BLAS_FOUND - set to true if a library implementing the BLAS interface
#     is found
#   BLAS_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#     and -L).
#   BLAS_LIBRARIES - uncached list of libraries (using full path name) to
#     link against to use BLAS
#   BLAS95_LIBRARIES - uncached list of libraries (using full path name)
#     to link against to use BLAS95 interface
#   BLAS95_FOUND - set to true if a library implementing the BLAS f95 interface
#     is found
#   BLA_STATIC  if set on this determines what kind of linkage we do (static)
#   BLA_VENDOR  if set checks only the specified vendor, if not set checks
#      all the possibilities
#   BLA_F95     if set on tries to find the f95 interfaces for BLAS/LAPACK
#
# ######### ## List of vendors (BLA_VENDOR) valid in this module #
# Goto,ATLAS PhiPACK,CXML,DXML,SunPerf,SCSL,SGIMATH,IBMESSL,Intel10_32
# (intel mkl v10 32 bit),Intel10_64lp (intel mkl v10 64 bit,lp thread
# model, lp64 model), # Intel10_64lp_seq (intel mkl v10 64
# bit,sequential code, lp64 model), # Intel( older versions of mkl 32
# and 64 bit), ACML,ACML_MP,ACML_GPU,Apple, NAS, Generic C/CXX should be
# enabled to use Intel mkl

#=============================================================================
# Copyright 2007-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

include(CheckFunctionExists)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_QUIET ${BLAS_FIND_QUIETLY})

set(_blas_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

macro(Check_Fortran_Libraries LIBRARIES _prefix _name _flags _list _thread)
# This macro checks for the existence of the combination of fortran libraries
# given by _list.  If the combination is found, this macro checks (using the
# Check_Fortran_Function_Exists macro) whether can link against that library
# combination using the name of a routine given by _name using the linker
# flags given by _flags.  If the combination of libraries is found and passes
# the link test, LIBRARIES is set to the list of complete library paths that
# have been found.  Otherwise, LIBRARIES is set to FALSE.

# N.B. _prefix is the prefix applied to the names of all cached variables that
# are generated internally and marked advanced by this macro.

  set(_libdir ${ARGN})

  set(_libraries_work TRUE)
  set(${LIBRARIES})
  set(_combined_name)
  if (NOT _libdir)
    if (WIN32)
      set(_libdir ENV LIB)
    elseif (APPLE)
      set(_libdir ENV DYLD_LIBRARY_PATH)
    else ()
      set(_libdir ENV LD_LIBRARY_PATH)
    endif ()
  endif ()

  foreach(_library ${_list})
    set(_combined_name ${_combined_name}_${_library})

    if(_libraries_work)
      if (BLA_STATIC)
        if (WIN32)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
        endif ()
        if (APPLE)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
        else ()
          set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
        endif ()
      else ()
        if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
          # for ubuntu's libblas3gf and liblapack3gf packages
          set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
        endif ()
      endif ()
      find_library(${_prefix}_${_library}_LIBRARY
        NAMES ${_library}
        PATHS ${_libdir}
        )
      mark_as_advanced(${_prefix}_${_library}_LIBRARY)
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
    endif()
  endforeach()

  if(_libraries_work)
    # Test this combination of libraries.
    set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_thread})
#    message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
    check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    set(CMAKE_REQUIRED_LIBRARIES)
    mark_as_advanced(${_prefix}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  endif()
  if(NOT _libraries_work)
    set(${LIBRARIES} FALSE)
  endif()
endmacro()

set(BLAS_LINKER_FLAGS)
set(BLAS_LIBRARIES)
set(BLAS95_LIBRARIES)
if (NOT $ENV{BLA_VENDOR} STREQUAL "")
  set(BLA_VENDOR $ENV{BLA_VENDOR})
else ()
  if(NOT BLA_VENDOR)
    set(BLA_VENDOR "All")
  endif()
endif ()

#BLAS in intel mkl 10 library? (em64t 64bit)
if (NOT WIN32)
  set(LM "-lm")
endif ()

if(BLAS_FIND_QUIETLY OR NOT BLAS_FIND_REQUIRED)
  find_package(Threads)
else()
  find_package(Threads REQUIRED)
endif()

set(BLAS_SEARCH_LIBS "")

    set(BLAS_mkl_SEARCH_SYMBOL sgemm)
    set(_LIBRARIES BLAS_LIBRARIES)
    if (WIN32)
      if (BLA_STATIC)
        set(BLAS_mkl_DLL_SUFFIX "")
      else()
        set(BLAS_mkl_DLL_SUFFIX "_dll")
      endif()

      # Find the main file (32-bit or 64-bit)
      set(BLAS_SEARCH_LIBS_WIN_MAIN "")
      list(APPEND BLAS_SEARCH_LIBS_WIN_MAIN
        "mkl_intel_c${BLAS_mkl_DLL_SUFFIX}")
      list(APPEND BLAS_SEARCH_LIBS_WIN_MAIN
        "mkl_intel_ilp64${BLAS_mkl_DLL_SUFFIX}")

      # Add threading/sequential libs
      set(BLAS_SEARCH_LIBS_WIN_THREAD "")
      if (NOT BLA_VENDOR STREQUAL "*_seq" OR BLA_VENDOR STREQUAL "All")
        # old version
        list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
          "libguide40 mkl_intel_thread${BLAS_mkl_DLL_SUFFIX}")
        # mkl >= 10.3
        list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
          "libiomp5md mkl_intel_thread${BLAS_mkl_DLL_SUFFIX}")
      endif()
      if (BLA_VENDOR STREQUAL "*_seq" OR BLA_VENDOR STREQUAL "All")
        list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
          "mkl_sequential${BLAS_mkl_DLL_SUFFIX}")
      endif()

      # Cartesian product of the above
      foreach (MAIN ${BLAS_SEARCH_LIBS_WIN_MAIN})
        foreach (THREAD ${BLAS_SEARCH_LIBS_WIN_THREAD})
          list(APPEND BLAS_SEARCH_LIBS
            "${MAIN} ${THREAD} mkl_core${BLAS_mkl_DLL_SUFFIX}")
        endforeach()
      endforeach()
    endif()

  if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    list(APPEND BLAS_SEARCH_LIBS 
     "mkl_intel_ilp64 mkl_intel_thread mkl_core iomp5")
  else ()
    list(APPEND BLAS_SEARCH_LIBS
     "mkl_intel_ilp64 mkl_gnu_thread mkl_core gomp")
  endif ()

  foreach (IT ${BLAS_SEARCH_LIBS})
    string(REPLACE " " ";" SEARCH_LIBS ${IT})
    if (${_LIBRARIES})
    else ()
      check_fortran_libraries(
        ${_LIBRARIES}
        BLAS
        ${BLAS_mkl_SEARCH_SYMBOL}
        ""
        "${SEARCH_LIBS}"
        "${CMAKE_THREAD_LIBS_INIT};${LM}"
        "${CMAKE_PREFIX_PATH}/lib/intel64")
    endif ()
  endforeach ()

  set(BLAS_FOUND FALSE)
  if(BLAS_LIBRARIES)
    find_path(BLAS_INCLUDE_DIR mkl.h)
    if (BLAS_INCLUDE_DIR)
      set(BLAS_FOUND TRUE)
    endif()
  endif()

  if(NOT BLAS_FIND_QUIETLY)
    if(BLAS_FOUND)
      message(STATUS "A library with BLAS API found.")
    else()
      if(BLAS_FIND_REQUIRED)
        message(FATAL_ERROR
        "A required library with BLAS API not found. Please specify library location."
        )
      else()
        message(STATUS
        "A library with BLAS API not found. Please specify library location."
        )
      endif()
    endif()
  endif()

cmake_pop_check_state()
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_blas_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
