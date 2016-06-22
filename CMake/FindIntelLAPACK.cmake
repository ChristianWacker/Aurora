# - Find LAPACK library
# This module finds an installed fortran library that implements the LAPACK
# linear-algebra interface (see http://www.netlib.org/lapack/).
#
# The approach follows that taken for the autoconf macro file, acx_lapack.m4
# (distributed at http://ac-archive.sourceforge.net/ac-archive/acx_lapack.html).
#
# This module sets the following variables:
#  LAPACK_FOUND - set to true if a library implementing the LAPACK interface
#    is found
#  LAPACK_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#    and -L).
#  LAPACK_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use LAPACK
#  BLA_STATIC  if set on this determines what kind of linkage we do (static)
#  BLA_VENDOR  if set checks only the specified vendor, if not set checks
#     all the possibilities

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

set(_lapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
include(CheckFortranFunctionExists)

set(LAPACK_FOUND FALSE)

macro(Check_Lapack_Libraries LIBRARIES _prefix _name _flags _list _blas _threads)
# This macro checks for the existence of the combination of fortran libraries
# given by _list.  If the combination is found, this macro checks (using the
# Check_Fortran_Function_Exists macro) whether can link against that library
# combination using the name of a routine given by _name using the linker
# flags given by _flags.  If the combination of libraries is found and passes
# the link test, LIBRARIES is set to the list of complete library paths that
# have been found.  Otherwise, LIBRARIES is set to FALSE.

# N.B. _prefix is the prefix applied to the names of all cached variables that
# are generated internally and marked advanced by this macro.

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
    if(UNIX AND BLA_STATIC)
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} "-Wl,--start-group" ${${LIBRARIES}} ${_blas} "-Wl,--end-group" ${_threads})
    else()
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blas} ${_threads})
    endif()
    check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    set(CMAKE_REQUIRED_LIBRARIES)
    mark_as_advanced(${_prefix}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  endif()

  if(_libraries_work)
    set(${LIBRARIES} ${${LIBRARIES}} ${_blas} ${_threads})
  else()
    set(${LIBRARIES} FALSE)
  endif()

endmacro()

set(LAPACK_LINKER_FLAGS)
set(LAPACK_LIBRARIES)

if (NOT BLAS_FOUND)
  if(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
    find_package(IntelBLAS)
  else()
    find_package(IntelBLAS REQUIRED)
  endif()
endif()

if(BLAS_FOUND)
  set(LAPACK_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
  if ($ENV{BLA_VENDOR} MATCHES ".+")
    set(BLA_VENDOR $ENV{BLA_VENDOR})
  else ()
    if(NOT BLA_VENDOR)
      set(BLA_VENDOR "All")
    endif()
  endif ()

  #intel lapack
  if (NOT WIN32)
    set(LM "-lm")
  endif()
  if (_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)
    if(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
      find_PACKAGE(Threads)
    else()
      find_package(Threads REQUIRED)
    endif()

    if(NOT LAPACK_LIBRARIES)
        # old
      check_lapack_libraries(
        LAPACK_LIBRARIES
        LAPACK
        cheev
        ""
        "mkl_lapack"
        "${BLAS_LIBRARIES}"
        "${CMAKE_THREAD_LIBS_INIT};${LM}"
        )
    endif()
    if(NOT LAPACK_LIBRARIES)
      # new >= 10.3
      check_lapack_libraries(
        LAPACK_LIBRARIES
        LAPACK
        cheev
        ""
        "mkl_gf_lp64"
        "${BLAS_LIBRARIES}"
        "${CMAKE_THREAD_LIBS_INIT};${LM}"
        )
    endif()
  endif()
else()
  message(STATUS "LAPACK requires BLAS")
endif()

if(LAPACK_LIBRARIES)
  set(LAPACK_FOUND TRUE)
else()
  set(LAPACK_FOUND FALSE)
endif()

if(NOT LAPACK_FIND_QUIETLY)
  if(LAPACK_FOUND)
    message(STATUS "A library with LAPACK API found.")
  else()
    if(LAPACK_FIND_REQUIRED)
      message(FATAL_ERROR
      "A required library with LAPACK API not found. Please specify library location."
      )
    else()
      message(STATUS
      "A library with LAPACK API not found. Please specify library location."
      )
    endif()
  endif()
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_lapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
