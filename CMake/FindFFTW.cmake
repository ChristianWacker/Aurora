##### Aurora/CMake/FindFFTW.cmake ##############################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

# Find FFTW
# Once done this will define
#  FFTW_FOUND         - system has FFTW
#  FFTW_INCLUDE_DIRS  - the FFTW include directory
#  FFTW_LIBRARIES     - link these to use FFTW

FIND_PATH(FFTW_INCLUDE_DIR fftw3.h)

FIND_LIBRARY(FFTW_LIBRARY_DOUBLE NAMES fftw3)
FIND_LIBRARY(FFTW_LIBRARY_FLOAT NAMES fftw3f)
FIND_LIBRARY(FFTW_LIBRARY_DOUBLE_OMP NAMES fftw3_omp)
FIND_LIBRARY(FFTW_LIBRARY_FLOAT_OMP NAMES fftw3f_omp)
FIND_LIBRARY(FFTW_LIBRARY_DOUBLE_THREADS NAMES fftw3_threads)
FIND_LIBRARY(FFTW_LIBRARY_FLOAT_THREADS NAMES fftw3f_threads)

SET(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
SET(FFTW_LIBRARIES ${FFTW_LIBRARY_DOUBLE_OMP}
                   ${FFTW_LIBRARY_FLOAT_OMP}
                   ${FFTW_LIBRARY_DOUBLE}
                   ${FFTW_LIBRARY_FLOAT})
#SET(FFTW_LIBRARIES ${FFTW_LIBRARY_DOUBLE_THREADS}
#                   ${FFTW_LIBRARY_FLOAT_THREADS}
#                   ${FFTW_LIBRARY_DOUBLE}
#                   ${FFTW_LIBRARY_FLOAT})

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if all
# listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW
                                  DEFAULT_MSG
                                  FFTW_LIBRARY_DOUBLE_OMP
                                  FFTW_LIBRARY_FLOAT_OMP
                                  FFTW_LIBRARY_DOUBLE
                                  FFTW_LIBRARY_FLOAT                                  
                                  FFTW_INCLUDE_DIRS)

MARK_AS_ADVANCED(FFTW_INCLUDE_DIRS
                 FFTW_LIBRARY_DOUBLE
                 FFTW_LIBRARY_FLOAT
                 FFTW_LIBRARY_DOUBLE_OMP
                 FFTW_LIBRARY_FLOAT_OMP
                 FFTW_LIBRARY_DOUBLE_THREADS
                 FFTW_LIBRARY_FLOAT_THREADS)
