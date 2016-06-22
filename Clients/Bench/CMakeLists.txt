##### Aurora/Clients/Bench/CMakeLists.txt ######################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(bench)

SET(HEADER_FILES)

SET(SOURCE_FILES
    Bench.cpp)

# FFTW
IF(AURORA_USE_FFTW)
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIRS})
ENDIF()

ADD_EXECUTABLE(bench ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(bench auroralib)

INSTALL(TARGETS bench
        RUNTIME DESTINATION bin)
