##### Aurora/UnitTest/CMakeLists.txt ###########################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(unitTest)

SET(HEADER_FILES
    RelativeError.hpp)

SET(SOURCE_FILES
    CBuffer2DTest.cpp
    CtfTest.cpp
    FftTest.cpp
    HistogramTest.cpp
    MathTest.cpp
    MemoryTest.cpp
    PhysicsTest.cpp
    QuadratureTest.cpp
    RBuffer2DTest.cpp
    SmallVectorTest.cpp
    StatisticsTest.cpp
    StringRefTest.cpp
    UnitTest.cpp
    UtilsTest.cpp)

INCLUDE_DIRECTORIES(SYSTEM ${CMAKE_SOURCE_DIR}/Libs/gtest/include)

# FFTW
IF(AURORA_USE_FFTW)
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIRS})
ENDIF()

IF(AURORA_UNITY_BUILD)
  AURORA_ENABLE_UNITY_BUILD(unitTest SOURCE_FILES)
ENDIF()

ADD_EXECUTABLE(unitTest ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(unitTest auroralib gtest)

IF(MSVC)
  SET_PROPERTY(TARGET unitTest
               APPEND PROPERTY LINK_FLAGS_RELEASE "/LTCG")
ENDIF()

ADD_CUSTOM_TARGET(check COMMAND unitTest --gtest_color=yes)

INSTALL(TARGETS unitTest
        RUNTIME DESTINATION bin)
