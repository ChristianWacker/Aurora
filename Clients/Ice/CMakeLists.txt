##### Aurora/Clients/Ice/CMakeLists.txt ########################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

Project(ice)

set(HEADER_FILES
    InelasticSimulation.hpp
    MutualCoherenceFunction.hpp)

set(SOURCE_FILES
    Ice.cpp
    InelasticSimulation.cpp
    MutualCoherenceFunction.cpp)

IF(AURORA_UNITY_BUILD)
  AURORA_ENABLE_UNITY_BUILD(ice SOURCE_FILES)
ENDIF()

# FFTW
IF(AURORA_USE_FFTW)
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIRS})
ENDIF()

ADD_EXECUTABLE(ice ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(ice auroralib)

IF(MSVC)
  SET_PROPERTY(TARGET ice
               APPEND PROPERTY LINK_FLAGS_RELEASE "/LTCG")
ENDIF()

INSTALL(TARGETS ice
        RUNTIME DESTINATION bin)
