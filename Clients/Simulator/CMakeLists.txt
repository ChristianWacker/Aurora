##### Aurora/Clients/Simulator/CMakeLists.txt ##################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(simulator)

SET(HEADER_FILES
    Simulator.hpp)

SET(SOURCE_FILES
    Simulator.cpp)

SET(FORM_FILES
    SimulatorMain.ui)

INCLUDE_DIRECTORIES(${QWT_INCLUDE_DIR})

ADD_EXECUTABLE(simulator ${HEADER_FILES} ${SOURCE_FILES} ${FORM_FILES})

TARGET_LINK_LIBRARIES(simulator auroralib qtsupport ${QWT_LIBRARY})

IF(MSVC)
  SET_PROPERTY(TARGET simulator
               APPEND PROPERTY LINK_FLAGS_RELEASE "/LTCG")
ENDIF ()

INSTALL(TARGETS simulator
        RUNTIME DESTINATION bin)
