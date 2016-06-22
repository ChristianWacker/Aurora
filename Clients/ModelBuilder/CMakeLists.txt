##### Aurora/Clients/ModelBuilder/CMakeLists.txt ###############################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(modelbuilder)

SET(HEADER_FILES
    ModelBuilder.hpp)

SET(SOURCE_FILES
    ModelBuilder.cpp)

SET(FORM_FILES
    ModelBuilderMain.ui)

ADD_EXECUTABLE(modelbuilder ${HEADER_FILES} ${SOURCE_FILES} ${FORM_FILES})

TARGET_LINK_LIBRARIES(modelbuilder auroralib qtsupport)

IF(MSVC)
  SET_PROPERTY(TARGET modelbuilder
               APPEND PROPERTY LINK_FLAGS_RELEASE "/LTCG")
ENDIF()

INSTALL(TARGETS modelbuilder
        RUNTIME DESTINATION bin)
