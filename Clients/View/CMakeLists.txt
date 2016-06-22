##### Aurora/Clients/View/CMakeLists.txt #######################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(view)

SET(HEADER_FILES
    View.hpp)

SET(SOURCE_FILES
    View.cpp)

SET(FORM_FILES
    MainView.ui)

INCLUDE_DIRECTORIES(${QWT_INCLUDE_DIR})

ADD_EXECUTABLE(view ${HEADER_FILES} ${SOURCE_FILES} ${FORM_FILES})

TARGET_LINK_LIBRARIES(view auroralib qtsupport ${QWT_LIBRARY})

IF(MSVC)
  SET_PROPERTY(TARGET view
               APPEND PROPERTY LINK_FLAGS_RELEASE "/LTCG")
ENDIF()

INSTALL(TARGETS view
        RUNTIME DESTINATION bin)
