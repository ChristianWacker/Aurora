##### Aurora/Clients/OpSys2/CMakeLists.txt #####################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(OpSys2)

SET(HEADER_FILES
    OpSys2.hpp)

SET(SOURCE_FILES
    OpSys2.cpp)

SET(FORM_FILES
    OpSys2Main.ui)

INCLUDE_DIRECTORIES(${QWT_INCLUDE_DIR})

ADD_EXECUTABLE(opsys2
               ${HEADER_FILES}
               ${SOURCE_FILES}
               ${FORM_FILES})

TARGET_LINK_LIBRARIES(opsys2 auroralib qtsupport ${QWT_LIBRARY})

IF(MSVC)
  set_property(TARGET opsys2
               APPEND PROPERTY LINK_FLAGS_RELEASE "/LTCG")
ENDIF()

INSTALL(TARGETS opsys2
        RUNTIME DESTINATION bin)
