##### Aurora/Clients/Subtract/CMakeLists.txt ###################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(Subtract)

SET(HEADER_FILES)

SET(SOURCE_FILES
     Subtract.cpp)

ADD_EXECUTABLE(subtract ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(subtract auroralib)

INSTALL(TARGETS subtract
        RUNTIME DESTINATION bin)
