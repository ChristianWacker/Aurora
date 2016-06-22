##### Aurora/Clients/EMule/CMakeLists.txt ######################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

Project(emule)

set(HEADER_FILES)

set(SOURCE_FILES
    EMule.cpp)

ADD_EXECUTABLE(emule ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(emule auroralib)

INSTALL(TARGETS emule
        RUNTIME DESTINATION bin)
