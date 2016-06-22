##### Aurora/Clients/Compare/CMakeLists.txt ####################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(compare)

SET(HEADER_FILES)

SET(SOURCE_FILES
    Compare.cpp)

ADD_EXECUTABLE(compare ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(compare auroralib)

INSTALL(TARGETS compare
        RUNTIME DESTINATION bin)
