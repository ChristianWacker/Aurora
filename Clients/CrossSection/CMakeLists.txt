##### Aurora/Clients/CrossSection/CMakeLists.txt ###############################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(crosssection)

SET(HEADER_FILES)

SET(SOURCE_FILES
    CrossSection.cpp)

ADD_EXECUTABLE(crosssection ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(crosssection auroralib)

INSTALL(TARGETS crosssection
        RUNTIME DESTINATION bin)
