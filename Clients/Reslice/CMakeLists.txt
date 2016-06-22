##### Aurora/Clients/Reslice/CMakeLists.txt ####################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(Reslice)

SET(HEADER_FILES)

SET(SOURCE_FILES
    Reslice.cpp)

ADD_EXECUTABLE(reslice ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(reslice auroralib)

INSTALL(TARGETS reslice
        RUNTIME DESTINATION bin)
