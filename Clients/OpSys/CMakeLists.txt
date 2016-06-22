##### Aurora/Clients/OpSys/CMakeLists.txt ######################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(opsys)

SET(HEADER_FILES)

SET(SOURCE_FILES
    OpSys.cpp)

ADD_EXECUTABLE(opsys ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(opsys auroralib)

INSTALL(TARGETS opsys
        RUNTIME DESTINATION bin)
