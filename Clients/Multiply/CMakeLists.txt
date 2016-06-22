##### Aurora/Clients/Multiply/CMakeLists.txt ###################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(multiply)

SET(HEADER_FILES)

SET(SOURCE_FILES
    Multiply.cpp)

ADD_EXECUTABLE(multiply ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(multiply auroralib)

INSTALL(TARGETS multiply
        RUNTIME DESTINATION bin)
