##### Aurora/Clients/Converter/CMakeLists.txt ##################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(converter)

SET(HEADER_FILES)

SET(SOURCE_FILES
    Converter.cpp)

ADD_EXECUTABLE(converter ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(converter auroralib)

INSTALL(TARGETS converter
        RUNTIME DESTINATION bin)
