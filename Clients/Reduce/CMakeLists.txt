##### Aurora/Clients/Reduce/CMakeLists.txt #####################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(reduce)

SET(HEADER_FILES)

SET(SOURCE_FILES
    Reduce.cpp)

ADD_EXECUTABLE(reduce ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(reduce auroralib)

INSTALL(TARGETS reduce
        RUNTIME DESTINATION bin)
