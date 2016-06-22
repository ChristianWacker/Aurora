##### Aurora/Clients/FormFactor/CMakeLists.txt #################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(formfactor)

SET(HEADER_FILES)

SET(SOURCE_FILES
    FormFactor.cpp)

ADD_EXECUTABLE(formfactor ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(formfactor auroralib)

INSTALL(TARGETS formfactor
        RUNTIME DESTINATION bin)
