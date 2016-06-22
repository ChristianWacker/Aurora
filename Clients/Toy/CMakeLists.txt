##### Aurora/Clients/Toy/CMakeLists.txt ########################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

Project(toy)

set(HEADER_FILES)

set(SOURCE_FILES
    Toy.cpp)

ADD_EXECUTABLE(toy ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(toy auroralib)

INSTALL(TARGETS toy
        RUNTIME DESTINATION bin)
