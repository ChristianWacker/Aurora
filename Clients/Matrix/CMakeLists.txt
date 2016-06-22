##### Aurora/Clients/Matrix/CMakeLists.txt #####################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(matrix)

SET(HEADER_FILE
    Matrix.hpp
    Pass.hpp
    Pass1.hpp
    Pass2.hpp
    Pass3.hpp
    PassBandlimit.hpp
    PassPower.hpp
    PassScanning.hpp
    PassTem.hpp)

SET(SOURCE_FILES
    Matrix.cpp
    Pass.cpp
    Pass1.cpp
    Pass2.cpp
    Pass3.cpp
    PassBandlimit.cpp
    PassPower.cpp
    PassScanning.cpp
    PassTem.cpp)

IF(AURORA_UNITY_BUILD)
  AURORA_ENABLE_UNITY_BUILD(matrix SOURCE_FILES)
ENDIF()

ADD_EXECUTABLE(matrix ${HEADER_FILES} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(matrix auroralib)

INSTALL(TARGETS matrix
        RUNTIME DESTINATION bin)
