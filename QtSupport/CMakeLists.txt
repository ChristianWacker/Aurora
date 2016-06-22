##### Aurora/QtSupport/CMakeLists.txt ##########################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(qtsupport)

# Run MOC when necessary
SET(CMAKE_AUTOMOC TRUE)

SET(QT_SUPPORT_HEADER_FILES
    QtSupport.hpp
    QView.hpp)

SET(QT_SUPPORT_SOURCE_FILES
    QtSupport.cpp
    QView.cpp)

IF(AURORA_UNITY_BUILD)
  AURORA_ENABLE_UNITY_BUILD(qtsupport QT_SUPPORT_SOURCE_FILES)
ENDIF()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ${QWT_INCLUDE_DIR})

ADD_LIBRARY(qtsupport SHARED
            ${QT_SUPPORT_HEADER_FILES}
            ${QT_SUPPORT_SOURCE_FILES})

TARGET_LINK_LIBRARIES(qtsupport auroralib
                                Qt5::Widgets
                                Qt5::PrintSupport
                                ${QWT_LIBRARY})

SET_PROPERTY(TARGET qtsupport
             APPEND PROPERTY COMPILE_DEFINITIONS AURORA_QT_SUPPORT_BUILD=1)

IF(MSVC)
  SET_PROPERTY(TARGET qtsupport
               APPEND PROPERTY LINK_FLAGS_RELEASE "/LTCG")
ENDIF()

INSTALL(TARGETS qtsupport
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
