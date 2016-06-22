##### Aurora/Clients/CMakeLists.txt ############################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(clients)

ADD_SUBDIRECTORY(Bench)
ADD_SUBDIRECTORY(Compare)
ADD_SUBDIRECTORY(Converter)
ADD_SUBDIRECTORY(CrossSection)
ADD_SUBDIRECTORY(EMule)
ADD_SUBDIRECTORY(FormFactor)
ADD_SUBDIRECTORY(Matrix)
ADD_SUBDIRECTORY(Multiply)
ADD_SUBDIRECTORY(OpSys)
ADD_SUBDIRECTORY(Reduce)
ADD_SUBDIRECTORY(Reslice)
ADD_SUBDIRECTORY(Subtract)
ADD_SUBDIRECTORY(Toy)

IF(FFTW_FOUND)
# ADD_SUBDIRECTORY(Ice)
ELSE()
  MESSAGE(STATUS "FFTW not found => Ice client will not be build.")
ENDIF()

IF(AURORA_BUILD_QT_CLIENTS)
  SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)
  # Run MOC when necessary
  SET(CMAKE_AUTOMOC TRUE)
  # Run UI compiler when necessary
  SET(CMAKE_AUTOUIC TRUE)

  ADD_SUBDIRECTORY(ModelBuilder)
  ADD_SUBDIRECTORY(OpSys2)
  ADD_SUBDIRECTORY(Simulator)
  ADD_SUBDIRECTORY(View)
ENDIF()
