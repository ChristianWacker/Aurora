##### Aurora/AuroraLib/CMakeLists.txt ##########################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(auroralib)

SET(AURORA_LIB_HEADER_FILES
    Atom.hpp
    Bandlimit.hpp
    BinaryIStream.hpp
    BinaryOStream.hpp
    Buffer2DInfo.hpp
    CBuffer2D.hpp
    CBuffer2DCodecTiff.hpp
    CommandLine.hpp
    CompilerAbstraction.hpp
    Complex.hpp
    Config.hpp.in
    Crystal.hpp
    Ctf.hpp
    CtfIncoherent.hpp
    Element.hpp
    ElementGaussian.hpp
    ElementHydrogen.hpp
    Energy.hpp
    Exception.hpp
    Factory.hpp
    Fft.hpp
    FftKiss.hpp
    Formatter.hpp
    FourierRingCorrelation.hpp
    Histogram.hpp
    Image.hpp
    ImageCodecTiff.hpp
    Lapack.hpp
    Laplace.hpp
    Loader.hpp
    MarginBuffer2D.hpp
    Math.hpp
    MatrixSupport.hpp
    Memory.hpp
    MemoryBlock.hpp
    MemoryTracker.hpp
    MicroscopeSimulation.hpp
    MultisliceParameters.hpp
    MultisliceOptions.hpp
    PhaseShift.hpp
    Physics.hpp
    Potential.hpp
    Prerequisites.hpp
    Progress.hpp
    Propagator.hpp
    Quadrature.hpp
    RBuffer2D.hpp
    Resource.hpp
    Sample.hpp
    SampleCodecAtom.hpp
    SampleCodecPdb.hpp
    Saver.hpp
    SmallVector.hpp
    Statistics.hpp
    StemDetector.hpp
    StemProbe.hpp
    StemSimulation.hpp
    StringRef.hpp
    TemSimulation.hpp
    TiffSupport.hpp
    Utils.hpp
    Version.hpp
    Windows.hpp)

SET(AURORA_LIB_SOURCE_FILES
    Bandlimit.cpp
    CBuffer2D.cpp
    CBuffer2DCodecTiff.cpp
    CommandLine.cpp
    Crystal.cpp
    Ctf.cpp
    CtfIncoherent.cpp
    Element.cpp
    ElementDoyle.cpp
    ElementKirkland.cpp
    ElementLobato.cpp
    ElementPeng.cpp
    ElementWeickenmeier.cpp
    ElementWentzel.cpp
    Exception.cpp
    Fft.cpp
    FftKiss.cpp
    Formatter.cpp
    FourierRingCorrelation.cpp
    Histogram.cpp
    Image.cpp
    ImageCodecTiff.cpp
    Laplace.cpp
    Math.cpp
    MatrixSupport.cpp
    MemoryBlock.cpp
    MicroscopeSimulation.cpp
    MultisliceParameters.cpp
    PhaseShift.cpp
    Potential.cpp
    Progress.cpp
    Propagator.cpp
    Quadrature.cpp
    RBuffer2D.cpp
    Sample.cpp
    SampleCodecAtom.cpp
    SampleCodecPdb.cpp
    StemDetector.cpp
    StemProbe.cpp
    StemSimulation.cpp
    StringRef.cpp
    TemSimulation.cpp
    Utils.cpp)

LIST(APPEND AURORA_LIB_HEADER_FILES ${CMAKE_BINARY_DIR}/AuroraLib/Config.hpp)

# TIFF library
LIST(APPEND AURORA_LIB_LIBRARIES ${TIFF_LIBRARIES})
INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIRS})

# FFTW
IF(AURORA_USE_FFTW)
  LIST(APPEND AURORA_LIB_HEADER_FILES FftFftw.hpp)
  LIST(APPEND AURORA_LIB_SOURCE_FILES FftFftw.cpp)
  LIST(APPEND AURORA_LIB_LIBRARIES ${FFTW_LIBRARIES})
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIRS})
ENDIF()

#LAPACK
IF(AURORA_USE_LAPACK)
  LIST(APPEND AURORA_LIB_LIBRARIES ${LAPACK_LIBRARIES})
ENDIF()

IF(AURORA_UNITY_BUILD)
  AURORA_ENABLE_UNITY_BUILD(AuroraLib AURORA_LIB_SOURCE_FILES)
ENDIF()

# Version number
# Create an empty Version.cpp to keep cmake happy
FILE(WRITE ${CMAKE_BINARY_DIR}/AuroraLib/Version.cpp)

LIST(APPEND AURORA_LIB_SOURCE_FILES ${CMAKE_BINARY_DIR}/AuroraLib/Version.cpp)

# Create the real Version.cpp
ADD_CUSTOM_TARGET(version
  ${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/AuroraLib/Version.cpp.in
                   -D DST=${CMAKE_BINARY_DIR}/AuroraLib/Version.cpp
                   -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                   -D WORKING_DIR=${CMAKE_SOURCE_DIR}
                   -P ${CMAKE_SOURCE_DIR}/CMake/Version.cmake)

ADD_LIBRARY(auroralib SHARED
            ${AURORA_LIB_HEADER_FILES}
            ${AURORA_LIB_SOURCE_FILES})

# auroralib depends on the target "version". Thus, before "auroralib" will be
# be build the version number needs to be updated.
ADD_DEPENDENCIES(auroralib version)

TARGET_LINK_LIBRARIES(auroralib ${AURORA_LIB_LIBRARIES})

SET_PROPERTY(TARGET auroralib
             APPEND PROPERTY COMPILE_DEFINITIONS AURORA_NONCLIENT_BUILD)

INSTALL(TARGETS auroralib
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
