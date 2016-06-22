##### Aurora/Docs/CMakeLists.txt ###############################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

PROJECT(Docs)

ADD_CUSTOM_TARGET(docs)

FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                 @ONLY)
  ADD_CUSTOM_TARGET(doxygen
                    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Building Doxygen documentation" VERBATIM)
  
  ADD_DEPENDENCIES(docs doxygen)
ENDIF()

FIND_PACKAGE(Sphinx)
IF(SPHINX_FOUND)
  SET(SPHINX_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")
  SET(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
  SET(SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  SET(SPHINX_OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/sphinx")

  CONFIGURE_FILE("${SPHINX_SOURCE_DIR}/conf.py.in"
                 "${SPHINX_BUILD_DIR}/conf.py"
                 @ONLY)

  ADD_CUSTOM_TARGET(sphinx
                    COMMAND ${SPHINX_EXECUTABLE}
                            -q -b html
                            -c "${SPHINX_BUILD_DIR}"
                            -d "${SPHINX_CACHE_DIR}"
                            "${SPHINX_SOURCE_DIR}"
                            "${SPHINX_OUT_DIR}"
                    COMMENT "Building Sphinx documentation")

  ADD_DEPENDENCIES(docs sphinx)
ENDIF()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/index.html
               ${CMAKE_CURRENT_BINARY_DIR}/index.html
               COPYONLY)
