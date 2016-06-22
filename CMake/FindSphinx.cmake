##### Aurora/CMake/FindSphinx.cmake ############################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

FIND_PROGRAM(SPHINX_EXECUTABLE NAMES sphinx-build
             HINTS $ENV{SPHINX_DIR} ${SPHINX_DIR}
             PATH_SUFFIXES bin
             DOC "Sphinx documentation generator")

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Sphinx DEFAULT_MSG SPHINX_EXECUTABLE)

MARK_AS_ADVANCED(SPHINX_EXECUTABLE)
