##### Aurora/CMake/Version.cmake ###############################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                WORKING_DIRECTORY "${WORKING_DIR}"
                OUTPUT_VARIABLE "AURORA_VERSION"
                OUTPUT_STRIP_TRAILING_WHITESPACE)
CONFIGURE_FILE(${SRC} ${DST} @ONLY)
