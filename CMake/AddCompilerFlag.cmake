##### Aurora/CMake/AddCompilerFlag.cmake #######################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

INCLUDE(CheckCXXCompilerFlag)

MACRO(AURORA_ADD_FLAG COMPILER_FLAGS FLAG)
  # Remove special characters from the flag
  STRING(REGEX REPLACE "-" "_" SANITIZED_FLAG ${FLAG})
  STRING(REGEX REPLACE "\\+" "p" SANITIZED_FLAG ${SANITIZED_FLAG})
  CHECK_CXX_COMPILER_FLAG("-${FLAG}" COMPILER_SUPPORTS_${SANITIZED_FLAG})
  IF(COMPILER_SUPPORTS_${SANITIZED_FLAG})
    SET(${COMPILER_FLAGS} "${${COMPILER_FLAGS}} -${FLAG}")
  ENDIF()
ENDMACRO()
