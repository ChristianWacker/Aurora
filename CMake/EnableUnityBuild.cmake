##### Aurora/CMake/UnityBuild.cmake ############################################
#
#    This file is distributed under the University of Illinois Open Source
#    License. See LICENSE.TXT for details.
#
################################################################################

FUNCTION(AURORA_ENABLE_UNITY_BUILD UNITY_SUFFIX SOURCE_VARIABLE_NAME)
  SET(source_files ${${SOURCE_VARIABLE_NAME}})

  # Exclude all source source files from the compilation process
  SET_SOURCE_FILES_PROPERTIES(${source_files} PROPERTIES HEADER_FILE_ONLY true)

  # Generate a unique filename for the unity build translation unit
  SET(unity_source_file
      ${CMAKE_CURRENT_BINARY_DIR}/UnityBuild${UNITY_SUFFIX}.cpp)

  # Create the unity source file
  FILE(WRITE ${unity_source_file}
       "// Automatically generated file for Unity Builds\n")

  # Add include statements for each source file
  FOREACH(source_file ${source_files})
    FILE(APPEND ${unity_source_file}
         "#include \"${CMAKE_CURRENT_SOURCE_DIR}/${source_file}\"\n")
  ENDFOREACH()

  # Add the unity source file to the list of source files
  SET(${SOURCE_VARIABLE_NAME} ${source_files}
                              ${unity_source_file} PARENT_SCOPE)

ENDFUNCTION()
