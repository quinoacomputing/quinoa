MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)

  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting ${PROJECT_NAME}_ENABLE_WrapExternal=OFF"
      " because ${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES or"
      " ${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES is ON!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_WrapExternal OFF)
  ENDIF()

  IF ("${PYTHON_EXECUTABLE}" STREQUAL "")
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting ${PROJECT_NAME}_ENABLE_WrapExternal=OFF"
      " because PYTHON_EXECUTABLE=''!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_WrapExternal OFF)
  ENDIF()

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)

  IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting ${PROJECT_NAME}_ENABLE_MixedLang=OFF"
      " because ${PROJECT_NAME}_ENABLE_Fortran='${${PROJECT_NAME}_ENABLE_Fortran}'!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_MixedLang OFF)
  ENDIF()

ENDMACRO()
