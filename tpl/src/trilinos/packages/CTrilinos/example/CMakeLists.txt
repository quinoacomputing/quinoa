
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  verySimple
  crsgraph
  )

IF (${PACKAGE_NAME}_ENABLE_Amesos)

  ADD_SUBDIRECTORIES(
    amesos
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_AztecOO)

  ADD_SUBDIRECTORIES(
    aztecoo
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  ADD_SUBDIRECTORIES(
    galeri
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_AztecOO)
    IF (${PACKAGE_NAME}_ENABLE_Amesos)

      ADD_SUBDIRECTORIES(
        ifpack
        )

    ENDIF()
  ENDIF()
ENDIF()
