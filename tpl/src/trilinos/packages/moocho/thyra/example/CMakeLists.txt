ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Stratimikos)
IF (${PACKAGE_NAME}_ENABLE_Stratimikos)

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
  IF (${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
    ADD_SUBDIRECTORY(NLPThyraEpetraModelEval4DOpt)
    ADD_SUBDIRECTORY(NLPThyraEpetraAdvDiffReactOpt)
    ADD_SUBDIRECTORY(DiagonalQuadraticResponseOnlyOpt)

    ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Rythmos)
    IF (${PACKAGE_NAME}_ENABLE_Rythmos)
      ADD_SUBDIRECTORY(DiagonalTransientInversion)
    ENDIF()

  ENDIF()

ENDIF()
