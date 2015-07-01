
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  verySimple
  )

IF(${PACKAGE_NAME}_ENABLE_OBJECT_ORIENTED)

  IF (${PACKAGE_NAME}_ENABLE_AztecOO)
    ADD_SUBDIRECTORIES(
      aztecoo
       )
  ENDIF()

  ADD_SUBDIRECTORIES(
   ADT_3D_Burgers_6th_Pade
   Abstract_Burgers_Central_Difference
   Concrete_Burgers_Central_Difference
   verySimpleObjectOriented
   FEpetra_power_method
    )
ENDIF()
