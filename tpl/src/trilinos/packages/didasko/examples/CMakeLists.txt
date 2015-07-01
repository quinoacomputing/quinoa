INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)

IF(HAVE_DIDASKO_AMESOS)
  ADD_SUBDIRECTORY(amesos)
ENDIF()

IF(HAVE_DIDASKO_AZTECOO)
  ADD_SUBDIRECTORY(aztecoo)
ENDIF()

IF(HAVE_DIDASKO_EPETRA)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

#Commented out because this relies on zoltan, even though no zoltan
#dependency is listed and there is actually a compile error when the
#proper defines are given to compile the test. The auto tools version
#of this test is not giving the proper defines and as such the code
#has been short circuited to not actually do anything.
#
#Newer comment:  the EpetraExt code does not work in serial, so we will
#build this test only if MPI is enabled.
IF(HAVE_DIDASKO_EPETRAEXT AND ${PACKAGE_NAME}_ENABLE_Zoltan AND TPL_ENABLE_MPI)
  ADD_SUBDIRECTORY(epetraext)
ENDIF()

IF(HAVE_DIDASKO_IFPACK)
  ADD_SUBDIRECTORY(ifpack)
ENDIF()

IF(HAVE_DIDASKO_ML)
  ADD_SUBDIRECTORY(ml)
ENDIF()

IF(HAVE_DIDASKO_TEUCHOS)
  ADD_SUBDIRECTORY(teuchos)
ENDIF()

IF(HAVE_DIDASKO_TRIUTILS)
  ADD_SUBDIRECTORY(triutils)
ENDIF()

IF(TPL_ENABLE_HYPRE)
  ADD_SUBDIRECTORY(hypre)
ENDIF()

IF(HAVE_EUCLID)
  ADD_SUBDIRECTORY(euclid)
ENDIF()
