## Add each new test subdirectory name to the list at the bottom and that's it.

INCLUDE(AddSubdirectories)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

######################################################################

FIND_PROGRAM(PERL_PATH perl)
IF(NOT PERL_PATH)
  MESSAGE("*** Perl not found. ForTrilinos unit tests will be skipped! ***")
  SET(FORTRILINOS_SKIP_UNITTESTS 1)
ENDIF()

######################################################################
######################################################################

IF(NOT FORTRILINOS_SKIP_UNITTESTS)

## EDIT BLOCK BELOW ONLY
## UNITTEST_* must be at the TOP of the list
ADD_SUBDIRECTORIES(
  UNITTEST_PARSER
  UNITTEST_LAUNCHER
  Utils
  Error
  Epetra_Comm
  Epetra_CrsMatrix
  Epetra_BlockMap
  Epetra_SerialComm
  Epetra_Map
  Epetra_MultiVector
  Epetra_Vector
  Epetra_Export
  Epetra_Import
  AztecOO_AztecOO
  )

IF(HAVE_MPI)
ADD_SUBDIRECTORIES(
  Epetra_MpiComm
)
ENDIF()

ENDIF()

######################################################################
######################################################################

