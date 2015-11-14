

#
# EpetraExtDiagonalQuadraticResponseOnlyOpt
#


SET(EPETRAEXT_DQROO_DIR
  ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalQuadraticResponseOnlyOpt)

SET(THYRA_DQROME_DIR
  thyra/core/test/nonlinear/models)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${EPETRAEXT_DQROO_DIR})
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_BINARY_DIR}/../${THYRA_DQROME_DIR})
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/../${THYRA_DQROME_DIR})


TRIBITS_ADD_EXECUTABLE_AND_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  SOURCES
    EpetraExtDiagonalQuadraticResponseOnlyOptMain.cpp
    ${EPETRAEXT_DQROO_DIR}/EpetraExt_DiagonalQuadraticResponseOnlyModelEvaluator.cpp
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  NAME_POSTFIX np1
  ARGS "--local-dim=100 --pt=0.0 --p0=2.0"
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  NAME_POSTFIX np2
  ARGS "--local-dim=100 --pt=0.0 --p0=2.0"
  COMM mpi
  NUM_MPI_PROCS 2
  )


TRIBITS_ADD_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  NAME_POSTFIX np3
  ARGS "--local-dim=100 --pt=0.0 --p0=2.0"
  COMM mpi
  NUM_MPI_PROCS 3
  )


TRIBITS_ADD_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  NAME_POSTFIX np4
  ARGS "--local-dim=100 --pt=0.0 --p0=2.0"
  COMM mpi
  NUM_MPI_PROCS 4
  )


#
# ThyraDiagonalQuadraticResponseOnlyOpt
#


TRIBITS_ADD_EXECUTABLE_AND_TEST( ThyraDiagonalQuadraticResponseOnlyOpt
  SOURCES
    ThyraDiagonalQuadraticResponseOnlyOptMain.cpp
  TESTONLYLIBS thyra_nonlinear_test_models
  NUM_MPI_PROCS 1
  )

# ToDo: Add some more tests ...
