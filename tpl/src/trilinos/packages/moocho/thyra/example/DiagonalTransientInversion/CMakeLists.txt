

SET(EPETRAEXT_DT_DIR
  ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalTransient)


TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${EPETRAEXT_DT_DIR})


TRIBITS_ADD_EXECUTABLE( DiagonalTransientInversion
  SOURCES
    DiagonalTransientInversionMain.cpp
    ${EPETRAEXT_DT_DIR}/EpetraExt_DiagonalTransientModel.cpp
  )


TRIBITS_ADD_TEST( DiagonalTransientInversion

  POSTFIX_AND_ARGS_0  basic_no1
    --moocho-options-file=MoochoDiagTransInvSmall.opt
    --params-file=diagTransInvSmall_1.xml
    --num-observations=1
    --final-time=1e-2 --scale-params-by=1.5
    --max-sens-error=0.07
    --max-p-inv-error=1e-12
    --no-print-valid-options --verb-level=VERB_EXTREME

  POSTFIX_AND_ARGS_1  inv_no1
    --moocho-options-file=MoochoDiagTransInvSmall.opt
    --params-file=diagTransInvSmall_1.xml
    --num-observations=1
    --final-time=1e-4 --scale-params-by=1.5
    --max-sens-error=1e-4
    --no-do-inv-prob
    --no-print-valid-options --verb-level=VERB_EXTREME

  NUM_MPI_PROCS 1

  )


TRIBITS_ADD_TEST( DiagonalTransientInversion

  CATEGORIES NIGHTLY

  POSTFIX_AND_ARGS_0  inv_no7
    --moocho-options-file=MoochoDiagTransInvMedium.opt
    --params-file=diagTransInvMedium.xml
    --num-observations=7
    --final-time=1e-2 --scale-params-by=1.5
    --max-sens-error=0.03
    --do-inv-prob --max-p-inv-error=1e-10
    --no-print-valid-options --verb-level=VERB_EXTREME

  POSTFIX_AND_ARGS_1  inv_no7_BDF
    --moocho-options-file=MoochoDiagTransInvBdfVariable.opt
    --params-file=diagTransInvBdfVariable.xml
    --use-BDF
    --num-observations=7
    --final-time=1e-2 --scale-params-by=1.5
    --max-sens-error=1e-3
    --do-inv-prob --max-p-inv-error=1e-10
    --no-print-valid-options --verb-level=VERB_EXTREME

  NUM_MPI_PROCS 1

  )


TRIBITS_COPY_FILES_TO_BINARY_DIR( DiagonalTransientInversion_CopyFiles
  DEST_FILES
    MoochoDiagTransInvSmall.opt
    diagTransInvSmall_1.xml
    MoochoDiagTransInvMedium.opt
    diagTransInvMedium.xml
    MoochoDiagTransInvBdfVariable.opt
    diagTransInvBdfVariable.xml
  SOURCE_PREFIX "_"
  )
