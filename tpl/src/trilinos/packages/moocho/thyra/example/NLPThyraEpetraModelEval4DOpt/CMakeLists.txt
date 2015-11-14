

SET(EPETRAEXT_4DOPT_DIR ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/4dopt)


TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${EPETRAEXT_4DOPT_DIR})


#
# NLPThyraEpetraModelEval4DOpt
#

MULTILINE_SET( XML_PARAMS_FWD_SOLVE
  "<ParameterList><Parameter name='Solve Mode' type='string' value='Forward Solve'/></ParameterList>"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NLPThyraEpetraModelEval4DOpt
  SOURCES
    NLPThyraEpetraModelEval4DOptMain.cpp
    ${EPETRAEXT_4DOPT_DIR}/EpetraModelEval4DOpt.cpp
  ARGS
    " "
    "--extra-xml-file=ExtraXML.txt"
    "--x00=0.9"
    "--x00=0.9 --extra-xml-file=ExtraXML.txt"
    "--xt0=0.9"
    "--xt0=0.9 --extra-xml-file=ExtraXML.txt"
    "--x00=0.5 --xU0=0.6"
    "--x01=0.5 --xU1=0.6"
    "--x00=2.0 --xL0=1.1"
    "--x01=2.0 --xL1=1.1"
    "--x00=2.0 --xL0=1.1 --x01=2.0 --xL1=1.1"
    "--x00=0.5 --xU0=0.6 --x01=0.5 --xU1=0.6"
    "--moocho-thyra-params-file=moochoThyraFdNand.xml --no-support-derivs --p00=1.0"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraNLPThyraEpetraModelEval4DOptCopyFiles
  DEST_FILES 
    moochoThyraFdNand.xml
    ExtraXML.txt
  SOURCE_PREFIX "_"
  EXEDEPS NLPThyraEpetraModelEval4DOpt
  )


#
# Advanded test for parameter lumping
# 


ASSERT_DEFINED(PYTHONINTERP_FOUND)
IF (PYTHONINTERP_FOUND)

  MULTILINE_SET( MOOCHO_OPTIONS 
    "MoochoSolver{test_nlp=false}:"
    "NLPSolverClientInterface{opt_tol=1e-10,feas_tol=1e-10,max_iter=7}"
    )
  
  MULTILINE_SET( MESSAGE_A
    "A) This is the original problem with no basis lumping which has a near zero\n"
    "objective value of 1e-17 and a parameter solution of {2.0,0.0}."
    )
  
  MULTILINE_SET( MESSAGE_B
    "B) This test shows that if you pick an inital guess for p where the\n"
    "single-column basis representation can indeed represent the solution, then\n"
    "you can in fact find the true optimal value for the full parameter solution.\n"
    "The true solution is {2.0,0.0}.  If you simply add a scalar multiple of the\n"
    "single basis column B={0.707107,0.707107}, then you will get the optimal\n"
    "solution of {2.0,0.0} with the reduced parameter of p=x(var_indep)={-1.0}\n"
    "which is what we see in the MoochoJournal.out file!.\n"
    )
  
  MULTILINE_SET( MESSAGE_C
    "C) This test shows that if you pick an inital guess for p where the\n"
    "single-column basis representation can not represent the solution, then you\n"
    "do not get the same answer to the fully parameterized optimization problem.\n"
    "In this case, the solution gives the final objective function value of 1.2\n"
    "and not about 1e-17!\n"
    )
  
  MULTILINE_SET( MESSAGE_D
    "D) This test shows that with a full-rank basis, you can solve the\n"
    "optimization problem for the the same as without the parameter lumping.\n"
    "In this problem we compute the original parameters p_orig={2.0,0.0} using\n"
    "the affine approximation parameters p=x(var_indep)={2.12132,-4.94975}.\n"
    "This is good confirmation that this capability is working correctly!\n"
    )
  
  SET( THYRA_COMPARE_VECTORS
    "${PROJECT_SOURCE_DIR}/packages/thyra/core/test/operator_vector/compare-thyra-spmd-multivector-files.py"
    )
  
  TRIBITS_ADD_ADVANCED_TEST(
    NLPThyraEpetraModelEval4DOpt_lumped_param
    COMM serial mpi
    OVERALL_NUM_MPI_PROCS 1
    OVERALL_WORKING_DIRECTORY TEST_NAME
    TEST_0 EXEC NLPThyraEpetraModelEval4DOpt
      MESSAGE "${MESSAGE_A}"
      WORKING_DIRECTORY nonlumped
      ARGS
        --p00=4.0 --p01=-5.0
        --moocho-extra-options="${MOOCHO_OPTIONS}"
        --moocho-thyra-params-file="../../lumped_param_nonlumped.xml"
        --echo-command-line
    TEST_1 EXEC NLPThyraEpetraModelEval4DOpt
      MESSAGE "${MESSAGE_B}"
      WORKING_DIRECTORY param1opt
      ARGS
        --p00=2.707107 --p01=0.707107
        --moocho-extra-options="${MOOCHO_OPTIONS}"
        --moocho-thyra-params-file="../../lumped_param_param1opt.xml"
        --echo-command-line
    TEST_2 CMND ${PYTHON_EXECUTABLE}
      ARGS
        "${THYRA_COMPARE_VECTORS}"
        --file1=nonlumped/p.out.01.00 --file2=param1opt/p.out.01.00 --norm-type=inf
        --min-norm-err=0.0 --max-norm-err=1e-10 --dump
    TEST_3 EXEC NLPThyraEpetraModelEval4DOpt
      MESSAGE "${MESSAGE_C}"
      WORKING_DIRECTORY param1nonopt
      ARGS
        --p00=4.0 --p01=0.0
        --moocho-extra-options="${MOOCHO_OPTIONS}"
        --moocho-thyra-params-file="../../lumped_param_param1opt.xml"
        --echo-command-line
    TEST_4 CMND ${PYTHON_EXECUTABLE}
      ARGS
        "${THYRA_COMPARE_VECTORS}"
        --file1=nonlumped/p.out.01.00 --file2=param1nonopt/p.out.01.00 --norm-type=inf
        --min-norm-err=1.0 --max-norm-err=1.1 --dump
    TEST_5 EXEC NLPThyraEpetraModelEval4DOpt
      MESSAGE "${MESSAGE_D}"
      WORKING_DIRECTORY param2opt
      ARGS
        --p00=4.0 --p01=-5.0
        --moocho-extra-options="${MOOCHO_OPTIONS}"
        --moocho-thyra-params-file="../../lumped_param_param2opt.xml"
        --echo-command-line
    TEST_6 CMND ${PYTHON_EXECUTABLE}
      ARGS
        "${THYRA_COMPARE_VECTORS}"
        --file1=nonlumped/p.out.01.00 --file2=param2opt/p.out.01.00 --norm-type=inf
        --min-norm-err=0.0 --max-norm-err=1e-10 --dump
    )
  
  
  TRIBITS_COPY_FILES_TO_BINARY_DIR( NLPThyraEpetraModelEval4DOptCopyFiles_lumped_param
    DEST_FILES
      lumped_param_nonlumped.xml
      lumped_param_param1opt.xml
      lumped_param_param2opt.xml
    )

ENDIF(PYTHONINTERP_FOUND)


#
# Multiperiod test
#


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_MPI)
IF (${PACKAGE_NAME}_ENABLE_MPI)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    NLPThyraEpetraMultiPointModelEval4DOptMain
    SOURCES
      NLPThyraEpetraMultiPointModelEval4DOptMain.cpp
      ${EPETRAEXT_4DOPT_DIR}/EpetraMultiPointModelEval4DOpt.cpp
    ARGS " "
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
  
ENDIF()
