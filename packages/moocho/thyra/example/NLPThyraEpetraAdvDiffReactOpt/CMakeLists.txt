

SET(EPETRAEXT_GLPAPP_DIR ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/GLpApp)


TRIBITS_INCLUDE_DIRECTORIES(
  REQUIRED_DURING_INSTALLATION_TESTING 
  ${PACKAGE_SOURCE_DIR}/../thyra/example/operator_vector
  ${EPETRAEXT_GLPAPP_DIR}
  )


SET( COMMON_SOURCES
  ${EPETRAEXT_GLPAPP_DIR}/GLpApp_AdvDiffReactOptModel.cpp
  ${EPETRAEXT_GLPAPP_DIR}/GLpApp_AdvDiffReactOptModelCreator.cpp
  ${EPETRAEXT_GLPAPP_DIR}/GLpApp_GLpYUEpetraDataPool.cpp
  ${EPETRAEXT_GLPAPP_DIR}/GenSQP_YUEpetraVector.cpp
  ${EPETRAEXT_GLPAPP_DIR}/rect2DMeshGenerator.cpp
  )


#
# NLPThyraEpetraAdvDiffReactOpt
#


TRIBITS_ADD_EXECUTABLE_AND_TEST( NLPThyraEpetraAdvDiffReactOpt
  SOURCES
    NLPThyraEpetraAdvDiffReactOptMain.cpp
    ${COMMON_SOURCES}
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST( NLPThyraEpetraAdvDiffReactOpt

  NAME_POSTFIX fwd

  POSTFIX_AND_ARGS_0 square1
    --moocho-thyra-params-file=moochoThyra.fwd.xml
    --geom-file-base=square.1
    --moocho-extra-options="NLPSolverClientInterface{max_iter=8,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_1 square1_rr5
    --moocho-thyra-params-file=moochoThyra.fwd.xml
    --geom-file-base=square.1 --x0=1.0 --p0=0.0 --reaction-rate=1e-5
    --moocho-extra-options="NLPSolverClientInterface{max_iter=6,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_2 square2 
    --moocho-thyra-params-file=moochoThyra.fwd.xml
    --geom-file-base=square.2
    --moocho-extra-options="NLPSolverClientInterface{max_iter=7,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_3 square2_rr5 
    --moocho-thyra-params-file=moochoThyra.fwd.xml
    --geom-file-base=square.2 --x0=1.0 --p0=0.0 --reaction-rate=1e-5
    --moocho-extra-options="NLPSolverClientInterface{max_iter=6,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  NUM_MPI_PROCS 1

  )


TRIBITS_ADD_TEST( NLPThyraEpetraAdvDiffReactOpt

  NAME_POSTFIX fwd

  POSTFIX_AND_ARGS_0 aztecoo_9x9.2
    --moocho-thyra-params-file=moochoThyra.fwd.xml
    --geom-file-base=square.9x9.2
    --moocho-extra-options="MoochoSolver{test_nlp=false}:NLPSolverClientInterface{opt_tol=1e-6,feas_tol=1e-8,max_iter=11,journal_output_level=PRINT_NOTHING}"
    --linear-solver-params-file=linearSolver.aztecoo.xml
    --echo-command-line

  COMM mpi
  NUM_MPI_PROCS 2

  )


TRIBITS_ADD_TEST( NLPThyraEpetraAdvDiffReactOpt
  NAME_POSTFIX fwd
  POSTFIX_AND_ARGS_0 aztecoo_genmsh_rr3
    --moocho-thyra-params-file=moochoThyra.fwd.xml
    --len-x=1.0 --len-y=10.0 --local-nx=3 --local-ny=4 --reaction-rate=1e-3
    --moocho-extra-options="MoochoSolver{test_nlp=false}:NLPSolverClientInterface{opt_tol=1e-3,feas_tol=1e-8,max_iter=11,journal_output_level=PRINT_NOTHING}"
    --linear-solver-params-file=linearSolver.aztecoo.xml
    --echo-command-line
  COMM mpi
  NUM_MPI_PROCS 1-4
  )


TRIBITS_ADD_TEST( NLPThyraEpetraAdvDiffReactOpt

  NAME_POSTFIX opt

  POSTFIX_AND_ARGS_0 firstOrder_square1
    --moocho-thyra-params-file=moochoThyra.opt.firstOrder.xml
    --geom-file-base=square.1
    --moocho-extra-options="NLPSolverClientInterface{max_iter=16,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_1 direct_square1
    --moocho-thyra-params-file=moochoThyra.opt.direct.xml
    --geom-file-base=square.1
    --moocho-extra-options="NLPSolverClientInterface{max_iter=16,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_2 firstOrder_square1_aztecoo
    --moocho-thyra-params-file=moochoThyra.opt.firstOrder.xml
    --geom-file-base=square.1
    --moocho-extra-options="NLPSolverClientInterface{max_iter=16,journal_output_level=PRINT_NOTHING}"
    --linear-solver-params-file=linearSolver.aztecoo.xml
    --echo-command-line

  POSTFIX_AND_ARGS_3 direct_square1_aztecoo
    --moocho-thyra-params-file=moochoThyra.opt.direct.xml
    --geom-file-base=square.1
    --moocho-extra-options="NLPSolverClientInterface{max_iter=16,journal_output_level=PRINT_NOTHING}"
    --linear-solver-params-file=linearSolver.aztecoo.xml
    --echo-command-line

  POSTFIX_AND_ARGS_4 square2_np1
    --geom-file-base=square.2 --np=1
    --moocho-extra-options="NLPSolverClientInterface{max_iter=16,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_5 square2_np4
    --geom-file-base=square.2 --np=4
    --moocho-extra-options="NLPSolverClientInterface{max_iter=39,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_6 square2_direct_beta3
    --moocho-thyra-params-file=moochoThyra.opt.direct.xml
    --geom-file-base=square.2 --beta=1e+3
    --moocho-extra-options="NLPSolverClientInterface{opt_tol=0.2,max_iter=25,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_7 square2_rr5
    --geom-file-base=square.2 --reaction-rate=1e+5
    --moocho-extra-options="NLPSolverClientInterface{max_iter=15,journal_output_level=PRINT_NOTHING}"
    --echo-command-line

  POSTFIX_AND_ARGS_8 square2_orthogonal
    --geom-file-base=square.2
    --moocho-extra-options="NLPSolverClientInterface{opt_tol=1e-4,feas_tol=1e-10,max_iter=52,journal_output_level=PRINT_ALGORITHM_STEPS}:DecompositionSystemStateStepBuilderStd{range_space_matrix=ORTHOGONAL}"
    --echo-command-line

  NUM_MPI_PROCS 1

  )


TRIBITS_ADD_TEST( NLPThyraEpetraAdvDiffReactOpt
  NAME_POSTFIX opt
  POSTFIX_AND_ARGS_0 genmesh_50_5_4
    --local-nx=50 --local-ny=5 --np=4
    --moocho-extra-options="NLPSolverClientInterface{opt_tol=1e-4,feas_tol=1e-10,max_iter=40,journal_output_level=PRINT_ALGORITHM_STEPS}:DecompositionSystemStateStepBuilderStd{range_space_matrix=ORTHOGONAL}"
    --echo-command-line
  NUM_MPI_PROCS 2
  )


#
# Run varied tests
#

IF (PERL_FOUND)


  SET(EXE_PATH
    "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_NLPThyraEpetraAdvDiffReactOpt.exe")


  # Helper function to add run-varied-tests cases.  This function is needed so
  # that we can use the PACKAGE_ADD_ADVANDED_TEST(...) function with
  # OVERALL_WORKING_DIRECTORY TESTNAME.  This is needed so that each test will
  # create its own working directory.  Before, each test would create and
  # write to the working directory run_tests and when running tests in
  # parallel, they would clober each other.  Now, they can all run
  # independently.

  FUNCTION(ADD_RUN_VARIED_TESTS_CASE  POSTFIX_NAME)
    TRIBITS_ADD_ADVANCED_TEST( run-varied-tests_fdobj_bb_direct_testing
      TEST_0
       CMND ${PERL_EXECUTABLE}
       ARGS "${CMAKE_CURRENT_BINARY_DIR}/run-varied-tests.pl" --exe="${EXE_PATH}"
         ${ARGN}
        PASS_REGULAR_EXPRESSION "FINAL END RESULT: TEST PASSED"
      OVERALL_WORKING_DIRECTORY TEST_NAME
      COMM serial
      )
  ENDFUNCTION()


  # Test cases


  ADD_RUN_VARIED_TESTS_CASE( fdobj_bb_direct
    --beta=0.0 --reaction-rate=1e-2
    --local-nx=5 --global-ny=5 --starting-num-procs=1 --max-num-procs=1
    --starting-model-np=3 --use-amesos --solve-fwd-prob --x0=0.001
    --p0=1.0 --p-solu-scale=20.0 --do-analytic
    --do-finite-diff --use-finite-diff-obj
    --fd-step-len=1e-8 --do-black-box --bb-fd-step-len=1e-7
    --use-direct --fwd-newton-tol=1e-10
    --extra-args="--linear-solver-params-used-file=linear-solver.used.xml"
    )

  ADD_RUN_VARIED_TESTS_CASE( fdobjcon_direct
    --beta=0.0 --reaction-rate=1e-2
    --local-nx=5 --global-ny=5 --starting-num-procs=1 --max-num-procs=1
    --starting-model-np=3 --use-amesos --solve-fwd-prob --x0=0.001
    --p0=1.0 --p-solu-scale=20.0 --do-analytic
    --do-finite-diff --use-finite-diff-obj --use-finite-diff-con
    --fd-step-len=1e-8
    --use-direct --fwd-newton-tol=1e-10
    --extra-args="--linear-solver-params-used-file=linear-solver.used.xml"
    )

  ADD_RUN_VARIED_TESTS_CASE( fdobjcon_direct_defaultobj
    --beta=0.0 --reaction-rate=1e-2
    --local-nx=5 --global-ny=5 --starting-num-procs=1 --max-num-procs=1
    --starting-model-np=3 --use-amesos --solve-fwd-prob --x0=0.001
    --p0=1.0 --p-solu-scale=20.0 --do-analytic --use-default-inv-obj
    --do-finite-diff --use-finite-diff-obj --use-finite-diff-con
    --fd-step-len=1e-8
    --use-direct --fwd-newton-tol=1e-10
    --extra-args="--linear-solver-params-used-file=linear-solver.used.xml"
    )

  ADD_RUN_VARIED_TESTS_CASE( fdobj_direct
    --beta=0.0 --reaction-rate=1e-2
    --local-nx=5 --global-ny=5 --starting-num-procs=1
    --max-num-procs=1 --starting-model-np=3 --use-amesos --solve-fwd-prob --x0=0.001
    --p0=1.0 --p-solu-scale=20.0 --do-analytic
    --do-finite-diff --use-finite-diff-obj
    --fd-step-len=1e-8
    --use-direct --fwd-newton-tol=1e-10
    --extra-args="--linear-solver-params-used-file=linear-solver.used.xml"
    )

  ADD_RUN_VARIED_TESTS_CASE( fdobj_direct_defaultobj
    --beta=0.0 --reaction-rate=1e-2
    --local-nx=5 --global-ny=5 --starting-num-procs=1
    --max-num-procs=1 --starting-model-np=3 --use-amesos --solve-fwd-prob --x0=0.001
    --p0=1.0 --p-solu-scale=20.0 --do-analytic --use-default-inv-obj
    --do-finite-diff --use-finite-diff-obj
    --fd-step-len=1e-8
    --use-direct --fwd-newton-tol=1e-10
    --extra-args="--linear-solver-params-used-file=linear-solver.used.xml"
    )

  ADD_RUN_VARIED_TESTS_CASE( fdcon_direct
    --beta=0.0 --reaction-rate=1e-2
    --local-nx=5 --global-ny=5 --starting-num-procs=1 --max-num-procs=1
    --starting-model-np=3 --use-amesos --solve-fwd-prob --x0=0.001
    --p0=1.0 --p-solu-scale=20.0 --do-analytic
    --do-finite-diff --use-finite-diff-con
    --fd-step-len=1e-8
    --use-direct --fwd-newton-tol=1e-10
    --extra-args="--linear-solver-params-used-file=linear-solver.used.xml"
    )

  ADD_RUN_VARIED_TESTS_CASE( fdcon_direct_defaultobj
    --beta=0.0 --reaction-rate=1e-2
    --local-nx=5 --global-ny=5 --starting-num-procs=1 --max-num-procs=1
    --starting-model-np=3 --use-amesos --solve-fwd-prob --x0=0.001
    --p0=1.0 --p-solu-scale=20.0 --do-analytic --use-default-inv-obj
    --do-finite-diff --use-finite-diff-con
    --fd-step-len=1e-8
    --use-direct --fwd-newton-tol=1e-10
    --extra-args="--linear-solver-params-used-file=linear-solver.used.xml"
    )


ENDIF (PERL_FOUND)


#
# MultiPeriodNLPThyraEpetraAdvDiffReactOpt
#


TRIBITS_ADD_EXECUTABLE_AND_TEST( MultiPeriodNLPThyraEpetraAdvDiffReactOpt
  SOURCES
    MultiPeriodNLPThyraEpetraAdvDiffReactOptMain.cpp
    ${COMMON_SOURCES}
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST( MultiPeriodNLPThyraEpetraAdvDiffReactOpt

  NAME_POSTFIX np1

  POSTFIX_AND_ARGS_0  nppc1_nomp_inner
    --np=1 --num-periods-per-cluster=1 --p-perturb-scaling=2.0 --beta=0.0
    --reaction-rate=1e-2 --state-inv-error=1e-6 --param-inv-error=1e-6
    --no-multi-period --use-inner-inverse
    --linear-solver-params-file=aztecooParams.xml

  POSTFIX_AND_ARGS_1  nppc1_inner
    --np=1 --num-periods-per-cluster=1 --p-perturb-scaling=2.0 --beta=0.0
    --reaction-rate=1e-2 --state-inv-error=1e-6 --param-inv-error=1e-6
    --multi-period --use-inner-inverse
    --linear-solver-params-file=aztecooParams.xml

  POSTFIX_AND_ARGS_2  nppc1_outer
    --np=1 --num-periods-per-cluster=1 --p-perturb-scaling=2.0 --beta=0.0
    --reaction-rate=1e-2 --state-inv-error=1e-6 --param-inv-error=1e-6
    --multi-period --use-outer-inverse
    --linear-solver-params-file=aztecooParams.xml

  POSTFIX_AND_ARGS_3  nppc4_outer_initaao
    --np=1 --num-periods-per-cluster=4 --p-perturb-scaling=2.0 --period-param-scale=0.80
    --beta=0.0 --reaction-rate=1e-2 --state-inv-error=1e-6 --param-inv-error=1e-6
    --multi-period --use-outer-inverse --init-solve-all-at-once
    --linear-solver-params-file=aztecooParams.xml

  POSTFIX_AND_ARGS_4  nppc4_outer_initcont
    --np=1 --num-periods-per-cluster=4 --p-perturb-scaling=2.0 --period-param-scale=0.80
    --beta=0.0 --reaction-rate=1e-2 --state-inv-error=1e-6 --param-inv-error=1e-6
    --multi-period --use-outer-inverse --init-solve-continuation
    --linear-solver-params-file=aztecooParams.xml

  POSTFIX_AND_ARGS_5  nppc4_outer_initcont_pps2
    --np=1 --num-periods-per-cluster=4 --p-perturb-scaling=2.0 --period-param-scale=2.0
    --beta=0.0 --reaction-rate=1e-2 --state-inv-error=1e-6 --param-inv-error=1e-6
    --multi-period --use-outer-inverse --init-solve-continuation
    --linear-solver-params-file=aztecooParams.xml

  POSTFIX_AND_ARGS_6  nppc4_outer_initcont_stateless
    --np=1 --num-periods-per-cluster=4 --p-perturb-scaling=2.0 --period-param-scale=0.80
    --beta=0.0 --reaction-rate=1e-2 --state-inv-error=1e-6 --param-inv-error=1e-6
    --multi-period --use-outer-inverse --init-solve-continuation
    --use-stateless-period-model
    --linear-solver-params-file=aztecooParams.xml

  NUM_MPI_PROCS 1-4

  STANDARD_PASS_OUTPUT

  )


TRIBITS_ADD_TEST( MultiPeriodNLPThyraEpetraAdvDiffReactOpt

  NAME_POSTFIX np4

  CATEGORIES NIGHTLY

  POSTFIX_AND_ARGS_0  nppc8
    --np=4 --num-periods-per-cluster=8 --p-perturb-scaling=2.0
    --period-param-scale=1.5 --beta=0.0
    --moocho-extra-options="MoochoSolver{test_nlp=false}:NLPSolverClientInterface{opt_tol=1e-10,feas_tol=1e-10,max_iter=40}"
    --state-inv-error=1e-6 --param-inv-error=1e-6 --multi-period --use-outer-inverse
    --init-solve-continuation --use-stateless-period-model
 
  NUM_MPI_PROCS 1

  STANDARD_PASS_OUTPUT

  )


TRIBITS_ADD_TEST( MultiPeriodNLPThyraEpetraAdvDiffReactOpt

  NAME_POSTFIX np4

  POSTFIX_AND_ARGS_0  nppc8_rr2
    --np=4 --num-periods-per-cluster=8 --p-perturb-scaling=2.0
    --period-param-scale=1.5 --beta=0.0 --reaction-rate=1e-2
    --moocho-extra-options="MoochoSolver{test_nlp=false}:NLPSolverClientInterface{opt_tol=1e-10,feas_tol=1e-10,max_iter=50}"
    --state-inv-error=1e-6 --param-inv-error=1e-6 --multi-period --use-outer-inverse
    --init-solve-continuation --use-stateless-period-model
    --linear-solver-params-file=aztecooParams.xml
 
  NUM_MPI_PROCS 1-4

  STANDARD_PASS_OUTPUT

  )


TRIBITS_ADD_TEST( MultiPeriodNLPThyraEpetraAdvDiffReactOpt

  NAME_POSTFIX dcspmdv

  POSTFIX_AND_ARGS_0  nppc1
    --num-procs-per-cluster=1 --skip-solve

  POSTFIX_AND_ARGS_1  nppc2
    --num-procs-per-cluster=2 --skip-solve

  NUM_MPI_PROCS 4  # Should be 2,4,6,8

  STANDARD_PASS_OUTPUT

  )


#
# Copy files
#

TRIBITS_COPY_FILES_TO_BINARY_DIR( NLPThyraEpetraAdvDiffReactOpt_CopyFiles1
  DEST_FILES
    moochoThyra.fwd.xml
    moochoThyra.opt.firstOrder.xml
    moochoThyra.opt.direct.xml
    linearSolver.aztecoo.xml
    run-varied-tests.pl
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR( NLPThyraEpetraAdvDiffReactOpt_CopyFiles2
  DEST_FILES aztecooParams.xml
  SOURCE_PREFIX "_"
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR( NLPThyraEpetraAdvDiffReactOpt_CopyDataFiles
  DEST_FILES square.1.000 square.2.000 square.9x9.2.000 square.9x9.2.001
  SOURCE_DIR ${EPETRAEXT_GLPAPP_DIR}/data/square
  )
