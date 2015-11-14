

TRIBITS_ADD_EXECUTABLE_AND_TEST( ExampleNLPBanded
  SOURCES
    ExampleNLPBandedMain.cpp
    NLPInterfacePack_ExampleNLPBanded.cpp
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST( ExampleNLPBanded

  POSTFIX_AND_ARGS_0  nD2_nI1_bw1_mI0_nsb
    --nD=2 --nI=1 --bw=1 --mI=0 --xo=0.1 --diag-scal=1e+1 --nlp-selects-basis

  POSTFIX_AND_ARGS_1  nD2_nI1_bw1_mI0
    --nD=2 --nI=2 --bw=2 --mI=0 --xo=0.1 --diag-scal=1e+1

  POSTFIX_AND_ARGS_2  nD1000_nI8_bw5_mI0_nsb_ds3
    --nD=1000 --nI=8 --bw=5 --mI=0 --nlp-selects-basis --xo=0.1 --diag-scal=1e+3

  POSTFIX_AND_ARGS_3  nD1000_nI300_bw20_mI0_nsb_ds2
    --nD=1000 --nI=300 --bw=20 --mI=0 --nlp-selects-basis --xo=0.1 --diag-scal=1e+2

  POSTFIX_AND_ARGS_4  nD8_nI2_bw2_mI4_ds1
    --nD=8 --nI=2 --bw=2 --mI=4 --xo=0.1 --diag-scal=1e+1

  POSTFIX_AND_ARGS_5  nD8_nI2_bw2_mI4_nsb_ds1
    --nD=100 --nI=50 --bw=100 --mI=100 --nlp-selects-basis --xo=0.1 --diag-scal=1e+1

  POSTFIX_AND_ARGS_6  nD100_nI4_bw20_mI0_ds3
    --nD=100 --nI=4 --bw=20 --mI=0 --xo=0.1 --diag-scal=1e+3

  POSTFIX_AND_ARGS_7  nD4_nI4_bw2_mI0_ds3_xIl
    --nD=4 --nI=4 --bw=2 --mI=0 --xo=0.1 --diag-scal=1e+3 --xIl=0.01

  POSTFIX_AND_ARGS_8  nD4_nI4_bw2_mI0_ds3_xIl_qpschur
    --nD=4 --nI=4 --bw=2 --mI=0 --xo=0.1 --diag-scal=1e+3 --xIl=0.01
    --moocho-extra-options="NLPAlgoConfigMamaJama{qp_solver=QPSCHUR}"

  POSTFIX_AND_ARGS_9  nD2_nI2_bw2_ds3_ic_xIl
    --nD=2 --nI=2 --ignore-constraints --xo=0.1 --xIl=0.01

  POSTFIX_AND_ARGS_10  nD2_nI2_bw2_ds3_ic_xIl_qpschur
    --nD=2 --nI=2 --ignore-constraints --xo=0.1 --xIl=0.01
    --moocho-extra-options="NLPAlgoConfigMamaJama{qp_solver=QPSCHUR}"
 
  POSTFIX_AND_ARGS_11  nD100_nI100_ic_xIl
    --nD=100 --nI=100 --ignore-constraints --xo=0.1 --xIl=0.01

  POSTFIX_AND_ARGS_12  nD100_nI100_ic_xIl_qpschur
    --nD=100 --nI=100 --ignore-constraints --xo=0.1 --xIl=0.01
    --moocho-extra-options="NLPAlgoConfigMamaJama{qp_solver=QPSCHUR}"

  POSTFIX_AND_ARGS_13  nD4_nI4_bw2_mI0_ds3_xIl_qpschur_check
    --nD=4 --nI=4 --bw=2 --mI=0 --xo=0.1 --diag-scal=1e+3 --xIl=0.01
    --moocho-extra-options="NLPSolverClientInterface{check_results=true}:NLPAlgoConfigMamaJama{qp_solver=QPSCHUR}"

  NUM_MPI_PROCS 1

  )


