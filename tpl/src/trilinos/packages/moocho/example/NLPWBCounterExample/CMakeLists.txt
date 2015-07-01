INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST( NLPWBCounterExample
  SOURCES
    NLPWBCounterExampleMain.cpp
    NLPInterfacePack_NLPWBCounterExample.cpp
  ARGS
    "--a=0.0 --b=1.0 --x1-init=0.0 --x2-init=0.0 --x3-init=0.0 --nlp-selects-basis --linear-obj"
    "--a=0.0 --b=1.0 --x1-init=0.0 --x2-init=0.0 --x3-init=0.0 --nlp-selects-basis --nonlinear-obj"
    "--a=0.0 --b=1.0 --x1-init=0.0 --x2-init=0.0 --x3-init=0.0 --no-nlp-selects-basis --nonlinear-obj"
    "--a=0.0 --b=0.01 --x1-init=0.0 --x2-init=0.0 --x3-init=0.0 --no-nlp-selects-basis --nonlinear-obj"
    "--a=0.0 --b=0.01 --x1-init=0.0 --x2-init=0.0 --x3-init=0.0 --no-nlp-selects-basis --linear-obj"
  NUM_MPI_PROCS 1
  )
