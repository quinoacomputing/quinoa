

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ../ExampleNLPDirect)


TRIBITS_ADD_EXECUTABLE_AND_TEST( ExampleNLPFirstOrder
  SOURCES
    ExampleNLPFirstOrderMain.cpp
    NLPInterfacePack_ExampleBasisSystem.cpp
    NLPInterfacePack_ExampleNLPFirstOrder.cpp
    NLPInterfacePack_ExampleNLPFirstOrderRun.cpp
    ../ExampleNLPDirect/AbstractLinAlgPack_exampleNLPDiagSetup.cpp
    ../ExampleNLPDirect/ExampleNLPDirectRTOps.c
    ../ExampleNLPDirect/NLPInterfacePack_ExampleNLPDirect.cpp
    ../ExampleNLPDirect/NLPInterfacePack_ExampleNLPObjGrad.cpp
  ARGS
    "--n=4 --xo=0.1"
    "--n=4 --xo=0.2"
    "--n=1000 --xo=0.2"
  NUM_MPI_PROCS 1
  )
