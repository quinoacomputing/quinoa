/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */

#ifdef HAVE_CTRILINOS_EXPERIMENTAL
// LOCA Objects
#include "LOCA.H"
#include "LOCA_Epetra.H"

// Trilinos Objects
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// User's application specific files
#include "CNOX_Interface.hpp" // Interface file to NOX

// Required for reading and writing parameter lists from xml format
// Configure Trilinos with --enable-teuchos-extended
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

using namespace std;

// Define some variables that are Global to the file that
// allow NOX to be called in 3 stages (init, solve, finalize)
static Teuchos::RCP<NOX::Solver::Generic> solver;
static Teuchos::RCP<NOX::Epetra::Vector> initialGuess;
static Teuchos::RCP<CNOX_Interface> interface;
static Teuchos::RCP<LOCA::GlobalData> globalData;

extern "C" {

//prototype
void printNewtonKrylovStats(Teuchos::ParameterList& nlParams);


void cnoxinit(int* nelems, double* statevector, int* mpi_comm_ignored,
              void* blackbox_res, void* blackbox_prec,
              void (*residualFunction)(double *, double *, int, void *),
              void (*precFunction)(double *, double *, int, double*, void *))
{
  try {
    TEUCHOS_TEST_FOR_EXCEPTION(!is_null(solver), logic_error,
         "Exception: cnoxinit() called with solver!=null: "
      << "did not cnoxfinish() since last cnoxinit() call!!");
 
    // Create a communicator for Epetra objects
    Epetra_MpiComm Comm( MPI_COMM_WORLD );

    // Get the process ID and the total number of processors
    int MyPID = Comm.MyPID();

    // Begin LOCA Solver ************************************
    //
    // Create and initialize the continuation/bifurcation parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("ContinuationParam", 0.0);

    // Create parameter (options) list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Read in the parameter list from a file
    cout << "Reading parameter list from \"input.xml\"" << endl;
    Teuchos::updateParametersFromXmlFile("input.xml", paramList.get());

    // Set some default parameters
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");
    Teuchos::ParameterList& locaStepperList = locaParamsList.sublist("Stepper");
    locaStepperList.set("Continuation Parameter", "ContinuationParam");

    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlDir = nlParams.sublist("Direction");
    Teuchos::ParameterList& nlNewton = nlDir.sublist("Newton");
    Teuchos::ParameterList& lsParams = nlNewton.sublist("Linear Solver");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("MyPID", MyPID);
     if (!nlPrintParams.isParameter("Output Information"))
        nlPrintParams.set("Output Information",
                          NOX::Utils::OuterIteration  +
//                        NOX::Utils::OuterIterationStatusTest +
//                        NOX::Utils::InnerIteration +
//                        NOX::Utils::Details +
                          NOX::Utils::LinearSolverDetails
//                        NOX::Utils::Warning +
//                        NOX::Utils::StepperIteration +
//                        NOX::Utils::StepperDetails +
//                        NOX::Utils::StepperParameters
                                );


    // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
    interface = Teuchos::rcp(new CNOX_Interface(nelems, statevector, pVector, Comm,
                                                blackbox_res, blackbox_prec,
                                                residualFunction, precFunction));
    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;

    // Interface is inherited from these 2 classes as well for user prec
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;
    Teuchos::RCP<Epetra_Operator> precOperator = interface;

    // Create the Epetra_Vector for the  state vector
    Teuchos::RCP<Epetra_Vector> soln = interface->getVector();

// The next 2 statements need to be changed if your interface can
// supply a Jacobian matrix!!! See nox examples epetra .....
    Teuchos::RCP<NOX::Epetra::MatrixFree> FD =
      Teuchos::rcp(new NOX::Epetra::MatrixFree(nlPrintParams, interface, soln));

    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = FD;

    // Create the linear systems
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams,
                                         iJac, FD, iPrec, precOperator, soln));

/* With a Jacobian, and no user-define preconditioner the previous 3 calls should be:

  Teuchos::RCP<Epetra_RowMatrix> Amat = Teuchos::rcp(interface->getJacobian(),false)
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  // Create the linear systems
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams,
                                                      iReq, iJac, Amat, soln))
*/

    // Create the loca vector
    initialGuess = Teuchos::RCP<NOX::Epetra::Vector>(new NOX::Epetra::Vector(soln));

    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create global data object
    globalData = LOCA::createGlobalData(paramList, epetraFactory);

    // Create the Group
    Teuchos::RCP<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams,
                                           iReq, *initialGuess, linsys,
                                           pVector));
    grp->computeF();

    // Create the Solver convergence test
    Teuchos::RCP<NOX::StatusTest::NormF> wrms =
      Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Convergence Tolerance",1.0e-8)));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(10));
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(wrms);
    combo->addStatusTest(maxiters);

    Teuchos::RCP<Teuchos::ParameterList> nlParamsRCP = Teuchos::rcp(&nlParams, false);
    solver = NOX::Solver::buildSolver(grp, combo, nlParamsRCP);

  } //end try block

  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }
}

void cnoxsolve(int* nelems, double* statevector,
              void* blackbox_res, void* blackbox_prec)
{

  try {
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(solver), logic_error,
         "Exception: cnoxsolve called with solver=null: "
      << "either did call cnoxinit first, or called cnoxfinish already");

    // reset solver with new initial Guess
    Epetra_Vector& initialGuessEV = initialGuess->getEpetraVector();
    for (int i=0; i<*nelems; i++) initialGuessEV[i] = statevector[i];
    solver->reset(*initialGuess);

    // reset interface with new blackbox
    interface->resetBlackbox(blackbox_res, blackbox_prec);

    // reset counter for iteration statistics
    Teuchos::ParameterList& nlParams =
         const_cast<Teuchos::ParameterList&>(solver->getList());
    nlParams.sublist("Direction").sublist("Newton").
                        sublist("Linear Solver").sublist("Output").
                        set("Total Number of Linear Iterations", 0);

    cout << "NOX Solving" << endl;
    NOX::StatusTest::StatusType status = solver->solve();

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out()
        << std::endl << "Final Parameters" << std::endl
        << "****************" << std::endl;
      solver->getList().print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    printNewtonKrylovStats(nlParams);

    // copy out final solution
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>
          (solver->getSolutionGroup().getX())).getEpetraVector();
    for (int i=0; i<*nelems; i++) statevector[i] = finalSolution[i];
  }
  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }
}

void cnoxfinish()
{
  LOCA::destroyGlobalData(globalData);
  initialGuess.release();
  interface.release();
  solver.release();
}

void printNewtonKrylovStats(Teuchos::ParameterList& nlParams)
{
  static int totalNewtonIters=0;
  static int totalKrylovIters=0;
  static int stepNum=0;
  int NewtonIters = nlParams.sublist("Output").get("Nonlinear Iterations", -1000);
  int KrylovIters = nlParams.sublist("Direction").sublist("Newton").
                    sublist("Linear Solver").sublist("Output").
                    get("Total Number of Linear Iterations", -1000);
  totalNewtonIters += NewtonIters;
  totalKrylovIters += KrylovIters;
  stepNum++;

  cout << "Convergence Stats: for step  #" << stepNum
       << " : Newton, Krylov, Kr/Ne: "
       << NewtonIters << "  " << KrylovIters << "  "
       << (double) KrylovIters / (double) (NewtonIters+1.0e-14) << endl;
  if (stepNum > 1)
    cout << "Convergence Stats: running total: Newton, Krylov, Kr/Ne, Kr/Step: "
       << totalNewtonIters << "  " << totalKrylovIters << "  "
       << (double) totalKrylovIters / (double) totalNewtonIters
       << "  " << (double) totalKrylovIters / (double) stepNum << endl;
}

} //extern "C"

#endif //HAVE_CTRILINOS_EXPERIMENTAL
