/*
// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "EpetraExt_DiagonalQuadraticResponseOnlyModelEvaluator.hpp"
#include "MoochoPack_MoochoThyraSolver.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif


int main( int argc, char* argv[] )
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::OSTab;
  using MoochoPack::MoochoSolver;
  using MoochoPack::MoochoThyraSolver;
  using Teuchos::CommandLineProcessor;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Teuchos::Time timer("");
  
  bool dummySuccess = true;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {
 
    //
    // A) Create the solver objects that will insert their command-line
    // options
    //

    MoochoThyraSolver solver;

    //
    // B) Get options from the command line
    //

    int localDim = 4;
    double pt = 0.0;
    double p0 = 0.1;
    double scale = 0.1;

    CommandLineProcessor clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    clp.setOption("local-dim", &localDim);
    clp.setOption("pt", &pt);
    clp.setOption("p0", &p0);
    clp.setOption("scale", &scale);

    solver.setupCLP(&clp);

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    solver.readParameters( out.get() );

    Teuchos::RCP<Epetra_Comm> comm = Teuchos::null;

#ifdef HAVE_MPI
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    comm = Teuchos::rcp(new Epetra_MpiComm(mpiComm));
#else
    comm = Teuchos::rcp(new Epetra_SerialComm());
#endif
 
    //
    // C) Create the Thyra::ModelEvaluator object
    //
 
    *out << "\nCreate EpetraExt::DiagonalQuadraticResponseOnlyModelEvaluator object ...\n";
 
    Teuchos::RCP<EpetraExt::DiagonalQuadraticResponseOnlyModelEvaluator>
      epetraModel = EpetraExt::diagonalQuadraticResponseOnlyModelEvaluator(
        comm, localDim, pt, p0, scale);
      
    *out << "\nCreate the Thyra::EpetraModelEvaluator wrapper object ...\n";
    
    Teuchos::RCP<Thyra::EpetraModelEvaluator>
      epetraThyraModel(new Thyra::EpetraModelEvaluator()); // Sets default options!
    epetraThyraModel->initialize(epetraModel, Teuchos::null);
    
    //
    // D) Solve the NLP
    //
    
    // Set the model
    solver.setModel(epetraThyraModel);
    
    // Read the initial guess if one exists
    solver.readInitialGuess(out.get());
    
    // Solve the NLP
    const MoochoSolver::ESolutionStatus	solution_status = solver.solve();
    
    // Write the solution to file
    solver.writeFinalSolution(out.get());
    
    // Write the final parameters to file
    solver.writeParamsFile();
    
    //
    // E) Return the solution status (0 if successful)
    //
    
    if(solution_status == MoochoSolver::SOLVE_RETURN_SOLVED)
      *out << "\nEnd Result: TEST PASSED\n";
    else
      *out << "\nEnd Result: TEST FAILED\n";
    
    return solution_status;
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, dummySuccess)
    
  return MoochoSolver::SOLVE_RETURN_EXCEPTION;

}
