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

#include <iostream>

#include "NLPInterfacePack_NLPWBCounterExample.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main( int argc, char* argv[] )
{
  using MoochoPack::MoochoSolver;
  using NLPInterfacePack::NLPWBCounterExample;
  using Teuchos::CommandLineProcessor;
  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {

    MoochoSolver  solver;
  
    //
    // Get options from the command line
    //
    
    double   xinit[3]          = { 0.0, 0.0, 0.0 };
    double   a                 = 0.0;
    double   b                 = 1.0;
    bool     nlp_selects_basis = true;
    bool     linear_obj        = true;

    CommandLineProcessor  clp(false); // don't throw exceptions
    solver.setup_commandline_processor(&clp);
    clp.setOption( "x1-init",  &xinit[0], "Initail guess for x(1)" );
    clp.setOption( "x2-init",  &xinit[1], "Initail guess for x(2)" );
    clp.setOption( "x3-init",  &xinit[2], "Initail guess for x(3)" );
    clp.setOption( "a",  &a, "Constant for c(1)" );
    clp.setOption( "b",  &b, "Constant for c(2)" );
    clp.setOption(
      "nlp-selects-basis", "no-nlp-selects-basis", &nlp_selects_basis
      ,"Determine if NLP will select basis" );
    clp.setOption(
      "linear-obj", "nonlinear-obj", &linear_obj
      ,"Determine if objective is linear" );
    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;
    
    //
    // Create the NLP and solve it
    //

    // Create the NLP
    NLPWBCounterExample
      nlp(xinit,a,b,nlp_selects_basis,linear_obj);
    // Set the NLP
    solver.set_nlp( Teuchos::rcp(&nlp,false) );
    // Solve the NLP
    const MoochoSolver::ESolutionStatus
      solution_status = solver.solve_nlp();
    
    //
    // Return the solution status (0 if sucessfull)
    //

    return solution_status;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,success)

  return MoochoSolver::SOLVE_RETURN_EXCEPTION;

}
