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

#include "NLPInterfacePack_ExampleNLPBanded.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"

int main( int argc, char* argv[] )
{
  namespace mp  = MoochoPack;
  namespace nlpip = NLPInterfacePack;
  using mp::MoochoSolver;
  using nlpip::NLP;
  using nlpip::ExampleNLPBanded;
  typedef nlpip::size_type  size_type;
  typedef nlpip::value_type value_type;
  using Teuchos::CommandLineProcessor;
  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {

    MoochoSolver solver;
  
    //
    // Get options from the command line
    //
    
    bool     show_options = false;
    int      nD = 1;
    int      nI = 1;
    int      bw = 1;
    int      mI = 0;
    double   xo = 1;
    bool     nlp_selects_basis = true;
    double   xDl = -NLP::infinite_bound();
    double   xDu = +NLP::infinite_bound(); 
    double   xIl = -NLP::infinite_bound(); 
    double   xIu = +NLP::infinite_bound();
    /*double   xDl = -1000.;
    double   xDu = +1000.;
    double   xIl = -1000.;
    double   xIu = +1000.;*/

    int      mU = 0;
    double   hl = -NLP::infinite_bound();
    double   hu = +NLP::infinite_bound();
    double   diag_scal = 1.0;
    double   diag_vary = 1.0;
    bool     sym_basis = false;
    double   f_offset  = 0.0;
    double   co        = 0.0;
    bool     ignore_constraints = false;
    
    CommandLineProcessor clp(false); // don't throw exceptions

    solver.setup_commandline_processor(&clp);

    clp.setOption( "show-options", "no-show-options", &show_options, "Show the commandline options or not." );
    clp.setOption( "nD",  &nD, "Number of dependent variables" );
    clp.setOption( "nI",  &nI, "Number of independent variables" );
    clp.setOption( "bw",  &bw, "Band width of the basis matrix" );
    clp.setOption( "mI",  &mI, "Number of general inequality constriants" );
    clp.setOption( "xo",  &xo, "Initial guess for x" );
    clp.setOption( "xDl", &xDl, "Lower bounds on xD" );
    clp.setOption( "xDu", &xDu, "Upper bounds on xD" );
    clp.setOption( "xIl", &xIl, "Lower bounds on xI" );
    clp.setOption( "xIu", &xIu, "Upper bounds on xI" );
//		clp.setOption( "mU",  &mU,  "Number of dependent equality constriants" );
    clp.setOption( "hl", &hl, "Lower bounds on general inequalities" );
    clp.setOption( "hu", &hu, "Upper bounds on general inequalities" );
    clp.setOption( "diag-scal", &diag_scal, "Scaling of the basis diagonal" );
    clp.setOption( "diag-vary", &diag_vary, "Variation of the basis diagonal scaling" );
    clp.setOption(
      "nlp-selects-basis", "no-nlp-selects-basis", &nlp_selects_basis
      ,"Determine if the NLP will select basis" );
    clp.setOption(
      "sym-basis", "unsym-basis", &sym_basis
      ,"Determine if the basis is symmetric" );
    clp.setOption( "f_offset", &f_offset, "Constant offset for objective function" );
    clp.setOption( "co", &co, "Constant term in general equalities" );
    clp.setOption(
      "ignore-constraints", "no-ignore-constraints", &ignore_constraints
      ,"Determine if constraints are ignored or not" );
  
    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;
    
    if(show_options) {
      std::cout << "\nPrinting commandline options used (options used shown as (default: \?\?\?) ...\n\n";
      clp.printHelpMessage(argv[0],std::cout);
    }
    
    //
    // Create and solve the NLP
    //

    ExampleNLPBanded
      nlp(nD,nI,bw,mU,mI,xo,xDl,xDu,xIl,xIu,hl,hu
        ,nlp_selects_basis,diag_scal,diag_vary
        ,sym_basis,f_offset,co,ignore_constraints
        );

    solver.set_nlp( Teuchos::rcp(&nlp,false) );

    const MoochoSolver::ESolutionStatus
      solution_status = solver.solve_nlp();
    
    return solution_status;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,success)

  return MoochoSolver::SOLVE_RETURN_EXCEPTION;
}
