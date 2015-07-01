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
//

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "NLPInterfacePack_ExampleNLPDirectRun.hpp"
#include "NLPInterfacePack_ExampleNLPDirect.hpp"
#include "MoochoPack_NLPAlgoConfigMamaJama.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"

MoochoPack::MoochoSolver::ESolutionStatus
NLPInterfacePack::ExampleNLPDirectRun(
  const VectorSpace&   vec_space
  ,value_type          xo
  ,bool                has_bounds
  ,bool                dep_bounded
  ,std::ostream*       console_out
  ,std::ostream*       error_out
  ,bool                throw_solve_exception
  ,std::ostream*       algo_out
  ,std::ostream*       summary_out
  ,std::ostream*       journal_out
  )
{
  using std::endl;
  using std::setw;
  namespace rcp = MemMngPack;
  using Teuchos::RCP;
  namespace ofsp = OptionsFromStreamPack;
  using ofsp::OptionsFromStream;
  namespace rsqp = MoochoPack;
  using rsqp::MoochoSolver;
  using rsqp::NLPAlgoConfigMamaJama;

  MoochoSolver::ESolutionStatus
    solve_return = MoochoSolver::SOLVE_RETURN_EXCEPTION;

  int err = 0;
  
  int w = 15;
  int prec = 8;

  if(console_out)
    *console_out
      << std::setprecision(prec)
      << std::scientific
      << "***************************************************\n"
      << "*** Running Tests on ExampleNLPDirect ***\n"
      << "***************************************************\n"
      << "\nUsing a vector space of type \'" << typeName(vec_space) << "\'"
      << "\nwith a dimension of vec_space.dim() = " << vec_space.dim()
      << std::endl;

  // Create the nlp
  ExampleNLPDirect
    nlp(VectorSpace::space_ptr_t(&vec_space,false),xo,has_bounds,dep_bounded);

  // Create the solver object and set it up
  MoochoSolver solver;
  solver.set_nlp(Teuchos::rcp(&nlp,false));                  // Set the NLP!
  solver.set_error_handling(                             // set up outputting
    throw_solve_exception
    ,Teuchos::rcp(error_out,false)
    );
  solver.set_console_out(Teuchos::rcp(console_out,false));
  solver.set_summary_out(Teuchos::rcp(summary_out,false));
  solver.set_journal_out(Teuchos::rcp(journal_out,false));
  solver.set_algo_out(   Teuchos::rcp(algo_out,false)   );

  // Run MOOCHO using the MamaJama configuration
  solve_return = solver.solve_nlp();

  return solve_return;
}
