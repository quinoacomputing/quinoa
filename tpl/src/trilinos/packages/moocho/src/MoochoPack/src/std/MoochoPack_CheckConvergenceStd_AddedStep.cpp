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

//#include <assert.h>

#include <ostream>

#include "MoochoPack_CheckConvergenceStd_AddedStep.hpp"
#include "MoochoPack_NLPAlgoContainer.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"

namespace MoochoPack {

CheckConvergenceStd_AddedStep::CheckConvergenceStd_AddedStep(
  Teuchos::RCP<CheckConvergence_Strategy> convergence_strategy
  )
  :
  convergence_strategy_(convergence_strategy)
  {}

bool CheckConvergenceStd_AddedStep::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
  {
  
    TEUCHOS_TEST_FOR_EXCEPTION(!convergence_strategy_.get(),
          std::logic_error,
          "Don't have a valid convergence_strategy in CheckConvergenceStd_AddedStep\n"
    );

  NLPAlgo	&algo	  = rsqp_algo(_algo);

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  const bool found_solution = convergence_strategy_->Converged(_algo);

  if( found_solution )
  {
    if( static_cast<int>(olevel) > static_cast<int>(PRINT_NOTHING) )
      out	<< "\nJackpot!  Found the solution!!!!!! (k = " << algo.state().k() << ")\n";
    algo.terminate(true);	// found min
    return false; // skip the other steps and terminate
  }
  
  if( static_cast<int>(olevel) > static_cast<int>(PRINT_NOTHING) )
    out	<< "\nHave not found the solution yet, have to keep going (k = " << algo.state().k() << ") :-(\n";
  
  // We are not at the solution so keep going
  return true;
  }

void CheckConvergenceStd_AddedStep::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
  {

    TEUCHOS_TEST_FOR_EXCEPTION(!convergence_strategy_.get(),
          std::logic_error,
          "Don't have a valid convergence_strategy in CheckConvergenceStd_AddedStep\n"
    );
  
  convergence_strategy_->print_step(algo, out, L);
  }

}	// end namespace MoochoPack
