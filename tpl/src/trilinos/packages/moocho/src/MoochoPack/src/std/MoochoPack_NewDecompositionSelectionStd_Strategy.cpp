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

#include "MoochoPack_NewDecompositionSelectionStd_Strategy.hpp"
#include "MoochoPack_MoochoAlgorithmStepNames.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_NLPAlgoState.hpp"

namespace MoochoPack {

NewDecompositionSelectionStd_Strategy::NewDecompositionSelectionStd_Strategy(
  const decomp_sys_handler_ptr_t   &decomp_sys_handler
  )
  :decomp_sys_handler_(decomp_sys_handler)
{}

bool NewDecompositionSelectionStd_Strategy::new_decomposition(
  NLPAlgo& algo, Algorithm::poss_type step_poss
  ,IterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
  )
{
  NLPAlgoState               &s     = algo.rsqp_state();
  EJournalOutputLevel     olevel = algo.algo_cntr().journal_output_level();
  std::ostream&           out    = algo.track().journal_out();

  // Check to see if we have a decomposition system set
  if( !get_decomp_sys_handler().get() ) {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
      out	<< "\nWe are asked to select a new basis but there is no\n"
          "decomposition system set so we have no choice but to terminiate\n"
          "the algorithm"
        << " (k = " << algo.state().k() << ")\n";
    }
    algo.terminate(false);
    return false;
  }

  // We may get an infinite loop here so make sure we are under the max
  // number of iterations.
  if( s.k() >= algo.algo_cntr().max_iter() ) {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
      out	<< "\nThe maximum number of rSQP iterations\n"
        << " have been exceeded so quit "
        << " (k = " << algo.state().k() << ")\n";
    }
    algo.terminate(false);
    return false;
  }

  // Select a new decomposition
  decomp_sys_handler().select_new_decomposition(true);
  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    out	<< "x_kp1 = x_k\n"
      << "k=k+1\n"
      << "goto EvalNewPoint\n";
  }
  s.x().set_k(1) = s.x().get_k(0);
  s.alpha().set_k(0) = 0.0;	// Show that no step was taken.
  algo.track().output_iteration( algo );
  s.next_iteration();
  algo.do_step_next( EvalNewPoint_name );
  return false;
}

void NewDecompositionSelectionStd_Strategy::print_new_decomposition(
  const NLPAlgo& algo, Algorithm::poss_type step_poss
  ,IterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
  ,std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "if k > max_iter then\n"
    << L << "  terminate the algorithm\n"
    << L << "end\n"		
    << L << "Select a new basis at current point\n"
    << L << "x_kp1 = x_k\n"
    << L << "alpha_k = 0\n"
    << L << "k=k+1\n"
    << L << "goto EvalNewPoint\n";
}

}	// end namespace MoochoPack 
