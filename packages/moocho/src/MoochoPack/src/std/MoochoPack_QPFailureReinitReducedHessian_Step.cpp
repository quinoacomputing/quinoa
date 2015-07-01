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

#include <ostream>
#include <typeinfo>

#include "MoochoPack_QPFailureReinitReducedHessian_Step.hpp"
#include "MoochoPack_MoochoAlgorithmStepNames.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"

namespace MoochoPack {

QPFailureReinitReducedHessian_Step::QPFailureReinitReducedHessian_Step(
  const null_space_step_ptr_t& null_space_step
  )
  :null_space_step_(null_space_step)
  ,last_qp_failure_k_(-100) // has not failed yet.
{}

bool QPFailureReinitReducedHessian_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  try {
    return null_space_step().do_step(_algo,step_poss,type,assoc_step_poss);
  }
  catch(const QPFailure& qp_excpt) {
    NLPAlgo              &algo   = rsqp_algo(_algo);
    NLPAlgoState         &s      = algo.rsqp_state();
    EJournalOutputLevel  olevel  = algo.algo_cntr().journal_output_level();
    std::ostream         &out    = algo.track().journal_out();

    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
      out	<< "\nQP failed! "
        << " (k = " << algo.state().k() << ")\n"
        << "QPFailure description: " << qp_excpt.what() << "\n";
    }
    if( s.k() >= algo.algo_cntr().max_iter() ) {
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
        out	<< "\nThe maximum number of rSQP iterations\n"
          << " have been exceeded so quit "
          << " (k = " << algo.state().k() << ")\n";
      }
      algo.terminate(false);
      return false;
    }
    if( last_qp_failure_k_ == s.k() ) {
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
        out	<< "\nThe QP failed again even with a new reduced Hessian rHL_k!"
          << " (k = " << algo.state().k() << ")\n"
          << "We quit!\n";
      }
      throw;
    }
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
      out	<< "\nWiping out all memory for rHL and going back to reinitalize it ..."
        << " (k = " << algo.state().k() << ")\n";
    }
    last_qp_failure_k_ = s.k(); // Remember this for later!
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "Wipe out all update rHL_{k} for all k\n"
        << "goto ReducedHessian\n";
    }
    s.rHL().set_all_not_updated();
    algo.do_step_next( ReducedHessian_name );
    return false;
  }
  return false;	// will never be executed.
}

void QPFailureReinitReducedHessian_Step::print_step(
  const Algorithm& algo
  ,poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  ,std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "do null space step : " << typeName(null_space_step()) << std::endl;
  null_space_step().print_step(algo,step_poss,type,assoc_step_poss,out,L+"  ");
  out
    << L << "end null space step\n"
    << L << "if QPFailure was thrown then\n"
    << L << "  if QP failed already then\n"
    << L << "    rethrow QPFailure\n"
    << L << "  end\n"
    << L << "  if k > max_iter then\n"
    << L << "    terminate the algorithm!\n"
    << L << "  end\n"
    << L << "  set all rHL_{k} to not updated\n"
    << L << "  goto ReducedHessian\n"
    << L << "end\n"
    ;
}

} // end namespace MoochoPack
