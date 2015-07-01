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
#include <iostream>
#include <math.h>

#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "NLPInterfacePack_NLPBarrier.hpp"
#include "MoochoPack_PostProcessBarrierLineSearch_Step.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

#define min(a,b) ( (a < b) ? a : b )
#define max(a,b) ( (a > b) ? a : b )

namespace MoochoPack {

PostProcessBarrierLineSearch_Step::PostProcessBarrierLineSearch_Step(
  Teuchos::RCP<NLPInterfacePack::NLPBarrier> barrier_nlp
  )
  :
  barrier_nlp_(barrier_nlp)
  {
  TEUCHOS_TEST_FOR_EXCEPTION(
    !barrier_nlp_.get(),
    std::logic_error,
    "PostProcessBarrierLineSearch_Step given NULL NLPBarrier."
    );
  }
  

bool PostProcessBarrierLineSearch_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
  {
  using Teuchos::dyn_cast;
  using IterationPack::print_algorithm_step;
  using AbstractLinAlgPack::Vp_StV;

  NLPAlgo            &algo   = dyn_cast<NLPAlgo>(_algo);
  IpState             &s      = dyn_cast<IpState>(_algo.state());

     EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();
  
  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
    }

  // Get iteration quantities...
  value_type& f_kp1 = s.f().set_k(+1);
  f_kp1 = barrier_nlp_->objective_term();

  VectorMutable& Gf_kp1 = s.Gf().set_k(+1);
  Gf_kp1 = *(barrier_nlp_->grad_objective_term());
  
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    {
    out << "\nf = " << f_kp1
      << "\nbarrier_term = " << barrier_nlp_->barrier_term() << std::endl;
    
    }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
    {
    out << "Gf = \n" << Gf_kp1
      << "\ngrad_barrier_term  = \n" << *(barrier_nlp_->grad_barrier_term());
    
    }
  return true;
  }


void PostProcessBarrierLineSearch_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
  {
  //const NLPAlgo   &algo = rsqp_algo(_algo);
  //const NLPAlgoState  &s    = algo.rsqp_state();
  out << L << "# Process out the temporary barrier term used for line search\n"
    << L << "f_k -= barrier_term_k\n";
  }
} // end namespace MoochoPack 
