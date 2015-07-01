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

#include <limits>
#include <ostream>

#include "MoochoPack_MeritFunc_DummyUpdate_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

bool MeritFunc_DummyUpdate_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  )
{
  using Teuchos::dyn_cast;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  IterQuantityAccess<MeritFuncNLP>
    &merit_func_nlp_iq = s.merit_func_nlp();

  if(!merit_func_nlp_iq.updated_k(0)) {
    const int last_updated_k = merit_func_nlp_iq.last_updated();
    MeritFuncNLP
      &merit_func_nlp_k =  ( last_updated_k != IterQuantity::NONE_UPDATED
                   ? merit_func_nlp_iq.set_k(0,last_updated_k)
                   : merit_func_nlp_iq.set_k(0) );
    MeritFuncNLPDirecDeriv
      &direc_deriv = dyn_cast<MeritFuncNLPDirecDeriv>(merit_func_nlp_k);

    direc_deriv.calc_deriv(
      s.Gf().get_k(0)
      ,NULL
      ,NULL
      ,NULL
      ,NULL
      ,s.d().get_k(0)
      );

  }

  return true;
}

void MeritFunc_DummyUpdate_Step::print_step(
  const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Simply sets the current merit function value for unconstrained linesearch\n"
    << L << "if merit_func_nlp_k not updated set merit_func_nlp_k = merit_func_nlp_k(last_updated)\n";
}

} // end namespace MoochoPack
