#if 0

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

#include "MoochoPack_MeritFunc_PenaltyParamUpdateWithMult_AddedStep.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"

namespace MoochoPack {

MeritFunc_PenaltyParamUpdateWithMult_AddedStep::MeritFunc_PenaltyParamUpdateWithMult_AddedStep(
      const merit_func_ptr_t& merit_func, value_type small_mu
    , value_type mult_factor, value_type kkt_near_sol )
  : MeritFunc_PenaltyParamUpdateGuts_AddedStep(merit_func,small_mu,mult_factor,kkt_near_sol)
{}

// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

bool MeritFunc_PenaltyParamUpdateWithMult_AddedStep::min_mu(
  NLPAlgoState& s, value_type* min_mu ) const
{
  if ( s.lambda().updated_k(0) ) {
    *min_mu = s.lambda().get_k(0).norm_inf();
    return true;
  }
  return false;
}

void MeritFunc_PenaltyParamUpdateWithMult_AddedStep::print_min_mu_step(
  std::ostream& out, const std::string& L ) const
{
  out
    << L << "if lambda_k is updated then\n"
    << L << "    min_mu = norm( lambda_k, inf )\n"
    << L << "    update_mu = true\n"
    << L << "else\n"
    << L << "    update_mu = false\n"
    << L << "endif\n"
    ;
}

}	// end namespace MoochoPack

#endif // 0
