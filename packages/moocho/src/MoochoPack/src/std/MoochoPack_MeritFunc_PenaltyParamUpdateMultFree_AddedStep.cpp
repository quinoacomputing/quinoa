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

#include <math.h>

#include <ostream>
#include <typeinfo>

#include "MoochoPack_MeritFunc_PenaltyParamUpdateMultFree_AddedStep.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"

namespace MoochoPack {

MeritFunc_PenaltyParamUpdateMultFree_AddedStep::MeritFunc_PenaltyParamUpdateMultFree_AddedStep(
  value_type    small_mu
  ,value_type   mult_factor
  ,value_type   kkt_near_sol
  )
  :MeritFunc_PenaltyParamUpdateGuts_AddedStep(small_mu,mult_factor,kkt_near_sol)
{}

// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

bool MeritFunc_PenaltyParamUpdateMultFree_AddedStep::min_mu(
  NLPAlgoState& s, value_type* min_mu
  ) const
{
  using AbstractLinAlgPack::dot;

  IterQuantityAccess<VectorMutable>
    &Gf_iq    = s.Gf(),
    &nu_iq    = s.nu(),
    &Ypy_iq   = s.Ypy(),
    &c_iq     = s.c();
  if ( Gf_iq.updated_k(0) && nu_iq.updated_k(0) && Ypy_iq.updated_k(0) && c_iq.updated_k(0) ) {
    // min_mu = abs((Gf_k+nu_k)'*Ypy_k) / norm(c_k,1)
    const value_type
      dot_Gf_Ypy = dot( Gf_iq.get_k(0), Ypy_iq.get_k(0) ),
      dot_nu_Ypy = dot( nu_iq.get_k(0), Ypy_iq.get_k(0) ),
      nrm_c      = c_iq.get_k(0).norm_1(),
      small_num  = std::numeric_limits<value_type>::min();
    *min_mu = ::fabs( dot_Gf_Ypy + dot_nu_Ypy ) / ( nrm_c + small_num );
    return true;
  }
  return false;
}

void MeritFunc_PenaltyParamUpdateMultFree_AddedStep::print_min_mu_step(
  std::ostream& out, const std::string& L ) const
{
  out
    << L << "if Gf_k, nu_k, Ypy_k and c_k are updated then\n"
    << L << "   min_mu = abs((Gf_k+nu_k)'*Ypy_k) / ( norm(c_k,1) + small_num )\n"
    << L << "   update_mu = true\n"
    << L << "else\n"
    << L << "   update_mu = false\n"
    << L << "endif\n"
    ;
}

}	// end namespace MoochoPack
