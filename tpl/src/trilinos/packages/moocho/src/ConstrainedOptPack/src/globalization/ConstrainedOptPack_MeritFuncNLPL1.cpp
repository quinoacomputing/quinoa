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

#include "ConstrainedOptPack_MeritFuncNLPL1.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace ConstrainedOptPack {

MeritFuncNLPL1::MeritFuncNLPL1()
  : deriv_(0.0), mu_(0.0)
{}

// Overridden from MeritFuncNLP

MeritFuncNLP& MeritFuncNLPL1::operator=(const MeritFuncNLP& merit_func)
{
  using Teuchos::dyn_cast;
  const MeritFuncNLPL1 &merit_func_l1 = dyn_cast<const MeritFuncNLPL1>(merit_func);
  if(this == &merit_func_l1)
    return *this; // assignment to self
  this->deriv_   = merit_func_l1.deriv_;
  this->mu_      = merit_func_l1.mu_;
  return *this;
}

value_type MeritFuncNLPL1::value(
  value_type       f
  ,const Vector    *c
  ,const Vector    *h
  ,const Vector    *hl
  ,const Vector    *hu
  ) const
{
  value_type phi_val_h = 0.0;
  if(h) {
    size_type   max_viol_i;
    value_type  max_viol, h_i, hlu_i;
    int         bnd_type;
    AbstractLinAlgPack::max_inequ_viol(*h,*hl,*hu,&max_viol_i,&max_viol,&h_i,&bnd_type,&hlu_i);
    if(max_viol_i) phi_val_h += mu_ * fabs(h_i - hlu_i);
  }
  return f + ( c ? mu_ * c->norm_1() : 0.0) + phi_val_h;
}

value_type MeritFuncNLPL1::deriv() const
{
  return deriv_;
}

void MeritFuncNLPL1::print_merit_func(std::ostream& out
  , const std::string& L ) const
{
  out
    << L << "*** Define L1 merit funciton (assumes Gc_k'*d_k + c_k = 0):\n"
    << L << "phi(f,c) = f + mu_k*( norm(c,1) + max_viol( hl <= h <= hu ) )\n"
    << L << "Dphi(x_k,d_k) = Gf_k' * d_k - mu*( norm(c_k,1)  + max_viol( hl <= h_k <= hu ) )\n";
}

// Overridden from MeritFuncNLPDirecDeriv

value_type MeritFuncNLPL1::calc_deriv(
  const Vector    &Gf_k
  ,const Vector   *c_k
  ,const Vector   *h_k
  ,const Vector   *hl
  ,const Vector   *hu
  ,const Vector   &d_k
  )
{
  using AbstractLinAlgPack::dot;
  TEUCHOS_TEST_FOR_EXCEPTION(
    h_k || hl || hu, std::logic_error
    ,"MeritFuncNLPL1::value(...) : Error! general inequalities are not supported yet" );
  return deriv_ = dot( Gf_k, d_k ) - ( c_k ? mu_ * c_k->norm_1() : 0.0 );
}

// Overridden from MeritFuncPenaltyParam

void MeritFuncNLPL1::mu(value_type mu)
{
  mu_ = mu;
}

value_type MeritFuncNLPL1::mu() const
{
  return mu_;
}

}	// end namespace ConstrainedOptPack 
