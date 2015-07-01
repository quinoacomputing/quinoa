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
#include <iostream>
#include <limits>

#include "NLPInterfacePack_NLPBarrier.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "Teuchos_Assert.hpp"

namespace NLPInterfacePack {

NLPBarrier::NLPBarrier()
  :
  barrier_term_(0.0),
  objective_term_(0.0),
  nlp_(Teuchos::null)
  {
  }


void NLPBarrier::InitializeFromNLP(
  Teuchos::RCP<NLP> original_nlp
  )
  {
  TEUCHOS_TEST_FOR_EXCEPTION(
    !original_nlp.get(),
    std::logic_error,
    "null nlp passed to NLPBarrier decorator"
    );

  nlp_ = Teuchos::rcp_dynamic_cast<NLPObjGrad>(original_nlp);

  TEUCHOS_TEST_FOR_EXCEPTION(
    !nlp_.get(),
    std::logic_error,
    "non NLPObjGrad NLP passed to NLPBarrier decorator"
    );
  }

void NLPBarrier::mu(const value_type mu)
  {
  mu_ = mu;
  }

value_type NLPBarrier::barrier_term() const
  {
  return barrier_term_;
  }

value_type NLPBarrier::objective_term() const
  {
  return objective_term_;
  }

const Teuchos::RCP<Vector> NLPBarrier::grad_barrier_term() const
  {
  return grad_barrier_term_;
  }

const Teuchos::RCP<Vector>  NLPBarrier::grad_objective_term() const
  {
  return grad_objective_term_;
  }


void NLPBarrier::calc_f(const Vector& x, bool newx) const
  {
  nlp_->calc_f(x, newx);
  value_type* f_val = nlp_->get_f();

  objective_term_ = *f_val;
  barrier_term_   = CalculateBarrierTerm(x);

  (*f_val) += barrier_term_;
  }

void NLPBarrier::calc_Gf(const Vector& x, bool newx) const
  {
  using AbstractLinAlgPack::inv_of_difference;

     nlp_->calc_Gf(x, newx);
  grad_objective_term_ = nlp_->get_Gf()->clone();

  //std::cout << "grad_objective_term=\n";
  //grad_objective_term_->output(std::cout);

  if (!grad_barrier_term_temp_.get())
    { grad_barrier_term_temp_ = grad_objective_term_->space().create_member(); }

  if (!grad_barrier_term_.get())
    { grad_barrier_term_ = grad_objective_term_->space().create_member(); }	

  *grad_barrier_term_temp_ = 0.0;
  *grad_barrier_term_ = 0.0;	

  inv_of_difference(mu_, nlp_->xu(), x, grad_barrier_term_.get());
   //std::cout << "mu*invXU=\n";
  //grad_barrier_term_->output(std::cout);

  inv_of_difference(mu_, x, nlp_->xl(), grad_barrier_term_temp_.get());
   //std::cout << "mu*invXL=\n";
  //grad_barrier_term_temp_->output(std::cout);

  grad_barrier_term_->axpy(-1.0, *grad_barrier_term_temp_);

  nlp_->get_Gf()->axpy(1.0, *grad_barrier_term_);

  //std::cout << "grad_objective_term with barrier=\n";
  //nlp_->get_Gf()->output(std::cout);
  }

void NLPBarrier::imp_calc_f(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
  {
  TEUCHOS_TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }

void NLPBarrier::imp_calc_c(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
  {
  TEUCHOS_TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }

void NLPBarrier::imp_calc_c_breve(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info_breve
  ) const
  {	
  TEUCHOS_TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }

void NLPBarrier::imp_calc_h_breve(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info_breve
  ) const
  {	
  TEUCHOS_TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }

void NLPBarrier::imp_calc_Gf(
  const Vector& x,
  bool newx, 
  const ObjGradInfo& obj_grad_info
  ) const
  {
  TEUCHOS_TEST_FOR_EXCEPT( !( false && !"This should never get called." ) );
  }


value_type NLPBarrier::CalculateBarrierTerm(const Vector& x) const
  {
  using AbstractLinAlgPack::log_bound_barrier;
  barrier_term_ = log_bound_barrier(x, xl(), xu());
//	std::cerr << "NLPBarrier::CalculateBarrierTerm(x) : (1) barrier_term_ = " << barrier_term_ << std::endl;
  barrier_term_ *= -mu_;
//	std::cerr << "NLPBarrier::CalculateBarrierTerm(x) : (2) barrier_term_ = " << barrier_term_ << std::endl;
  return barrier_term_;
  }

}	// end namespace NLPInterfacePack
