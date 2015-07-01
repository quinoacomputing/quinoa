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

#include "NLPInterfacePack_NLPSecondOrder.hpp"
#include "Teuchos_Assert.hpp"

namespace {
  const char name_HL[] = "HL";
}

namespace NLPInterfacePack {

// constructors

NLPSecondOrder::NLPSecondOrder()
  : HL_(NULL)
{}


void NLPSecondOrder::initialize(bool test_setup) {
  num_HL_evals_ = 0;
  NLPFirstOrder::initialize(test_setup);
}

// <<std aggr>> members for HL

void NLPSecondOrder::set_HL(MatrixSymOp* HL)
{
  HL_ = HL;
}

MatrixSymOp* NLPSecondOrder::get_HL()
{
  return StandardCompositionRelationshipsPack::get_role_name(HL_, false, name_HL);
}

MatrixSymOp& NLPSecondOrder::HL()
{
  return StandardCompositionRelationshipsPack::role_name(HL_, false, name_HL);
}

const MatrixSymOp& NLPSecondOrder::HL() const
{
  return StandardCompositionRelationshipsPack::role_name(HL_, false, name_HL);
}

void NLPSecondOrder::unset_quantities()
{
  NLPFirstOrder::unset_quantities();
  HL_ = NULL;
}

// calculations

void NLPSecondOrder::calc_HL(
  const Vector& x, const Vector* lambda, bool newpoint
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( lambda  && this->m()  == 0, std::logic_error, "" );
#endif
  StandardCompositionRelationshipsPack::assert_role_name_set(HL_, "NLP::calc_HL()", name_HL);
  imp_calc_HL(x,lambda,newpoint,second_order_info());
}

size_type NLPSecondOrder::num_HL_evals() const
{
  return num_HL_evals_;
}

} // namespace NLPInterfacePack
