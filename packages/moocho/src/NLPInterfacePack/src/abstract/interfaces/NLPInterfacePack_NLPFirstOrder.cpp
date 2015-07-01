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

#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "Teuchos_Assert.hpp"

namespace {
  const char name_Gc[] = "Gc";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPFirstOrder::NLPFirstOrder()
  : Gc_(NULL)
{}

void NLPFirstOrder::initialize(bool test_setup) {
  num_Gc_evals_ = 0;
  NLPObjGrad::initialize(test_setup);
}

// BasisSystem

const NLPFirstOrder::basis_sys_ptr_t
NLPFirstOrder::basis_sys() const
{
  return Teuchos::null;
}

// <<std aggr>> members for Gc

void NLPFirstOrder::set_Gc(MatrixOp* Gc)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
  Gc_ = Gc;
}

MatrixOp* NLPFirstOrder::get_Gc()
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
  return StandardCompositionRelationshipsPack::get_role_name(Gc_, false, name_Gc);
}

MatrixOp& NLPFirstOrder::Gc()
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
  return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

const MatrixOp& NLPFirstOrder::Gc() const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
  return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

void NLPFirstOrder::unset_quantities()
{
  NLPObjGrad::unset_quantities();
  Gc_ = NULL;
}

// calculations

void NLPFirstOrder::calc_Gc(const Vector& x, bool newx) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
  StandardCompositionRelationshipsPack::assert_role_name_set(Gc_, "NLP::calc_Gc()", name_Gc);
  imp_calc_Gc(x,newx,first_order_info());
  num_Gc_evals_++;
}

size_type NLPFirstOrder::num_Gc_evals() const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
  return num_Gc_evals_;
}

}	// end namespace NLPInterfacePack 
