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

#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "RTOp_TOp_assign_scalar.h"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_TOp_scale_vector.h"
#include "RTOpPack_RTOpC.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace {

// vector scalar assignment operator
RTOpPack::RTOpC& assign_scalar_op()
{
  static RTOpPack::RTOpC          assign_scalar_op_;
  return(assign_scalar_op_);
}
// vector assignment operator
static RTOpPack::RTOpC          assign_vec_op;
// scale vector
static RTOpPack::RTOpC          scale_vector_op;

// Simple class for an object that will initialize the operator objects
class init_rtop_server_t {
public:
  init_rtop_server_t() {
    // Vector scalar assignment operator
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_construct(0.0,&assign_scalar_op().op()));
    // Vector assignment operator
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_vectors_construct(&assign_vec_op.op()));
    // Operator scale_vector
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_scale_vector_construct(0.0,&scale_vector_op.op()));
  }
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace AbstractLinAlgPack {

// Clone

MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutable::mv_clone()
{
  multi_vec_mut_ptr_t
    new_mv = this->space_cols().create_members(this->cols());
  const MultiVector*  multi_vecs[]      = { this };
  MultiVectorMutable* targ_multi_vecs[] = { new_mv.get() };
  AbstractLinAlgPack::apply_op(APPLY_BY_COL,assign_vec_op,1,multi_vecs,1,targ_multi_vecs,NULL);
  return new_mv;
}

// Sub-view methods

MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutable::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng)
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: return a MultiVectorMutableSubView object.
  // Note that the MultiVectorMutableSubView class should derive from
  // MultiVectorSubView.
  return Teuchos::null;
}

// Overridden from MatrixOp

void MultiVectorMutable::zero_out()
{
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_set_alpha(0.0,&assign_scalar_op().op()));
  MultiVectorMutable* targ_multi_vecs[] = { this };
  AbstractLinAlgPack::apply_op(APPLY_BY_COL,assign_scalar_op(),0,NULL,1,targ_multi_vecs,NULL);
}

void MultiVectorMutable::Mt_S( value_type alpha )
{
  if( alpha == 0.0 ) {
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_set_alpha(alpha,&assign_scalar_op().op()));
    MultiVectorMutable* targ_multi_vecs[] = { this };
    AbstractLinAlgPack::apply_op(APPLY_BY_COL,assign_scalar_op(),0,NULL,1,targ_multi_vecs,NULL);
  }
  else if( alpha != 1.0 ) {
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_scale_vector_set_alpha(alpha,&scale_vector_op.op()));
    MultiVectorMutable* targ_multi_vecs[] = { this };
    AbstractLinAlgPack::apply_op(APPLY_BY_COL,scale_vector_op,0,NULL,1,targ_multi_vecs,NULL);
  }
}

MatrixOp& MultiVectorMutable::operator=(const MatrixOp& mwo_rhs)
{
  const MultiVector *mv_rhs = dynamic_cast<const MultiVector*>(&mwo_rhs);
  if(mv_rhs) {
    const MultiVector*  multi_vecs[]      = { mv_rhs };
    MultiVectorMutable* targ_multi_vecs[] = { this };
    AbstractLinAlgPack::apply_op(APPLY_BY_COL,assign_vec_op,1,multi_vecs,1,targ_multi_vecs,NULL);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Get column by column or row by row
  }
  return *this;
}

MatrixOp::mat_mut_ptr_t
MultiVectorMutable::clone()
{
  return this->mv_clone();
}

bool MultiVectorMutable::Mp_StM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs
  ) const
{
  return false; // ToDo: Specialize!
}

bool MultiVectorMutable::Mp_StM(
  value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs
  )
{
  return false; // ToDo: Specialize!
}

// Overridden form MultiVector

MultiVector::multi_vec_ptr_t MultiVectorMutable::mv_clone() const
{
  return const_cast<MultiVectorMutable*>(this)->mv_clone();
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::col(index_type j) const
{
  return const_cast<MultiVectorMutable*>(this)->col(j);
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::row(index_type i) const
{
  return const_cast<MultiVectorMutable*>(this)->row(i);
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::diag(int k) const
{
  return const_cast<MultiVectorMutable*>(this)->diag(k);
}

MultiVectorMutable::multi_vec_ptr_t
MultiVectorMutable::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  return const_cast<MultiVectorMutable*>(this)->mv_sub_view(row_rng,col_rng);
}

} // end namespace AbstractLinAlgPack
