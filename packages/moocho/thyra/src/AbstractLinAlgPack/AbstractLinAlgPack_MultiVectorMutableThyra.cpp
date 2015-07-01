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

#include <assert.h>

#include "AbstractLinAlgPack_MultiVectorMutableThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

MultiVectorMutableThyra::MultiVectorMutableThyra()
{}

MultiVectorMutableThyra::MultiVectorMutableThyra(
  const Teuchos::RCP<Thyra::MultiVectorBase<value_type> >& thyra_multi_vec
  )
{
  this->initialize(thyra_multi_vec);
}

void MultiVectorMutableThyra::initialize(
  const Teuchos::RCP<Thyra::MultiVectorBase<value_type> >& thyra_multi_vec
  )
{
  namespace mmp = MemMngPack;
  TEUCHOS_TEST_FOR_EXCEPTION(
    thyra_multi_vec.get()==NULL, std::invalid_argument
    ,"MultiVectorMutableThyra::initialize(thyra_multi_vec): Error!"
    );
  MatrixOpThyra::initialize(thyra_multi_vec);
}

Teuchos::RCP<Thyra::MultiVectorBase<value_type> > 
MultiVectorMutableThyra::set_uninitialized()
{
  Teuchos::RCP<Thyra::MultiVectorBase<value_type> >
    tmp_thyra_multi_vec = cast_thyra_multi_vec();
  MatrixOpThyra::set_uninitialized();
  return tmp_thyra_multi_vec;
}

Teuchos::RCP<const Thyra::MultiVectorBase<value_type> >
MultiVectorMutableThyra::thyra_multi_vec() const
{
  return Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<value_type> >(this->thyra_linear_op());
}

// Overridden from MatrixOpThyra

void MultiVectorMutableThyra::initialize(
  const Teuchos::RCP<const Thyra::LinearOpBase<value_type> >& thyra_linear_op
  )
{
  namespace mmp = MemMngPack;
  this->initialize(
    Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<value_type> >(
      Teuchos::rcp_const_cast<Thyra::LinearOpBase<value_type> >(thyra_linear_op)
      )
    );
}

// Overridden from MatrixOp

MatrixOp::mat_mut_ptr_t
MultiVectorMutableThyra::clone()
{
  return this->MatrixOpThyra::clone();
}

MatrixOp& MultiVectorMutableThyra::operator=(const MatrixOp& mwo_rhs)
{
  return this->MultiVectorMutable::operator=(mwo_rhs);
}

void MultiVectorMutableThyra::Vp_StMtV(
  VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  ,const Vector& v_rhs2, value_type beta
  ) const
{
  this->MatrixOpThyra::Vp_StMtV(v_lhs,alpha,trans_rhs1,v_rhs2,beta);
}

bool MultiVectorMutableThyra::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ,value_type beta
  ) const
{
  return this->MatrixOpThyra::Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta);
}

// Overridden from MultiVector

MultiVector::access_by_t
MultiVectorMutableThyra::access_by() const
{
  return COL_ACCESS;
}

void MultiVectorMutableThyra::apply_op(
  EApplyBy apply_by, const RTOpPack::RTOp& primary_op
  ,const size_t num_multi_vecs,      const MultiVector*   multi_vecs[]
  ,const size_t num_targ_multi_vecs, MultiVectorMutable*  targ_multi_vecs[]
  ,RTOpPack::ReductTarget* reduct_objs[]
  ,const index_type primary_first_ele,   const index_type primary_sub_dim, const index_type primary_global_offset
  ,const index_type secondary_first_ele, const index_type secondary_sub_dim
  ) const
{
  MultiVector::apply_op(
    apply_by,primary_op,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
    ,reduct_objs
    ,primary_first_ele,primary_sub_dim,primary_global_offset
    ,secondary_first_ele,secondary_sub_dim
    ); // ToDo: Specialize!
}

void MultiVectorMutableThyra::apply_op(
  EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
  ,const size_t num_multi_vecs,      const MultiVector*   multi_vecs[]
  ,const size_t num_targ_multi_vecs, MultiVectorMutable*  targ_multi_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type primary_first_ele, const index_type primary_sub_dim, const index_type primary_global_offset
  ,const index_type secondary_first_ele, const index_type secondary_sub_dim
  ) const
{
  MultiVector::apply_op(
    apply_by,primary_op,secondary_op,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
    ,reduct_obj
    ,primary_first_ele,primary_sub_dim,primary_global_offset
    ,secondary_first_ele,secondary_sub_dim
    ); // ToDo: Specialize!
}

// Overridden from MultiVectorMutable

MultiVectorMutable::vec_mut_ptr_t
MultiVectorMutableThyra::col(index_type j)
{
  return Teuchos::rcp(new VectorMutableThyra(cast_thyra_multi_vec()->col(j-1)));
}

MultiVectorMutable::vec_mut_ptr_t
MultiVectorMutableThyra::row(index_type i)
{
  return Teuchos::null;
}

MultiVectorMutable::vec_mut_ptr_t
MultiVectorMutableThyra::diag(int k)
{
  return Teuchos::null;
}

MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutableThyra::mv_sub_view(const Range1D& row_rng_in, const Range1D& col_rng_in)
{
  const index_type  this_rows = this->rows();
  const Range1D     row_rng = RangePack::full_range(row_rng_in,1,this->rows());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(row_rng.lbound()==1 && row_rng.ubound()==this_rows), std::invalid_argument
    ,"MultiVectorMutableThyra::mv_sub_view(thyra_multi_vec): Error, can not handle subviews of the"
    " elements in a row yet!"
    );
  return Teuchos::rcp(new MultiVectorMutableThyra(cast_thyra_multi_vec()->subView(convert(col_rng_in))));
}

// private

Teuchos::RCP<Thyra::MultiVectorBase<value_type> >
MultiVectorMutableThyra::cast_thyra_multi_vec()
{
  namespace mmp = MemMngPack;
  return Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<value_type> >(
    Teuchos::rcp_const_cast<Thyra::LinearOpBase<value_type> >(this->thyra_linear_op())
    );
}

} // end namespace AbstractLinAlgPack
