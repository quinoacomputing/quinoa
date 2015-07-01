/*
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
*/

// /////////////////////////////////////////////////////////////////////////
// MultiVector.cpp

#include <assert.h>

#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_MatrixSymDiag.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_AssertOp.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"

namespace {

// Map from EApplyBy to Transp
inline
BLAS_Cpp::Transp
to_trans(AbstractLinAlgPack::EApplyBy apply_by)
{
  return ( apply_by == AbstractLinAlgPack::APPLY_BY_ROW
      ? BLAS_Cpp::no_trans
      : BLAS_Cpp::trans
      );
}

// Return a row or a column vector from a multi-vector

inline 
AbstractLinAlgPack::MultiVector::vec_ptr_t
vec(
  const AbstractLinAlgPack::MultiVector&      multi_vec
  ,const AbstractLinAlgPack::size_type        k
  ,AbstractLinAlgPack::EApplyBy               apply_by
  )
{
  return ( apply_by == AbstractLinAlgPack::APPLY_BY_ROW
      ? multi_vec.row(k)
      : multi_vec.col(k)
      );
}

inline 
AbstractLinAlgPack::MultiVectorMutable::vec_mut_ptr_t
vec(
  AbstractLinAlgPack::MultiVectorMutable*         multi_vec
  ,const AbstractLinAlgPack::size_type            k
  ,AbstractLinAlgPack::EApplyBy                   apply_by
  )
{
  return ( apply_by == AbstractLinAlgPack::APPLY_BY_ROW
      ? multi_vec->row(k)
      : multi_vec->col(k)
      );
}

// Implement a matrix-matrix multiplication with a diagonal matrix.
//
// op(C) = b*op(C) + a*D*op(B)
//
bool mat_vec(
  const AbstractLinAlgPack::value_type        &a
  ,const AbstractLinAlgPack::MatrixOp         &D_mwo  // Diagonal matrix?
  ,BLAS_Cpp::Transp                           D_trans
  ,const AbstractLinAlgPack::MultiVector      &B
  ,BLAS_Cpp::Transp                           B_trans
  ,const AbstractLinAlgPack::value_type       &b
  ,BLAS_Cpp::Transp                           C_trans
  ,AbstractLinAlgPack::MatrixOp               *C
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;

  typedef AbstractLinAlgPack::MultiVector          MV;
  typedef AbstractLinAlgPack::MultiVectorMutable   MVM;
  using AbstractLinAlgPack::size_type;
  using AbstractLinAlgPack::Vector;
  using AbstractLinAlgPack::MatrixOp;
  using AbstractLinAlgPack::MultiVectorMutable;
  using AbstractLinAlgPack::MatrixSymDiag;
  using AbstractLinAlgPack::ele_wise_prod;
  using LinAlgOpPack::Vt_S;
  
  AbstractLinAlgPack::Mp_MtM_assert_compatibility(C,C_trans,D_mwo,D_trans,B,B_trans);

  MultiVectorMutable
    *Cmv = dynamic_cast<MultiVectorMutable*>(C);
  const MatrixSymDiag
    *D = dynamic_cast<const MatrixSymDiag*>(&D_mwo);
  if( !Cmv || !D || !(Cmv->access_by() & ( C_trans == no_trans ? MV::COL_ACCESS : MV::ROW_ACCESS ))
    || !(B.access_by() & ( B_trans == no_trans ? MV::COL_ACCESS : MV::ROW_ACCESS ))
    )
  {
    return false;
  }
  //
  // op(C).col(j) = b*op(C).col(j) + a*ele_wise_prod(D_diag,op(B).col(j)), for j = 1...op(C).cols()
  //
  const Vector  &D_diag = D->diag();
  const size_type
    opC_cols = BLAS_Cpp::cols( Cmv->rows(), Cmv->cols(), C_trans );
  for( size_type j = 1; j <= opC_cols; ++j ) {
    MV::vec_ptr_t
      opB_col_j = ( B_trans == no_trans ? B.col(j)    : B.row(j) );
    MVM::vec_mut_ptr_t
      opC_col_j = ( C_trans == no_trans ? Cmv->col(j) : Cmv->row(j) );
    Vt_S( opC_col_j.get(), b );
    ele_wise_prod( a, D_diag, *opB_col_j, opC_col_j.get() );
  }	
  return true;
}

} // end namespace

namespace AbstractLinAlgPack {

MultiVector::multi_vec_ptr_t
MultiVector::mv_clone() const
{
  return Teuchos::null;
}

MultiVector::multi_vec_ptr_t
MultiVector::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: return a MultiVectorSubView object.
  // Note that the MultiVectorSubView class should derive from MatrixOpSubView
  // so that a client can rely on the MatrixOpSubView interface.
  return Teuchos::null;
}

void MultiVector::apply_op(
  EApplyBy apply_by, const RTOpPack::RTOp& prim_op
  ,const size_t num_multi_vecs,      const MultiVector*   multi_vecs[]
  ,const size_t num_targ_multi_vecs, MultiVectorMutable*  targ_multi_vecs[]
  ,RTOpPack::ReductTarget* reduct_objs[]
  ,const index_type prim_first_ele_in, const index_type prim_sub_dim_in, const index_type prim_global_offset_in
  ,const index_type sec_first_ele_in, const index_type sec_sub_dim_in
  ) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  // ToDo: Validate the input!

  // Get the primary and secondary dimmensions.
  const index_type
    prim_dim     = ( apply_by == APPLY_BY_ROW ? rows()          : cols()   ),
    sec_dim      = ( apply_by == APPLY_BY_ROW ? cols()          : rows()   ),
    prim_sub_dim = ( prim_sub_dim_in != 0     ? prim_sub_dim_in : prim_dim ),
    sec_sub_dim  = ( sec_sub_dim_in != 0      ? sec_sub_dim_in  : sec_dim  );
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 < prim_sub_dim && prim_sub_dim <= prim_dim  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 < sec_sub_dim  && sec_sub_dim  <= sec_dim  ) );

  //
  // Apply the reduction/transformation operator and trnasform the target
  // vectors and reduce each of the reduction objects.
  //

  Workspace<MultiVector::vec_ptr_t>             vecs_s(wss,num_multi_vecs);
  Workspace<const Vector*>                      vecs(wss,num_multi_vecs);
  Workspace<MultiVectorMutable::vec_mut_ptr_t>  targ_vecs_s(wss,num_targ_multi_vecs);
  Workspace<VectorMutable*>                     targ_vecs(wss,num_targ_multi_vecs);

  {for(size_type j = sec_first_ele_in; j <= sec_first_ele_in - 1 + sec_sub_dim; ++j) {
    // Fill the arrays of vector arguments 
    {for(size_type k = 0; k < num_multi_vecs; ++k) {
      vecs_s[k] = vec( *multi_vecs[k], j, apply_by );
      vecs[k] = vecs_s[k].get();
    }}
    {for(size_type k = 0; k < num_targ_multi_vecs; ++k) {
      targ_vecs_s[k] = vec( targ_multi_vecs[k], j, apply_by );
      targ_vecs[k] = targ_vecs_s[k].get();
    }}
    // Apply the reduction/transformation operator
    AbstractLinAlgPack::apply_op(
      prim_op
      ,num_multi_vecs,      num_multi_vecs      ? &vecs[0]      : NULL
      ,num_targ_multi_vecs, num_targ_multi_vecs ? &targ_vecs[0] : NULL
      ,reduct_objs ? reduct_objs[j-1] : NULL
      ,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
      );
  }}

  // At this point all of the designated targ vectors in the target multi-vectors have
  // been transformed and all the reduction objects in reduct_obj[] have accumulated
  // the reductions.

}

void MultiVector::apply_op(
  EApplyBy apply_by, const RTOpPack::RTOp& prim_op, const RTOpPack::RTOp& sec_op
  ,const size_t num_multi_vecs,      const MultiVector*   multi_vecs[]
  ,const size_t num_targ_multi_vecs, MultiVectorMutable*  targ_multi_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type prim_first_ele_in, const index_type prim_sub_dim_in, const index_type prim_global_offset_in
  ,const index_type sec_first_ele_in, const index_type sec_sub_dim_in
  ) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  // ToDo: Validate the input!

  // Get the primary and secondary dimmensions.
  const index_type
    prim_dim    = ( apply_by == APPLY_BY_ROW ? rows()         : cols()  ),
    sec_dim     = ( apply_by == APPLY_BY_ROW ? cols()         : rows()  ),
    sec_sub_dim = ( sec_sub_dim_in != 0      ? sec_sub_dim_in : sec_dim );
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 < sec_sub_dim && sec_sub_dim <= sec_dim  ) );

  // Create a temporary buffer for the reduction objects of the primary reduction
  // so that we can call the companion version of this method.
  Workspace<Teuchos::RCP<RTOpPack::ReductTarget> >   rcp_reduct_objs(wss,sec_sub_dim);
  Workspace<RTOpPack::ReductTarget*>                         reduct_objs(wss,sec_sub_dim);
  for(index_type k = 0; k < sec_sub_dim; ++k) {
    rcp_reduct_objs[k] = prim_op.reduct_obj_create();
    reduct_objs[k] = &*rcp_reduct_objs[k];
  }
  
  // Call the campanion version that accepts an array of reduction objects
  this->apply_op(
    apply_by, prim_op
    ,num_multi_vecs,       multi_vecs
    ,num_targ_multi_vecs,  targ_multi_vecs
    ,&reduct_objs[0]
    ,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
    ,sec_first_ele_in,  sec_sub_dim_in
    );

  // Reduce all the reduction objects using the secondary reduction operator
  // into one reduction object and free the intermedate reduction objects.
  for(index_type k = 0; k < sec_sub_dim; ++k) {
    sec_op.reduce_reduct_objs( *reduct_objs[k], Teuchos::ptr(reduct_obj) );
  }
}

// Overridden form MatrixOp

MatrixOp::mat_ptr_t
MultiVector::clone() const
{
  return this->mv_clone();
}

MatrixOp::mat_ptr_t
MultiVector::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  return mv_sub_view(row_rng,col_rng);
}

bool MultiVector::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ,value_type beta
  ) const
{
  return mat_vec(
    alpha
    ,mwo_rhs1,trans_rhs1
    ,*this,trans_rhs2
    ,beta,BLAS_Cpp::no_trans,mwo_lhs
    );
}

bool MultiVector::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ,value_type beta
  ) const
{
  return mat_vec(
    alpha
    ,mwo_rhs2,BLAS_Cpp::trans_not(trans_rhs2)
    ,*this,BLAS_Cpp::trans_not(trans_rhs1)
    ,beta,BLAS_Cpp::trans,mwo_lhs
    );
}

} // end namespace AbstractLinAlgPack

// nonmembers

void AbstractLinAlgPack::apply_op(
  EApplyBy                        apply_by
  ,const RTOpPack::RTOp           &primary_op
  ,const size_t                   num_multi_vecs
  ,const MultiVector*             multi_vecs[]
  ,const size_t                   num_targ_multi_vecs
  ,MultiVectorMutable*            targ_multi_vecs[]
  ,RTOpPack::ReductTarget*        reduct_objs[]
  ,const index_type               primary_first_ele
  ,const index_type               primary_sub_dim
  ,const index_type               primary_global_offset
  ,const index_type               secondary_first_ele
  ,const index_type               secondary_sub_dim
  )
{
  if(num_multi_vecs)
    multi_vecs[0]->apply_op(
      apply_by,primary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_objs,primary_first_ele,primary_sub_dim,primary_global_offset
      ,secondary_first_ele,secondary_sub_dim
      );
  else if(num_targ_multi_vecs)
    targ_multi_vecs[0]->apply_op(
      apply_by,primary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_objs,primary_first_ele,primary_sub_dim,primary_global_offset
      ,secondary_first_ele,secondary_sub_dim
      );
}

void AbstractLinAlgPack::apply_op(
  EApplyBy                        apply_by
  ,const RTOpPack::RTOp           &primary_op
  ,const RTOpPack::RTOp           &secondary_op
  ,const size_t                   num_multi_vecs
  ,const MultiVector*             multi_vecs[]
  ,const size_t                   num_targ_multi_vecs
  ,MultiVectorMutable*            targ_multi_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
  ,const index_type               primary_first_ele
  ,const index_type               primary_sub_dim
  ,const index_type               primary_global_offset
  ,const index_type               secondary_first_ele
  ,const index_type               secondary_sub_dim
  )
{
  if(num_multi_vecs)
    multi_vecs[0]->apply_op(
      apply_by,primary_op,secondary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_obj,primary_first_ele,primary_sub_dim,primary_global_offset
      ,secondary_first_ele,secondary_sub_dim
      );
  else if(num_targ_multi_vecs)
    targ_multi_vecs[0]->apply_op(
      apply_by,primary_op,secondary_op
      ,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
      ,reduct_obj,primary_first_ele,primary_sub_dim,primary_global_offset
      ,secondary_first_ele,secondary_sub_dim
      );
}
