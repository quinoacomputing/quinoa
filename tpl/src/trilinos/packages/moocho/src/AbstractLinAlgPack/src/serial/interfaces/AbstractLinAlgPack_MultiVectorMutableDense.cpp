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

#include "AbstractLinAlgPack_MultiVectorMutableDense.hpp"
#include "AbstractLinAlgPack_VectorMutableDense.hpp"
#include "AbstractLinAlgPack_MatrixSymOpGetGMSSymMutable.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

MultiVectorMutableDense::MultiVectorMutableDense(
  const size_type                    rows
  ,const size_type                   cols
  )
{
  this->initialize(rows,cols);
}

MultiVectorMutableDense::MultiVectorMutableDense(
  DMatrixSlice                     gms
  ,BLAS_Cpp::Transp                  gms_trans
  ,const release_resource_ptr_t&     gms_release
  )
{
  this->initialize(gms,gms_trans,gms_release);
}

void MultiVectorMutableDense::initialize(
  const size_type                    rows
  ,const size_type                   cols
  )
{
  namespace rcp = MemMngPack;
  namespace rmp = MemMngPack;
  typedef Teuchos::RCP<DMatrix> vec_ptr_t;
  vec_ptr_t gms_ptr = Teuchos::rcp(new DMatrix(rows,cols));
  this->initialize(
    (*gms_ptr)()
    ,BLAS_Cpp::no_trans
    ,Teuchos::rcp(
      new rmp::ReleaseResource_ref_count_ptr<DMatrix>(
        gms_ptr
        )
      )
    );
}

void MultiVectorMutableDense::initialize(
  DMatrixSlice                     gms
  ,BLAS_Cpp::Transp                  gms_trans
  ,const release_resource_ptr_t&     gms_release
  )
{
  gms_.bind(gms);
  gms_trans_   = gms_trans;
  gms_release_ = gms_release;
}

// Overridden from MatrixOpGetGMS

const DMatrixSlice MultiVectorMutableDense::get_gms_view() const
{
  if(gms_trans_ == BLAS_Cpp::trans) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: We need to create a copy and transpose it!
  }
  return get_gms(); // No memory to allocate!
}

void MultiVectorMutableDense::free_gms_view(const DMatrixSlice* gms_view) const
{
  if(gms_trans_ == BLAS_Cpp::trans) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: We need to free the copy that we created in get_gms_view()
  }
  else {
    // Nothing to free!
  }
}

// Overridden from MatrixOpGetGMSMutable

DMatrixSlice MultiVectorMutableDense::get_gms_view()
{
  if(gms_trans_ == BLAS_Cpp::trans) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: We need to create a copy and transpose it!
  }
  return set_gms(); // No memory to allocate!
}

void MultiVectorMutableDense::commit_gms_view(DMatrixSlice* gms_view)
{
  if(gms_trans_ == BLAS_Cpp::trans) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: We need to free the copy that we created in get_gms_view()
  }
  else {
    // Nothing to free!
  }
}

// Overridden from MultiVector

MultiVectorMutableDense::access_by_t
MultiVectorMutableDense::access_by() const
{
  return ROW_ACCESS | COL_ACCESS | DIAG_ACCESS; // row, column and diagonal access is available!
}

// Overridden from MultiVectorMutable

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::col(index_type j)
{
  namespace rcp = MemMngPack;
  return Teuchos::rcp(
    new VectorMutableDense(
      DenseLinAlgPack::col( set_gms(), gms_trans(), j )
      ,Teuchos::null ) );
}

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::row(index_type i)
{
  namespace rcp = MemMngPack;
  return Teuchos::rcp(
    new VectorMutableDense(
      DenseLinAlgPack::row( set_gms(), gms_trans(), i )
      ,Teuchos::null ) );
}

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::diag(int k)
{
  namespace rcp = MemMngPack;
  return Teuchos::rcp(
    new VectorMutableDense(
      gms_.diag( gms_trans() == BLAS_Cpp::no_trans ? k : -k )
      ,Teuchos::null ) );
}

MultiVectorMutableDense::multi_vec_mut_ptr_t
MultiVectorMutableDense::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng)
{
  namespace rcp = MemMngPack;
  return Teuchos::rcp(
    new MultiVectorMutableDense(
      gms_(
        gms_trans() == BLAS_Cpp::no_trans   ? row_rng : col_rng
        ,gms_trans() == BLAS_Cpp::no_trans  ? col_rng : row_rng )
      ,gms_trans()
      ,Teuchos::null ) );
}

// Overridden from MatrixBase

size_type MultiVectorMutableDense::rows() const
{
  return BLAS_Cpp::rows( get_gms().rows(), get_gms().cols(), gms_trans() );
}

size_type MultiVectorMutableDense::cols() const
{
  return BLAS_Cpp::cols( get_gms().rows(), get_gms().cols(), gms_trans() );
}

// Overridden from MatrixOp

void MultiVectorMutableDense::zero_out()
{
  gms_ = 0.0;
}

void MultiVectorMutableDense::Mt_S( value_type alpha )
{
  DenseLinAlgPack::Mt_S(&gms_,alpha);
}

MatrixOp& MultiVectorMutableDense::operator=(const MatrixOp& mwo_rhs)
{
  DenseLinAlgPack::assign( &set_gms(), MatrixDenseEncap(mwo_rhs)(), gms_trans() );
  return *this;
}

std::ostream& MultiVectorMutableDense::output(std::ostream& out) const
{
  if(gms_trans() == BLAS_Cpp::no_trans)
    return out << gms_;
  return MatrixOpSerial::output(out);
}

bool MultiVectorMutableDense::Mp_StM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs
  ) const
{
  return MultiVectorMutable::Mp_StM(mwo_lhs,alpha,trans_rhs);
}

bool MultiVectorMutableDense::Mp_StM(
  value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs
  )
{
  return MultiVectorMutable::Mp_StM(alpha,M_rhs,trans_rhs);
}

bool MultiVectorMutableDense::syrk(
  BLAS_Cpp::Transp M_trans, value_type alpha
  ,value_type beta, MatrixSymOp* sym_lhs
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    sym_lhs == NULL, std::invalid_argument
    ,"MultiVectorMutableDense::syrk(...) : Error!" );
#endif
  MatrixSymOpGetGMSSymMutable
    *sym_get_lhs = dynamic_cast<MatrixSymOpGetGMSSymMutable*>(sym_lhs);
  if(!sym_get_lhs)
    return false;
  MatrixDenseSymMutableEncap  sym_gms_lhs(sym_get_lhs);
  DenseLinAlgPack::syrk( M_trans, alpha, get_gms(), beta, &sym_gms_lhs() );
  return true;
}

bool MultiVectorMutableDense::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ,value_type beta ) const
{
  if(MultiVector::Mp_StMtM(mwo_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2,beta))
    return true;
  return MatrixOpSerial::Mp_StMtM(mwo_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2,beta);
}

bool MultiVectorMutableDense::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ,value_type beta ) const
{
  if(MultiVector::Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta))
    return true;
  return MatrixOpSerial::Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta);
}

// Overridden from MatrixOpSerial

void MultiVectorMutableDense::Vp_StMtV(
  DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& vs_rhs2, value_type beta) const
{
  DenseLinAlgPack::Vp_StMtV(
    vs_lhs,alpha,gms_,BLAS_Cpp::trans_trans(gms_trans(),trans_rhs1),vs_rhs2,beta);
}

void MultiVectorMutableDense::Vp_StMtV(
  DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2, value_type beta) const
{
  AbstractLinAlgPack::Vp_StMtV(
    vs_lhs,alpha,gms_,BLAS_Cpp::trans_trans(gms_trans(),trans_rhs1),sv_rhs2,beta);
}

// protected

// Overridden from MultiVector

void MultiVectorMutableDense::apply_op(
  EApplyBy apply_by, const RTOpPack::RTOp& primary_op
  ,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
  ,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
  ,RTOpPack::ReductTarget* reduct_objs[]
  ,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
  ,const index_type secondary_first_ele , const index_type secondary_sub_dim 
  ) const
{
  MultiVector::apply_op(
    apply_by, primary_op, num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs
    ,reduct_objs
    ,primary_first_ele, primary_sub_dim, primary_global_offset, secondary_first_ele, secondary_sub_dim
    ); // ToDo: Provide Specialized implementation if needed?
}

void MultiVectorMutableDense::apply_op(
  EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
  ,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
  ,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
  ,const index_type secondary_first_ele , const index_type secondary_sub_dim 
  ) const
{
  MultiVector::apply_op(
    apply_by, primary_op, secondary_op, num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs
    ,reduct_obj
    ,primary_first_ele, primary_sub_dim, primary_global_offset, secondary_first_ele, secondary_sub_dim
    ); // ToDo: Provide Specialized implementation if needed?
}

} // end namespace AbstractLinAlgPack
