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

#include <limits>
#include <vector>

#include "AbstractLinAlgPack_MatrixSymPosDefCholFactor.hpp"
#include "AbstractLinAlgPack_BFGS_helpers.hpp"
#include "AbstractLinAlgPack_rank_2_chol_update.hpp"
#include "AbstractLinAlgPack_VectorMutableDense.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack_MatrixOpGetGMS.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_MatrixSymDiag.hpp"
#include "AbstractLinAlgPack_VectorSpaceFactory.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgLAPack.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_delete_row_col.hpp"
#include "DenseLinAlgPack_assert_print_nan_inf.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "ProfileHackPack_profile_hack.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"

#ifdef HAVE_MOOCHO_FORTRAN
#  define ALAP_DCHUD_DECL FORTRAN_FUNC_DECL_UL( void, DCHUD, dchud )
#  define ALAP_DCHUD_CALL FORTRAN_FUNC_CALL_UL( DCHUD, dchud )
#else
#  define ALAP_DCHUD_DECL void dchud_c
#  define ALAP_DCHUD_CALL dchud_c
#endif

// Helper functions
extern "C" {
  ALAP_DCHUD_DECL ( FortranTypes::f_dbl_prec* R
    , const FortranTypes::f_int& LDR, const FortranTypes::f_int& P
    , FortranTypes::f_dbl_prec* X, FortranTypes::f_dbl_prec* Z
    , const FortranTypes::f_int& LDZ, const FortranTypes::f_int& NZ
    , FortranTypes::f_dbl_prec* Y, FortranTypes::f_dbl_prec*  RHO
    , FortranTypes::f_dbl_prec* C, FortranTypes::f_dbl_prec* S );
} // end extern "C"

namespace {
//
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
//
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }
} // end namespace

namespace AbstractLinAlgPack {

// Constructors/initalizers

MatrixSymPosDefCholFactor::MatrixSymPosDefCholFactor()
  : maintain_original_(true), maintain_factor_(false)
  , factor_is_updated_(false), allocates_storage_(true)
  , max_size_(0), M_size_(0), M_l_r_(1), M_l_c_(1)
  , U_l_r_(1), U_l_c_(1), scale_(0.0), is_diagonal_(false)
  , pivot_tols_(0.0,0.0,0.0) // This is what DPOTRF(...) uses!
{}

MatrixSymPosDefCholFactor::MatrixSymPosDefCholFactor(
  DMatrixSlice                      *MU_store
  ,const release_resource_ptr_t&    release_resource_ptr
  ,size_type                        max_size
  ,bool                             maintain_original
  ,bool                             maintain_factor
  ,bool                             allow_factor
  ,bool                             set_full_view
  ,value_type                       scale
  )
  : pivot_tols_(0.0,0.0,0.0) // This is what DPOTRF(...) uses!
{
  init_setup(MU_store,release_resource_ptr,max_size,maintain_original
         ,maintain_factor,allow_factor,set_full_view,scale);
}

void MatrixSymPosDefCholFactor::init_setup(
  DMatrixSlice                      *MU_store
  ,const release_resource_ptr_t&    release_resource_ptr
  ,size_type                        max_size
  ,bool                             maintain_original
  ,bool                             maintain_factor
  ,bool                             allow_factor
  ,bool                             set_full_view
  ,value_type                       scale
  )
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  maintain_original || maintain_factor  ) );
  if( MU_store == NULL ) {
    maintain_original_ = maintain_original;
    maintain_factor_   = maintain_factor;
    factor_is_updated_ = maintain_factor;
    allocates_storage_ = true; // We will be able to allocate our own storage!
    release_resource_ptr_ = Teuchos::null; // Free any bound resource
    MU_store_.bind( DMatrixSlice(NULL,0,0,0,0) ); // Unbind this!
    max_size_ = max_size;
    M_size_ = 0;
    M_l_r_ = M_l_c_ = 1;
    if( !maintain_factor && !allow_factor )
      U_l_r_ = 0;  // Do not allow the factor to be computed
    else
      U_l_r_ = 1;  // Allow the factor to be computed
    scale_ = +1.0;
    is_diagonal_ = false;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT( !(  MU_store->rows()  ) );
    allocates_storage_ = false; // The client allocated the storage!
    MU_store_.bind(*MU_store);
    release_resource_ptr_ = release_resource_ptr;
    max_size_ = my_min( MU_store->rows(), MU_store->cols() ) - 1;
    if( set_full_view ) {
      TEUCHOS_TEST_FOR_EXCEPT( !(  scale != 0.0  ) );
      this->set_view(
        max_size_
        ,scale,maintain_original,1,1
        ,maintain_factor, allow_factor ? 1: 0, allow_factor ? 1 : 0
        );
    }
    else {
      this->set_view(
        0
        ,0.0,maintain_original,1,1
        ,maintain_factor, allow_factor ? 1: 0, allow_factor ? 1 : 0
        );
    }
  }
}

void MatrixSymPosDefCholFactor::set_view(
  size_t            M_size
  ,value_type       scale
  ,bool             maintain_original
  ,size_t           M_l_r
  ,size_t           M_l_c
  ,bool             maintain_factor
  ,size_t           U_l_r
  ,size_t           U_l_c
  )
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  maintain_original || maintain_factor  ) );
  if( max_size_ )
    allocate_storage(max_size_);
  else
    allocate_storage( my_max( M_l_r + M_size, M_l_c + M_size ) - 1 );
  // Check the preconditions
  if( maintain_original ) {
    TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= M_l_r && M_l_r <= M_l_c  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  M_l_r+M_size <= MU_store_.rows()  ) );
  }
  if( maintain_factor ) {
    TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= U_l_r && U_l_r >= U_l_c  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  U_l_c+M_size <= MU_store_.cols()  ) );
  }
  // Set the members
  maintain_original_    = maintain_original;
  maintain_factor_      = maintain_factor;
  is_diagonal_          = false;
  if( M_size ) {
    max_size_             = my_min( MU_store_.rows() - U_l_r, MU_store_.cols() - U_l_c );
    M_size_               = M_size;
    scale_                = scale;
    M_l_r_                = M_l_r;
    M_l_c_                = M_l_c;
    U_l_r_                = U_l_r;
    U_l_c_                = U_l_c;
    factor_is_updated_    = maintain_factor;
  }
  else {
    max_size_             = 0;
    M_size_               = 0;
    scale_                = 0.0;
    M_l_r_                = 0;
    M_l_c_                = 0;
    U_l_r_                = U_l_r;
    U_l_c_                = U_l_r;
    factor_is_updated_    = maintain_factor;
  }
}

void MatrixSymPosDefCholFactor::pivot_tols( PivotTolerances pivot_tols )
{
  pivot_tols_ = pivot_tols;
}

MatrixSymAddDelUpdateable::PivotTolerances MatrixSymPosDefCholFactor::pivot_tols() const
{
  return pivot_tols_;
}

// Overridden from MatrixBase

size_type MatrixSymPosDefCholFactor::rows() const
{
  return M_size_;
}

// Overridden from MatrixOp

void MatrixSymPosDefCholFactor::zero_out()
{
  this->init_identity(this->space_cols(),0.0);
}

std::ostream& MatrixSymPosDefCholFactor::output(std::ostream& out_arg) const
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  Teuchos::OSTab tab(out);
  if( M_size_ ) {
    if( maintain_original_ ) {
      *out
        << "Unfactored symmetric matrix stored as lower triangle (ignore upper nonzeros):\n"
        << M().gms();
    }
    if( factor_is_updated_ ) {
      *out
        << "Matrix scaling M = scale*U'*U, scale = " << scale_ << std::endl
        << "Upper cholesky factor U (ignore lower nonzeros):\n"
        << U().gms();
    }
  }
  else {
    *out << "0 0\n";
  }
  return out_arg;
}

bool MatrixSymPosDefCholFactor::Mp_StM(
  MatrixOp* m_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs
  ) const
{
  MatrixSymOpGetGMSSymMutable
    *symwo_gms_lhs = dynamic_cast<MatrixSymOpGetGMSSymMutable*>(m_lhs);
  if(!symwo_gms_lhs)
    return false;
  MatrixDenseSymMutableEncap sym_lhs(symwo_gms_lhs);
  const DMatrixSliceSym M       = this->M();
  DenseLinAlgPack::Mp_StM(
    &DMatrixSliceTriEle(sym_lhs().gms(),sym_lhs().uplo())
    ,alpha
    ,DMatrixSliceTriEle(M.gms(),M.uplo())
    );

  return true;
}

bool MatrixSymPosDefCholFactor::Mp_StM(
  value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !maintain_original_, std::logic_error
    ,"MatrixSymPosDefCholFactor::Mp_StM(alpha,M_rhs,trans_rhs): Error, Current implementation "
    "can not perform this operation unless the original matrix is being maintained." );
  // Perform the operation
  bool did_op  = false;
  bool diag_op = false;
  if(const MatrixSymOpGetGMSSym *symwo_gms_rhs = dynamic_cast<const MatrixSymOpGetGMSSym*>(&M_rhs)) {
    DMatrixSliceSym               M = this->M();
    MatrixDenseSymEncap   sym_rhs(*symwo_gms_rhs);
    DenseLinAlgPack::Mp_StM(
      &DMatrixSliceTriEle(M.gms(),M.uplo())
      ,alpha
      ,DMatrixSliceTriEle(sym_rhs().gms(),sym_rhs().uplo())
      );
    did_op  = true;
    diag_op = false;
  }
  else if(const MatrixSymDiag *symwo_diag_rhs = dynamic_cast<const MatrixSymDiag*>(&M_rhs)) {
    DMatrixSliceSym            M = this->M();
    VectorDenseEncap   sym_rhs_diag(symwo_diag_rhs->diag());
    LinAlgOpPack::Vp_StV( &M.gms().diag(), alpha, sym_rhs_diag() );
    did_op  = true;
    diag_op = true;
  }
  else if(const MatrixSymOp *symwo_rhs = dynamic_cast<const MatrixSymOp*>(&M_rhs)) {
    // ToDo: Implement this!
  }
  // Properly update the state of *this.
  // If only the original is updated
  if(did_op) {
    if( diag_op && is_diagonal_ )
      this->init_diagonal(VectorMutableDense(this->M().gms().diag(),Teuchos::null));
    else
      this->initialize(this->M());
    return true;
  }
  return false;
}

bool MatrixSymPosDefCholFactor::syrk(
  const MatrixOp      &mwo_rhs
  ,BLAS_Cpp::Transp   M_trans
  ,value_type         alpha
  ,value_type         beta
  )
{
  MatrixDenseSymMutableEncap  sym_gms_lhs(this);
  const MatrixOpGetGMS *mwo_rhs_gms = dynamic_cast<const MatrixOpGetGMS*>(&mwo_rhs);
  if(mwo_rhs_gms) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement
    return true;
  }
  else {
    // Here we will give up on symmetry and just compute the whole product S = op(mwo_rhs)*op(mwo_rhs')
    DenseLinAlgPack::DMatrixSliceTriEle tri_ele_gms_lhs = DenseLinAlgPack::tri_ele(sym_gms_lhs().gms(),sym_gms_lhs().uplo());
    if(beta==0.0)        DenseLinAlgPack::assign( &tri_ele_gms_lhs, 0.0 );
    else if(beta!=1.0)   DenseLinAlgPack::Mt_S( &tri_ele_gms_lhs, beta );
    const VectorSpace                 &spc = mwo_rhs.space_rows();
    const index_type                  m    = spc.dim();
    VectorSpace::multi_vec_mut_ptr_t  S    = spc.create_members(m);
    S->zero_out();
    LinAlgOpPack::M_MtM( S.get(), mwo_rhs, M_trans, mwo_rhs, BLAS_Cpp::trans_not(M_trans) );
    // Copy S into sym_ghs_lhs
    if( sym_gms_lhs().uplo() == BLAS_Cpp::lower ) {
      for( index_type j = 1; j <= m; ++j ) {
        LinAlgOpPack::Vp_StV( &sym_gms_lhs().gms().col(j)(j,m), alpha, VectorDenseEncap(*S->col(j)->sub_view(j,m))() );
      }
    }
    else {
      for( index_type j = 1; j <= m; ++j ) {
        LinAlgOpPack::Vp_StV( &sym_gms_lhs().gms().col(j)(1,j), alpha, VectorDenseEncap(*S->col(j)->sub_view(1,j))() );
      }
    }
    return true;
  }
  return false;
}

// Overridden from MatrixOpSerial

void MatrixSymPosDefCholFactor::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const DVectorSlice& x, value_type b
  ) const
{
  using Teuchos::implicit_ref_cast;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::Vp_StMtV(...DVectorSlice...)" );
#endif
  assert_initialized();

  DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(),
    implicit_ref_cast<const MatrixBase>(*this).rows(), cols(), no_trans, x.dim() );

  if( maintain_original_ ) {
    //
    // M = symmetric
    //
    // y = b*y + a*M*x
    //
    DenseLinAlgPack::Vp_StMtV( y, a, M(), no_trans, x, b );
  }
  else {
    //
    // M = scale*U'*U
    //
    // y = b*y + a*op(M)*x
    //   = b*y = scale*a*U'*U*x
    //
    const DMatrixSliceTri
      U = this->U();

    if( b == 0.0 ) {
      // No temporary needed
      //
      // y = U*x
      DenseLinAlgPack::V_MtV( y, U, no_trans, x ); 
      // y = U'*y
      DenseLinAlgPack::V_MtV( y, U, trans, *y ); 
      // y *= scale*a
      if( a != 1.0 || scale_ != 1.0 )
        DenseLinAlgPack::Vt_S( y, scale_*a );
    }
    else {
      // We need a temporary
      DVector t;
      // t = U*x
      DenseLinAlgPack::V_MtV( &t, U, no_trans, x );
      // t = U'*t
      DenseLinAlgPack::V_MtV( &t(), U, trans, t() );
      // y *= b
      if(b != 1.0)
        DenseLinAlgPack::Vt_S( y, b );
      // y += scale*a*t
      DenseLinAlgPack::Vp_StV( y, scale_*a, t() );
    }
  }
}

void MatrixSymPosDefCholFactor::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const SpVectorSlice& x, value_type b
  ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::Vp_StMtV(...SpVectorSlice...)" );
#endif
  assert_initialized();
  if( maintain_original_ ) {
    const DMatrixSlice M = this->M().gms(); // This is lower triangular!
    const size_type n = M.rows();
    DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(), n, n, no_trans, x.dim() );
    DenseLinAlgPack::Vt_S(y,b); // y = b*y
    //
    // Compute product column by column corresponding to x_itr->index() + x.offset()
    //
    // y += a * M * e(i) * x(i)
    //
    // [ y(1:i-1) ] += a * x(i) * [ ...  M(1:i-1,i) ... ]  stored as M(i,1:i-1) 
    // [ y(i:n)   ]               [ ...  M(i:n,i)   ... ]  stored as M(i:n,i)
    //
    for( SpVectorSlice::const_iterator x_itr = x.begin(); x_itr != x.end(); ++x_itr ) {
      const size_type i = x_itr->index() + x.offset();
      if( i > 1 )
        DenseLinAlgPack::Vp_StV( &(*y)(1,i-1), a * x_itr->value(), M.row(i)(1,i-1) );
      DenseLinAlgPack::Vp_StV( &(*y)(i,n), a * x_itr->value(), M.col(i)(i,n) );
    }
  }
  else {
    MatrixOpSerial::Vp_StMtV(y,a,M_trans,x,b); // ToDo: Specialize when needed!
  }
}

void MatrixSymPosDefCholFactor::Vp_StPtMtV(
  DVectorSlice* y, value_type a, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp H_trans, const DVectorSlice& x, value_type b
  ) const
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::Vp_StPtMtV(...DVectorSlice...)" );
#endif
  assert_initialized();
  MatrixOpSerial::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // ToDo: Specialize when needed!
}

void MatrixSymPosDefCholFactor::Vp_StPtMtV(
  DVectorSlice* y, value_type a, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp H_trans, const SpVectorSlice& x, value_type b
  ) const
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::Vp_StPtMtV(...SpVectorSlice...)" );
#endif
  assert_initialized();
  if( maintain_original_ ) {
    DenseLinAlgPack::Vt_S(y,b); // y = b*y
    const DMatrixSlice M = this->M().gms(); // This is lower triangular!
    // Compute product column by corresponding to x_itr->index() + x.offset()
    /*
    if( P.is_identity() ) {
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement
    }
    else {
      for( SpVectorSlice::const_iterator x_itr = x.begin(); x_itr != x.end(); ++x_itr ) {
        const size_type i = x_itr->index() + x.offset();
        

      }
    }
    */
    MatrixOpSerial::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // ToDo: Specialize when needed!
  }
  else {
    MatrixOpSerial::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // ToDo: Specialize when needed!
  }
}

// Overridden from MatrixSymOpSerial

void MatrixSymPosDefCholFactor::Mp_StPtMtP(
  DMatrixSliceSym* S, value_type a
  , EMatRhsPlaceHolder dummy_place_holder
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , value_type b ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::Mp_StPtMtP(...)" );
#endif
  assert_initialized();
  if( !maintain_original_ ) {
    MatrixSymOpSerial::Mp_StPtMtP(S,a,dummy_place_holder,P,P_trans,b);
  }
  else {
    MatrixSymOpSerial::Mp_StPtMtP(S,a,dummy_place_holder,P,P_trans,b); // ToDo: Override when needed!
  }
}

// Overridden from MatrixNonsingSerial

void MatrixSymPosDefCholFactor::V_InvMtV(
  DVectorSlice* y, BLAS_Cpp::Transp M_trans, const DVectorSlice& x
  ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::V_InvMtV(...DVectorSlice...)" );
#endif
  assert_initialized();

  //
  // M = scale*U'*U
  //
  // y = inv(op(M))*x
  // =>
  // op(M)*y = x
  // =>
  // scale*U'*U*y = x
  // =>
  // y = (1/scale)*inv(U)*inv(U')*x
  //
  update_factorization();
  const DMatrixSliceTri
    U = this->U();
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(), U.rows(), U.cols(), no_trans, x.dim() );
  // y = inv(U')*x
  DenseLinAlgPack::V_InvMtV( y, U, trans, x );
  // y = inv(U)*y
  DenseLinAlgPack::V_InvMtV( y, U, no_trans, *y );
  // y *= (1/scale)
  if( scale_ != 1.0 )
    DenseLinAlgPack::Vt_S( y, 1.0/scale_ );
}

void MatrixSymPosDefCholFactor::V_InvMtV(
  DVectorSlice* y, BLAS_Cpp::Transp M_trans, const SpVectorSlice& x
  ) const
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::V_InvMtV(...SpVectorSlice...)" );
#endif
  assert_initialized();
  MatrixNonsingSerial::V_InvMtV(y,M_trans,x);
}

// Overridden from MatrixSymNonsingSerial

void MatrixSymPosDefCholFactor::M_StMtInvMtM(
    DMatrixSliceSym* S, value_type a, const MatrixOpSerial& B
  , BLAS_Cpp::Transp B_trans, EMatrixDummyArg dummy_arg
  ) const
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::M_StMtInvMtM(...)" );
#endif

//	// Uncomment to use the defalut implementation (for debugging)
//	MatrixSymFactorized::M_StMtInvMtM(S,a,B,B_trans,dummy_arg); return;

  using BLAS_Cpp::trans;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans_not;
  using BLAS_Cpp::upper;
  using BLAS_Cpp::nonunit;
  using DenseLinAlgPack::tri;
  using DenseLinAlgPack::syrk;
  using DenseLinAlgPack::M_StInvMtM;
  using LinAlgOpPack::assign;

  assert_initialized();
  update_factorization();

  //
  // S = a * op(B) * inv(M) * op(B)'
  // 
  // M = scale*U'*U
  // =>
  // inv(M) = scale*inv(U'*U) = scale*inv(U)*inv(U')
  // =>
  // S = scale*a * op(B) * inv(U) * inv(U') * op(B)'
  // 
  // T = op(B)'
  // 
  // T = inv(U') * T (inplace with BLAS)
  // 
  // S = scale*a * T' * T
  // 
  DenseLinAlgPack::MtM_assert_sizes( 
    rows(), cols(), no_trans
    ,B.rows(), B.cols(), trans_not(B_trans)
    );
  DenseLinAlgPack::Mp_MtM_assert_sizes(
    S->rows(), S->cols(), no_trans
    ,B.rows(), B.cols(), B_trans
    ,B.rows(), B.cols(), trans_not(B_trans)
    );
  // T = op(B)'
  DMatrix T;
  assign( &T, B, trans_not(B_trans) );
  // T = inv(U') * T (inplace with BLAS)
  M_StInvMtM( &T(), 1.0, this->U(), trans, T(), no_trans );
  // S = scale*a * T' * T
  syrk( trans, scale_*a, T(), 0.0, S );
}

// Overridden from MatrixSymDenseInitialize

void MatrixSymPosDefCholFactor::initialize( const DMatrixSliceSym& M )
{
  // Initialize without knowing the inertia but is must be p.d.
  this->initialize(
    M, M.rows(), maintain_factor_
    , MatrixSymAddDelUpdateable::Inertia()
    , MatrixSymAddDelUpdateable::PivotTolerances()
    );
}

// Overridden from MatrixSymOpGetGMSSym

const DenseLinAlgPack::DMatrixSliceSym MatrixSymPosDefCholFactor::get_sym_gms_view() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !maintain_original_, std::logic_error
    ,"MatrixSymPosDefCholFactor::get_sym_gms_view(): Error, maintain_original must be "
    "true in order to call this method!" );
  return this->M();
}

void MatrixSymPosDefCholFactor::free_sym_gms_view(const DenseLinAlgPack::DMatrixSliceSym* sym_gms_view) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !maintain_original_, std::logic_error
    ,"MatrixSymPosDefCholFactor::free_sym_gms_view(...): Error, maintain_original must be "
    "true in order to call this method!" );
  // Nothing todo
}

// Overridden from MatrixSymOpGetGMSSymMutable

DenseLinAlgPack::DMatrixSliceSym MatrixSymPosDefCholFactor::get_sym_gms_view()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !maintain_original_, std::logic_error
    ,"MatrixSymPosDefCholFactor::get_sym_gms_view(): Error, maintain_original must be "
    "true in order to call this method!" );
  return this->M();
}

void MatrixSymPosDefCholFactor::commit_sym_gms_view(DenseLinAlgPack::DMatrixSliceSym* sym_gms_view)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !maintain_original_, std::logic_error
    ,"MatrixSymPosDefCholFactor::commit_sym_gms_view(...): Error, maintain_original must be "
    "true in order to call this method!" );
  this->initialize(*sym_gms_view);
}

// Overridden from MatrixExtractInvCholFactor

void MatrixSymPosDefCholFactor::extract_inv_chol( DMatrixSliceTriEle* InvChol ) const
{
  assert_initialized();
  update_factorization();
  //
  // The matrix is represented as the upper cholesky factor:
  //   M = scale * U' * U
  //
  //   inv(M) = inv(scale*U*U') = (1/sqrt(scale))*inv(U)*(1/sqrt(scale))*inv(U')
  //          = UInv * UInv'
  // =>
  //   UInv = (1/sqrt(scale))*inv(U)
  //
   // Here scale > 0 or an exception will be thrown.
  //
  // Compute the inverse cholesky factor as:
  //
  // Upper cholesky:
  //    sqrt(scale) * U * UInv = I => InvChol = UInv = (1/sqrt(scale))*inv(U) * I
  //    
  // Lower cholesky:
  //    sqrt(scale) * L * LInv = I => InvChol = LInv = (1/sqrt(scale))*inv(U) * inv(U') * I
  //
  TEUCHOS_TEST_FOR_EXCEPTION(
    scale_ < 0.0, std::logic_error
    ,"MatrixSymPosDefCholFactor::extract_inv_chol(...) : "
    "Error, we can not compute the inverse cholesky factor "
    "af a negative definite matrix." );
  DenseLinAlgPack::assign( &InvChol->gms(), 0.0 );  // Set InvChol to identity first.
  InvChol->gms().diag() = 1.0;
  DenseLinAlgPack::M_StInvMtM(                      // Comput InvChol using Level-3 BLAS
      &InvChol->gms(), 1.0 / std::sqrt(scale_), U()
    , InvChol->uplo() == BLAS_Cpp::upper ? BLAS_Cpp::no_trans : BLAS_Cpp::trans
    , InvChol->gms(), BLAS_Cpp::no_trans );
}

// Overridden from MatrixSymSecantUpdateble

void MatrixSymPosDefCholFactor::init_identity( const VectorSpace& space_diag, value_type alpha )
{
  const size_type n = space_diag.dim();
  allocate_storage( max_size_ ? max_size_ : n );
  //
  // B = alpha*I = = alpha*I*I = scale*U'*U
  // =>
  // U = sqrt(|alpha|) * I
  //
  // Here we will set scale = sign(alpha)
  //
  const value_type scale = alpha > 0.0 ? +1.0: -1.0; // Explicitly set the scale
  resize_and_zero_off_diagonal(n,scale);
  if( maintain_original_ ) {
    M().gms().diag()(1,n) = alpha;
  }
  if( maintain_factor_ ) {
    U().gms().diag()(1,n) = std::sqrt(std::fabs(alpha));
    factor_is_updated_ = true;
  }
  is_diagonal_ = true;
}

void MatrixSymPosDefCholFactor::init_diagonal( const Vector& diag_in )
{
  VectorDenseEncap diag_encap(diag_in);
  const DVectorSlice diag = diag_encap(); // When diag_encap is destroyed, bye-bye view!

  allocate_storage( max_size_ ? max_size_ : diag.dim() );
  //
  // M = scale * U' * U = scale * (1/scale)*diag^(1/2) * (1/scale)*diag^(1/2)
  //
  // Here we will set scale = sign(diag(1)) and validate the rest
  //
  if( diag.dim() == 0 ) {
    M_size_ = 0;
    return; // We are unsizing this thing
  }
  const value_type scale = diag(1) > 0.0 ? +1.0: -1.0;
  resize_and_zero_off_diagonal(diag.dim(),scale);
  if( maintain_original_ ) {
    M().gms().diag() = diag;
    // ToDo: validate that scale*diag > 0
  }
  if( maintain_factor_ ) {
    DVectorSlice U_diag = U().gms().diag();
    U_diag = diag;
    if( scale_ != 1.0 )
      DenseLinAlgPack::Vt_S( &U_diag, 1.0/scale_ );
    DenseLinAlgPack::sqrt( &U_diag, U_diag );
    DenseLinAlgPack::assert_print_nan_inf( U_diag, "(1/scale)*diag", true, NULL );
    factor_is_updated_ = true;
  }
  is_diagonal_ = true;
}

void MatrixSymPosDefCholFactor::secant_update(
  VectorMutable* s_in, VectorMutable* y_in, VectorMutable* Bs_in
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using DenseLinAlgPack::dot;
  using DenseLinAlgPack::norm_2;
  using DenseLinAlgPack::norm_inf;
  namespace rcp = MemMngPack;

  assert_initialized();

  // Validate the input
  TEUCHOS_TEST_FOR_EXCEPT( !(  s_in && y_in  ) );
  DenseLinAlgPack::Vp_MtV_assert_sizes( y_in->dim(), M_size_, M_size_, no_trans, s_in->dim() );

  // Get the serial vectors
  VectorDenseMutableEncap s_encap(*s_in);
  VectorDenseMutableEncap y_encap(*y_in);
  VectorDenseMutableEncap Bs_encap( Bs_in ? *Bs_in : *y_in); // Must pass something on
  DVectorSlice
    *s  = &s_encap(),   // When s_encap, y_encap and Bs_encap are destroyed
    *y  = &y_encap(),   // these views go bye-bye!
    *Bs = ( Bs_in ? &Bs_encap() : NULL );

  // Check skipping the BFGS update
  const value_type
    sTy	      = dot(*s,*y),
    sTy_scale = sTy*scale_;
  std::ostringstream omsg;
  if( !BFGS_sTy_suff_p_d(
      VectorMutableDense((*s)(),Teuchos::null)
      ,VectorMutableDense((*y)(),Teuchos::null)
      ,&sTy_scale,&omsg,"\nMatrixSymPosDefCholFactor::secant_update(...)"
      )
    )
  {
    throw UpdateSkippedException( omsg.str() );	
  }
  // Compute Bs if it was not passed in
  DVector Bs_store;
  DVectorSlice Bs_view;
  if( !Bs ) {
    LinAlgOpPack::V_MtV( &Bs_store, *this, no_trans, *s );
    Bs_view.bind( Bs_store() );
    Bs = &Bs_view;
  }
  // Check that s'*Bs is positive and if not then throw exception
  const value_type sTBs = dot(*s,*Bs);
  TEUCHOS_TEST_FOR_EXCEPTION(
    scale_*sTBs <= 0.0 && scale_ > 0.0, UpdateFailedException
    ,"MatrixSymPosDefCholFactor::secant_update(...) : "
    "Error, B can't be positive definite if s'*Bs <= 0.0" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    scale_*sTBs <= 0.0 && scale_ <= 0.0, UpdateFailedException
    ,"MatrixSymPosDefCholFactor::secant_update(...) : "
    "Error, B can't be negative definite if s'*Bs >= 0.0" );
  if( maintain_original_ ) {
    //
    // Compute the BFGS update of the original, nonfactored matrix
    //
    // Just preform two symmetric updates.
    //
    // B = B + (-1/s'*Bs) * Bs*Bs' + (1/s'*y) * y*y'
    //
    DMatrixSliceSym M = this->M();
    DenseLinAlgPack::syr( -1.0/sTBs, *Bs, &M );
    DenseLinAlgPack::syr( 1.0/sTy, *y, &M );
  }
  if( maintain_factor_ ) {
    //
    // Compute the BFGS update for the cholesky factor
    //
    // If this implementation is based on the one in Section 9.2, page 198-201 of:
    //
    // Dennis, J.E., R.B. Schnabel
    // "Numerical Methods for Unconstrained Optimization"
    //
    // Given that we have B = scale*U'*U and the BFGS update:
    //
    // B_new = B + y*y'/(y'*s) - (B*s)*(B*s)'/(s'*B*s)
    //
    // We can rewrite it in the following form:
    //
    // B_new = scale*(U' + a*u*v')*(U + a*v*u') = scale*U_new'*U_new
    //
    // where:
    //     v = sqrt(y'*s/(s'*B*s))*U*s
    //     u = y - U'*v
    //     a = 1/(v'*v)
    //
    DMatrixSliceTri U = this->U();
    // v = sqrt(y'*s/(s'*B*s))*U*s
    DVectorSlice v = *s; // Reuse s as storage for v
    DenseLinAlgPack::V_MtV( &v, U, no_trans, v ); // Direct call to xSYMV(...)
    DenseLinAlgPack::Vt_S( &v, std::sqrt( sTy / sTBs ) );
    // u = (y - U'*v)
    DVectorSlice u = *y;  // Reuse y as storage for u
    DenseLinAlgPack::Vp_StMtV( &u, -1.0, U, trans, v );
    // a = 1/(v'*v)
    const value_type  a = 1.0/dot(v,v);
    // Perform Givens rotations to make Q*(U' + a*u*v') -> U_new upper triangular:
    //
    // B_new = scale*(U' + a*u*v')*Q'*Q*(U + a*v*u') = scale*U_new'*U_new
    rank_2_chol_update(
      a, &v, u, v.dim() > 1 ? &U.gms().diag(-1) : NULL
      , &DMatrixSliceTriEle(U.gms(),BLAS_Cpp::upper), no_trans );
  }
  else {
    factor_is_updated_ = false;
  }
  is_diagonal_ = false;
}

// Overridden from MatrixSymAddDelUpdateble

void MatrixSymPosDefCholFactor::initialize(
  value_type         alpha
  ,size_type         max_size
  )
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::initialize(alpha,max_size)" );
#endif
  allocate_storage(max_size);

  if( alpha == 0.0 ) 
    throw SingularUpdateException(
      "MatrixSymPosDefCholFactor::initialize(...): "
      "Error, alpha == 0.0, matrix is singular.", 0.0 );

  // Resize the view
  if( maintain_original_ ) {
    M_l_r_ = 1;
    M_l_c_ = 1;
  }
  if( U_l_r_ ) {
    U_l_r_ = 1;
    U_l_c_ = 1;
  }
  M_size_ = 1;
  max_size_ = my_min( MU_store_.rows(), MU_store_.cols() ) - 1;
  scale_ = alpha > 0.0 ? +1.0 : -1.0;
  // Update the matrix
  if( maintain_original_ ) {
    M().gms()(1,1) = alpha;
  }
  if( U_l_r_ ) {
    U().gms()(1,1) = std::sqrt( scale_ * alpha );
    factor_is_updated_ = true;
  }
  is_diagonal_ = false;
}

void MatrixSymPosDefCholFactor::initialize(
  const DMatrixSliceSym      &A
  ,size_type         max_size
  ,bool              force_factorization
  ,Inertia           inertia
  ,PivotTolerances   pivot_tols
  )
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::initialize(A,max_size...)" );
#endif
  typedef MatrixSymAddDelUpdateable::Inertia   Inertia;

  allocate_storage(max_size);

  const size_type
    n = A.rows();

  // Validate proper usage of inertia parameter
  TEUCHOS_TEST_FOR_EXCEPT( !(  inertia.zero_eigens == Inertia::UNKNOWN || inertia.zero_eigens == 0  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  (inertia.neg_eigens == Inertia::UNKNOWN && inertia.pos_eigens == Inertia::UNKNOWN ) )
      || ( inertia.neg_eigens == n && inertia.pos_eigens == 0 )
      || ( inertia.neg_eigens == 0 && inertia.pos_eigens == n )
    );

  // We can infer if the matrix is p.d. or n.d. by the sign of the diagonal
  // elements.  If a matrix is s.p.d. (s.n.d) then A(i,i) > 0 (A(i,i) < 0)
  // for all i = 1...n so we will just check the first i.
  const value_type
    a_11 = A.gms()(1,1);
  const int sign_a_11 = ( a_11 == 0.0 ? 0 : ( a_11 > 0.0 ? +1 : -1 ) );
  if( sign_a_11 == 0.0 )
    std::logic_error(
      "MatrixSymPosDefCholFactor::initialize(...) : "
      "Error, A can not be positive definite or negative definete "
      "if A(1,1) == 0.0" );
  if( inertia.pos_eigens == n && sign_a_11 < 0 )
    std::logic_error(
      "MatrixSymPosDefCholFactor::initialize(...) : "
      "Error, A can not be positive definite "
      "if A(1,1) < 0.0" );
  if( inertia.neg_eigens == n && sign_a_11 > 0 )
    std::logic_error(
      "MatrixSymPosDefCholFactor::initialize(...) : "
      "Error, A can not be negative definite "
      "if A(1,1) > 0.0" );
  // Now we have got the sign
  const value_type
    scale = (value_type)sign_a_11;
  // Setup the view
  set_view(
    n,scale,maintain_original_,1,1
    ,maintain_factor_, U_l_r_ ? 1 : 0, U_l_r_ ? 1 : 0
    );
  // Now set the matrix and update the factors
  if( maintain_original_ ) {
    // Set M = S
    DenseLinAlgPack::assign( 
      &DMatrixSliceTriEle( M().gms(), BLAS_Cpp::lower )
      ,DMatrixSliceTriEle( A.gms(), A.uplo() )
      );
  }
  if( maintain_factor_ || force_factorization ) {
    // Copy S into U for an inplace factorization.
    DMatrixSliceTriEle U_ele = DMatrixSliceTriEle( U().gms(), BLAS_Cpp::upper );
    DenseLinAlgPack::assign( &U_ele, DMatrixSliceTriEle( A.gms(), A.uplo() ) );
    if( sign_a_11 < 0 )
      DenseLinAlgPack::Mt_S( &U_ele, -1.0 );
    try {
      DenseLinAlgLAPack::potrf( &U_ele );
      factor_is_updated_ = true;
    }
    catch( const DenseLinAlgLAPack::FactorizationException& excpt ) {
      M_size_ = 0;  // set unsized
      throw SingularUpdateException( excpt.what(), 0.0 );
    }
    catch(...) {
      M_size_ = 0;
      throw;
    }
    // Validate that the tolerances are met and throw appropriate
    // exceptions.  We already know that the matrix is technically
    // p.d. or n.d.  Now we must determine gamma = (min(|diag|)/max(|diag|))^2
    value_type
      min_diag = std::numeric_limits<value_type>::max(),
      max_diag = 0.0;
    DVectorSlice::iterator
      U_itr = U_ele.gms().diag().begin(),
      U_end = U_ele.gms().diag().end();
    while( U_itr != U_end ) {
      const value_type U_abs = std::abs(*U_itr++);
      if(U_abs < min_diag) min_diag = U_abs;
      if(U_abs > max_diag) max_diag = U_abs;
    }
    const value_type gamma = (min_diag*min_diag)/(max_diag*max_diag);
    // Validate gamma
    PivotTolerances use_pivot_tols = pivot_tols_;
    if( pivot_tols.warning_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.warning_tol = pivot_tols.warning_tol;
    if( pivot_tols.singular_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.singular_tol = pivot_tols.singular_tol;
    if( pivot_tols.wrong_inertia_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.wrong_inertia_tol = pivot_tols.wrong_inertia_tol;
    const bool throw_exception = (
      gamma == 0.0
      || gamma <= use_pivot_tols.warning_tol
      || gamma <= use_pivot_tols.singular_tol
      );
    // Create message and throw exceptions
    std::ostringstream onum_msg, omsg;
    if(throw_exception) {
      onum_msg
        << "gamma = (min(|diag(U)(k)|)/|max(|diag(U)(k)|))^2 = (" << min_diag <<"/"
        << max_diag << ")^2 = " << gamma;
      omsg
        << "MatrixSymPosDefCholFactor::initialize(...): ";
      if( gamma <= use_pivot_tols.singular_tol ) {
        M_size_ = 0;  // The initialization failed!
        omsg
          << "Singular update!\n" << onum_msg.str() << " <= singular_tol = "
          << use_pivot_tols.singular_tol;
        throw SingularUpdateException( omsg.str(), gamma );
      }
      else if( gamma <= use_pivot_tols.warning_tol ) {
        omsg
          << "Warning, near singular update!\nsingular_tol = " << use_pivot_tols.singular_tol
          << " < " << onum_msg.str() << " <= warning_tol = "
          << use_pivot_tols.warning_tol;
        throw WarnNearSingularUpdateException( omsg.str(), gamma ); // The initialization still succeeded through
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPT(true); // Only local programming error?
      }
    }
  }
  else {
    factor_is_updated_ = false; // The factor is not updated!
  }
}

size_type MatrixSymPosDefCholFactor::max_size() const
{
  return max_size_;
}

MatrixSymAddDelUpdateable::Inertia
MatrixSymPosDefCholFactor::inertia() const
{
  typedef MatrixSymAddDelUpdateable MSADU;
  typedef MSADU::Inertia Inertia;
  return ( M_size_
       ? ( scale_ > 0.0
         ? Inertia(0,0,M_size_)
         : Inertia(M_size_,0,0) )
       : Inertia(0,0,0) );
}

void MatrixSymPosDefCholFactor::set_uninitialized()
{
  M_size_ = 0;
}

void MatrixSymPosDefCholFactor::augment_update(
  const DVectorSlice  *t
  ,value_type        alpha
  ,bool              force_refactorization
  ,EEigenValType     add_eigen_val
  ,PivotTolerances   pivot_tols
  )
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::augment_udpate(...)" );
#endif
  using Teuchos::implicit_ref_cast;
  using DenseLinAlgPack::dot;
  using DenseLinAlgPack::norm_inf;
  typedef MatrixSymAddDelUpdateable  MSADU;

  assert_initialized();

  // Validate the input
  TEUCHOS_TEST_FOR_EXCEPTION(
    implicit_ref_cast<const MatrixBase>(*this).rows() >= max_size(),
    MaxSizeExceededException,
    "MatrixSymPosDefCholFactor::augment_update(...) : "
    "Error, the maximum size would be exceeded." );
  TEUCHOS_TEST_FOR_EXCEPTION(
    t && t->dim() != M_size_, std::length_error
    ,"MatrixSymPosDefCholFactor::augment_update(...): "
    "Error, t.dim() must be equal to this->rows()." );
  if( !(add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN || add_eigen_val != MSADU::EIGEN_VAL_ZERO ) )
    throw SingularUpdateException(
      "MatrixSymPosDefCholFactor::augment_update(...): "
      "Error, the client has specified a singular update in add_eigen_val.", -1.0 );
  if( alpha == 0.0 )
    throw SingularUpdateException(
      "MatrixSymPosDefCholFactor::augment_update(...): "
      "Error, alpha == 0.0 so the matrix is not positive definite "
      "or negative definite.", -1.0 );
  if( scale_ > 0.0 && alpha < 0.0 )
    throw WrongInertiaUpdateException(
      "MatrixSymPosDefCholFactor::augment_update(...): "
      "Error, alpha < 0.0 so the matrix is not postivie definite ", -1.0 );
  if( scale_ < 0.0 && alpha > 0.0 )
    throw WrongInertiaUpdateException(
      "MatrixSymPosDefCholFactor::augment_update(...): "
      "Error, alpha > 0.0 so the matrix is not negative definite ", -1.0 );
  if( scale_ > 0.0 && !( add_eigen_val == MSADU::EIGEN_VAL_POS || add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN ) )
    throw WrongInertiaUpdateException(
      "MatrixSymPosDefCholFactor::augment_update(...): "
      "Error, alpha > 0.0 but the user specified a non-positive definite new matrix.", -1.0 );
  if( scale_ < 0.0 && !( add_eigen_val == MSADU::EIGEN_VAL_NEG || add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN ) )
    throw WrongInertiaUpdateException(
      "MatrixSymPosDefCholFactor::augment_update(...): "
      "Error, alpha > 0.0 but the user specified a non-positive definite new matrix.", -1.0 );
  // First try to augment the factor to verify that the matrix is still p.d. or n.d.
  bool throw_exception = false; // If true then throw exception
  std::ostringstream omsg;      // Will be set if an exception has to be thrown.
  value_type gamma;             // ...
  if( maintain_factor_ ) {
    //
    // The update is:
    //
    // B_new = [  B     t   ]
    //         [  t'  alpha ]
    //
    //       = scale * [     U'*U         (1/scale)*t   ]
    //                 [ (1/scale)*t'   (1/scale)*alpha ]
    //
    // We seek to find a new cholesky factor of the form:
    //
    // U_new = [ U11  u12  ]
    //         [  0   u22  ]
    //
    // B_new = scale*U_new'*U_new
    //
    //       = scale * [ U11'   0  ] * [ U11   u12 ]
    //                 [ u12'  u22 ]   [  0    u22 ]
    //
    //       = scale * [  U11'*U11          U11'*u12  ]
    //                 [  u12'*U11   u12'*u12 + u22^2 ]
    //
    // From the above we see that:
    // =>
    //   U11 = U
    //   u12 = inv(U') * (1/scale) * t
    //   u22 = sqrt( (1/scale)*alpha - u12'*u12 );
    //
    // We must compute gamma relative to the LU factorization
    // so we must square ||U11.diag()||inf.

    // Get references to the storage for the to-be-updated parts for the new factor.
    DVectorSlice u12 = MU_store_.col(U_l_c_+M_size_+1)(U_l_r_,U_l_r_+M_size_-1);
    value_type &u22 = MU_store_(U_l_r_+M_size_,U_l_c_+M_size_+1);
    // u12 = inv(U') * (1/scale) * t
    if(t) {
      DenseLinAlgPack::V_InvMtV( &u12, U(), BLAS_Cpp::trans, *t );
      if( scale_ != 1.0 ) DenseLinAlgPack::Vt_S( &u12, 1.0/scale_ );
    }
    else {
      u12 = 0.0;
    }
    // u22^2 = (1/scale)*alpha - u12'*u12;
    const value_type
      u22sqr          = (1/scale_) * alpha - ( t ? dot( u12, u12 ) : 0.0 ),
      u22sqrabs       = std::abs(u22sqr),
      nrm_U_diag      = norm_inf(U().gms().diag()),
      sqr_nrm_U_diag  = nrm_U_diag * nrm_U_diag;
    // Calculate gamma in proper context
    gamma = u22sqrabs / sqr_nrm_U_diag;
    // Check gamma
    const bool
      correct_inertia = u22sqr > 0.0;
    PivotTolerances use_pivot_tols = pivot_tols_;
    if( pivot_tols.warning_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.warning_tol = pivot_tols.warning_tol;
    if( pivot_tols.singular_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.singular_tol = pivot_tols.singular_tol;
    if( pivot_tols.wrong_inertia_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.wrong_inertia_tol = pivot_tols.wrong_inertia_tol;
    throw_exception = (
      gamma == 0.0
      || correct_inertia && gamma <= use_pivot_tols.warning_tol
      || correct_inertia && gamma <= use_pivot_tols.singular_tol
      || !correct_inertia
      );
    // Create message and throw exceptions
    std::ostringstream onum_msg;
    if(throw_exception) {
      onum_msg
        << "gamma = |u22^2|/(||diag(U11)||inf)^2 = |" << u22sqr <<"|/("
        << nrm_U_diag << ")^2 = " << gamma;
      omsg
        << "MatrixSymPosDefCholFactor::augment_update(...): ";
      if( correct_inertia && gamma <= use_pivot_tols.singular_tol ) {
        omsg
          << "Singular update!\n" << onum_msg.str() << " <= singular_tol = "
          << use_pivot_tols.singular_tol;
        throw SingularUpdateException( omsg.str(), gamma );
      }
      else if( !correct_inertia && gamma <= use_pivot_tols.singular_tol ) {
        omsg
          << "Singular update!\n" << onum_msg.str() << " <= wrong_inertia_tol = "
          << use_pivot_tols.wrong_inertia_tol;
        throw SingularUpdateException( omsg.str(), gamma );
      }
      
      else if( !correct_inertia ) {
        omsg
          << "Indefinite update!\n" << onum_msg.str() << " >= wrong_inertia_tol = "
          << use_pivot_tols.wrong_inertia_tol;
        throw WrongInertiaUpdateException( omsg.str(), gamma );
      }
      else if( correct_inertia && gamma <= use_pivot_tols.warning_tol ) {
        omsg
          << "Warning, near singular update!\nsingular_tol = " << use_pivot_tols.singular_tol
          << " < " << onum_msg.str() << " <= warning_tol = "
          << use_pivot_tols.warning_tol;
        // Don't throw the exception till the very end!
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPT(true); // Only local programming error?
      }
    }
    // u22 = sqrt(u22^2)
    u22 = std::sqrt(u22sqrabs);
  }
  else {
    factor_is_updated_ = false;
  }
  // Now augment the original
  if( maintain_original_ ) {
    //
    // M_new = [  M    t    ]
    //         [  t'  alpha ]
    //
    DVectorSlice M12 = MU_store_.row(M_l_r_+M_size_+1)(M_l_c_,M_l_c_+M_size_-1);
    if(t)
      M12 = *t;
    else
      M12 = 0.0;
    MU_store_(M_l_r_+M_size_+1,M_l_c_+M_size_) = alpha;
  }
  ++M_size_; // Enlarge the matrix by one
  is_diagonal_ = false;
  if( throw_exception )
    throw WarnNearSingularUpdateException(omsg.str(),gamma);
}

void MatrixSymPosDefCholFactor::delete_update(
  size_type        jd
  ,bool            force_refactorization
  ,EEigenValType   drop_eigen_val
  ,PivotTolerances pivot_tols
  )
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::delete_udpate(...)" );
#endif
  typedef MatrixSymAddDelUpdateable MSADU;

  TEUCHOS_TEST_FOR_EXCEPTION(
    jd < 1 || M_size_ < jd, std::out_of_range
    ,"MatrixSymPosDefCholFactor::delete_update(jd,...): "
    "Error, the indice jd must be 1 <= jd <= rows()" );

  TEUCHOS_TEST_FOR_EXCEPT( !( drop_eigen_val == MSADU::EIGEN_VAL_UNKNOWN
    || (scale_ > 0.0 && drop_eigen_val == MSADU::EIGEN_VAL_POS)
    || (scale_ < 0.0 && drop_eigen_val == MSADU::EIGEN_VAL_NEG)
    ) );

  if( maintain_original_ ) {
    //
    // Here we have the lower portion of M partitioned as:
    //
    //  1 |\
    //    |  \
    //    |    \
    //    | M11  \
    //    |________\ _
    // jd |_________|_|
    //    |         | |\
    //    |         | |  \
    //    |         | |    \
    //    |   M31   | | M33  \
    //  n |         | |        \
    //    ----------------------
    //    1         jd         n
    //
    // We want to move M31 up one row and M33 up and to the left.
    //
    DenseLinAlgPack::delete_row_col( jd, &DenseLinAlgPack::nonconst_tri_ele( M().gms(), BLAS_Cpp::lower ) );
  }
  if( maintain_factor_ ) {
    //
    // Here we have U partitioned as:
    //
    // 1              jd      n        
    // ------------------------- 1
    // \             | |       |
    //   \    U11    | |  U13  |
    //     \         | |       |
    //       \  u12->| |  u23' |
    //         \     | |   |   |
    //           \   | |  \./  |
    //             \ |_|_______|
    //        u22 -> |_|_______| jd  
    //                 \       |
    //                   \ U33 |
    //                     \   |
    //                       \ | n
    //
    // To preform the update we need to first update U33 as:
    //   U33'*U33 + u23'*u23 ==> U33'*U33
    // Then we need to move U12 (one column at a time) to the
    // left one column and overwrite u12.
    // Then we will move the updated U33 (one column at a time)
    // up and to the left one position to overwrite u22 and
    // the left part of u23'.  We then decrease the dimension
    // of U by one and we are finished updating the factorization.
    //
    // See RAB notes from 1/21/99 and 1/26/99 for details on this update.
    //
    const size_type n = M_size_;
    // Resize workspace if it has not been done yet.
    size_type work_size = 3 * max_size_;
    if(work_.dim() < work_size)
      work_.resize(work_size);
    // Update the factors
    {
      DMatrixSlice U = this->U().gms();
      // Update U33 where it sits.
      if(jd < n) {
        size_type _size = n-jd;	// Set storage for u23, c and s
        Range1D rng(1,_size);
        DVectorSlice
          u23 = work_(rng),
          c   = work_(rng+_size),
          s   = work_(rng+2*_size);
        Range1D U_rng(jd+1,n);  // Set U33 and u23
        DMatrixSlice U33 = U(U_rng,U_rng);
        u23 = U.row(jd)(U_rng);
        // Update U33
        value_type dummy;
        ALAP_DCHUD_CALL (
          U33.col_ptr(1), U33.max_rows()
          ,U_rng.size(), u23.start_ptr(), &dummy, 1, 0, &dummy
          ,&dummy, c.start_ptr(), s.start_ptr() );
      }
      // Move U13 and U33 to delete row and column jd
      DenseLinAlgPack::delete_row_col( jd, &DenseLinAlgPack::nonconst_tri_ele( U, BLAS_Cpp::upper ) );
    }
  }
  else {
    factor_is_updated_ = false;
  }
  // Strink the size of M and U
  --M_size_;
}

// Overridden from Serializable

// ToDo: Refactor this code and create an external utility matrix
// serialization class that will convert from matrix type to matrix
// type.

void MatrixSymPosDefCholFactor::serialize( std::ostream &out ) const
{
  // Write key words on top line
  out << build_serialization_string() << std::endl;
  // Set the precision (very important!)
  out.precision(std::numeric_limits<value_type>::digits10+4);
  // Write the dimmension
  out << M_size_ << std::endl;
  if(M_size_) {
    // Write the matrix values
    if( maintain_original_ ) {
      const DMatrixSliceSym M = this->M();
      write_matrix( M.gms(), M.uplo(), out );
    }
    else {
      const DMatrixSliceTri U = this->U();
      write_matrix( U.gms(), U.uplo(), out );
    }
  }
  // ToDo: You need to write both M and U if both are computed!
}

void MatrixSymPosDefCholFactor::unserialize( std::istream &in )
{
  // Get the keywords for the matrix type
  std::string keywords;
  std::getline( in, keywords, '\n' );
  // For now make sure the types are exactly the same!
  const std::string this_keywords = build_serialization_string();
  TEUCHOS_TEST_FOR_EXCEPTION(
    this_keywords != keywords, std::logic_error
    ,"MatrixSymPosDefCholFactor::unserialize(...): Error, the matrix type being read in from file of "
    "\'"<<keywords<<"\' does not equal the type expected of \'"<<this_keywords<<"\'!"
    );
  // Read in the dimension of the matrix
  in >> M_size_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    M_size_ < 0, std::logic_error
    ,"MatrixSymPosDefCholFactor::unserialize(...): Error, read in a size of M_size = "<<M_size_<<" < 0!"
    );
  allocate_storage(M_size_);
  // Read in the matrix into storage
  if(maintain_original_) {
    M_l_r_ = M_l_c_ = 1;
    DMatrixSliceSym M = this->M();
    read_matrix( in, M.uplo(), &M.gms() );
  }
  else {
    U_l_r_ = U_l_c_ = 1;
    DMatrixSliceTri U = this->U();
    read_matrix( in, U.uplo(), &U.gms() );
  }
}

// Private member functions

void MatrixSymPosDefCholFactor::assert_storage() const
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  MU_store_.rows()  ) );
}

void MatrixSymPosDefCholFactor::allocate_storage(size_type max_size) const
{
  namespace rcp = MemMngPack;
  if( allocates_storage_ && MU_store_.rows() < max_size + 1 ) {
    // We have the right to allocate storage so lets just do it.
    Teuchos::RCP<DMatrix>
      MU_store = Teuchos::rcp(new DMatrix( max_size + 1, max_size + 1 ));
    typedef MemMngPack::ReleaseResource_ref_count_ptr<DMatrix> ptr_t;
    const_cast<MatrixSymPosDefCholFactor*>(this)->release_resource_ptr_ = Teuchos::rcp(new ptr_t(MU_store));
    const_cast<MatrixSymPosDefCholFactor*>(this)->MU_store_.bind( (*MU_store)() );
    const_cast<MatrixSymPosDefCholFactor*>(this)->max_size_ = max_size;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT( !(  MU_store_.rows()  >= max_size + 1  ) );
  }
}

void MatrixSymPosDefCholFactor::assert_initialized() const
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  M_size_  ) );
}

void MatrixSymPosDefCholFactor::resize_and_zero_off_diagonal(size_type n, value_type scale)
{
  using DenseLinAlgPack::nonconst_tri_ele;
  using DenseLinAlgPack::assign;

  TEUCHOS_TEST_FOR_EXCEPT( !(  n <= my_min( MU_store_.rows(), MU_store_.cols() ) - 1  ) );

  // Resize the views
  set_view(
    n, scale, maintain_original_, 1, 1, maintain_factor_
    ,U_l_r_ ? 1 : 0, U_l_r_ ? 1 : 0 );
  if( maintain_original_ ) {
    if(!is_diagonal_ && n > 1)
      assign( &nonconst_tri_ele( M().gms()(2,n,1,n-1), BLAS_Cpp::lower ), 0.0 );
  }
  if( U_l_r_ ) {
    if(!is_diagonal_ && n > 1)
      assign( &nonconst_tri_ele( U().gms()(1,n-1,2,n), BLAS_Cpp::upper ), 0.0 );
  }
}

void MatrixSymPosDefCholFactor::update_factorization() const
{
  if( factor_is_updated_ ) return; // The factor should already be updated.
  TEUCHOS_TEST_FOR_EXCEPTION(
    U_l_r_ == 0, std::logic_error
    ,"MatrixSymPosDefCholFactor::update_factorization() : "
    "Error, U_l_r == 0 was set in MatrixSymPosDefCholFactor::set_view(...) "
    "and therefore we can not update the factorization in the provided storage space." );
  // Update the factorization from scratch!
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::update_factorization(...) ... update" );
#endif
  MatrixSymPosDefCholFactor
    *nc_this = const_cast<MatrixSymPosDefCholFactor*>(this);
  DMatrixSliceTriEle U = DenseLinAlgPack::nonconst_tri_ele( nc_this->U().gms(), BLAS_Cpp::upper );
  DenseLinAlgPack::assign( &U, DenseLinAlgPack::tri_ele( M().gms(), BLAS_Cpp::lower ) );  // Copy in the original
  {
#ifdef PROFILE_HACK_ENABLED
    ProfileHackPack::ProfileTiming profile_timing( "MatrixSymPosDefCholFactor::update_factorization(...) ... potrf" );
#endif
    DenseLinAlgLAPack::potrf( &U );
  }
  nc_this->factor_is_updated_ = true;
}

std::string MatrixSymPosDefCholFactor::build_serialization_string() const
{
  std::string str = "SYMMETRIC POS_DEF";
  if( !maintain_original_ )
    str.append(" CHOL_FACTOR");
  if( maintain_original_ )
    str.append(" LOWER");
  else
    str.append(" UPPER");
  return str;
}

void MatrixSymPosDefCholFactor::write_matrix( const DMatrixSlice &Q, BLAS_Cpp::Uplo Q_uplo, std::ostream &out )
{
  const int Q_dim = Q.rows();
  if( Q_uplo == BLAS_Cpp::lower ) {
    for( int i = 1; i <= Q_dim; ++i ) {
      for( int j = 1; j <= i; ++j ) {
        out << " " << Q(i,j);
      }
      out << std::endl;
    }
  }
  else {
    for( int i = 1; i <= Q_dim; ++i ) {
      for( int j = i; j <= Q_dim; ++j ) {
        out << " " << Q(i,j);
      }
      out << std::endl;
    }
  }
}

void MatrixSymPosDefCholFactor::read_matrix(  std::istream &in, BLAS_Cpp::Uplo Q_uplo, DMatrixSlice *Q_out )
{
  DMatrixSlice &Q = *Q_out;
  const int Q_dim = Q.rows();
  if( Q_uplo == BLAS_Cpp::lower ) {
    for( int i = 1; i <= Q_dim; ++i ) {
      for( int j = 1; j <= i; ++j ) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(
          in.eof(), std::logic_error
          ,"MatrixSymPosDefCholFactor::read_matrix(in,lower,Q_out): Error, not finished reading in matrix yet (i="<<i<<",j="<<j<<")!"
          );
#endif
        in >> Q(i,j);
      }
    }
  }
  else {
    for( int i = 1; i <= Q_dim; ++i ) {
      for( int j = i; j <= Q_dim; ++j ) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(
          in.eof(), std::logic_error
          ,"MatrixSymPosDefCholFactor::read_matrix(in,upper,Q_out): Error, not finished reading in matrix yet (i="<<i<<",j="<<j<<")!"
          );
#endif
        in >> Q(i,j);
      }
    }
  }
}

} // end namespace AbstractLinAlgPack
