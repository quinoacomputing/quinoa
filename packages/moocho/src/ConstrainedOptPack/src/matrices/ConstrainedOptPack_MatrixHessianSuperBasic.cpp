#if 0

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

#include "ConstrainedOptPack_MatrixHessianSuperBasic.hpp"
#include "ConstrainedOptPack_initialize_Q_R_Q_X.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOpOut.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptPack {

MatrixHessianSuperBasic::MatrixHessianSuperBasic()
  : n_(0)
{}

void MatrixHessianSuperBasic::initialize(
  size_type            n
  ,size_type           n_R
  ,const size_type     i_x_free[]
  ,const size_type     i_x_fixed[]
  ,const EBounds       bnd_fixed[]
  ,const B_RR_ptr_t&   B_RR_ptr
  ,const B_RX_ptr_t&   B_RX_ptr
  ,BLAS_Cpp::Transp    B_RX_trans
  ,const B_XX_ptr_t&   B_XX_ptr
  )
{
  using DenseLinAlgPack::Mp_M_assert_sizes;
  using BLAS_Cpp::no_trans;

  const size_type
    n_X = n - n_R;

    // Validate input arguments

  // i_x_free
  if( 0 < n_R && n_R < n && i_x_free == NULL ) {
    throw std::invalid_argument(
      "MatrixHessianSuperBasic::initialize(...) : Error, "
      "i_x_free can not be NULL when 0 < n_R < n" );
  }
  // i_x_fixed
  if( 0 < n_X && n_X < n && i_x_fixed == NULL ) {
    throw std::invalid_argument(
      "MatrixHessianSuperBasic::initialize(...) : Error, "
      "i_x_fixed can not be NULL when 0 < n-n_R < n" );
  }
  // bnd_fixed
  if( 0 < n_X && bnd_fixed == NULL ) {
    throw std::invalid_argument(
      "MatrixHessianSuperBasic::initialize(...) : Error, "
      "bnd_fixed can not be NULL when 0 < n-n_R" );
  }
  // B_RR
  if(n_R > 0 ) {
    if( !B_RR_ptr.get() )
      throw std::invalid_argument(
        "MatrixHessianSuperBasic::initialize(...) : Error, "
        "B_RR_ptr.get() can not be NULL when n_R > 0" );
    Mp_M_assert_sizes( n_R, n_R, no_trans, B_RR_ptr->rows(), B_RR_ptr->cols(), no_trans );
  }
  // op(B_RX)
  if( n_R < n ) {
    if( B_RX_ptr.get() ) {
      Mp_M_assert_sizes( n_R, n_X, no_trans, B_RX_ptr->rows(), B_RX_ptr->cols(), B_RX_trans );
    }
  }
  // B_XX
  if( n_R < n ) {
    if( !B_XX_ptr.get() )
      throw std::invalid_argument(
        "MatrixHessianSuperBasic::initialize(...) : Error, "
        "B_XX_ptr.get() can not be NULL if n_R < n" );
    Mp_M_assert_sizes( n_X, n_X, no_trans, B_XX_ptr->rows(), B_XX_ptr->cols(), no_trans );
  }

  // Setup Q_R and Q_X and validate i_x_free[] and i_x_fixed[]
  const bool Q_R_is_idenity = (n_R == n && i_x_fixed == NULL );
  if( Q_R_is_idenity ) {
    Q_R_row_i_.resize(0);
    Q_R_col_j_.resize(0);
  }
  else {
    Q_R_row_i_.resize(n_R);
    Q_R_col_j_.resize(n_R);
  }
  Q_X_row_i_.resize(n_X);
  Q_X_col_j_.resize(n_X);
  bool test_setup = true;  // ToDo: Make this an input parameter!
  initialize_Q_R_Q_X(
    n_R,n_X,i_x_free,i_x_fixed,test_setup
    ,!Q_R_is_idenity ? &Q_R_row_i_[0] : NULL
    ,!Q_R_is_idenity ? &Q_R_col_j_[0] : NULL
    ,&Q_R_
    ,n_X ? &Q_X_row_i_[0] : NULL
    ,n_X ? &Q_X_col_j_[0] : NULL
    ,&Q_X_
  );

  // Setup bnd_fixed
  bnd_fixed_.resize(n_X);
  {for(size_type i = 0; i < n_X; ++i) bnd_fixed_[i] = bnd_fixed[i]; }

  // Set the rest of the arguments
  n_           = n;
  n_R_         = n_R;
  B_RR_ptr_    = B_RR_ptr;
  B_RX_ptr_    = B_RX_ptr;
  B_RX_trans_  = B_RX_trans;
  B_XX_ptr_    = B_XX_ptr;

}

// Overridden from Matrix

size_type MatrixHessianSuperBasic::rows() const
{
  return n_;
}

// Overridden from MatrixOp

void MatrixHessianSuperBasic::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp B_trans
  , const DVectorSlice& x, value_type b
  ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;
  using AbstractLinAlgPack::V_MtV;
  using LinAlgOpPack::V_MtV;
  assert_initialized();
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->size(), n_, n_, B_trans, x.size() );
  if( n_ == n_R_ ) {
    //
    // B = Q_R*B_RR*Q_R'
    //
    // y = b*y + a*Q_R*B_RR*Q_R'*x
    //
    if( Q_R().is_identity() ) {
      AbstractLinAlgPack::Vp_StMtV(y,a,*this->B_RR_ptr(),no_trans,x,b);
    }
    else {
      DVector Q_R_x;
      V_MtV( &Q_R_x, Q_R(), trans, x );
      AbstractLinAlgPack::Vp_StPtMtV(y,a,Q_R(),no_trans,*this->B_RR_ptr(),no_trans,Q_R_x(),b);
    }
  }
  else if( n_R_ == 0 ) {
    //
    // B = Q_X *B_XX * Q_X'
    //
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement this!
  }
  else {
    //
    // B = [ Q_R  Q_X  ] * [   B_RR      op(B_RX) ] * [ Q_R' ]
    //                     [ op(B_RX')      B_XX  ]   [ Q_X' ]
    //
    // y = b*y + a*op(B)*x
    //
    // y = b*y + a * [ Q_R  Q_X ] * [   B_RR      op(B_RX) ] * [ Q_R' ] * x
    //                              [ op(B_RX')      B_XX  ]   [ Q_X' ]
    //
    // y = b*y + a*Q_R*B_RR*x_R      + a*Q_R*op(B_RX)*x_X
    //         + a*Q_X*op(B_RX')*x_R + a*Q_X*B_XX*x_X
    // where:
    //     x_R = Q_R'*x
    //     x_X = Q_X'*x
    //
    SpVector
      x_R,
      x_X;
    // x_R = Q_R'*x
    V_MtV( &x_R, Q_R(), trans, x );
    // x_X = Q_X'*x
    V_MtV( &x_X, Q_X(), trans, x );
    // y = b*y + a*Q_R*B_RR*x_R
    AbstractLinAlgPack::Vp_StPtMtV(
      y, a, Q_R(), no_trans, *B_RR_ptr(), no_trans, x_R(), b );
    // y += a*Q_R*op(B_RX)*x_X + a*Q_X*op(B_RX')*x_R
    if( B_RX_ptr().get() ) {
      AbstractLinAlgPack::Vp_StPtMtV(
        y, a, Q_R(), no_trans, *B_RX_ptr(), B_RX_trans(), x_X() );
      AbstractLinAlgPack::Vp_StPtMtV(
        y, a, Q_X(), no_trans, *B_RX_ptr(), trans_not(B_RX_trans()), x_R() );
    }
    // y += a*Q_X*B_XX*x_X
    AbstractLinAlgPack::Vp_StPtMtV(
      y, a, Q_X(), no_trans, *B_XX_ptr(), no_trans, x_X() );
  }
}

void MatrixHessianSuperBasic::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp B_trans
  , const SpVectorSlice& x, value_type b
  ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;
  using AbstractLinAlgPack::V_MtV;
  using LinAlgOpPack::V_MtV;
  assert_initialized();
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->size(), n_, n_, B_trans, x.size() );
  if( n_ == n_R_ ) {
    //
    // B = Q_R*B_RR*Q_R'
    //
    // y = b*y + a*Q_R*B_RR*Q_R'*x
    //
    if( Q_R().is_identity() ) {
      AbstractLinAlgPack::Vp_StMtV(y,a,*this->B_RR_ptr(),no_trans,x,b);
    }
    else {
      SpVector Q_R_x;
      AbstractLinAlgPack::V_MtV( &Q_R_x, Q_R(), trans, x );
      AbstractLinAlgPack::Vp_StPtMtV(y,a,Q_R(),no_trans,*this->B_RR_ptr(),no_trans,Q_R_x(),b);
    }
  }
  else if( n_R_ == 0 ) {
    //
    // B = Q_X *B_XX * Q_X'
    //
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement this!
  }
  else {
    //
    // B = [ Q_R  Q_X  ] * [   B_RR      op(B_RX) ] * [ Q_R' ]
    //                     [ op(B_RX')      B_XX  ]   [ Q_X' ]
    //
    // y = b*y + a*op(B)*x
    //
    // y = b*y + a * [ Q_R  Q_X ] * [   B_RR      op(B_RX) ] * [ Q_R' ] * x
    //                              [ op(B_RX')      B_XX  ]   [ Q_X' ]
    //
    // y = b*y + a*Q_R*B_RR*x_R      + a*Q_R*op(B_RX)*x_X
    //         + a*Q_X*op(B_RX')*x_R + a*Q_X*B_XX*x_X
    // where:
    //     x_R = Q_R'*x
    //     x_X = Q_X'*x
    //
    SpVector
      x_R,
      x_X;
    // x_R = Q_R'*x
    V_MtV( &x_R, Q_R(), trans, x );
    // x_X = Q_X'*x
    V_MtV( &x_X, Q_X(), trans, x );
    // y = b*y + a*Q_R*B_RR*x_R
    AbstractLinAlgPack::Vp_StPtMtV(
      y, a, Q_R(), no_trans, *B_RR_ptr(), no_trans, x_R(), b );
    // y += a*Q_R*op(B_RX)*x_X + a*Q_X*op(B_RX')*x_R
    if( B_RX_ptr().get() ) {
      AbstractLinAlgPack::Vp_StPtMtV(
        y, a, Q_R(), no_trans, *B_RX_ptr(), B_RX_trans(), x_X() );
      AbstractLinAlgPack::Vp_StPtMtV(
        y, a, Q_X(), no_trans, *B_RX_ptr(), trans_not(B_RX_trans()), x_R() );
    }
    // y += a*Q_X*B_XX*x_X
    AbstractLinAlgPack::Vp_StPtMtV(
      y, a, Q_X(), no_trans, *B_XX_ptr(), no_trans, x_X() );
  }
}

void MatrixHessianSuperBasic::Vp_StPtMtV(
  DVectorSlice* y, value_type a
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp M_trans
  , const DVectorSlice& x, value_type b ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;
  using AbstractLinAlgPack::V_MtV;
  using DenseLinAlgPack::Vt_S;
  using LinAlgOpPack::V_MtV;
  namespace slap = AbstractLinAlgPack;

  assert_initialized();

  //
  // y = b*y + a * op(P) * B * x
  //
  // =>
  //
  // y = b*y + a * op(P)*(Q_R*B_RR*Q_R' + Q_R*op(B_RX)*Q_X' + Q_X*op(B_RX')*Q_R + Q_X*B_XX*Q_X')*x
  //   
  //   = b*y + a*op(P)*Q_R*B_RR*Q_R'*x     + a*op(P)*Q_R*op(B_RX)*Q_X'*x
  //         + a*op(P)*Q_X*op(B_RX')*Q_R*x + a*op(P)*Q_X*B_XX*Q_X'*x
  //
  // In order to implement the above as efficiently as possible we need to minimize the
  // computations with the constituent matrices.  First off we will compute
  // Q_RT_x = Q_R'*x (O(n_R)) and Q_XT_x = Q_X'*x (O(n_R)) neglect any terms where
  // Q_RT_x.nz() == 0 or Q_XT_x.nz() == 0.  We will also determine if op(P)*Q_R == 0 (O(n_R))
  // or op(P)*Q_X == 0 (O(n_X)) and neglect these terms if the are zero.
  // Hopefully this work will allow us to skip as many computations as possible.
  //
  LinAlgOpPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans
    , BLAS_Cpp::rows( rows(), cols(), M_trans) );
  LinAlgOpPack::Vp_MtV_assert_sizes( BLAS_Cpp::cols( P.rows(), P.cols(), P_trans)
    ,rows(),cols(),M_trans,x.size());
  // Q_R'*x
  SpVector Q_RT_x;
  if(n_R_) {
    slap::V_MtV( &Q_RT_x, Q_R(), trans, x );
  }
  // Q_X'*x
  SpVector Q_XT_x;
  if(n_ > n_R_) {
    slap::V_MtV( &Q_XT_x, Q_X(), trans, x );
  }
  // op(P)*Q_R overlap
  size_type P_Q_R_nz = 0;
  AbstractLinAlgPack::intersection( P, P_trans, Q_R(), no_trans, &P_Q_R_nz );
  // op(P)*Q_X overlap
  size_type P_Q_X_nz = 0;
  AbstractLinAlgPack::intersection( P, P_trans, Q_X(), no_trans, &P_Q_X_nz );
  // y = b*y
  if(b==0.0)      *y = 0.0;
  else if(b!=1.0) Vt_S(y,b);
  // 
  DVector t; // ToDo: use workspace
  // y += a*op(P)*Q_R*B_RR*Q_R'*x
  if( P_Q_R_nz && Q_RT_x.nz() ) {
    t.resize(n_);
    slap::Vp_StPtMtV( &t(), 1.0, Q_R(), no_trans, *B_RR_ptr(), no_trans, Q_RT_x() );
    slap::Vp_StMtV( y, a, P, P_trans, t() );
  }
  // y += a*op(P)*Q_R*op(B_RX)*Q_X'*x
  if( P_Q_R_nz && B_RX_ptr().get() && Q_XT_x.nz() ) {
    t.resize(n_);
    slap::Vp_StPtMtV( &t(), 1.0, Q_R(), no_trans, *B_RX_ptr(), B_RX_trans(), Q_XT_x() );
    slap::Vp_StMtV( y, a, P, P_trans, t() );
  }
  // y += a*op(P)*Q_X*op(B_RX')*Q_R*x
  if( P_Q_X_nz && B_RX_ptr().get() && Q_RT_x.nz() ) {
    t.resize(n_);
    slap::Vp_StPtMtV( &t(), 1.0, Q_X(), no_trans, *B_RX_ptr(), trans_not(B_RX_trans()), Q_RT_x() );
    slap::Vp_StMtV( y, a, P, P_trans, t() );
  }
  // y += a*op(P)*Q_X*B_XX*Q_X'*x
  if( P_Q_X_nz && Q_XT_x.nz() ) {
    t.resize(n_);
    slap::Vp_StPtMtV( &t(), 1.0, Q_X(), no_trans, *B_XX_ptr(), no_trans, Q_XT_x() );
    slap::Vp_StMtV( y, a, P, P_trans, t() );
  }
}

value_type MatrixHessianSuperBasic::transVtMtV(
  const SpVectorSlice& x1, BLAS_Cpp::Transp B_trans
  , const SpVectorSlice& x2 ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  assert_initialized();
  DenseLinAlgPack::Vp_MtV_assert_sizes( x1.size(), rows(), cols(), B_trans, x1.size() );
  if( n_ == n_R_ ) {
    //
    // B = Q_R*B_RR*Q_R'
    //
    // a = x1'*Q_R*B_RR*Q_R'*x2
    //
    if( Q_R().is_identity() ) {
      return AbstractLinAlgPack::transVtMtV( x1, *B_RR_ptr(), no_trans, x2 );
    }
    else {
      if( x1.overlap(x2) == DenseLinAlgPack::SAME_MEM ) {
        SpVector Q_RT_x2;
        AbstractLinAlgPack::V_MtV( &Q_RT_x2, Q_R(), trans, x2 );
        SpVectorSlice Q_RT_x2_slc = Q_RT_x2();
        return AbstractLinAlgPack::transVtMtV(
          Q_RT_x2_slc, *B_RR_ptr(), no_trans, Q_RT_x2_slc );
       }
      else {
        SpVector Q_RT_x2;
        AbstractLinAlgPack::V_MtV( &Q_RT_x2, Q_R(), trans, x2 );
        SpVector Q_RT_x1;
        AbstractLinAlgPack::V_MtV( &Q_RT_x1, Q_R(), trans, x1 );
        return AbstractLinAlgPack::transVtMtV(
          Q_RT_x1(), *B_RR_ptr(), no_trans, Q_RT_x2() );
      }
    }
  }
  else if( n_R_ == 0 ) {
    //
    // B = Q_X *B_XX * Q_X'
    //
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement this!
  }
  else {
    //
    // B = [ Q_R  Q_X  ] * [   B_RR      op(B_RX) ] * [ Q_R' ]
    //                     [ op(B_RX')      B_XX  ]   [ Q_X' ]
    //
    //
    // a = x1'*B*x2
    // =>
    // a = x1' * [ Q_R  Q_X  ] * [   B_RR      op(B_RX) ] * [ Q_R' ] * x2
    //                           [ op(B_RX')      B_XX  ]   [ Q_X' ]
    //
    // a = x1'*Q_R*B_RR*Q_R'*x2 + 2*x1'*Q_R*op(B_RX)*Q_X'*x2 + x1'*Q_X*B_XX*Q_X'*x2
    //
    if( x1.overlap(x2) == DenseLinAlgPack::SAME_MEM ) {
      // a = x1'*Q_R*B_RR*Q_R'*x1 + 2*x1'*Q_R*op(B_RX)*Q_X'*x1 + x1'*Q_X*B_XX*Q_X'*x1
      SpVector Q_RT_x1;
      if( Q_R().nz() )
        AbstractLinAlgPack::V_MtV( &Q_RT_x1, Q_R(), trans, x1 );
      SpVector Q_XT_x1;
      if( Q_X().nz() )
        AbstractLinAlgPack::V_MtV( &Q_XT_x1, Q_X(), trans, x1 );
      SpVectorSlice Q_RT_x1_slc = Q_RT_x1();
      SpVectorSlice Q_XT_x1_slc = Q_XT_x1();
      return
        ( Q_R().nz()
          ? AbstractLinAlgPack::transVtMtV(
            Q_RT_x1_slc, *B_RR_ptr(), no_trans, Q_RT_x1_slc )
          : 0.0
          )
        + 2*(  B_RX_ptr().get() && Q_R().nz() && Q_X().nz()
             ? AbstractLinAlgPack::transVtMtV(
               Q_RT_x1_slc, *B_RX_ptr(), B_RX_trans(), Q_XT_x1_slc )
             : 0.0
          )
        + ( Q_X().nz()
          ? AbstractLinAlgPack::transVtMtV(
            Q_XT_x1_slc, *B_XX_ptr(), no_trans, Q_XT_x1_slc )
          : 0.0
          );
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement this!
    }
  }
  return 0.0; // Will never be executed!
}

// Private

void MatrixHessianSuperBasic::assert_initialized() const
{
  if( !n_ )
    throw std::logic_error(
      "MatrixHessianSuperBasic::assert_initialized() : Error, "
      "The matrix is not initialized yet" );
}

} // end namespace ConstrainedOptPack

#endif // 0
