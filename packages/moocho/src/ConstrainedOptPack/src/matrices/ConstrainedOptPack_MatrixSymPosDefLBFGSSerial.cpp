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
//
// Let's define a compact representation for the matrix B^{k} and
// its inverse H^{k} = inv(B^{k}).
//
// Bk = (1/gk)*I - [ (1/gk)*S  Y ] * inv( [ (1/gk)*S'S   L ]     [ (1/gk)*S' ]
//                                        [    L'       -D ] ) * [    Y'     ]
//                                        \________________/
//                                                Q
//
// Hk = gk*I + [ S  gk*Y ] * [ inv(R')*(D+gk*Y'Y)*inv(R)     -inv(R') ] * [   S'  ]
//                           [            -inv(R)                0    ]   [ gk*Y' ]
//
// where:
//
// gk = gamma_k <: R
//
//			[	s^{1}'*y^{1}	s^{1}'*y^{2}	...		s^{1}'*y^{m}	]
//	S'Y =	[	s^{2}'*y^{1}	s^{2}'*y^{2}	...		s^{2}'*y^{m}	] <: R^(m x m)
//			[	.				.						.				]
//			[	s^{m}'*y^{1}	s^{m}'*y^{2}	...		s^{m}'*y^{m}	]
// 
//			[	s^{1}'*y^{1}	0				...		0				]
//	D =		[	0				s^{2}'*y^{2}	...		0				] <: R^(m x m)
//			[	.				.						.				]
//			[	0				0				...		s^{m}'*y^{m}	]
//
//	R = upper triangular part of S'Y
// 
//	L =	lower tirangular part of S'Y with zeros on the diagonal
//

#include <assert.h>

#include "ConstrainedOptPack_MatrixSymPosDefLBFGS.hpp"
#include "ConstrainedOptPack/src/AbstractLinAlgPack_BFGS_helpers.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgLAPack.hpp"

namespace {

  using DenseLinAlgPack::DVectorSlice;
  using DenseLinAlgPack::DMatrixSlice;

  /** \brief Compute Cb = Lb * inv(Db) * Lb' (see update_Q()).
    *
    * Here:
    *		Lb is lower triangular.
    *		Cb is upper triangular.
    *		Db_diag is the diagonal of Db
    */
  void comp_Cb( const DMatrixSlice& Lb, const DVectorSlice& Db_diag
    , DMatrixSlice* Cb );

}	// end namespace

namespace ConstrainedOptPack {

// /////////////////////////////////
// Inline private member functions

inline
const DMatrixSliceTri MatrixSymPosDefLBFGS::R() const
{
  return DenseLinAlgPack::tri( STY_(1,m_bar_,1,m_bar_), BLAS_Cpp::upper, BLAS_Cpp::nonunit );
}

inline
const DMatrixSliceTri MatrixSymPosDefLBFGS::Lb() const
{
  return DenseLinAlgPack::tri( STY_(2,m_bar_,1,m_bar_-1), BLAS_Cpp::lower, BLAS_Cpp::nonunit );
}

inline
DMatrixSlice MatrixSymPosDefLBFGS::STY()
{
  return STY_(1,m_bar_,1,m_bar_);
}

inline
const DMatrixSlice MatrixSymPosDefLBFGS::STY() const
{
  return STY_(1,m_bar_,1,m_bar_);
}

inline
DMatrixSliceSym MatrixSymPosDefLBFGS::STS()
{
  return DenseLinAlgPack::nonconst_sym( STSYTY_(2,m_bar_+1,1,m_bar_),BLAS_Cpp::lower );
}

inline
const DMatrixSliceSym MatrixSymPosDefLBFGS::STS() const
{
  return DenseLinAlgPack::sym( STSYTY_(2,m_bar_+1,1,m_bar_),BLAS_Cpp::lower );
}

inline
DMatrixSliceSym MatrixSymPosDefLBFGS::YTY()
{
  return DenseLinAlgPack::nonconst_sym( STSYTY_(1,m_bar_,2,m_bar_+1),BLAS_Cpp::upper );
}

inline
const DMatrixSliceSym MatrixSymPosDefLBFGS::YTY() const
{
  return DenseLinAlgPack::sym( STSYTY_(1,m_bar_,2,m_bar_+1),BLAS_Cpp::upper );
}

// ///////////////////////
// Nonlinined functions

MatrixSymPosDefLBFGS::MatrixSymPosDefLBFGS(
  size_type   max_size
    ,size_type  m
  ,bool       maintain_original
  ,bool       maintain_inverse
  ,bool       auto_rescaling
  )
{
  initial_setup(max_size,m,maintain_original,maintain_inverse,auto_rescaling);
}

void MatrixSymPosDefLBFGS::initial_setup(
  size_type   max_size
    ,size_type  m
  ,bool       maintain_original
  ,bool       maintain_inverse
  ,bool       auto_rescaling
  )
{
  // Validate input
  if( !maintain_original && !maintain_inverse )
    throw std::invalid_argument(
      "MatrixSymPosDefLBFGS::initial_setup(...) : "
      "Error, both maintain_original and maintain_inverse can not both be false!" );
  if( m < 1 )
    throw std::invalid_argument(
      "MatrixSymPosDefLBFGS::set_num_updates_stored(m) : "
      "Error, the number of storage locations must be > 0" );
  maintain_original_ = maintain_original;
  maintain_inverse_  = maintain_inverse;
  m_                 = m;
  n_                 = 0; // make uninitialized
  n_max_             = max_size;
  num_secant_updates_= 0;
}

// Overridden from Matrix

size_type MatrixSymPosDefLBFGS::rows() const
{
  return n_;
}

// Overridden from MatrixOp

std::ostream& MatrixSymPosDefLBFGS::output(std::ostream& out) const
{
  assert_initialized();
  out << "*** Limited Memory BFGS matrix.\n"
    << "Conversion to dense =\n";
  MatrixOp::output(out);
  out << "\n*** Stored quantities\n"
    << "\ngamma_k = " << gamma_k_ << std::endl;
  if( m_bar_ ) {
    out	<< "\nS =\n" << S()
      << "\nY =\n" << Y()
      << "\nS'Y =\n" << STY_(1,m_bar_,1,m_bar_)
      << "\nlower(S'S) \\ zero diagonal \\ upper(Y'Y) =\n"
        << STSYTY_(1,m_bar_+1,1,m_bar_+1)
      << "\nQ updated? = " << Q_updated_ << std::endl
      << "\nCholesky of schur complement of Q, QJ =\n" << QJ_(1,m_bar_,1,m_bar_);
  }
  return out;
}

MatrixOp& MatrixSymPosDefLBFGS::operator=(const MatrixOp& m)
{	
  const MatrixSymPosDefLBFGS *p_m = dynamic_cast<const MatrixSymPosDefLBFGS*>(&m);
  if(p_m) {
    if( p_m == this ) return *this;	// assignment to self
    // Important: Assign all members here.
    auto_rescaling_      = p_m->auto_rescaling_;
    maintain_original_   = p_m->maintain_original_;
    original_is_updated_ = p_m->original_is_updated_;
    maintain_inverse_    = p_m->maintain_inverse_;
    inverse_is_updated_  = p_m->inverse_is_updated_;
    n_max_ 		         = p_m->n_max_;
    n_	 		         = p_m->n_;
    m_			         = p_m->m_;
    m_bar_		         = p_m->m_bar_;
    k_bar_		         = p_m->k_bar_;
    num_secant_updates_  = p_m->num_secant_updates_;
    gamma_k_	         = p_m->gamma_k_;
    S_			         = p_m->S_;
    Y_			         = p_m->Y_;
    STY_		         = p_m->STY_;
    STSYTY_		         = p_m->STSYTY_;
    Q_updated_	         = p_m->Q_updated_;
    QJ_			         = p_m->QJ_;
  }
  else {
    throw std::invalid_argument("MatrixSymPosDefLBFGS::operator=(const MatrixOp& m)"
      " : The concrete type of m is not a subclass of MatrixSymPosDefLBFGS as expected" );
  }
  return *this;
}

// Level-2 BLAS

void MatrixSymPosDefLBFGS::Vp_StMtV(
    DVectorSlice* y, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& x, value_type beta) const
{
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::V_MtV;

  using DenseLinAlgPack::Vt_S;
  using DenseLinAlgPack::Vp_StV;
  using DenseLinAlgPack::Vp_StMtV;

  assert_initialized();

  TEUCHOS_TEST_FOR_EXCEPT( !(  original_is_updated_  ) ); // For now just always update

  // y = b*y + Bk * x
  //
  // y = b*y + (1/gk)*x - [ (1/gk)*S  Y ] * inv(Q) * [ (1/gk)*S' ] * x
  //                                                 [     Y'    ]
  // Perform the following operations (in order):
  //
  // y = b*y
  //
  // y += (1/gk)*x
  //
  // t1 = [ (1/gk)*S'*x ]		<: R^(2*m)
  //		[      Y'*x   ]
  //
  // t2 =	inv(Q) * t1			<: R^(2*m)
  //
  // y += -(1/gk) * S * t2(1:m)
  //
  // y += -1.0 * Y * t2(m+1,2m)

  const value_type
    invgk = 1.0 / gamma_k_;

  // y = b*y

  if( beta == 0.0 )
    *y = beta;
  else
    Vt_S( y, beta );

  // y += (1/gk)*x

  Vp_StV( y, invgk, x );

  if( !m_bar_ )
    return;	// No updates have been added yet.

  // Get workspace

  if( work_.size() < 4 * m_ )
    work_.resize( 4 * m_ );

  const size_type
    mb = m_bar_;

  const size_type
    t1s = 1,
    t1n = 2*mb,
    t2s = t1s+t1n,
    t2n = 2*mb;

  DVectorSlice
    t1 = work_(	t1s,	t1s + t1n - 1	),
    t2 = work_(	t2s,	t2s + t2n - 1 	);

  const DMatrixSlice
    &S = this->S(),
    &Y = this->Y();

  // t1 = [ (1/gk)*S'*x ]
  //		[      Y'*x   ]

  V_StMtV( &t1(1,mb), invgk, S, BLAS_Cpp::trans, x );
  V_MtV( &t1(mb+1,2*mb), Y, BLAS_Cpp::trans, x );

  // t2 =	inv(Q) * t1

  V_invQtV( &t2, t1 );

  // y += -(1/gk) * S * t2(1:m)

  Vp_StMtV( y, -invgk, S, BLAS_Cpp::no_trans, t2(1,mb) );

  // y += -1.0 * Y * t2(m+1,2m)

  Vp_StMtV( y, -1.0, Y, BLAS_Cpp::no_trans, t2(mb+1,2*mb) );

}

// Overridden from MatrixWithOpFactorized

// Level-2 BLAS

void MatrixSymPosDefLBFGS::V_InvMtV( DVectorSlice* y, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& x ) const
{
  using DenseLinAlgPack::V_mV;
  using DenseLinAlgPack::V_StV;
  using DenseLinAlgPack::V_InvMtV;

  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::Vp_MtV;
  using LinAlgOpPack::Vp_StMtV;

  assert_initialized();

  TEUCHOS_TEST_FOR_EXCEPT( !(  inverse_is_updated_  ) ); // For now just always update

  // y = inv(Bk) * x = Hk * x
  //
  // = gk*x + [S gk*Y] * [ inv(R')*(D+gk*Y'Y)*inv(R)     -inv(R') ] * [   S'  ] * x
  //                     [            -inv(R)                0    ]   [ gk*Y' ]
  //
  // Perform in the following (in order):
  //
  // y = gk*x
  //
  // t1 = [   S'*x  ]					<: R^(2*m)
  //      [ gk*Y'*x ]
  //
  // t2 = inv(R) * t1(1:m)			<: R^(m)
  //
  // t3 = - inv(R') * t1(m+1,2*m)		<: R^(m)
  //
  // t4 = gk * Y'Y * t2				<: R^(m)
  //
  // t4 += D*t2
  //
  // t5 = inv(R') * t4				<: R^(m)
  //
  // t5 += t3
  //
  // y += S*t5
  //
  // y += -gk*Y*t2

  // y = gk*x
  V_StV( y, gamma_k_, x );

  const size_type
    mb = m_bar_;	
  
  if( !mb )
    return;	// No updates have been performed.

  // Get workspace

  if( work_.size() < 6*m_ )
    work_.resize( 6*m_ );

  const size_type
    t1s		= 1,
    t1n		= 2*mb,
    t2s		= t1s + t1n,
    t2n		= mb,
    t3s		= t2s + t2n,
    t3n		= mb,
    t4s		= t3s + t3n,
    t4n		= mb,
    t5s		= t4s + t4n,
    t5n		= mb;

  DVectorSlice
    t1	= work_( t1s, t1s + t1n - 1 ),
    t2	= work_( t2s, t2s + t2n - 1 ),
    t3	= work_( t3s, t3s + t3n - 1 ),
    t4	= work_( t4s, t4s + t4n - 1 ),
    t5	= work_( t5s, t5s + t5n - 1 );

  const DMatrixSlice
    &S = this->S(),
    &Y = this->Y();

  const DMatrixSliceTri
    &R = this->R();

  const DMatrixSliceSym
    &YTY = this->YTY();

  // t1 = [   S'*x  ]
  //      [ gk*Y'*x ]
  V_MtV( &t1(1,mb), S, BLAS_Cpp::trans, x );
  V_StMtV( &t1(mb+1,2*mb), gamma_k_, Y, BLAS_Cpp::trans, x );

  // t2 = inv(R) * t1(1:m)
  V_InvMtV( &t2, R, BLAS_Cpp::no_trans, t1(1,mb) );

  // t3 = - inv(R') * t1(m+1,2*m)
  V_mV( &t3, t1(mb+1,2*mb) );
  V_InvMtV( &t3, R, BLAS_Cpp::trans, t3 );

  // t4 = gk * Y'Y * t2
  V_StMtV( &t4, gamma_k_, YTY, BLAS_Cpp::no_trans, t2 );

  // t4 += D*t2
  Vp_DtV( &t4, t2 );

  // t5 = inv(R') * t4
  V_InvMtV( &t5, R, BLAS_Cpp::trans, t4 );

  // t5 += t3
  Vp_V( &t5, t3 );

  // y += S*t5
  Vp_MtV( y, S, BLAS_Cpp::no_trans, t5 );

  // y += -gk*Y*t2
  Vp_StMtV( y, -gamma_k_, Y, BLAS_Cpp::no_trans, t2 );

}

// Overridden from MatrixSymSecant

void MatrixSymPosDefLBFGS::init_identity(size_type n, value_type alpha)
{
  // Validate input
  if( alpha <= 0.0 ) {
    std::ostringstream omsg;
    omsg
      << "MatrixSymPosDefLBFGS::init_identity(n,alpha) : Error, "
      << "alpha = " << alpha << " <= 0 is not allowed!";
    throw std::invalid_argument( omsg.str() );
  }
  if( n_max_ == 0 ) {
    n_max_ = n;
  }
  else if( n > n_max_ ) {
    std::ostringstream omsg;
    omsg
      << "MatrixSymPosDefLBFGS::init_identity(n,alpha) : Error, "
      << "n = " << n << " > max_size = " << n_max_;
    throw std::invalid_argument( omsg.str() );
  }

  // Resize storage
  S_.resize( n_max_, m_ );
  Y_.resize( n_max_, m_ );
  STY_.resize( m_, m_ );
  STSYTY_.resize( m_+1, m_+1 );
  STSYTY_.diag(0) = 0.0;

  gamma_k_ = 1.0/alpha;

  // Initialize counters
  k_bar_	= 0;
  m_bar_	= 0;

  n_ = n;	 // initialized;
  original_is_updated_ = true; // This will never change for now
  inverse_is_updated_  = true; // This will never change for now
  num_secant_updates_  = 0;
}

void MatrixSymPosDefLBFGS::init_diagonal(const DVectorSlice& diag)
{
  using DenseLinAlgPack::norm_inf;
  init_identity( diag.size(), norm_inf(diag) );
}

void MatrixSymPosDefLBFGS::secant_update(
  DVectorSlice* s, DVectorSlice* y, DVectorSlice* Bs)
{
  using DenseLinAlgPack::dot;
  using DenseLinAlgPack::norm_2;

  using LinAlgOpPack::V_MtV;

  assert_initialized();

  // Check skipping the BFGS update
  const value_type
    sTy	      = dot(*s,*y);
  std::ostringstream omsg;
  if( !BFGS_sTy_suff_p_d(*s,*y,&sTy,&omsg,"MatrixSymPosDefLBFGS::secant_update(...)" ) ) {
    throw UpdateSkippedException( omsg.str() );	
  }

  try {

  // Update counters
  if( k_bar_ == m_ ) {
//		// We are at the end storage so loop back around again
//		k_bar_ = 1;
    // We are at the end of the storage so remove the oldest stored update
    // and move updates to make room for the new update.  This has to be done for the
    // the matrix to behave properly
    {for( size_type k = 1; k <= m_-1; ++k ) {
      S_.col(k) = S_.col(k+1);                              // Shift S.col() to the left
      Y_.col(k) = Y_.col(k+1);                              // Shift Y.col() to the left
      STY_.col(k)(1,m_-1) = STY_.col(k+1)(2,m_);            // Move submatrix STY(2,m-1,2,m-1) up and left
      STSYTY_.col(k)(k+1,m_) = STSYTY_.col(k+1)(k+2,m_+1);  // Move triangular submatrix STS(2,m-1,2,m-1) up and left
      STSYTY_.col(k+1)(1,k) = STSYTY_.col(k+2)(2,k+1);      // Move triangular submatrix YTY(2,m-1,2,m-1) up and left
    }}
  }
  else {
    k_bar_++;
  }
  if( m_bar_ < m_ ) {
    // This is the first few updates where we have not maxed out the storage.
    m_bar_++;
  }

  // Set the update vectors
  S_.col(k_bar_)(1,n_) = *s;
  Y_.col(k_bar_)(1,n_) = *y;

  // /////////////////////////////////////////////////////////////////////////////////////
  // Update S'Y
  //
  // Update the row and column k_bar
  //
  //	S'Y = 
  //
  //	[	s(1)'*y(1)		...		s(1)'*y(k_bar)		...		s(1)'*y(m_bar)		]
  //	[	.						.							.					] row
  //	[	s(k_bar)'*y(1)	...		s(k_bar)'*y(k_bar)	...		s(k_bar)'*y(m_bar)	] k_bar
  //	[	.						.							.					]
  //	[	s(m_bar)'*y(1)	...		s(m_bar)'*y(k_bar)	...		s(m_bar)'*y(m_bar)	]
  //
  //								col k_bar
  //
  // Therefore we set:
  //	(S'Y)(:,k_bar) =  S'*y(k_bar)
  //	(S'Y)(k_bar,:) =  s(k_bar)'*Y

  const DMatrixSlice
    &S = this->S(),
    &Y = this->Y();

  //	(S'Y)(:,k_bar) =  S'*y(k_bar)
  V_MtV( &STY_.col(k_bar_)(1,m_bar_), S, BLAS_Cpp::trans, Y.col(k_bar_) );

  //	(S'Y)(k_bar,:)' =  Y'*s(k_bar)
  V_MtV( &STY_.row(k_bar_)(1,m_bar_), Y, BLAS_Cpp::trans, S.col(k_bar_) );

  // /////////////////////////////////////////////////////////////////
  // Update S'S
  //
  //	S'S = 
  //
  //	[	s(1)'*s(1)		...		symmetric					symmetric			]
  //	[	.						.							.					] row
  //	[	s(k_bar)'*s(1)	...		s(k_bar)'*s(k_bar)	...		symmetric			] k_bar
  //	[	.						.							.					]
  //	[	s(m_bar)'*s(1)	...		s(m_bar)'*s(k_bar)	...		s(m_bar)'*s(m_bar)	]
  //
  //								col k_bar
  //
  // Here we will update the lower triangular part of S'S.  To do this we
  // only need to compute:
  //		t = S'*s(k_bar) = { s(k_bar)' * [ s(1),..,s(k_bar),..,s(m_bar) ]  }'
  // then set the appropriate rows and columns of S'S.

  if( work_.size() < m_ )
    work_.resize(m_);

  // work = S'*s(k_bar)
  V_MtV( &work_(1,m_bar_), S, BLAS_Cpp::trans, S.col(k_bar_) );

  // Set row elements
  STSYTY_.row(k_bar_+1)(1,k_bar_) = work_(1,k_bar_);
  // Set column elements
  STSYTY_.col(k_bar_)(k_bar_+1,m_bar_+1) = work_(k_bar_,m_bar_);

  // /////////////////////////////////////////////////////////////////////////////////////
  // Update Y'Y
  //
  // Update the row and column k_bar
  //
  //	Y'Y = 
  //
  //	[	y(1)'*y(1)		...		y(1)'*y(k_bar)		...		y(1)'*y(m_bar)		]
  //	[	.						.							.					] row
  //	[	symmetric		...		y(k_bar)'*y(k_bar)	...		y(k_bar)'*y(m_bar)	] k_bar
  //	[	.						.							.					]
  //	[	symmetric		...		symmetric			...		y(m_bar)'*y(m_bar)	]
  //
  //								col k_bar
  //
  // Here we will update the upper triangular part of Y'Y.  To do this we
  // only need to compute:
  //		t = Y'*y(k_bar) = { y(k_bar)' * [ y(1),..,y(k_bar),..,y(m_bar) ]  }'
  // then set the appropriate rows and columns of Y'Y.

  // work = Y'*y(k_bar)
  V_MtV( &work_(1,m_bar_), Y, BLAS_Cpp::trans, Y.col(k_bar_) );

  // Set row elements
  STSYTY_.col(k_bar_+1)(1,k_bar_) = work_(1,k_bar_);
  // Set column elements
  STSYTY_.row(k_bar_)(k_bar_+1,m_bar_+1) = work_(k_bar_,m_bar_);

  // /////////////////////////////
  // Update gamma_k

  // gamma_k = s'*y / y'*y
  if(auto_rescaling_)
    gamma_k_ = STY_(k_bar_,k_bar_) / STSYTY_(k_bar_,k_bar_+1);

  // We do not initially update Q unless we have to form a matrix-vector
  // product later.
  
  Q_updated_ = false;
  num_secant_updates_++;

  }	//	end try
  catch(...) {
    // If we throw any exception the we should make the matrix uninitialized
    // so that we do not leave this object in an inconsistant state.
    n_ = 0;
    throw;
  }
}

// Overridden from MatrixSymAddDelUpdateble

void MatrixSymPosDefLBFGS::initialize(
  value_type         alpha
  ,size_type         max_size
  )
{
  // Validate input
  if( alpha <= 0.0 ) {
    std::ostringstream omsg;
    omsg
      << "MatrixSymPosDefLBFGS::initialize(alpha,max_size) : Error, "
      << "alpha = " << alpha << " <= 0 is not allowed!";
    throw std::invalid_argument( omsg.str() );
  }
  n_max_ = max_size;
  this->init_identity(1,alpha);
}

void MatrixSymPosDefLBFGS::initialize(
  const DMatrixSliceSym      &A
  ,size_type         max_size
  ,bool              force_factorization
  ,Inertia           inertia
  ,PivotTolerances   pivot_tols
  )
{
  throw std::runtime_error(
    "MatrixSymPosDefLBFGS::initialize(A,max_size,force_refactorization,inertia) : Error, "
    "This function is undefined for this subclass.  I am so sorry for this terrible hack!" );
}

size_type MatrixSymPosDefLBFGS::max_size() const
{
  return n_max_;
}

MatrixSymAddDelUpdateable::Inertia MatrixSymPosDefLBFGS::inertia() const
{
  return Inertia(0,0,n_);
}

void MatrixSymPosDefLBFGS::set_uninitialized()
{
  n_ = 0;
}

void MatrixSymPosDefLBFGS::augment_update(
  const DVectorSlice  *t
  ,value_type        alpha
  ,bool              force_refactorization
  ,EEigenValType     add_eigen_val
  ,PivotTolerances   pivot_tols
  )
{
  assert_initialized();
  if( n_ == n_max_ ) {
    std::ostringstream omsg;
    omsg
      << "MatrixSymPosDefLBFGS::augment_update(...) : Error, "
      << "this->rows() = " << n_ << " == this->max_size() = " << n_max_
      << " so we can't allow the matrix to grow!";
    throw std::invalid_argument( omsg.str() );
  }
  if( t ) {
    throw std::invalid_argument(		
      "MatrixSymPosDefLBFGS::augment_update(...) : Error, "
      "t must be NULL in this implemention.  Sorry for this hack" );
  }
  if( alpha <= 0.0 ) {
    std::ostringstream omsg;
    omsg
      << "MatrixSymPosDefLBFGS::augment_update(...) : Error, "
      << "alpha = " << alpha << " <= 0 is not allowed!";
    throw std::invalid_argument( omsg.str() );
  }
  if( add_eigen_val == MatrixSymAddDelUpdateable::EIGEN_VAL_NEG ) {
    std::ostringstream omsg;
    omsg
      << "MatrixSymPosDefLBFGS::augment_update(...) : Error, "
      << "add_eigen_val == EIGEN_VAL_NEG is not allowed!";
    throw std::invalid_argument( omsg.str() );
  }
  //
  // Here we will do the simplest thing possible.  We will just  set:
  //
  // [ S ] -> S       [ Y ] -> Y
  // [ 0 ]            [ 0 ]
  //
  // and let the new matrix be:
  //
  // [ B      0     ] -> B
  // [ 0  1/gamma_k ]
  //
  // Nothing else, not even Q, needs to be updated!
  //
  S_.row(n_+1)(1,m_bar_) = 0.0;
  Y_.row(n_+1)(1,m_bar_) = 0.0;
  ++n_;
}

void MatrixSymPosDefLBFGS::delete_update(
  size_type          jd
  ,bool              force_refactorization
  ,EEigenValType     drop_eigen_val
  ,PivotTolerances   pivot_tols
  )
{
  assert_initialized();
  //
  // Removing a symmetric row and column jd is the same a removing row
  // S(jd,:) from S and row Y(jd,:) from Y.  At the same time we must
  // update S'*Y, S'*S and Y'*Y.  To see how to update these matrices
  // not that we can represent each column of S and Y as:
  //
  //           [ S(1:jd-1,k) ]                 [ Y(1:jd-1,k) ]
  // S(:,k) =  [ S(jd,k)     ]     , Y(:,k) =  [ Y(jd,k)     ]  , k = 1...m_bar
  //           [ S(jd+1:n,k) ]                 [ Y(jd+1:n,k) ]
  //
  // Using the above, we can write:
  //
  // (S'*Y)(p,q) = S(1:jd-1,p)'*Y(1:jd-1,q) + S(jd,p)*Y(jd,q) + S(jd+1:n,p)'*Y(jd+1:n,q)
  //     , for p = 1...m_bar, q = 1...m_bar
  //
  // We see that the new S'*Y and the old differ by only the term S(jd,p)*Y(jd,q).  Therefore, we
  // only need to subtract off this term in each of the elements in order to update S'*Y for the
  // deletion of this element jd.  To see how to do this with BLAS, first consider subtracting
  // of the terms by column as:
  //
  // (S'*Y)(:,q) <- (S'*Y)(:,q) - S(jd,:)'*Y(jd,q)
  //     , for q = 1...m_bar
  // 
  // Then, if we put all of the columns together we get:
  //
  // (S'*Y)(:,:) <- (S'*Y)(:,:) - S(jd,:)'*Y(jd,:)
  // =>
  // (S'*Y) <- (S'*Y) - S.row(jd)*Y.row(jd)'
  //
  // In otherwords the above update operation is just an unsymmetric rank-1 update
  //
  // Similar updates for S'*S and Y'*Y are derived by just substituting matrices
  // in to the above update for S'*Y:
  //
  // (S'*S) <- (S'*S) - S.row(jd)*S.row(jd)'
  // 
  // (Y'*Y) <- (Y'*Y) - Y.row(jd)*Y.row(jd)'
  //
  // These updates are symmetric rank-1 updates.
  //
  DMatrixSlice S = this->S();
  DMatrixSlice Y = this->Y();
  DMatrixSlice STY = this->STY();
  DMatrixSliceSym        STS = this->STS();
  DMatrixSliceSym        YTY = this->YTY();
  // (S'*Y) <- (S'*Y) - S.row(jd)*Y.row(jd)'
  DenseLinAlgPack::ger( -1.0, S.row(jd), Y.row(jd), &STY );
  // (S'*S) <- (S'*S) - S.row(jd)*S.row(jd)'
  DenseLinAlgPack::syr( -1.0, S.row(jd), &STS );
  // (Y'*Y) <- (Y'*Y) - Y.row(jd)*Y.row(jd)'
  DenseLinAlgPack::syr( -1.0, Y.row(jd), &YTY );
  // Remove row jd from S and Y one column at a time
  // (one element at a time!)
  if( jd < n_ ) {
    {for( size_type k = 1; k <= m_bar_; ++k ) {
      value_type *ptr = S.col_ptr(k);
      std::copy( ptr + jd, ptr + n_, ptr + jd - 1 );
    }}
    {for( size_type k = 1; k <= m_bar_; ++k ) {
      value_type *ptr = Y.col_ptr(k);
      std::copy( ptr + jd, ptr + n_, ptr + jd - 1 );
    }}
  }
  // Update the size
  --n_;
  Q_updated_ = false;
}

// Private member functions

void MatrixSymPosDefLBFGS::Vp_DtV( DVectorSlice* y, const DVectorSlice& x ) const
{
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->size(), m_bar_, m_bar_
    , BLAS_Cpp::no_trans, x.size() );

  DVectorSlice::const_iterator
    d_itr	= STY_.diag(0).begin(),
    x_itr	= x.begin();
  DVectorSlice::iterator
    y_itr	= y->begin();

  while( y_itr != y->end() )
    *y_itr++ += (*d_itr++) * (*x_itr++);		
}

//
// We need to update the factorizations to solve for:
//
// x = inv(Q) * y   =>   Q * x = y
//
//	[ (1/gk)*S'S	 L	] * [ x1 ] = [ y1 ]
//	[      L'		-D	]   [ x2 ]   [ y2 ]
//
// We will solve the above system using a schur complement:
//
// C = (1/gk)*S'S + L*inv(D)*L'
//
// According to the referenced paper, C is p.d. so:
//
// C = J*J'
//
// We then compute the solution as:
//
// x1 = inv(C) * ( y1 + L*inv(D)*y2 )
// x2 = - inv(D) * ( y2 - L'*x1 )
//
// Therefore we will just update the factorization C = J*J'
// where the factor J is stored in QJ_.
//

void MatrixSymPosDefLBFGS::update_Q() const
{
  using DenseLinAlgPack::tri;
  using DenseLinAlgPack::tri_ele;
  using DenseLinAlgPack::Mp_StM;

  //
  // We need update the factorizations to solve for:
  //
  // x = inv(Q) * y
  //
  //	[ y1 ]	=	[ (1/gk)*S'S	 L	] * [ x1 ]
  //	[ y2 ]		[      L'		-D	]   [ x2 ]
  //
  // We will solve the above system using the schur complement:
  //
  // C = (1/gk)*S'S + L*inv(D)*L'
  //
  // According to the referenced paper, C is p.d. so:
  //
  // C = J*J'
  //
  // We then compute the solution as:
  //
  // x1 = inv(C) * ( y1 + L*inv(D)*y2 )
  // x2 = - inv(D) * ( y2 - L'*x1 )
  //
  // Therefore we will just update the factorization C = J*J'
  //

  // Form the upper triangular part of C which will become J
  // which we are using storage of QJ

  if( QJ_.rows() < m_ )
    QJ_.resize( m_, m_ );

  const size_type
    mb = m_bar_;

  DMatrixSlice
    C = QJ_(1,mb,1,mb);

  // C = L * inv(D) * L'
  //
  // Here L is a strictly lower triangular (zero diagonal) matrix where:
  //
  // L = [ 0  0 ]
  //     [ Lb 0 ]
  //
  // Lb is lower triangular (nonzero diagonal)
  //
  // Therefore we compute L*inv(D)*L' as:
  //
  // C = [ 0	0 ] * [ Db  0  ] * [ 0  Lb' ]
  //	   [ Lb 0 ]   [ 0   db ]   [ 0   0  ]
  //
  //   = [ 0  0  ] = [ 0      0     ]
  //     [ 0  Cb ]   [ 0  Lb*Db*Lb' ]
  //
  // We need to compute the upper triangular part of Cb.

  C.row(1) = 0.0;
  if( mb > 1 )
    comp_Cb( STY_(2,mb,1,mb-1), STY_.diag(0)(1,mb-1), &C(2,mb,2,mb) );

  // C += (1/gk)*S'S

  const DMatrixSliceSym &STS = this->STS();
  Mp_StM( &C, (1/gamma_k_), tri( STS.gms(), STS.uplo(), BLAS_Cpp::nonunit )
    , BLAS_Cpp::trans );

  // Now perform a cholesky factorization of C
  // After this factorization the upper triangular part of QJ
  // (through C) will contain the cholesky factor.

  DMatrixSliceTriEle C_upper = tri_ele( C, BLAS_Cpp::upper );
  DenseLinAlgLAPack::potrf( &C_upper );

  Q_updated_ = true;
}

void MatrixSymPosDefLBFGS::V_invQtV( DVectorSlice* x, const DVectorSlice& y ) const
{
  using DenseLinAlgPack::sym;
  using DenseLinAlgPack::tri;
  using DenseLinAlgPack::Vp_StV;
  using DenseLinAlgPack::V_InvMtV;

  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_MtV;


  // Solve the system 
  //
  // Q * x = y
  //
  // Using the schur complement factorization as described above.

  const size_type
    mb = m_bar_;

  if(!Q_updated_) {
    update_Q();
  }

  DVectorSlice
    x1 = (*x)(1,mb),
    x2 = (*x)(mb+1,2*mb);

  const DVectorSlice
    y1 = y(1,mb),
    y2 = y(mb+1,2*mb);

  // //////////////////////////////////////
  // x1 = inv(C) * ( y1 + L*inv(D)*y2 )
  //     = inv(J'*J) * r
  //     = inv(J) * inv(J') * r

  {	// x1 = inv(D) * y2
    DVectorSlice::const_iterator
      d_itr = STY_.diag(0).begin(),
      y2_itr = y2.begin();
    DVectorSlice::iterator
      x1_itr = x1.begin();
    while( x1_itr != x1.end() )
      *x1_itr++ = *y2_itr++ / *d_itr++;
  }

  // x1 = L * x1
  //
  //    = [ 0  0 ] * [ x1(1:mb-1) ]
  //      [ Lb 0 ]   [ x1(mb)     ]
  //
  //    = [ 0             ]
  //      [ Lb*x1(1:mb-1) ]
  //
  if( mb > 1 ) {
    // x1(2,mb) = x1(1,mb-1) ( copy from mb-1 to mb then mb-2 to mb-1
    // etc. so that we don't overwrite the elements we need to copy ).
    DVectorSlice
      x1a = x1(1,mb-1),
      x1b = x1(2,mb);
    std::copy( x1a.rbegin(), x1a.rend(), x1b.rbegin() );
    V_MtV( &x1b, Lb(), BLAS_Cpp::no_trans, x1b );
  }
  x1(1) = 0.0;

  // r = x1 += y1
  Vp_V( &x1, y1 );

  // x1 = inv(J') * r
  const DMatrixSliceTri J = tri( QJ_(1,mb,1,mb), BLAS_Cpp::upper, BLAS_Cpp::nonunit );
  V_InvMtV( &x1, J, BLAS_Cpp::trans, x1 );

  // x1 = inv(J) * x1
  V_InvMtV( &x1, J, BLAS_Cpp::no_trans, x1 );

  // /////////////////////////////////////
  // x2 = inv(D) * ( - y2 + L'*x1 )

  // x2 = L'*x1
  //
  //    = [ 0  Lb' ] * [ x1(1)    ]
  //      [ 0  0   ]   [ x1(2,mb) ]
  //
  //    = [ Lb'*x1(2,mb) ]
  //      [      0       ]
  //
  if( mb > 1 ) {
    V_MtV( &x2(1,mb-1), Lb(), BLAS_Cpp::trans, x1(2,mb) );
  }
  x2(mb) = 0.0;

  // x2 += -y2
  Vp_StV( &x2, -1.0, y2 );

  // x2 = inv(D) * x2
  {
    DVectorSlice::const_iterator
      d_itr = STY_.diag(0).begin();
    DVectorSlice::iterator
      x2_itr = x2.begin();
    while( x2_itr != x2.end() )
      *x2_itr++ /= *d_itr++;
  }
}

void MatrixSymPosDefLBFGS::assert_initialized() const
{
  if(!n_)
    throw std::logic_error( "MatrixSymPosDefLBFGS::assert_initialized() : "
      "Error, matrix not initialized" );
}

}	// end namespace ConstrainedOptPack 

namespace {

void comp_Cb( const DMatrixSlice& Lb, const DVectorSlice& Db_diag
  , DMatrixSlice* Cb )
{
  // Lb * inv(Db) * Lb =
  //
  // [ l(1,1)						]   [ dd(1)					]   [ l(1,1)	l(2,1)	...	l(p,1)	]
  // [ l(2,1)	l(2,2)				]   [		dd(2)			]   [			l(2,2)	...	l(p,2)	]
  // [ .		.        .			] * [			.			] * [					.	.		]
  // [ l(p,1)	l(p,2)	...	l(p,p)	]   [				dd(p)	]   [						l(p,p)	]
  //
  //
  // [ l(1,1)*dd(1)*l(1,1)	l(1,1)*dd(1)*l(2,1)							...		l(1,1)*dd(1)*l(p,1)							]
  // [ symmetric				l(2,1)*dd(1)*l(2,1) + l(2,2)*dd(2)*l(2,2)	...		l(2,1)*dd(1)*l(p,1) + l(2,2)*dd(2)*l(p,2)	]		
  // [ .						.											...
  // [ symmetric				symmetric									...		sum( l(p,i)*dd(i)*l(p,i), i=1,..,p )		]
  //
  // Therefore we can express the upper triangular elemetns of Cb as:
  //
  // Cb(i,j) = sum( l(i,k)*dd(k)*l(j,k), k = 1,..,i )

  typedef DenseLinAlgPack::size_type size_type;
  typedef DenseLinAlgPack::value_type value_type;

  TEUCHOS_TEST_FOR_EXCEPT( !(  Lb.rows() == Cb->rows() && Cb->rows() == Db_diag.size()  ) );

  const size_type p = Db_diag.size();

  for( size_type i = 1; i <= p; ++i ) {
    for( size_type j = i; j <= p; ++j ) {
      value_type &c = (*Cb)(i,j) = 0.0;
      for( size_type k = 1; k <= i; ++k ) {
        c += Lb(i,k) * Lb(j,k) / Db_diag(k);
      }
    }
  }

  // ToDo: Make the above operation more efficent if needed!
}

}	// end namespace

#endif // 0
