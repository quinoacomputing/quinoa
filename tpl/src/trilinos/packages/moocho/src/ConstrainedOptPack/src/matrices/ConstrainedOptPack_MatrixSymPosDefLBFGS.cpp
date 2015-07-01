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

#include <typeinfo>

#include "ConstrainedOptPack_MatrixSymPosDefLBFGS.hpp"
#include "AbstractLinAlgPack_BFGS_helpers.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgLAPack.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"

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
    size_type   m
  ,bool       maintain_original
  ,bool       maintain_inverse
  ,bool       auto_rescaling
  )
{
  initial_setup(m,maintain_original,maintain_inverse,auto_rescaling);
}

void MatrixSymPosDefLBFGS::initial_setup(
  size_type m,
  bool maintain_original,
  bool maintain_inverse,
  bool auto_rescaling
  )
{
  // Validate input
  TEUCHOS_TEST_FOR_EXCEPTION(
    !maintain_original && !maintain_inverse, std::invalid_argument
    ,"MatrixSymPosDefLBFGS::initial_setup(...) : "
    "Error, both maintain_original and maintain_inverse can not both be false!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    m < 1, std::invalid_argument
    ,"MatrixSymPosDefLBFGS::set_num_updates_stored(m) : "
    "Error, the number of storage locations must be > 0" );
  vec_spc_           = Teuchos::null;
  maintain_original_ = maintain_original;
  maintain_inverse_  = maintain_inverse;
  m_                 = m;
  n_                 = 0; // make uninitialized
  num_secant_updates_= 0;
  auto_rescaling_    = auto_rescaling;
}

// Overridden from MatrixOp

const VectorSpace& MatrixSymPosDefLBFGS::space_cols() const
{
  return *vec_spc_;
}

std::ostream& MatrixSymPosDefLBFGS::output(std::ostream& out) const
{
  assert_initialized();
  out << "\n*** Limited Memory BFGS matrix\n";
  out << "\nConversion to dense =\n";
  MatrixOp::output(out);
  out << "\nStored quantities\n"
    << "\nn       = " << n_
    << "\nm       = " << m_
    << "\nm_bar   = " << m_bar_
    << "\ngamma_k = " << gamma_k_ << std::endl;
  if( m_bar_ ) {
    out	<< "\nS =\n" << *S()
      << "\nY =\n" << *Y()
      << "\nS'Y =\n" << STY_(1,m_bar_,1,m_bar_)
      << "\nlower(S'S) \\ zero diagonal \\ upper(Y'Y) =\n"
        << STSYTY_(1,m_bar_+1,1,m_bar_+1)
      << "\nQ updated? = " << Q_updated_ << std::endl;
    if(Q_updated_)
      out << "\nCholesky of schur complement of Q, QJ =\n" << QJ_(1,m_bar_,1,m_bar_);
  }
  return out;
}

MatrixOp& MatrixSymPosDefLBFGS::operator=(const MatrixOp& mwo)
{	
  const MatrixSymPosDefLBFGS *p_m = dynamic_cast<const MatrixSymPosDefLBFGS*>(&mwo);
  if(p_m) {
    if( p_m == this ) return *this;	// assignment to self
    // Important: Assign all members here.
    auto_rescaling_      = p_m->auto_rescaling_;
    maintain_original_   = p_m->maintain_original_;
    original_is_updated_ = p_m->original_is_updated_;
    maintain_inverse_    = p_m->maintain_inverse_;
    inverse_is_updated_  = p_m->inverse_is_updated_;
    vec_spc_             = p_m->vec_spc_.get() ? p_m->vec_spc_->clone() : Teuchos::null;
    n_	 		         = p_m->n_;
    m_			         = p_m->m_;
    m_bar_		         = p_m->m_bar_;
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      true,std::invalid_argument
      ,"MatrixSymPosDefLBFGS::operator=(const MatrixOp& mwo) : Error, "
      "The concrete type of mwo \'" << typeName(mwo) << "\' is not "
      "as subclass of MatrixSymPosDefLBFGS as required" );
  }
  return *this;
}

// Level-2 BLAS

void MatrixSymPosDefLBFGS::Vp_StMtV(
    VectorMutable* y, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const Vector& x, value_type beta
  ) const
{
  using AbstractLinAlgPack::Vt_S;
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::Vp_StMtV;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::V_MtV;
  typedef VectorDenseEncap         vde;
  typedef VectorDenseMutableEncap  vdme;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

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
  Vt_S( y, beta );

  // y += (1/gk)*x
  Vp_StV( y, invgk, x );

  if( !m_bar_ )
    return;	// No updates have been added yet.

  const multi_vec_ptr_t
    S = this->S(),
    Y = this->Y();

  // Get workspace

  const size_type
    mb = m_bar_;

  Workspace<value_type>  t1_ws(wss,2*mb);
  DVectorSlice                 t1(&t1_ws[0],t1_ws.size());
  Workspace<value_type>  t2_ws(wss,2*mb);
  DVectorSlice                 t2(&t2_ws[0],t2_ws.size());

  VectorSpace::vec_mut_ptr_t
    t = S->space_rows().create_member();

  // t1 = [ (1/gk)*S'*x ]
  //		[      Y'*x   ]

  V_StMtV( t.get(), invgk, *S, BLAS_Cpp::trans, x );
  t1(1,mb) = vde(*t)();
  V_MtV( t.get(), *Y, BLAS_Cpp::trans, x );
  t1(mb+1,2*mb) = vde(*t)();

  // t2 =	inv(Q) * t1
  V_invQtV( &t2, t1 );

  // y += -(1/gk) * S * t2(1:m)
  (vdme(*t)() = t2(1,mb));
  Vp_StMtV( y, -invgk, *S, BLAS_Cpp::no_trans, *t );

  // y += -1.0 * Y * t2(m+1,2m
  (vdme(*t)() = t2(mb+1,2*mb));
  Vp_StMtV( y, -1.0, *Y, BLAS_Cpp::no_trans, *t );

}

// Overridden from MatrixOpNonsing

void MatrixSymPosDefLBFGS::V_InvMtV(
  VectorMutable* y, BLAS_Cpp::Transp trans_rhs1
  , const Vector& x
  ) const
{
  using AbstractLinAlgPack::Vp_StMtV;
  using DenseLinAlgPack::V_InvMtV;
  using LinAlgOpPack::V_mV;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::Vp_MtV;
  using DenseLinAlgPack::Vp_StMtV;
  typedef VectorDenseEncap         vde;
  typedef VectorDenseMutableEncap  vdme;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

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

  const multi_vec_ptr_t
    S = this->S(),
    Y = this->Y();

  // Get workspace

  Workspace<value_type>    t1_ws(wss,2*mb);
  DVectorSlice                   t1(&t1_ws[0],t1_ws.size());
  Workspace<value_type>    t2_ws(wss,mb);
  DVectorSlice                   t2(&t2_ws[0],t2_ws.size());
  Workspace<value_type>    t3_ws(wss,mb);
  DVectorSlice                   t3(&t3_ws[0],t3_ws.size());
  Workspace<value_type>    t4_ws(wss,mb);
  DVectorSlice                   t4(&t4_ws[0],t4_ws.size());
  Workspace<value_type>    t5_ws(wss,mb);
  DVectorSlice                   t5(&t5_ws[0],t5_ws.size());

  VectorSpace::vec_mut_ptr_t
    t = S->space_rows().create_member();

  const DMatrixSliceTri
    &R = this->R();

  const DMatrixSliceSym
    &YTY = this->YTY();

  // t1 = [   S'*x  ]
  //      [ gk*Y'*x ]
  V_MtV( t.get(), *S, BLAS_Cpp::trans, x );
  t1(1,mb) = vde(*t)();
  V_StMtV( t.get(), gamma_k_, *Y, BLAS_Cpp::trans, x );
  t1(mb+1,2*mb) = vde(*t)();

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
  (vdme(*t)() = t5);
  Vp_MtV( y, *S, BLAS_Cpp::no_trans, *t );

  // y += -gk*Y*t2
  (vdme(*t)() = t2);
  Vp_StMtV( y, -gamma_k_, *Y, BLAS_Cpp::no_trans, *t );

}

// Overridden from MatrixSymSecant

void MatrixSymPosDefLBFGS::init_identity( const VectorSpace& space_diag, value_type alpha )
{
  // Validate input
  TEUCHOS_TEST_FOR_EXCEPTION(
    alpha <= 0.0, std::invalid_argument
    ,"MatrixSymPosDefLBFGS::init_identity(n,alpha) : Error, "
    "alpha = " << alpha << " <= 0 is not allowed!" );
  
  // Set the vector space
  vec_spc_ = space_diag.clone();
  vec_spc_.get();

  // Set storage
  S_ = vec_spc_->create_members(m_);
  Y_ = vec_spc_->create_members(m_);
  TEUCHOS_TEST_FOR_EXCEPT( !( S_.get() ) );
  TEUCHOS_TEST_FOR_EXCEPT( !( Y_.get() ) );
  STY_.resize( m_, m_ );
  STSYTY_.resize( m_+1, m_+1 );
  STSYTY_.diag(0) = 0.0;

  gamma_k_ = 1.0/alpha;

  // Initialize counters
  m_bar_	= 0;

  n_ = vec_spc_->dim();        // initialized;
  original_is_updated_ = true; // This will never change for now
  inverse_is_updated_  = true; // This will never change for now
  num_secant_updates_  = 0;    // reset this to zero
}

void MatrixSymPosDefLBFGS::init_diagonal(const Vector& diag)
{
  init_identity( diag.space(), diag.norm_inf() );
}

void MatrixSymPosDefLBFGS::secant_update(
  VectorMutable* s, VectorMutable* y, VectorMutable* Bs
  )
{
  using AbstractLinAlgPack::BFGS_sTy_suff_p_d;
  using AbstractLinAlgPack::dot;
  using LinAlgOpPack::V_MtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

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
  if( m_bar_ == m_ ) {
    // We are at the end of the storage so remove the oldest stored update
    // and move updates to make room for the new update.  This has to be done for the
    // the matrix to behave properly
    {for( size_type k = 1; k <= m_-1; ++k ) {
      S_->col(k) = S_->col(k+1);                            // Shift S.col() to the left
      Y_->col(k) = Y_->col(k+1);                            // Shift Y.col() to the left
      STY_.col(k)(1,m_-1) = STY_.col(k+1)(2,m_);            // Move submatrix STY(2,m-1,2,m-1) up and left
      STSYTY_.col(k)(k+1,m_) = STSYTY_.col(k+1)(k+2,m_+1);  // Move triangular submatrix STS(2,m-1,2,m-1) up and left
      STSYTY_.col(k+1)(1,k) = STSYTY_.col(k+2)(2,k+1);      // Move triangular submatrix YTY(2,m-1,2,m-1) up and left
    }}
    // ToDo: Create an abstract interface, call it MultiVectorShiftVecs, to rearrange S and Y all at once.
    // This will be important for parallel performance.
  }
  else {
    m_bar_++;
  }
  // Set the update vectors
  *S_->col(m_bar_) = *s;
  *Y_->col(m_bar_) = *y;

  // /////////////////////////////////////////////////////////////////////////////////////
  // Update S'Y
  //
  // Update the row and column m_bar
  //
  //	S'Y = 
  //
  //	[	s(1)'*y(1)		...		s(1)'*y(m_bar)		...		s(1)'*y(m_bar)		]
  //	[	.						.							.					] row
  //	[	s(m_bar)'*y(1)	...		s(m_bar)'*y(m_bar)	...		s(m_bar)'*y(m_bar)	] m_bar
  //	[	.						.							.					]
  //	[	s(m_bar)'*y(1)	...		s(m_bar)'*y(m_bar)	...		s(m_bar)'*y(m_bar)	]
  //
  //								col m_bar
  //
  // Therefore we set:
  //	(S'Y)(:,m_bar) =  S'*y(m_bar)
  //	(S'Y)(m_bar,:) =  s(m_bar)'*Y

  const multi_vec_ptr_t
    S = this->S(),
    Y = this->Y();

  VectorSpace::vec_mut_ptr_t
    t = S->space_rows().create_member();  // temporary workspace

  //	(S'Y)(:,m_bar) =  S'*y(m_bar)
  V_MtV( t.get(), *S, BLAS_Cpp::trans, *y );
  STY_.col(m_bar_)(1,m_bar_) = VectorDenseEncap(*t)();

  //	(S'Y)(m_bar,:)' =  Y'*s(m_bar)
  V_MtV( t.get(), *Y, BLAS_Cpp::trans, *s );
  STY_.row(m_bar_)(1,m_bar_) = VectorDenseEncap(*t)();

  // /////////////////////////////////////////////////////////////////
  // Update S'S
  //
  //	S'S = 
  //
  //	[	s(1)'*s(1)		...		symmetric					symmetric			]
  //	[	.						.							.					] row
  //	[	s(m_bar)'*s(1)	...		s(m_bar)'*s(m_bar)	...		symmetric			] m_bar
  //	[	.						.							.					]
  //	[	s(m_bar)'*s(1)	...		s(m_bar)'*s(m_bar)	...		s(m_bar)'*s(m_bar)	]
  //
  //								col m_bar
  //
  // Here we will update the lower triangular part of S'S.  To do this we
  // only need to compute:
  //		t = S'*s(m_bar) = { s(m_bar)' * [ s(1),..,s(m_bar),..,s(m_bar) ]  }'
  // then set the appropriate rows and columns of S'S.

  Workspace<value_type>   work_ws(wss,m_bar_);
  DVectorSlice                  work(&work_ws[0],work_ws.size());

  // work = S'*s(m_bar)
  V_MtV( t.get(), *S, BLAS_Cpp::trans, *s );
  work = VectorDenseEncap(*t)();

  // Set row elements
  STSYTY_.row(m_bar_+1)(1,m_bar_) = work;
  // Set column elements
  STSYTY_.col(m_bar_)(m_bar_+1,m_bar_+1) = work(m_bar_,m_bar_);

  // /////////////////////////////////////////////////////////////////////////////////////
  // Update Y'Y
  //
  // Update the row and column m_bar
  //
  //	Y'Y = 
  //
  //	[	y(1)'*y(1)		...		y(1)'*y(m_bar)		...		y(1)'*y(m_bar)		]
  //	[	.						.							.					] row
  //	[	symmetric		...		y(m_bar)'*y(m_bar)	...		y(m_bar)'*y(m_bar)	] m_bar
  //	[	.						.							.					]
  //	[	symmetric		...		symmetric			...		y(m_bar)'*y(m_bar)	]
  //
  //								col m_bar
  //
  // Here we will update the upper triangular part of Y'Y.  To do this we
  // only need to compute:
  //		t = Y'*y(m_bar) = { y(m_bar)' * [ y(1),..,y(m_bar),..,y(m_bar) ]  }'
  // then set the appropriate rows and columns of Y'Y.

  // work = Y'*y(m_bar)
  V_MtV( t.get(), *Y, BLAS_Cpp::trans, *y );
  work = VectorDenseEncap(*t)();

  // Set row elements
  STSYTY_.col(m_bar_+1)(1,m_bar_) = work;
  // Set column elements
  STSYTY_.row(m_bar_)(m_bar_+1,m_bar_+1) = work(m_bar_,m_bar_);

  // /////////////////////////////
  // Update gamma_k

  // gamma_k = s'*y / y'*y
  if(auto_rescaling_)
    gamma_k_ = STY_(m_bar_,m_bar_) / STSYTY_(m_bar_,m_bar_+1);

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

// Private member functions

void MatrixSymPosDefLBFGS::Vp_DtV( DVectorSlice* y, const DVectorSlice& x ) const
{
  DenseLinAlgPack::Vp_MtV_assert_sizes(
    y->dim(), m_bar_, m_bar_, BLAS_Cpp::no_trans, x.dim() );

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
  try {
    DenseLinAlgLAPack::potrf( &C_upper );
  }
  catch( const DenseLinAlgLAPack::FactorizationException &fe ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, UpdateFailedException
      ,"Error, the factorization of Q which should be s.p.d. failed with"
      " the error message: {" << fe.what() << "}";
      );
  }

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

void comp_Cb(
  const DMatrixSlice& Lb, const DVectorSlice& Db_diag
  ,DMatrixSlice* Cb
  )
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

  TEUCHOS_TEST_FOR_EXCEPT( !(  Lb.rows() == Cb->rows() && Cb->rows() == Db_diag.dim()  ) ); // only a local error!

  const size_type p = Db_diag.dim();

  for( size_type i = 1; i <= p; ++i ) {
    for( size_type j = i; j <= p; ++j ) {
      value_type &c = (*Cb)(i,j) = 0.0;
      for( size_type k = 1; k <= i; ++k ) {
        c += Lb(i,k) * Lb(j,k) / Db_diag(k);
      }
    }
  }

  // ToDo: Make the above operation more efficent if needed! (i.e. write
  // it in fortran or something?).
}

}	// end namespace
