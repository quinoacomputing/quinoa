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


#include "Moocho_ConfigDefs.hpp"


#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPKWIK


#include <assert.h>

#include <vector>

#include "ConstrainedOptPack_QPSolverRelaxedQPKWIK.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp"
#include "AbstractLinAlgPack_SortByDescendingAbsValue.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_sparse_bounds.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

namespace QPKWIKNEW_CppDecl {

// Declarations that will link to the fortran object file.
// These may change for different platforms

using FortranTypes::f_int;			// INTEGER
using FortranTypes::f_real;			// REAL
using FortranTypes::f_dbl_prec;		// DOUBLE PRECISION
using FortranTypes::f_logical;		// LOGICAL

// ////////////////////////////////////////////////////
// Declarations to link with Fortran QPKWIK procedures

namespace Fortran {
extern "C" {

FORTRAN_FUNC_DECL_UL(void,QPKWIKNEW,qpkwiknew) (
  const f_int& N, const f_int& M1, const f_int& M2, const f_int& M3
  ,const f_dbl_prec GRAD[], f_dbl_prec UINV[], const f_int& LDUINV
  ,const f_int IBND[], const f_dbl_prec BL[], const f_dbl_prec BU[]
  ,const f_dbl_prec A[], const f_int& LDA, const f_dbl_prec YPY[]
  ,const f_int& IYPY, const f_int& WARM, f_dbl_prec NUMPARAM[], const f_int& MAX_ITER
  ,f_dbl_prec X[], f_int* NACTSTORE, f_int IACTSTORE[], f_int* INF
  ,f_int* NACT, f_int IACT[], f_dbl_prec UR[], f_dbl_prec* EXTRA
  ,f_int* ITER, f_int* NUM_ADDS, f_int* NUM_DROPS
  ,f_int ISTATE[], const f_int& LRW, f_dbl_prec RW[]
  );

FORTRAN_FUNC_DECL_UL_(f_int,QPKWIKNEW_LISTATE,qpkwiknew_listate) (
  const f_int& n, const f_int& m1, const f_int& m2, const f_int& m3);

FORTRAN_FUNC_DECL_UL_(f_int,QPKWIKNEW_LRW,qpkwiknew_lrw) (
  const f_int& n, const f_int& m1, const f_int& m2, const f_int& m3);

} // end extern "C"
} // end namespace Fortran

// //////////////////////////////////
// QPKWIK interface functions

// Solve a QP using QPKWIK.
//
// See the Fortran file for documentation.  C++ programs should use this interface.
inline
void qpkwiknew ( 
  const f_int& n, const f_int& m1, const f_int& m2, const f_int& m3
  ,const f_dbl_prec grad[], f_dbl_prec uinv[], const f_int& lduinv
  ,const f_int ibnd[], const f_dbl_prec bl[], const f_dbl_prec bu[]
  ,const f_dbl_prec a[], const f_int& lda, const f_dbl_prec ypy[]
  ,const f_int& iypy, const f_int& warm, f_dbl_prec numparam[], const f_int& max_iter
  ,f_dbl_prec x[], f_int* nactstore, f_int iactstore[], f_int* inf
  ,f_int* nact, f_int iact[], f_dbl_prec ur[], f_dbl_prec* extra
  ,f_int* iter, f_int* num_adds, f_int* num_drops
  ,f_int istate[], const f_int& lrw, f_dbl_prec rw[]
  )
{
  Fortran::FORTRAN_FUNC_CALL_UL(QPKWIKNEW,qpkwiknew) (
    n, m1, m2, m3, grad, uinv, lduinv
    , ibnd, bl, bu, a, lda, ypy, iypy, warm, numparam, max_iter, x, nactstore
    , iactstore, inf, nact, iact, ur, extra, iter, num_adds, num_drops, istate
    , lrw, rw
    );
}

// Get the length of the integer state variables
inline
f_int qpkwiknew_listate(const f_int& n, const f_int& m1, const f_int& m2
            , const f_int& m3)
{
  return Fortran::FORTRAN_FUNC_CALL_UL_(QPKWIKNEW_LISTATE,qpkwiknew_listate) (n, m1, m2, m3);
}

// Get the length of the real (double precision) workspace
inline
f_int qpkwiknew_lrw(const f_int& n, const f_int& m1, const f_int& m2
          , const f_int& m3)
{
  return Fortran::FORTRAN_FUNC_CALL_UL_(QPKWIKNEW_LRW,qpkwiknew_lrw) (n, m1, m2, m3);
}

} // end namespace QPKWIKNEW_CppDecl

// /////////////////////////////////////
// Local helpers

namespace {

template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }

using FortranTypes::f_int;
typedef DenseLinAlgPack::value_type value_type;

enum EConstraintType { NU_L, NU_U, GAMA_L, GAMA_U, LAMBDA, RELAXATION };
char constraint_type_name[6][15] = { "NU_L", "NU_U", "GAMA_L", "GAMA_U", "LAMBDA", "RELAXATION" };

EConstraintType constraint_type( const f_int m1, const f_int m2, const f_int m3, const f_int j )
{
  if     (1 <= j			 && j <= m1			 ) return NU_L;
  else if(m1+1 <= j		 && j <= m1+m2		 ) return GAMA_L;
  else if(m1+m2+1 <= j	 && j <= 2*m1+m2	 ) return NU_U;
  else if(2*m1+m2+1 <= j	 && j <= 2*m1+2*m2	 ) return GAMA_U;
  else if(2*m1+2*m2+1 <= j && j <= 2*m1+2*m2+m3) return LAMBDA;
  else if( j == 2*m1+2*m2+m3 + 1				 ) return RELAXATION;
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return NU_L;	// should never be exectuted
}

f_int constraint_index( const f_int m1, const f_int m2, const f_int m3, const f_int ibnd[]
            , const EConstraintType type, const f_int j )
{
  switch(type) {
    case NU_L		: return ibnd[j-1];
    case GAMA_L		: return j-m1;
    case NU_U		: return ibnd[j-m1-m2-1];
    case GAMA_U		: return j-2*m1-m2;
    case LAMBDA		: return j-2*m1-2*m2;
    case RELAXATION	: return 0;
  }
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return 0;	// should never be exectuted
}

}	// end namespace

// ///////////////////////////////////////
// Members for QPSolverRelaxedQPKWIK

namespace ConstrainedOptPack {

QPSolverRelaxedQPKWIK::QPSolverRelaxedQPKWIK(
  value_type        max_qp_iter_frac
  ,value_type       infinite_bound
  )
  :max_qp_iter_frac_(max_qp_iter_frac)
  ,infinite_bound_(infinite_bound)
  ,N_(0)
  ,M1_(0)
  ,M2_(0)
  ,M3_(0)
{
  NUMPARAM_[0] = 1e-10;	// SMALL
  NUMPARAM_[1] = 1e-20;	// VSMALL
  NUMPARAM_[2] = 1e+20;	// VLARGE
}

QPSolverRelaxedQPKWIK::~QPSolverRelaxedQPKWIK()
{
  this->release_memory();
}

// Overridden from QPSolverRelaxed

QPSolverStats
QPSolverRelaxedQPKWIK::get_qp_stats() const
{
  return qp_stats_;
}

void QPSolverRelaxedQPKWIK::release_memory()
{
  // Todo: resize to zero all the workspace!
}

QPSolverStats::ESolutionType
QPSolverRelaxedQPKWIK::imp_solve_qp(
  std::ostream* out, EOutputLevel olevel, ERunTests test_what
  ,const Vector& g, const MatrixSymOp& G
  ,value_type etaL
  ,const Vector* dL, const Vector* dU
  ,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
  ,const Vector* eL, const Vector* eU
  ,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
  ,value_type* obj_d
  ,value_type* eta, VectorMutable* d
  ,VectorMutable* nu
  ,VectorMutable* mu, VectorMutable* Ed
  ,VectorMutable* lambda, VectorMutable* Fd
  )
{
  using Teuchos::dyn_cast;
  using DenseLinAlgPack::nonconst_tri_ele;
  using LinAlgOpPack::dot;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::assign;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_MtV;
  using AbstractLinAlgPack::EtaVector;
  using AbstractLinAlgPack::transVtMtV;
  using AbstractLinAlgPack::num_bounded;
  using ConstrainedOptPack::MatrixExtractInvCholFactor;

  // /////////////////////////
  // Map to QPKWIK input

  // Validate that rHL is of the proper type.
  const MatrixExtractInvCholFactor &cG
    = dyn_cast<const MatrixExtractInvCholFactor>(G);

  // Determine the number of sparse bounds on variables and inequalities.
  // By default set for the dense case
  const value_type inf = this->infinite_bound();
  const size_type
    nd              = d->dim(),
    m_in            = E  ? b->dim()                  : 0,
    m_eq            = F  ? f->dim()                  : 0,
    nvarbounds      = dL ? num_bounded(*dL,*dU,inf)  : 0,
    ninequbounds    = E  ? num_bounded(*eL,*eU,inf)  : 0,
    nequalities     = F  ? f->dim()                  : 0;

  // Determine if this is a QP with a structure different from the
  // one just solved.
  
  const bool same_qp_struct = (  N_ == nd && M1_ == nvarbounds && M2_ == ninequbounds && M3_ == nequalities );

  /////////////////////////////////////////////////////////////////
  // Set the input parameters to be sent to QPKWIKNEW

  // N
  N_ = nd;

  // M1
  M1_ = nvarbounds;

  // M2
  M2_ = ninequbounds;

  // M3
  M3_ = nequalities;

  // GRAD
  GRAD_ = VectorDenseEncap(g)();

  // UINV_AUG
  //
  // UINV_AUG = [ sqrt(bigM)  0  ]
  //            [ 0           L' ]
  //
  UINV_AUG_.resize(N_+1,N_+1);
  cG.extract_inv_chol( &nonconst_tri_ele( UINV_AUG_(2,N_+1,2,N_+1), BLAS_Cpp::upper ) );
  UINV_AUG_(1,1) = 1.0 / ::sqrt( NUMPARAM_[2] );
  UINV_AUG_.col(1)(2,N_+1) = 0.0;
  UINV_AUG_.row(1)(2,N_+1) = 0.0;

  // LDUINV_AUG
  LDUINV_AUG_ = UINV_AUG_.rows();

  // IBND, BL , BU, A, LDA, YPY

  IBND_INV_.resize( nd + m_in);
  std::fill( IBND_INV_.begin(), IBND_INV_.end(), 0 ); // Initialize the zero
  IBND_.resize( my_max( 1, M1_ + M2_ ) );
  BL_.resize( my_max( 1, M1_ + M2_ ) );
  BU_.resize( my_max( 1, M1_ + M2_ + M3_ ) );
  LDA_ = my_max( 1, M2_ + M3_ );
  A_.resize( LDA_, (  M2_ + M3_ > 0 ? N_ : 1 ) );
  YPY_.resize( my_max( 1, M1_ + M2_ ) );
  if(M1_)
    YPY_(1,M1_) = 0.0; // Must be for this QP interface

  // Initialize variable bound constraints
  if( dL ) {
    VectorDenseEncap dL_de(*dL);
    VectorDenseEncap dU_de(*dU);
    // read iterators
    AbstractLinAlgPack::sparse_bounds_itr
      dLU_itr( dL_de().begin(), dL_de().end()
          ,dU_de().begin(), dU_de().end()
          ,inf );
    // written iterators
    IBND_t::iterator
      IBND_itr = IBND_.begin(),
      IBND_end = IBND_.begin() + M1_;
    DVector::iterator
      BL_itr = BL_.begin(),
      BU_itr = BU_.begin(),
      YPY_itr = YPY_.begin();
    // Loop
    for( size_type ibnd_i = 1; IBND_itr != IBND_end; ++ibnd_i, ++dLU_itr ) {
      IBND_INV_[dLU_itr.index()-1] = ibnd_i;
      *IBND_itr++ = dLU_itr.index();
      *BL_itr++	= dLU_itr.lbound();
      *BU_itr++	= dLU_itr.ubound();
      *YPY_itr++	= 0.0; // Must be zero with this QP interface
    }
  }

  // Initialize inequality constraints
  
  if(M2_) {
    VectorDenseEncap eL_de(*eL);
    VectorDenseEncap eU_de(*eU);
    VectorDenseEncap b_de(*b);
    AbstractLinAlgPack::sparse_bounds_itr
      eLU_itr( eL_de().begin(), eL_de().end()
          ,eU_de().begin(), eU_de().end()
          ,inf );
    if( M2_ < m_in ) {
      // Initialize BL, BU, YPY and A for sparse bounds on general inequalities
      // written iterators
      DVector::iterator
        BL_itr		= BL_.begin() + M1_,
        BU_itr		= BU_.begin() + M1_,
        YPY_itr		= YPY_.begin() + M1_;
      IBND_t::iterator
        ibnds_itr	= IBND_.begin() + M1_;
      // loop
      for(size_type i = 1; i <= M2_; ++i, ++eLU_itr, ++ibnds_itr ) {
        TEUCHOS_TEST_FOR_EXCEPT( !( !eLU_itr.at_end() ) );
        const size_type k      = eLU_itr.index();
        *BL_itr++              = eLU_itr.lbound();
        *BU_itr++              = eLU_itr.ubound();
        *YPY_itr++             = b_de()(k);
        *ibnds_itr             = k;  // Only for my record, not used by QPKWIK
        IBND_INV_[nd+k-1]      = M1_ + i;
        // Add the corresponding row of op(E) to A
        // y == A.row(i)'
        // y' = e_k' * op(E) => y = op(E')*e_k
        DVectorSlice y = A_.row(i);
        EtaVector e_k(k,eL_de().dim());
        V_MtV( &y( 1, N_ ), *E, BLAS_Cpp::trans_not(trans_E), e_k() ); // op(E')*e_k
      }
    }
    else {
      //
      // Initialize BL, BU, YPY and A for dense bounds on general inequalities
      //
      // Initialize BL(M1+1:M1+M2), BU(M1+1:M1+M2)
      // and IBND(M1+1:M1+M2) = identity (only for my record, not used by QPKWIK)
      DVector::iterator
        BL_itr		= BL_.begin() + M1_,
        BU_itr		= BU_.begin() + M1_;
      IBND_t::iterator
        ibnds_itr	= IBND_.begin() + M1_;
      for(size_type i = 1; i <= m_in; ++i ) {
        if( !eLU_itr.at_end() && eLU_itr.index() == i ) {
          *BL_itr++ = eLU_itr.lbound();
          *BU_itr++ = eLU_itr.ubound();
          ++eLU_itr;
        }
        else {
          *BL_itr++ = -inf;
          *BU_itr++ = +inf;
        }
        *ibnds_itr++     = i;
        IBND_INV_[nd+i-1]= M1_ + i;
      }
      // A(1:M2,1:N) = op(E)
      assign( &A_(1,M2_,1,N_), *E, trans_E );
      // YPY
      YPY_(M1_+1,M1_+M2_) = b_de();
    }
  }

  // Initialize equalities

  if(M3_) {
    V_StV( &BU_( M1_ + M2_ + 1, M1_ + M2_ + M3_ ), -1.0, VectorDenseEncap(*f)() );
    assign( &A_( M2_ + 1, M2_ + M3_, 1, N_ ), *F, trans_F );
  }

  // IYPY
  IYPY_ = 1; // ???

  // WARM
  WARM_ = 0; // Cold start by default

  // MAX_ITER
  MAX_ITER_ = static_cast<f_int>(max_qp_iter_frac() * N_);

  // INF
  INF_ = ( same_qp_struct ? 1 : 0 );
  
  // Initilize output, internal state and workspace quantities.
  if(!same_qp_struct) {
    X_.resize(N_);
    NACTSTORE_ = 0;
    IACTSTORE_.resize(N_+1);
    IACT_.resize(N_+1);
    UR_.resize(N_+1);
    ISTATE_.resize( QPKWIKNEW_CppDecl::qpkwiknew_listate(N_,M1_,M2_,M3_) );
    LRW_ = QPKWIKNEW_CppDecl::qpkwiknew_lrw(N_,M1_,M2_,M3_);
    RW_.resize(LRW_);
  }

  // /////////////////////////////////////////////
  // Setup a warm start form the input arguments
  //
  // Interestingly enough, QPKWIK sorts all of the
  // constraints according to scaled multiplier values
  // and mixes equality with inequality constriants.
  // It seems to me that you should start with equality
  // constraints first.

  WARM_      = 0;
  NACTSTORE_ = 0;

  if( m_eq ) {
    // Add equality constraints first since we know these will
    // be part of the active set.
    for( size_type j = 1; j <= m_eq; ++j ) {
      IACTSTORE_[NACTSTORE_] = 2*M1_ + 2*M2_ + j;
      ++NACTSTORE_;
    }
  }
  if( ( nu && nu->nz() ) || ( mu && mu->nz() ) ) {
    // Add inequality constraints
    const size_type
      nu_nz = nu ? nu->nz() : 0,
      mu_nz = mu ? mu->nz() : 0;
    // Combine all the multipliers for the bound and general inequality
    // constraints and sort them from the largest to the smallest.  Hopefully
    // the constraints with the larger multiplier values will not be dropped
    // from the active set.
    SpVector gamma( nd + 1 + m_in , nu_nz + mu_nz );
    typedef SpVector::element_type ele_t;
    if(nu && nu_nz) {
      VectorDenseEncap nu_de(*nu);
      DVectorSlice::const_iterator
        nu_itr = nu_de().begin(),
        nu_end = nu_de().end();
      index_type i = 1;
      while( nu_itr != nu_end ) {
        for( ; *nu_itr == 0.0; ++nu_itr, ++i );
        gamma.add_element(ele_t(i,*nu_itr));
        ++nu_itr; ++i;
      }
    }
    if(mu && mu_nz) {
      VectorDenseEncap mu_de(*mu);
      DVectorSlice::const_iterator
        mu_itr = mu_de().begin(),
        mu_end = mu_de().end();
      index_type i = 1;
      while( mu_itr != mu_end ) {
        for( ; *mu_itr == 0.0; ++mu_itr, ++i );
        gamma.add_element(ele_t(i+nd,*mu_itr));
        ++mu_itr; ++i;
      }
    }
    std::sort( gamma.begin(), gamma.end()
      , AbstractLinAlgPack::SortByDescendingAbsValue() );
    // Now add the inequality constraints in decreasing order
    const SpVector::difference_type o = gamma.offset();
    for( SpVector::const_iterator itr = gamma.begin(); itr != gamma.end(); ++itr ) {
      const size_type  j   = itr->index() + o;
      const value_type val = itr->value();
      if( j <= nd ) { // Variable bound
        const size_type ibnd_i = IBND_INV_[j-1];
        TEUCHOS_TEST_FOR_EXCEPT( !( ibnd_i ) );
        IACTSTORE_[NACTSTORE_]
          = (val < 0.0
             ? ibnd_i               // lower bound (see IACT(*))
             : M1_ + M2_ + ibnd_i   // upper bound (see IACT(*))
            );
        ++NACTSTORE_;
      }
      else if( j <= nd + m_in ) { // General inequality constraint
        const size_type ibnd_i = IBND_INV_[j-1]; // offset into M1_ + ibnd_j
        TEUCHOS_TEST_FOR_EXCEPT( !( ibnd_i ) );
        IACTSTORE_[NACTSTORE_]
          = (val < 0.0
             ? ibnd_i               // lower bound (see IACT(*))
             : M1_ + M2_ + ibnd_i   // upper bound (see IACT(*))
            );
        ++NACTSTORE_;
      }
    }
  }
  if( NACTSTORE_ > 0 )
    WARM_ = 1;

  // /////////////////////////
  // Call QPKWIK

  if( out && olevel > PRINT_NONE ) {
    *out
      << "\nCalling QPKWIK to solve QP problem ...\n";
  }

  QPKWIKNEW_CppDecl::qpkwiknew(
    N_, M1_, M2_, M3_, &GRAD_(1), &UINV_AUG_(1,1), LDUINV_AUG_, &IBND_[0]
    ,&BL_(1), &BU_(1), &A_(1,1), LDA_, &YPY_(1), IYPY_, WARM_, NUMPARAM_, MAX_ITER_, &X_(1)
    ,&NACTSTORE_, &IACTSTORE_[0], &INF_, &NACT_, &IACT_[0], &UR_[0], &EXTRA_
    ,&ITER_, &NUM_ADDS_, &NUM_DROPS_, &ISTATE_[0], LRW_, &RW_[0]
    );

  // ////////////////////////
  // Map from QPKWIK output

  // eta
  *eta = EXTRA_;
  // d
  (VectorDenseMutableEncap(*d))() = X_();
  // nu (simple variable bounds) and mu (general inequalities)
  if(nu) *nu = 0.0;
  if(mu) *mu = 0.0;
  // ToDo: Create VectorDenseMutableEncap views for faster access!
  {for(size_type i = 1; i <= NACT_; ++i) {
    size_type j = IACT_[i-1];
    EConstraintType type = constraint_type(M1_,M2_,M3_,j);
    FortranTypes::f_int idc = constraint_index(M1_,M2_,M3_,&IBND_[0],type,j);
    switch(type) {
      case NU_L:
        nu->set_ele( idc , -UR_(i) );
        break;
      case GAMA_L:
        mu->set_ele( IBND_[ M1_ + idc - 1 ], -UR_(i) );
        break;
      case NU_U:
        nu->set_ele( idc, UR_(i)) ;
        break;
      case GAMA_U:
        mu->set_ele( IBND_[ M1_ + idc - 1 ], UR_(i) );
        break;
      case LAMBDA:
        lambda->set_ele( idc, UR_(i) );
        break;
    }
  }}
  // obj_d (This could be updated within QPKWIK in the future)
  if(obj_d) {
    // obj_d = g'*d + 1/2 * d' * G * g
    *obj_d = dot(g,*d) + 0.5 * transVtMtV(*d,G,BLAS_Cpp::no_trans,*d);
  }
  // Ed (This could be updated within QPKWIK in the future)
  if(Ed) {
    V_MtV( Ed, *E, trans_E, *d );
  }
  // Fd (This could be updated within QPKWIK in the future)
  if(Fd) {
    V_MtV( Fd, *F, trans_F, *d );
  }
  // Set the QP statistics
  QPSolverStats::ESolutionType solution_type;
  if( INF_ >= 0 ) {
    solution_type = QPSolverStats::OPTIMAL_SOLUTION;
  }
  else if( INF_ == -1 ) { // Infeasible constraints
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, QPSolverRelaxed::Infeasible
      ,"QPSolverRelaxedQPKWIK::solve_qp(...) : Error, QP is infeasible" );
  }
  else if( INF_ == -2 ) { // LRW too small
    TEUCHOS_TEST_FOR_EXCEPT( !( INF_ != -2 ) );  // Local programming error?
  }
  else if( INF_ == -3 ) { // Max iterations exceeded
    solution_type = QPSolverStats::DUAL_FEASIBLE_POINT;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true); // Unknown return value!
  }
  qp_stats_.set_stats(
    solution_type, QPSolverStats::CONVEX
    ,ITER_, NUM_ADDS_, NUM_DROPS_
    ,WARM_==1, *eta > 0.0 );

  return qp_stats_.solution_type();
}


}	// end namespace ConstrainedOptPack


#endif // CONSTRAINED_OPTIMIZATION_PACK_USE_QPKWIK

