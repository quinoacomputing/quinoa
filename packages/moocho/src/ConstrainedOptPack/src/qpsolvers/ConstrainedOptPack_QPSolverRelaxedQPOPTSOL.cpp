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

#include <vector>

#include "ConstrainedOptPack_QPSolverRelaxedQPOPTSOL.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_SortByDescendingAbsValue.hpp"
#include "AbstractLinAlgPack_sparse_bounds.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "ProfileHackPack_profile_hack.hpp"

// /////////////////////////////////////////////////////////////////
//
// This subclass uses a relaxation of the equality and inequality
// constraints.  The mapping to the arguments of QPOPT or QPSOL
// is done as follows.
//
//  QP formulation:
//  ---------------
//
//  min          g'*d + 1/2*d'*G*d + (eta + 1/2*eta^2)*M
//  d <: R^n
//         
//  s.t.
//               etaL <=  eta
//               dL   <=  d                       <= dU
//               eL   <=  op(E)*d - b*eta         <= eU
//                        op(F)*d + (1 - eta) * f  = 0
//
//  Rearranged to :
//  ---------------
//
//  min          [ g', M ] * [  d  ] + 1/2 * [ d', eta ] * [ G  0 ] * [  d  ]
//                           [ eta ]                       [ 0  M ]   [ eta ]
//
//  s.t.         [  bL  ]    [   I  ,  0 ]              [ dU  ]
//               [ etaL ] <= [   0  ,  1 ] * [  d  ] <= [ inf ]
//               [  eL  ]    [ op(E), -b ]   [ eta ]    [ eU  ]
//               [  -f  ]    [ op(F), -f ]              [ -f  ]
//
//  Which maps to the QPSOL interface which is:
//  -------------------------------------------
//
//  min           CVEC' * X + 1/2 * X'* H * X
//
//  s.t.          BL <= [    X   ] <= BU
//                      [  A * X ]
//
//  Which gives us:
//
//  X    = [ d; eta ] 
//  CVEC = [ g; M ]
//  H    = [ G, 0; 0, M ]
//  BL   = [ bL, etaL, eL, -f ]
//  BU   = [ bU, inf,  eU, -f ]
//  A    = [ op(E), -b; op(F), -f ]
//

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::Mp_StM;
  using AbstractLinAlgPack::Vp_StMtV;
}

// ///////////////////////////////////////
// Members for QPSolverRelaxedQPOPTSOL

namespace ConstrainedOptPack {

QPSolverRelaxedQPOPTSOL::QPSolverRelaxedQPOPTSOL()
  :N_(0)
  ,bigM_(1e+10)
  ,use_as_bigM_(1e+10)
  ,G_(NULL)
{}

QPSolverRelaxedQPOPTSOL::~QPSolverRelaxedQPOPTSOL()
{
  this->release_memory();
}

const MatrixOp* QPSolverRelaxedQPOPTSOL::G() const
{
  return G_;
}

value_type QPSolverRelaxedQPOPTSOL::use_as_bigM() const
{
  return use_as_bigM_;
}

// Overridden from QPSolverRelaxed

QPSolverStats
QPSolverRelaxedQPOPTSOL::get_qp_stats() const
{
  return qp_stats_;
}

void QPSolverRelaxedQPOPTSOL::release_memory()
{
  // Todo: resize to zero all the matrices and vectors
}

QPSolverStats::ESolutionType
QPSolverRelaxedQPOPTSOL::imp_solve_qp(
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

  using AbstractLinAlgPack::VectorDenseEncap;
  using AbstractLinAlgPack::VectorDenseMutableEncap;

#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "QPSolverRelaxedQPOPTSOL::imp_solve_qp(...)" );
#endif

  const size_type n = d->dim();
  const value_type inf = this->infinite_bound();
  
  //
  // Map to the input arguments for QPOPT or QPSOL
  //

  // N
  N_ = n + 1; // With relaxation

  // NCLIN
  n_inequ_bnds_ = ( E ? AbstractLinAlgPack::num_bounded(*eL,*eU,inf) : 0 );
  NCLIN_ = n_inequ_bnds_ + (F ? f->dim() : 0);

  // A, BL, BU
  A_.resize(NCLIN_,N_);
  BL_.resize(N_+NCLIN_);
  BU_.resize(N_+NCLIN_);
  if(dL) {
    VectorDenseEncap dL_de(*dL);
    BL_(1,n) = dL_de();
  }
  else {
    BL_(1,n) = -inf;
  }
  if(dU) {
    VectorDenseEncap dU_de(*dU);
    BU_(1,n) = dU_de();
  }
  else {
    BU_(1,n) = -inf;
  }
  BL_(N_) = etaL;
  BU_(N_) = +inf;
  TEUCHOS_TEST_FOR_EXCEPTION(
    E!=NULL, std::logic_error
    ,"Error, the QPOPT/QPSOL wrapper has not been updated for general inequalities yet!"
    );
/* ToDo: Update this when needed!
  if( E ) {
    i_inequ_bnds_.resize(n_inequ_bnds_);
    if( n_inequ_bnds_ < b->dim() ) {
      // Initialize BL, BU, and A for sparse bounds on general inequalities
      //
      // read iterators
      AbstractLinAlgPack::sparse_bounds_itr
        eLU_itr( eL->begin(), eL->end(), eL->offset()
             , eU->begin(), eU->end(), eU->offset(), inf );
      // written iterators
      DVector::iterator
        BL_itr		= BL_.begin() + N_,
        BU_itr		= BU_.begin() + N_;
      ibnds_t::iterator
        ibnds_itr	= i_inequ_bnds_.begin();
      // loop
      for(size_type i = 1; i <= n_inequ_bnds_; ++i, ++eLU_itr, ++ibnds_itr ) {
        TEUCHOS_TEST_FOR_EXCEPT( !( !eLU_itr.at_end() ) );
        const size_type k      = eLU_itr.indice();
        *BL_itr++              = eLU_itr.lbound();
        *BU_itr++              = eLU_itr.ubound();
        *ibnds_itr             = k;  // Only for my record
        // Add the corresponding row of [ op(E), -b ] to A
        // y == A.row(i)
        // y(1,n) = op(E')*e_k
        DVectorSlice y = A_.row(i);
        AbstractLinAlgPack::EtaVector e_k(k,eL->dim());
        LinAlgOpPack::V_MtV( &y(1,n), *E, BLAS_Cpp::trans_not(trans_E), e_k() ); // op(E')*e_k
        // y(n+1) = -b(k)
        y(n+1) = -(*b)(k);
      }
    }
    else {
      // Initialize BL, BU and A for dense bounds on general inequalities
      //
      // Initialize BL(N+1:N+n_inequ_bnds), BU(N+1:N+n_inequ_bnds)
      // and i_inequ_bnds_ = identity (only for my record, not used by QPKWIK)
      AbstractLinAlgPack::sparse_bounds_itr
        eLU_itr( eL->begin(), eL->end(), eL->offset()
             , eU->begin(), eU->end(), eU->offset(), inf );
      DVector::iterator
        BL_itr		= BL_.begin() + N_,
        BU_itr		= BU_.begin() + N_;
      ibnds_t::iterator
        ibnds_itr	= i_inequ_bnds_.begin();
      for(size_type i = 1; i <= n_inequ_bnds_; ++i, ++eLU_itr, ++ibnds_itr ) {
        TEUCHOS_TEST_FOR_EXCEPT( !( !eLU_itr.at_end() ) );
        const size_type k      = eLU_itr.indice();
        *BL_itr++              = eLU_itr.lbound();
        *BU_itr++              = eLU_itr.ubound();
        *ibnds_itr             = k;  // Only for my record
      }
      // A(1:n_inequ_bnds,1:n) = op(E)
      LinAlgOpPack::assign( &A_(1,n_inequ_bnds_,1,n), *E, trans_E );
      // A(1:n_inequ_bnds,n+1) = -b
      LinAlgOpPack::V_StV( &A_.col(n+1)(1,n_inequ_bnds_), -1.0, *b );
    }
  }
*/
  TEUCHOS_TEST_FOR_EXCEPTION(
    F!=NULL, std::logic_error
    ,"Error, the QPOPT/QPSOL wrapper has not been updated for general equalities yet!"
    );
/* ToDo: Update this when needed!
  if( F ) {
    // BL(N+n_inequ_bnds+1:N+NCLIN) = -f
    LinAlgOpPack::V_StV( &BL_(N_+n_inequ_bnds_+1,N_+NCLIN_), -1.0, *f );
    // BU(N+n_inequ_bnds+1:N+NCLIN) = -f
    LinAlgOpPack::V_StV( &BU_(N_+n_inequ_bnds_+1,N_+NCLIN_), -1.0, *f );
    // A(n_inequ_bnds+1:NCLIN,1:n) = op(F)
    LinAlgOpPack::assign( &A_(n_inequ_bnds_+1,NCLIN_,1,n), *F, trans_F );
    // A(n_inequ_bnds+1:NCLIN,n+1) = -f
    LinAlgOpPack::V_StV( &A_.col(n+1)(n_inequ_bnds_+1,NCLIN_), -1.0, *f );
  }
*/
  
  // CVEC
  CVEC_.resize(N_);
  CVEC_(1,n) = VectorDenseEncap(g)();
  CVEC_(n+1) = bigM_;

  // HESS
  G_ = &G; // That's all there is to it!

  // ISTATE
  ISTATE_.resize(N_+NCLIN_);
  std::fill( ISTATE_.begin(), ISTATE_.end(), 0 ); // cold start!
  ISTATE_[n] = 1; // Make eta >= etaL active

  // X
  X_.resize(N_);
  X_(1,n) = VectorDenseEncap(*d)();
  X_(n+1) = *eta;

  // AX
  // will be resized by QPOPT but not QPSOL

  // CLAMBDA
  CLAMDA_.resize(N_+NCLIN_);

  // LIWORK, IWORK
  LIWORK_ = liwork(N_,NCLIN_);
  if(static_cast<f_int>(IWORK_.size()) < LIWORK_)	IWORK_.resize(LIWORK_);

  // LWORK, WORK
  LWORK_ = lrwork(N_,NCLIN_);
  if(static_cast<f_int>(WORK_.size()) < LWORK_) WORK_.resize(LWORK_);

  // We need to initialize some warm start information if
  // it was given by the user!
  bool warm_start = false;
  if( (nu && nu->nz()) || (mu && mu->nz() ) ) {
    // Let's try a warm start
    if(nu) {
      VectorDenseEncap nu_de(*nu);
      for(int i = 1; i <= n; ++i ) {
        if( nu_de()(i) < 0.0 )
          ISTATE_[i-1] = 1; // Lower bound is active
        else if( nu_de()(i) > 0.0 )
          ISTATE_[i-1] = 2; // Upper bound is active
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      mu!=NULL, std::logic_error
      ,"Error, the QPOPT/QPSOL wrapper has not been updated for general inequalities yet!"
      );
/* ToDo: Update below when needed!
    if(mu) {
      const SpVectorSlice::difference_type o = mu->offset();
      for( SpVectorSlice::const_iterator itr = mu->begin(); itr != mu->end(); ++itr ) {
        if( itr->value() < 0.0 )
          ISTATE_[ itr->indice() + o + n ] = 1; // Lower bound is active
        else if( itr->value() > 0.0 )
          ISTATE_[ itr->indice() + o + n ] = 2; // Upper bound is active
      }
    }
*/
    warm_start = true;
  }

  //
  // Solve the QP using QPOPT or QPSOL
  //

  const EInform inform_return = call_qp_solver(warm_start);

  //
  // Map from the output from QPOPT or QPSOL
  //

  // d
  {
    VectorDenseMutableEncap d_de(*d);
    d_de() = X_(1,n);
  }
  
  // eta
  *eta = X_(n+1);

  // obj_d
  if(obj_d)
    *obj_d = OBJ_ - (*eta) * bigM_ - 0.5 * (*eta)*(*eta) * use_as_bigM_; 

  // nu
  if(nu) {
    VectorDenseMutableEncap nu_de(*nu);
    nu_de() = 0.0;
    ISTATE_t::const_iterator
      istate_itr = ISTATE_.begin();
    DVector::const_iterator
      clamda_itr = CLAMDA_.begin();
    for( size_type i = 1; i <= n; ++i, ++istate_itr, ++clamda_itr ) {
      const f_int state = *istate_itr;
      switch(state) {
        case -2: // The lower bound is violated by more than feas_tol
        case -1: // The upper bound is violated by more than feas_tol
          // What do we do?
          break;
        case 0: // Within bounds by more than feas_tol
          break;
        case 1: // lower bound is active
        case 2: // upper bound is active
        case 3: // the bounds are equal and are satisfied
          nu_de()(i) = -(*clamda_itr); // Different sign convention
          break;
        case 4: // Temporaraly fixed at current value
          // What do we do?
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here!
      }
    }
  }
  
  // mu
  TEUCHOS_TEST_FOR_EXCEPTION(
    n_inequ_bnds_!=0, std::logic_error
    ,"Error, the QPOPT/QPSOL wrapper has not been updated for general inequalities yet!"
    );
/* ToDo: Update below when needed!
  if( n_inequ_bnds_ ) {
    mu->resize(b->dim(),n_inequ_bnds_);
    typedef SpVector::element_type ele_t;
    ISTATE_t::const_iterator
      istate_itr = ISTATE_.begin() + N_;
    DVector::const_iterator
      clamda_itr = CLAMDA_.begin() + N_;
    ibnds_t::const_iterator
      bnd_itr = i_inequ_bnds_.begin();
    for( size_type k = 1; k <= n_inequ_bnds_; ++k, ++istate_itr, ++clamda_itr, ++bnd_itr )
    {
      const f_int state = *istate_itr;
      const size_type j = *bnd_itr;
      switch(state) {
          case -2: // The lower bound is violated by more than feas_tol
          case -1: // The upper bound is violated by more than feas_tol
          // What do we do?
          break;
          case 0: // Within bounds by more than feas_tol
          break;
          case 1: // lower bound is active
          case 2: // upper bound is active
          case 3: // the bounds are equal and are satisfied
          mu->add_element(ele_t(j,-(*clamda_itr))); // Different sign!
          break;
          case 4: // Temporaraly fixed at current value
          // What do we do?
          break;
          default:
          TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here!
      }
    }
    mu->assume_sorted(true);
  }
  else if(E) {
    mu->resize( eL->dim(), 0 );
  }
*/

  TEUCHOS_TEST_FOR_EXCEPTION(
    F!=NULL, std::logic_error
    ,"Error, the QPOPT/QPSOL wrapper has not been updated for general equalities yet!"
    );
/* ToDo: Update this when needed!

  // lambda
  if( F ) {
    LinAlgOpPack::V_StV( lambda, -1.0, CLAMDA_(N_+n_inequ_bnds_+1,N_+NCLIN_) );
    // Validate istate
    ISTATE_t::const_iterator
      istate_itr = ISTATE_.begin() + N_ + n_inequ_bnds_;
    for( size_type k = 1; k <= f->dim(); ++k, ++istate_itr ) {
      TEUCHOS_TEST_FOR_EXCEPT( !(  *istate_itr == 3  ) );
    }
  }

  // Ed, Fd
  if( E && AX_.size() && eL->dim() == n_inequ_bnds_ ) {
    if( Ed ) { // Ed = AX + b*eta
      *Ed = AX_(1,n_inequ_bnds_);
      if( *eta > 0.0 )
        LinAlgOpPack::Vp_StV( Ed, *eta, *b );
    }
    if( Fd ) { // Fd = AX + f*eta
      *Fd = AX_(n_inequ_bnds_+1,NCLIN_);
      if( *eta > 0.0 )
        LinAlgOpPack::Vp_StV( Fd, *eta, *f );
    }
  }
  else {
    if(Ed)
      LinAlgOpPack::V_MtV( Ed, *E, trans_E, *d );
    if(Fd)
      LinAlgOpPack::V_MtV( Fd, *F, trans_F, *d );
  }

*/

  //
  // Setup the QP statistics
  //

  QPSolverStats::ESolutionType solution_type  = QPSolverStats::SOLUTION_TYPE_NOT_KNOWN;
  QPSolverStats::EConvexity    convexity_type = QPSolverStats::CONVEXITY_NOT_KNOWN;
  switch(inform_return) {
      case STRONG_LOCAL_MIN:
      solution_type  = QPSolverStats::OPTIMAL_SOLUTION;
      convexity_type =  QPSolverStats::CONVEX;
      break;
      case WEAK_LOCAL_MIN:
      solution_type  = QPSolverStats::OPTIMAL_SOLUTION;
      convexity_type =  QPSolverStats::NONCONVEX;
      break;
      case MAX_ITER_EXCEEDED:
      solution_type  = QPSolverStats::PRIMAL_FEASIBLE_POINT;
      convexity_type =  QPSolverStats::CONVEXITY_NOT_KNOWN;
      break;
      case OTHER_ERROR:
      solution_type  = QPSolverStats::SUBOPTIMAL_POINT;
      convexity_type =  QPSolverStats::CONVEXITY_NOT_KNOWN;
      break;
  }
  qp_stats_.set_stats(
    solution_type, convexity_type
    ,ITER_, QPSolverStats::NOT_KNOWN, QPSolverStats::NOT_KNOWN
    ,warm_start, *eta > 0.0 );

  return qp_stats_.solution_type();
}

}	// end namespace ConstrainedOptPack
