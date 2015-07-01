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
// Here we map from the QPSolverRelaxed QP formulation to the LOQO QP formulation.
//
// QPSolverRelaxed QP formulation:
// ------------------------------
//
// min     g'*d + 1/2 * d'*G*d + (eta + 1/2*eta^2) * M
// s.t.    dL   <= d                   <= dU
//         etaL <= eta
//         eL   <= op(E)*d - b*eta     <= eU
//                 op(F)*d + (1-eta)*f  = 0
//
// LOQO QP formulation:
// -------------------
//
// min      c'*x + 1/2 * x'*Q*x
// s.t.     b <= A*x <= b + r
//          l <= x <= u
//
// Mapping =>
//
// LOQO   QPSolverRelaxed 
// ----   ---------------
// x      [ d; eta ]
// c      [ g; M ]
// Q      [ G, 0; 0, M ]
// A      [ op(E), -b; op(F), -f ]
// b      [ eL; -f ]
// r      [ eU-eL; 0 ]
// l      [ dL, etaL ]
// u      [ dU, +inf ]
//
// Above, in the LOQO formulation all singly bounded inequalities
// must be formulated as b(j) <= A(j,:)*x with r(j) = inf.  This
// will require some fudging since eL(j) == -inf may be true in some
// cases.  Here we will have to exchange eL(j) and eU(j) and use
// A(j,:) = -op(E)(j,:).
//

#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_LOQO

#include <assert.h>

#include <vector>

#include "ConstrainedOptPack_QPSolverRelaxedLOQO.hpp"
#include "ConstrainedOptPack/src/AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_SortByDescendingAbsValue.hpp"
#include "AbstractLinAlgPack_sparse_bounds.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_sparse_bounds.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Midynamic_cast_verbose.h"
#include "MiWorkspacePack.h"

extern "C" {
#include "loqo.h"     // -I$(LOQODIR)
#include "myalloc.h"  // -I$(LOQODIR)
} // end extern "C"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::Mp_StM;
  using AbstractLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptPack {

// ///////////////////////////////////////
// Members for QPSolverRelaxedLOQO::InitLOQOHessianJacobian

void QPSolverRelaxedLOQO::InitLOQOHessianJacobian::init_hess_jacob(
  const MatrixOp& G, const value_type bigM
  , const MatrixOp* E, BLAS_Cpp::Transp trans_E, const DVectorSlice* b
  , const int loqo_b_stat[], const size_type num_inequal
  , const MatrixOp* F, BLAS_Cpp::Transp trans_F, const DVectorSlice* f
  , void* _loqo_lp
  ) const
{

  LOQO* loqo_lp = (LOQO*)_loqo_lp;

  const size_type
    nd    = G.rows(),
    m_in  = E ? b->size() : 0,
    m_eq  = F ? f->size() : 0;

  TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp->n == nd + 1  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp->m == num_inequal + m_eq  ) );

  // This default implementation assumes G, E and F are completely dense!

  //
  // Setup Q
  //

  loqo_lp->qnz = nd*nd + 1;
  MALLOC( loqo_lp->Q, loqo_lp->qnz, double );
  MALLOC( loqo_lp->iQ, loqo_lp->qnz, int );
  MALLOC( loqo_lp->kQ, nd+2, int );
  // Setup kQ[] and iQ[]
  {for( size_type j = 1; j <= nd; ++j ) {
    loqo_lp->kQ[j-1] = nd*(j-1);
    for( size_type i = 1; i <= nd; ++i )
      loqo_lp->iQ[ loqo_lp->kQ[j-1] + (i-1) ] = i-1; // zero based in LOQO
  }}
  loqo_lp->kQ[nd]              = nd*nd;
  loqo_lp->iQ[loqo_lp->kQ[nd]] = nd; // zero based in LOQO
  loqo_lp->kQ[nd+1]            = nd*nd + 1;
  // Setup Q[]
  {
    DMatrixSlice Q( loqo_lp->Q, nd*nd, nd, nd, nd );
    LinAlgOpPack::assign( &Q, G, BLAS_Cpp::no_trans );
    loqo_lp->Q[nd*nd] = bigM;
  }

  //
  // Setup A
  //

  loqo_lp->nz = (num_inequal+m_eq) * (nd+1);
  MALLOC( loqo_lp->A,  loqo_lp->nz, double );
  MALLOC( loqo_lp->iA, loqo_lp->nz, int );
  MALLOC( loqo_lp->kA, nd+2, int );

  if( num_inequal == m_in ) {
    // All the inequalities have finite bounds
    // Setup kA[] and iA[]
    {for( size_type j = 1; j <= nd+1; ++j ) {
      loqo_lp->kA[j-1] = (m_in+m_eq)*(j-1);
      for( size_type i = 1; i <= m_in+m_eq; ++i )
        loqo_lp->iA[ loqo_lp->kA[j-1] + (i-1) ] = i-1; // zero based in LOQO
    }}
    loqo_lp->kA[nd+1] = (m_in+m_eq)*(nd+1);
    // Setup A[]
    DMatrixSlice A( loqo_lp->A, loqo_lp->nz, loqo_lp->m, loqo_lp->m, nd+1 );
    if(E) {
      LinAlgOpPack::assign( &A(1,m_in,1,nd), *E, trans_E );  // A(1:m_in,1:nd) = op(E)
      LinAlgOpPack::V_StV( &A.col(nd+1)(1,m_in), -1.0, *b ); // A(1:m_in,nd+1) = -b
    }
    if(F) {
      LinAlgOpPack::assign( &A(m_in+1,m_in+m_eq,1,nd), *F, trans_F );  // A(m_in+1:m_in+m_eq,1:nd) = op(F)
      LinAlgOpPack::V_StV( &A.col(nd+1)(m_in+1,m_in+m_eq), -1.0, *f ); // A(m_in+1:m_in+m_eq,nd+1) = -f
    }
  }
  else {
    // At least one of the inequality constriants has
    // both infinite upper and lower bounds.
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Finish this!
  }
  
  // Loop through and adjust A for absent lower bound and using upper bound
  if( num_inequal ) {
    DMatrixSlice A( loqo_lp->A, loqo_lp->nz, loqo_lp->m, loqo_lp->m, nd+1 );
    for(size_type k = 1; k <= num_inequal; ++k ) {
      const int j = loqo_b_stat[k-1];
      if( j < 0 )
        DenseLinAlgPack::Vt_S( &A.row(j), -1.0 );
    }
  }

}

// ///////////////////////////////////////
// Members for QPSolverRelaxedLOQO

QPSolverRelaxedLOQO::QPSolverRelaxedLOQO(
  const init_hess_jacob_ptr_t  init_hess_jacob
  ,value_type                  bigM
  ,value_type                  nonbinding_lag_mult
  )
  :init_hess_jacob_(init_hess_jacob)
  ,bigM_(bigM)
  ,nonbinding_lag_mult_(nonbinding_lag_mult)
{
//	bigM_ = 1.0; // Just test this!
  nonbinding_lag_mult_ = 1e-6;
}

QPSolverRelaxedLOQO::~QPSolverRelaxedLOQO()
{
  this->release_memory();
}

// Overridden from QPSolverRelaxed

QPSolverStats
QPSolverRelaxedLOQO::get_qp_stats() const
{
  return qp_stats_;
}

void QPSolverRelaxedLOQO::release_memory()
{
  // Todo: resize to zero all the workspace!
}

QPSolverStats::ESolutionType
QPSolverRelaxedLOQO::imp_solve_qp(
      std::ostream* out, EOutputLevel olevel, ERunTests test_what
    , const DVectorSlice& g, const MatrixOp& G
    , value_type etaL
    , const SpVectorSlice& dL, const SpVectorSlice& dU
    , const MatrixOp* E, BLAS_Cpp::Transp trans_E, const DVectorSlice* b
      , const SpVectorSlice* eL, const SpVectorSlice* eU
    , const MatrixOp* F, BLAS_Cpp::Transp trans_F, const DVectorSlice* f
    , value_type* obj_d
    , value_type* eta, DVectorSlice* d
    , SpVector* nu
    , SpVector* mu, DVectorSlice* Ed
    , DVectorSlice* lambda, DVectorSlice* Fd
  )
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = wsp::default_workspace_store.get();

  const value_type inf_bnd  = std::numeric_limits<value_type>::max();
//	const value_type real_big = 1e+20;
  const value_type real_big = HUGE_VAL;

  const size_type
    nd   = g.size(),
    m_in = E ? b->size() : 0,
    m_eq = F ? f->size() : 0;

  //
  // Create a LOQO QP definition struct
  //

  LOQO *loqo_lp = openlp();
  TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp  ) );

  //
  // Setup loqo_r and loqo_b and count the number of actual
  // constraints.
  //

  // LOQO's b vector storage
  MALLOC( loqo_lp->b, m_in+m_eq, double ); // May not use all of this storage
  DVectorSlice loqo_b( loqo_lp->b, m_in+m_eq );
  // LOQO's r vector storage
  MALLOC( loqo_lp->r, m_in+m_eq, double ); // May not use all of this storage
  DVectorSlice loqo_r( loqo_lp->r, m_in+m_eq );
  // Gives status of b.
  //                  /  j : if eL(j) > -inf_bnd
  // loqo_b_stat(k) = |
  //                  \ -j : if eL(j) <= -inf_bnd && eU(j) < +inf_bnd
  //
  // , for k = 1...num_inequal
  //
  Workspace<int>               loqo_b_stat_ws(wss,m_in); // May not use all of this
  DenseLinAlgPack::VectorSliceTmpl<int>  loqo_b_stat(&loqo_b_stat_ws[0],loqo_b_stat_ws.size());
  std::fill( loqo_b_stat.begin(), loqo_b_stat.end(), 0 ); // Initialize to zero

  // Fill up loqo_b, loqo_r and loqo_b_stat
  size_type num_inequal = 0; // The actual number of bouned general inequalities
  if(E) {
    // Read iterators
    AbstractLinAlgPack::sparse_bounds_itr
      eLU_itr( eL->begin(), eL->end(), eL->offset()
           , eU->begin(), eU->end(), eU->offset(), inf_bnd );
    // written iterators
    DVectorSlice::iterator
      b_itr		= loqo_b.begin(),
      r_itr		= loqo_r.begin();
    DenseLinAlgPack::VectorSliceTmpl<int>::iterator
      b_stat_itr  = loqo_b_stat.begin();
    // loop
    for( int k = 1; !eLU_itr.at_end(); ++k, ++eLU_itr, ++b_itr, ++r_itr, ++b_stat_itr, ++num_inequal )
    {
      const size_type j = eLU_itr.indice();
      if(eLU_itr.lbound() > -inf_bnd) {
        *b_itr = eLU_itr.lbound();
        *r_itr = eLU_itr.ubound() >= inf_bnd ? real_big : eLU_itr.ubound() - eLU_itr.lbound();
        *b_stat_itr = j; // We need to make A(k,:) = [ +op(E)(j,:), -b(j) ]
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPT( !( eLU_itr.ubound() < +inf_bnd ) );
        *b_itr = -eLU_itr.ubound();
        *r_itr = eLU_itr.lbound() <= -inf_bnd ? real_big : - eLU_itr.lbound() + eLU_itr.ubound();
        *b_stat_itr = -j; // We need to make A(k,:) = [ -op(E)(j,:), +b(j) ]
      }
    }
  }
  if(F) {
    LinAlgOpPack::V_StV( &loqo_b(num_inequal+1,num_inequal+m_eq), -1.0, *f );
    loqo_r(num_inequal+1,num_inequal+m_eq) = 0.0;
  }

  //
  // Setup the QP dimensions
  //

  loqo_lp->n = nd+1;
  loqo_lp->m = num_inequal + m_eq;

  //
  // Setup loqo_c, loqo_l and loqo_u
  //

  // LOQO's c vector storage
  MALLOC( loqo_lp->c, nd+1, double );
  DVectorSlice loqo_c( loqo_lp->c, nd+1 );
  loqo_c(1,nd) = g;
  loqo_c(nd+1) = bigM();

  // LOQO's l vector storage
  MALLOC( loqo_lp->l, nd+1, double );
  DVectorSlice loqo_l( loqo_lp->l, nd+1 );
  std::fill( loqo_l.begin(), loqo_l.end(), -real_big );
  {
    SpVectorSlice::const_iterator
      dL_itr = dL.begin(),
      dL_end = dL.end();
    for( ; dL_itr != dL_end; ++dL_itr )
      loqo_l( dL_itr->indice() + dL.offset() ) = dL_itr->value();
  }
  loqo_l(nd+1) = etaL;

  // LOQO's u vector storage
  MALLOC( loqo_lp->u, nd+1, double );
  DVectorSlice loqo_u( loqo_lp->u, nd+1 );
  std::fill( loqo_u.begin(), loqo_u.end(), +real_big );
  {
    SpVectorSlice::const_iterator
      dU_itr = dU.begin(),
      dU_end = dU.end();
    for( ; dU_itr != dU_end; ++dU_itr )
      loqo_u( dU_itr->indice() + dU.offset() ) = dU_itr->value();
  }
  loqo_u(nd+1) = +real_big;
  
  //
  // Setup the objective and constraint matrices (using strategy interface).
  //

  init_hess_jacob().init_hess_jacob(
    G,bigM(),E,trans_E,b,&loqo_b_stat[0],num_inequal,F,trans_F,f
    ,loqo_lp);

  //
  // Setup the starting point
  //

  MALLOC( loqo_lp->x, nd+1, double );
  DVectorSlice loqo_x( loqo_lp->x, nd+1 );
  loqo_x(1,nd) = *d;
  loqo_x(nd+1) = *eta;

  //
  // Set some control parameters
  //
  
//	strcpy( loqo_lp->name, "loqo_qp" );
  loqo_lp->quadratic = 1;
  loqo_lp->convex    = 1;
  switch( olevel ) {
    case PRINT_NONE:
      loqo_lp->verbose = 0;
      break;
    case PRINT_BASIC_INFO:
      loqo_lp->verbose = 1;
      break;
    case PRINT_ITER_SUMMARY:
      loqo_lp->verbose = 2;
      break;
    case PRINT_ITER_STEPS:
      loqo_lp->verbose = 3;
      break;
    case PRINT_ITER_ACT_SET:
      loqo_lp->verbose = 4;
      break;
    case PRINT_ITER_VECTORS:
      loqo_lp->verbose = 5;
      break;
    case PRINT_EVERY_THING:
      loqo_lp->verbose = 6;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  //
  // Solve the QP
  //

  if( out && olevel >= PRINT_BASIC_INFO ) {
    *out << "\nSolving QP using LOQO ...\n";
    out->flush();
  }
  
  const int loqo_status = solvelp(loqo_lp);

  if( out && olevel >= PRINT_BASIC_INFO ) {
    *out << "\nLOQO returned status = " << loqo_status << "\n";
  }

  //
  // Map the solution to the output arguments
  //

  TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp->x  ) );
  DVectorSlice loqo_x_sol( loqo_lp->x, nd+1 );

  // d
  *d    = loqo_x_sol(1,nd);

  // eta
  *eta  = loqo_x_sol(nd+1);

  // obj_d
  if(obj_d)
    *obj_d = loqo_lp->primal_obj - (*eta + 0.5 * (*eta)*(*eta)) * bigM();

  // nu
  if(nu) {
    nu->resize(nd,nd);
    TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp->z  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp->s  ) );
    const DVectorSlice
      loqo_z(loqo_lp->z,loqo_lp->n),   // Multipliers for l - x <= 0
      loqo_s(loqo_lp->s,loqo_lp->n);   // Multipliers for x - u <= 0
    DVectorSlice::const_iterator
      z_itr = loqo_z.begin(),
      s_itr = loqo_s.begin();
    typedef SpVector::element_type ele_t;
    for( size_type i = 1; i <= nd; ++i, ++z_itr, ++s_itr ) {
      if( *z_itr > *s_itr && *z_itr >= nonbinding_lag_mult() ) {
        // Lower bound is active
        nu->add_element(ele_t(i,-(*z_itr)));
      }
      else if( *s_itr > *z_itr && *s_itr >= nonbinding_lag_mult() ) {
        // Upper bound is active
        nu->add_element(ele_t(i,+(*s_itr)));
      }
    }
    // We could look at z(nd+1) and s(nd+1) for the value of kappa?
    nu->assume_sorted(true);
  }

  // mu
  if(mu) {
    mu->resize(m_in,num_inequal);
    DenseLinAlgPack::VectorSliceTmpl<int>::iterator
      b_stat_itr  = loqo_b_stat.begin();
    TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp->v  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp->q  ) );
    const DVectorSlice
      loqo_v(loqo_lp->v,loqo_lp->m),   // Multipliers for b <= A*x
      loqo_q(loqo_lp->q,loqo_lp->m);   // Multipliers for A*x <= b + r
    DVectorSlice::const_iterator
      v_itr = loqo_v.begin(),
      q_itr = loqo_q.begin();
    // loop
    typedef SpVector::element_type ele_t;
    for( size_type k = 1; k <= num_inequal; ++k, ++b_stat_itr, ++v_itr, ++q_itr ) {
      const int j = *b_stat_itr;
      if( *v_itr > *q_itr && *v_itr >= nonbinding_lag_mult() ) {
        // Lower bound is active
        if( j < 0 ) // We had to flip this since it was really and upper bound
          mu->add_element(ele_t(-j,+(*v_itr)));
        else // This really was a lower bound
          mu->add_element(ele_t(+j,-(*v_itr)));
      }
      else if( *q_itr > *v_itr && *q_itr >= nonbinding_lag_mult() ) {
        // Upper bound is active
        mu->add_element(ele_t(+j,+(*q_itr)));
      }
    }
  }

  // Ed
  if(Ed) {
    LinAlgOpPack::V_MtV( Ed, *E, trans_E, *d );
  }

  // lambda
  if(lambda) {
    TEUCHOS_TEST_FOR_EXCEPT( !(  loqo_lp->y  ) );
    const DVectorSlice
      loqo_y(loqo_lp->y,loqo_lp->m);         // Multipliers for equalities
    DVectorSlice::const_iterator
      y_itr = loqo_y.begin() + num_inequal;  // Get iterators to equalities
    DVectorSlice::iterator
      lambda_itr = lambda->begin();
    // loop
    for( size_type k = 1; k <= m_eq; ++k, ++y_itr, ++lambda_itr ) {
      *lambda_itr = -(*y_itr);
    }
  }

  // Fd
  if(Fd) {
    LinAlgOpPack::V_MtV( Fd, *F, trans_F, *d );
  }

  //
  // Setup the QP statistics
  //

  QPSolverStats::ESolutionType solution_type = QPSolverStats::OPTIMAL_SOLUTION; // Assume this?
  switch( loqo_status ) { // I had to find this out by trial and error!
      case 0:
      solution_type = QPSolverStats::OPTIMAL_SOLUTION;
      break;
    case 2:
      solution_type = QPSolverStats::DUAL_FEASIBLE_POINT;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  qp_stats_.set_stats(
    solution_type, QPSolverStats::CONVEX
    ,loqo_lp->iter, QPSolverStats::NOT_KNOWN, QPSolverStats::NOT_KNOWN
    ,false, *eta > 0.0 );

  //
  // Clean up dynamically allocated memory for LOQO
  //

  inv_clo();          // frees memory associated with matrix factorization
  closelp(loqo_lp);   // frees all allocated arrays with free(...).

  return qp_stats_.solution_type();

}

}	// end namespace ConstrainedOptPack

#endif // CONSTRAINED_OPTIMIZATION_PACK_USE_LOQO
