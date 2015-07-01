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

#include "ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorSpaceFactory.hpp"
#include "AbstractLinAlgPack_SortByDescendingAbsValue.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_VectorSpaceSerial.hpp"
#include "AbstractLinAlgPack_sparse_bounds.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "ProfileHackPack_profile_hack.hpp"

namespace ConstrainedOptPack {

QPSolverRelaxedQPSchur::QPSolverRelaxedQPSchur(
  const init_kkt_sys_ptr_t&  init_kkt_sys
  ,const constraints_ptr_t&  constraints
  ,value_type			max_qp_iter_frac
  ,value_type         max_real_runtime
  ,QPSchurPack::ConstraintsRelaxedStd::EInequalityPickPolicy
                      inequality_pick_policy
  ,ELocalOutputLevel	print_level
  ,value_type 		bounds_tol
  ,value_type 		inequality_tol
  ,value_type 		equality_tol
  ,value_type			loose_feas_tol
  ,value_type			dual_infeas_tol
  ,value_type			huge_primal_step
  ,value_type			huge_dual_step
  ,value_type			bigM
  ,value_type			warning_tol
  ,value_type			error_tol
  ,size_type          iter_refine_min_iter
  ,size_type          iter_refine_max_iter
  ,value_type         iter_refine_opt_tol
  ,value_type         iter_refine_feas_tol
  ,bool               iter_refine_at_solution
  ,value_type         pivot_warning_tol
  ,value_type         pivot_singular_tol
  ,value_type         pivot_wrong_inertia_tol
  ,bool               add_equalities_initially
  )
  :init_kkt_sys_(init_kkt_sys)
  ,constraints_(constraints)
  ,max_qp_iter_frac_(max_qp_iter_frac)
  ,max_real_runtime_(max_real_runtime)
  ,inequality_pick_policy_(inequality_pick_policy)
  ,print_level_(print_level)
  ,bounds_tol_(bounds_tol)
  ,inequality_tol_(inequality_tol)
  ,equality_tol_(equality_tol)
  ,loose_feas_tol_(loose_feas_tol)
  ,dual_infeas_tol_(dual_infeas_tol)
  ,huge_primal_step_(huge_primal_step)
  ,huge_dual_step_(huge_dual_step)
  ,bigM_(bigM)
  ,warning_tol_(warning_tol)
  ,error_tol_(error_tol)
  ,iter_refine_min_iter_(iter_refine_min_iter)
  ,iter_refine_max_iter_(iter_refine_max_iter)
  ,iter_refine_opt_tol_(iter_refine_opt_tol)
  ,iter_refine_feas_tol_(iter_refine_feas_tol)
  ,iter_refine_at_solution_(iter_refine_at_solution)
  ,pivot_warning_tol_(pivot_warning_tol)
  ,pivot_singular_tol_(pivot_singular_tol)
  ,pivot_wrong_inertia_tol_(pivot_wrong_inertia_tol)
  ,add_equalities_initially_(add_equalities_initially)
{}

QPSolverRelaxedQPSchur::~QPSolverRelaxedQPSchur()
{
  this->release_memory();
}

// Overridden from QPSolverRelaxed

QPSolverStats
QPSolverRelaxedQPSchur::get_qp_stats() const
{
  return qp_stats_;
}

void QPSolverRelaxedQPSchur::release_memory()
{
  // Nothing to release!
}

QPSolverStats::ESolutionType
QPSolverRelaxedQPSchur::imp_solve_qp(
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
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  using LinAlgOpPack::V_mV;
  typedef QPSchurPack::ConstraintsRelaxedStd constr_t;

#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "QPSolverRelaxedQPSchur::imp_solve_qp(...)" );
#endif

  const size_type
    nd   = g.dim(),
    m_in = E ? BLAS_Cpp::rows(E->rows(),E->cols(),trans_E) : 0,
    m_eq = F ? BLAS_Cpp::rows(F->rows(),F->cols(),trans_F) : 0;

  VectorDenseEncap  g_de(g);

  VectorSpace::space_ptr_t
    space_d_eta = d->space().small_vec_spc_fcty()->create_vec_spc(nd+1); // ToDo: Generalize!

  // ///////////////////////////////
  // Setup the initial KKT system

  InitKKTSystem::i_x_free_t     i_x_free;
  InitKKTSystem::i_x_fixed_t    i_x_fixed;
  InitKKTSystem::bnd_fixed_t    bnd_fixed;
  InitKKTSystem::j_f_decomp_t   j_f_decomp;
  size_type n_R_tmp;
  init_kkt_sys().initialize_kkt_system(
    g,G,etaL,dL,dU,F,trans_F,f,d,nu
    ,&n_R_tmp,&i_x_free,&i_x_fixed,&bnd_fixed,&j_f_decomp
    ,&b_X_,&Ko_,&fo_ );
  const size_type
    n_R = n_R_tmp,
    n_X = nd + 1 - n_R; // fixed variables in d and eta
  TEUCHOS_TEST_FOR_EXCEPT( !(  i_x_free.size() == 0 || i_x_free.size()  >= n_R  ) );  // Todo: Make an exception!
  TEUCHOS_TEST_FOR_EXCEPT( !(  i_x_fixed.size() >= n_X  ) );  // Todo: Make an exception!
  TEUCHOS_TEST_FOR_EXCEPT( !(  bnd_fixed.size() >= n_X  ) ); // Todo: Make and exception!

  // //////////////////////////////
  // Initialize constraints object

  // Setup j_f_undecomp
  const bool all_f_undecomp = F ? j_f_decomp.size() == 0 : true;
  const size_type
    m_undecomp = F ? f->dim() - j_f_decomp.size() : 0;
  typedef std::vector<size_type> j_f_undecomp_t;
  j_f_undecomp_t j_f_undecomp;
  if( m_undecomp && !all_f_undecomp ) {
    j_f_undecomp.resize(m_undecomp);
    // Create a full lookup array to determine if a constraint
    // is decomposed or not.  We need this to fill the array
    // j_f_undecomp[] (which is sorted).
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement this!
  }

  // initialize constraints object
  constraints_->initialize(
    space_d_eta
    ,etaL,dL,dU,E,trans_E,b,eL,eU,F,trans_F,f
    ,m_undecomp, m_undecomp && !all_f_undecomp ? &j_f_undecomp[0] : NULL
    ,Ed
    ,!add_equalities_initially()  // If we add equalities the the schur complement intially
                                  // then we don't need to check if they are violated.
    );
  // ToDo: Add j_f_decomp to the above constraints class!

  // ///////////////////////////
  // Initialize the QP object

  // g_relaxed_ = [ g; bigM ]
  g_relaxed_.resize(nd+1);
  g_relaxed_(1,nd) = g_de();
  g_relaxed_(nd+1) = bigM();

  // G_relaxed_ = [ G, zeros(...); zeros(...), bigM ]
  bigM_vec_.initialize(1); // dim == 1
  bigM_vec_ = bigM();	
  G_relaxed_.initialize(
    Teuchos::rcp(&dyn_cast<const MatrixSymOpNonsing>(G),false)
    ,Teuchos::rcp(&bigM_vec_,false)
    ,space_d_eta
    );
  
  qp_.initialize(
    g_relaxed_(),G_relaxed_,NULL
    ,n_R, i_x_free.size() ? &i_x_free[0] : NULL
    ,&i_x_fixed[0],&bnd_fixed[0]
    ,b_X_(),*Ko_,fo_(),constraints_.get()
    ,out,test_what==RUN_TESTS,warning_tol(),error_tol()
    ,int(olevel)>=int(PRINT_ITER_VECTORS)
    );

  // ///////////////////////////////////////////////////////
  // Setup for a warm start (changes to initial KKT system)

  typedef std::vector<int> 					ij_act_change_t;
  typedef std::vector<EBounds>				bnds_t;
  size_type			num_act_change = 0; // The default is a cold start
  const size_type     max_num_act_change = 2*nd;
  ij_act_change_t		ij_act_change(max_num_act_change);
  bnds_t				bnds(max_num_act_change);
  // Go ahead and add the equality constraints.  If these are linearly
  // dependent let's hope that QPSchur can handle this and still do a
  // good job of things.  This is a scary thing to do!
  if( m_eq && add_equalities_initially() ) {
    for( size_type j = 1; j <= m_eq; ++j ) {
      ij_act_change[num_act_change] = (nd + 1) + m_in + j;
      bnds[num_act_change]          = EQUALITY;
      ++num_act_change;
    }
  }
  // We will add fixed (EQUALITY) variable bounds to the initial active set
  // (if it is not already an intially fixed variable).  If fixing a variable
  // causes the KKT system to become singular then we are in real trouble!
  // We should add these eairly on since they will not be freed.
  if( dL ) {
    const QPSchurPack::QP::x_init_t &x_init = qp_.x_init();
    const value_type inf_bnd = this->infinite_bound();
    VectorDenseEncap dL_de(*dL);
    VectorDenseEncap dU_de(*dU);
    // read iterators
    AbstractLinAlgPack::sparse_bounds_itr
      dLU_itr( dL_de().begin(), dL_de().end()
          ,dU_de().begin(), dU_de().end()
          ,inf_bnd );
    for( ; !dLU_itr.at_end(); ++dLU_itr ) {
      if( dLU_itr.lbound() == dLU_itr.ubound() && x_init(dLU_itr.index()) == FREE ) {
        ij_act_change[num_act_change] = dLU_itr.index();
        bnds[num_act_change]          = EQUALITY;
        ++num_act_change;
      }
    }
  }
  // Add inequality constriants to the list from nu and mu
  if( ( nu && nu->nz() ) || ( m_in && mu->nz() ) ) {
    //
    // Setup num_act_change, ij_act_change, and bnds for a warm start!
    //
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
    // Now add the inequality constraints in decreasing order (if they are
    // not already initially fixed variables)
    const QPSchurPack::QP::x_init_t &x_init = qp_.x_init();
    if(gamma.nz()) {
      const SpVector::difference_type o = gamma.offset();
      for( SpVector::const_iterator itr = gamma.begin(); itr != gamma.end(); ++itr ) {
        const size_type i =  itr->index() + o;
        if( i <= nd && x_init(i) != FREE )
          continue; // This variable is already initially fixed
        // This is not an initially fixed variable so add it
        ij_act_change[num_act_change] = i;
        bnds[num_act_change]
          = itr->value() > 0.0 ? UPPER : LOWER;
        ++num_act_change;
      }
    }
  }
  // We need to loop through x_init() and nu() in order and look for variables
  // that are initially fixed in x_init() but are not present in nu().  For these
  // variables we need to free them in ij_act_change[].
  if( dL && nu->nz() ) {
    QPSchurPack::QP::x_init_t::const_iterator
      x_init_itr = qp_.x_init().begin();
    VectorDenseEncap nu_de(*nu);
    DVectorSlice::const_iterator
      nu_itr = nu_de().begin();
    for( size_type i = 1; i <= nd; ++i, ++x_init_itr, ++nu_itr ) {
      if( *x_init_itr != FREE && *x_init_itr != EQUALITY ) {
        // This is an initially fixed upper or lower bound.
        // Look for lagrange multiplier stating that it is
        // still fixed.
        if( *nu_itr != 0.0 ) {
          // This active bound is present but lets make sure
          // that it is still the same bound
          if( ( *x_init_itr == LOWER && *nu_itr > 0 )
            || ( *x_init_itr == UPPER && *nu_itr < 0 ) )
          {
            // The bound has changed from upper to lower or visa-versa!
            ij_act_change[num_act_change] = i;
            bnds[num_act_change]
              = *nu_itr > 0.0 ? UPPER : LOWER;
            ++num_act_change;
          }
        }
        else {
          // This initially fixed variable is not fixed in nu so lets free it!
          ij_act_change[num_act_change] = -i;
          bnds[num_act_change]          = FREE;
          ++num_act_change;
        }
      }
    }
  }
  // Consider the relaxation variable!
  if( *eta > etaL) {
    ij_act_change[num_act_change] = -int(nd+1);
    bnds[num_act_change]          = FREE;
    ++num_act_change;
  }		

  // Set the output level
  QPSchur::EOutputLevel qpschur_olevel;
  switch( print_level() ) {
    case USE_INPUT_ARG: {
      // Use the input print level
      switch( olevel ) {
        case PRINT_NONE:
          qpschur_olevel = QPSchur::NO_OUTPUT;
          break;
        case PRINT_BASIC_INFO:
          qpschur_olevel = QPSchur::OUTPUT_BASIC_INFO;
          break;
        case PRINT_ITER_SUMMARY:
          qpschur_olevel = QPSchur::OUTPUT_ITER_SUMMARY;
          break;
        case PRINT_ITER_STEPS:
          qpschur_olevel = QPSchur::OUTPUT_ITER_STEPS;
          break;
        case PRINT_ITER_ACT_SET:
        case PRINT_ITER_VECTORS:
          qpschur_olevel = QPSchur::OUTPUT_ACT_SET;
          break;
        case PRINT_EVERY_THING:
          qpschur_olevel = QPSchur::OUTPUT_ITER_QUANTITIES;
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true);
      }
      break;
    }
    case NO_OUTPUT:
      qpschur_olevel = QPSchur::NO_OUTPUT;
      break;
    case OUTPUT_BASIC_INFO:
      qpschur_olevel = QPSchur::OUTPUT_BASIC_INFO;
      break;
    case OUTPUT_ITER_SUMMARY:
      qpschur_olevel = QPSchur::OUTPUT_ITER_SUMMARY;
      break;
    case OUTPUT_ITER_STEPS:
      qpschur_olevel = QPSchur::OUTPUT_ITER_STEPS;
      break;
    case OUTPUT_ACT_SET:
      qpschur_olevel = QPSchur::OUTPUT_ACT_SET;
      break;
    case OUTPUT_ITER_QUANTITIES:
      qpschur_olevel = QPSchur::OUTPUT_ITER_QUANTITIES;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  //
  // Set options for ConstraintsRelaxedStd.
  // 
  if( bounds_tol() > 0.0 )
    constraints_->bounds_tol(bounds_tol());
  if( inequality_tol() > 0.0 )
    constraints_->inequality_tol(inequality_tol());
  if( equality_tol() > 0.0 )
    constraints_->equality_tol(equality_tol());
  constraints_->inequality_pick_policy(inequality_pick_policy());

  //
  // Set options for QPSchur.
  // 
  qp_solver_.max_iter(static_cast<index_type>(max_qp_iter_frac()*nd) );
  qp_solver_.max_real_runtime( max_real_runtime() );
  qp_solver_.feas_tol( constraints_->bounds_tol() );	// Let's assume the bound tolerance is the tightest
  if(loose_feas_tol() > 0.0)
    qp_solver_.loose_feas_tol( loose_feas_tol() );
  else
    qp_solver_.loose_feas_tol( 10.0 * qp_solver_.feas_tol() );
  if(dual_infeas_tol() > 0.0)
    qp_solver_.dual_infeas_tol( dual_infeas_tol() );
  if(huge_primal_step() > 0.0)
    qp_solver_.huge_primal_step( huge_primal_step() );
  if(huge_dual_step() > 0.0)
    qp_solver_.huge_dual_step( huge_dual_step() );
  qp_solver_.set_schur_comp( QPSchur::schur_comp_ptr_t( &schur_comp_, false ) );
  qp_solver_.warning_tol( warning_tol() );
  qp_solver_.error_tol( error_tol() );
  qp_solver_.iter_refine_min_iter( iter_refine_min_iter() );
  qp_solver_.iter_refine_max_iter( iter_refine_max_iter() );
  qp_solver_.iter_refine_opt_tol( iter_refine_opt_tol() );
  qp_solver_.iter_refine_feas_tol( iter_refine_feas_tol() );
  qp_solver_.iter_refine_at_solution( iter_refine_at_solution() );
  qp_solver_.pivot_tols(
    MatrixSymAddDelUpdateable::PivotTolerances(
      pivot_warning_tol(), pivot_singular_tol(), pivot_wrong_inertia_tol()
      ));
  
  //
  // Solve the QP with QPSchur
  // 
  DVector _x(nd+1);		// solution vector [ d; eta ]
  SpVector _mu;			// lagrange multipliers for variable bounds [ nu; kappa ]
  SpVector _lambda_breve;	// solution for extra general constraints [ mu; lambda ]
  size_type qp_iter = 0, num_adds = 0, num_drops = 0;
  QPSchur::ESolveReturn
    solve_returned
      = qp_solver_.solve_qp(
        qp_
        ,num_act_change, num_act_change ? &ij_act_change[0] : NULL
        ,num_act_change ? &bnds[0] : NULL
        ,out, qpschur_olevel
        ,test_what==RUN_TESTS ? QPSchur::RUN_TESTS : QPSchur::NO_TESTS
        ,&_x(), &_mu, NULL, &_lambda_breve
        ,&qp_iter, &num_adds, &num_drops
        );

  // Set the solution

  // d
  (VectorDenseMutableEncap(*d))() = _x(1,nd);
  // nu
  if( nu ) {
    *nu = 0.0;
    const SpVector::difference_type o = _mu.offset();
    if( _mu.nz() ) {
      for(SpVector::const_iterator _mu_itr = _mu.begin(); _mu_itr != _mu.end(); ++_mu_itr) {
        typedef SpVector::element_type ele_t;
        if( _mu_itr->index() + o <= nd ) // don't add multiplier for eta <= etaL
          nu->set_ele( _mu_itr->index() + o, _mu_itr->value() );
      }
    }
  }
  // mu, lambda	
  if( m_in || m_eq ) {
    *eta = _x(nd+1);	// must be non-null
    *mu = 0.0;
    const SpVector::difference_type o = _lambda_breve.offset();
    if(_lambda_breve.nz()) {
      for(SpVector::const_iterator itr = _lambda_breve.begin();
        itr != _lambda_breve.end();
        ++itr)
      {
        typedef SpVector::element_type ele_t;
        if( itr->index() + o <= m_in ) {
          mu->set_ele( itr->index() + o, itr->value() );
        }
        else {
          lambda->set_ele( itr->index() + o - m_in, itr->value() );
        }
      }
    }
  }
  // obj_d (This could be updated within QPSchur in the future)
  if(obj_d) {
    // obj_d = g'*d + 1/2 * d' * G * g
    *obj_d = AbstractLinAlgPack::dot(g,*d)
      + 0.5 * AbstractLinAlgPack::transVtMtV(*d,G,BLAS_Cpp::no_trans,*d);
  }
  // Ed
  if(Ed && E) {
    switch(constraints_->inequality_pick_policy()) {
      case constr_t::ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY:
        if(solve_returned == QPSchur::OPTIMAL_SOLUTION)
          break; // Ed already computed (see ConstraintsRelaxedStd::pick_violated())
      case constr_t::ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY:
        break; // Ed already computed (see ConstraintsRelaxedStd::pick_violated())
      default:
        // We need to compute Ed
        LinAlgOpPack::V_MtV( Ed, *E, trans_E, *d );
    }
  }
  // Fd (This could be updated within ConstraintsRelaxedStd in the future)
  if(Fd) {
    LinAlgOpPack::V_MtV( Fd, *F, trans_F, *d );
  }
  // Set the QP statistics
  QPSolverStats::ESolutionType solution_type;
  QPSolverStats::EConvexity convexity = QPSolverStats::CONVEX;
  switch( solve_returned ) {
    case QPSchur::OPTIMAL_SOLUTION:
      solution_type = QPSolverStats::OPTIMAL_SOLUTION;
      break;
    case QPSchur::MAX_ITER_EXCEEDED:
      solution_type = QPSolverStats::DUAL_FEASIBLE_POINT;
      break;
    case QPSchur::MAX_RUNTIME_EXEEDED_FAIL:
      solution_type = QPSolverStats::SUBOPTIMAL_POINT;
      break;
    case QPSchur::MAX_RUNTIME_EXEEDED_DUAL_FEAS:
      solution_type = QPSolverStats::DUAL_FEASIBLE_POINT;
      break;
    case QPSchur::MAX_ALLOWED_STORAGE_EXCEEDED:
      solution_type = QPSolverStats::DUAL_FEASIBLE_POINT;
      break;
    case QPSchur::INFEASIBLE_CONSTRAINTS:
    case QPSchur::NONCONVEX_QP:
      convexity = QPSolverStats::NONCONVEX;
    case QPSchur::DUAL_INFEASIBILITY:
    case QPSchur::SUBOPTIMAL_POINT:
      solution_type = QPSolverStats::SUBOPTIMAL_POINT;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  qp_stats_.set_stats(
    solution_type,convexity,qp_iter,num_adds,num_drops
    , num_act_change > 0 || n_X > 1, *eta > 0.0 );

  return qp_stats_.solution_type();
}

}	// end namespace ConstrainedOptPack
