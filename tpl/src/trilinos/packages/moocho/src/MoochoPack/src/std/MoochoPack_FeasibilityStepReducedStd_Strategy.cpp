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

#include "MoochoPack_FeasibilityStepReducedStd_Strategy.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

FeasibilityStepReducedStd_Strategy::FeasibilityStepReducedStd_Strategy(
  const quasi_range_space_step_ptr_t   &quasi_range_space_step
  ,const qp_solver_ptr_t               &qp_solver
  ,const qp_tester_ptr_t               &qp_tester
  ,EQPObjective                        qp_objective
  ,EQPTesting                          qp_testing
  )
  :quasi_range_space_step_(quasi_range_space_step)
  ,qp_solver_(qp_solver)
  ,qp_tester_(qp_tester)
  ,qp_objective_(qp_objective)
  ,qp_testing_(qp_testing)
  ,dl_iq_(dl_name)
  ,du_iq_(du_name)
  ,current_k_(-1)
{}

bool FeasibilityStepReducedStd_Strategy::compute_feasibility_step(
  std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
  ,const Vector& xo, const Vector& c_xo, VectorMutable* w
    )
{
  using Teuchos::dyn_cast;

/* Todo: UPdate below code!

  // problem dimensions
  const size_type
    n = algo->nlp().n(),
    m = algo->nlp().m(),
    r  = s->equ_decomp().size();

  // Compute the quasi-range space step Ywy
  Workspace<value_type> Ywy_ws(wss,xo.size());
  DVectorSlice                Ywy(&Ywy_ws[0],Ywy_ws.size());
  if(!quasi_range_space_step().solve_quasi_range_space_step(
    out,olevel,algo,s,xo,c_xo,&Ywy ))
    return false;

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\n||Ywy||2     = " << DenseLinAlgPack::norm_2(Ywy);
    out << std::endl;
  }
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nYwy = \n" << Ywy;
  }

  //
  // Set the bounds on the null space QP subproblem:
  //
  // d_bounds_k.l <= (xo - x_k) + (1-eta) * Ywy + Z*wz <= d_bounds_k.u
  // =>
  // bl <= Z*wz - eta * Ywy <= bu
  //
  // where:   bl = d_bounds_k.l - (xo - x_k) - Ywy
  //          bu = d_bounds_k.u - (xo - x_k) - Ywy
  //
  // Above, fix the variables that are at an active bound as equalities
  // to maintain the same active-set.
  //
  const SparseBounds
    &d_bounds = d_bounds_(*s).get_k(0);
  const SpVectorSlice
    &dl = d_bounds.l,
    &du = d_bounds.u;
  const DVector
    &x_k = s->x().get_k(0).v();
  const SpVector
    &nu_k = s->nu().get_k(0);
  TEUCHOS_TEST_FOR_EXCEPT( !( nu_k.is_sorted() ) );
  SpVector bl(n,n), bu(n,n);
  sparse_bounds_itr
    d_bnds_itr(dl.begin(),dl.end(),dl.offset(),du.begin(),du.end(),du.offset());
  SpVector::const_iterator
    nu_itr     = nu_k.begin(),
    nu_end     = nu_k.end();
  for( ; !d_bnds_itr.at_end(); ++d_bnds_itr ) {
    typedef SpVectorSlice::element_type ele_t;
    const size_type i = d_bnds_itr.indice();
    while( nu_itr != nu_end && nu_itr->indice() + nu_k.offset() < i )
      ++nu_itr;
    if( nu_itr != nu_end && nu_itr->indice() + nu_k.offset() == i ) {
      const value_type
        act_bnd = nu_itr->value() > 0.0 ? d_bnds_itr.ubound() : d_bnds_itr.lbound();
      bl.add_element(ele_t( i, act_bnd - xo(i) + x_k(i) - Ywy(i) ));
      bu.add_element(ele_t( i, act_bnd - xo(i) + x_k(i) - Ywy(i) ));
    }
    else {
      if( d_bnds_itr.lbound() != -d_bnds_itr.big_bnd() )
        bl.add_element(ele_t(i,d_bnds_itr.lbound()  - xo(i) + x_k(i) - Ywy(i) ));
      if( d_bnds_itr.ubound() != +d_bnds_itr.big_bnd() )
        bu.add_element(ele_t(i, d_bnds_itr.ubound() - xo(i) + x_k(i) - Ywy(i) ));
    }
  }
  bl.assume_sorted(true);
  bu.assume_sorted(true);
  //
  // Setup the objective function for the null space QP subproblem
  //
  // 
  // OBJ_MIN_FULL_STEP
  //    min 1/2 * (Y*wy + Z*wz)'*(Y*wy + Z*wz) = 1/2*wy'*Y'*Y*wy + (Z'*Y*wy)'*wz + 1/2*wz'*(Z'*Z)*wz
  //    => grad = (Z'*Y*wy), Hess = Z'*Z
  //
  // OBJ_MIN_WZ
  //    min 1/2 * wz'*wz => grad = 0, Hess = I
  //
  // OBJ_RSQP
  //    min qp_grad_k'*wz + 1/2 * wz'*rHL_k*wz
  //    => grad = qp_grad, Hess = rHL_k
  //
  const MatrixOp
    &Z_k = s->Z().get_k(0);
  if( current_k_ != s->k() ) {
    if( qp_objective() != OBJ_RSQP )
      grad_store_.resize(n-r);
    if( qp_objective() == OBJ_MIN_FULL_STEP )
      Hess_store_.resize(n-r+1,n-r+1);
  }
  DVectorSlice grad;
  switch(qp_objective())
  {
      case OBJ_MIN_FULL_STEP: // grad = (Z'*Ywy), Hess = Z'*Z
    {
      grad.bind( grad_store_() );
      if( current_k_ != s->k() ) {
        // grad = (Z'*Ywy)
        LinAlgOpPack::V_MtV( &grad, Z_k, BLAS_Cpp::trans, Ywy );
        // Hess = Z'*Z
        DMatrixSliceSym S(Hess_store_(2,n-r+1,1,n-r),BLAS_Cpp::lower); // Must be strictly lower triangular here!
        Z_k.syrk( BLAS_Cpp::trans, 1.0, 0.0, &S ); // S = 1.0*Z'*Z + 0.0*S
        MatrixSymPosDefCholFactor
          *H_ptr = NULL;
        if( Hess_ptr_.get() == NULL || dynamic_cast<const MatrixSymPosDefCholFactor*>(Hess_ptr_.get()) == NULL )
          Hess_ptr_ = new MatrixSymPosDefCholFactor;
        H_ptr = const_cast<MatrixSymPosDefCholFactor*>(dynamic_cast<const MatrixSymPosDefCholFactor*>(Hess_ptr_.get()));
        TEUCHOS_TEST_FOR_EXCEPT( !( H_ptr ) ); // Should not be null!
        H_ptr->init_setup(
          &Hess_store_()  // The original matrix is stored in the lower triangular part (below diagonal)!
          ,NULL           // Nothing to deallocate
          ,n-r
          ,true           // maintain the original factor
          ,false          // don't maintain the factor
          ,true           // allow the factor to be computed if needed
          ,true           // Set the view
          ,1.0            // Scale the matrix by one
          );
      }
      break;
    }
      case OBJ_MIN_NULL_SPACE_STEP: // grad = 0, Hess = I
    {
      grad.bind( grad_store_() );
      MatrixSymIdent
        *H_ptr = NULL;
      if( Hess_ptr_.get() == NULL || dynamic_cast<const MatrixSymIdent*>(Hess_ptr_.get()) == NULL )
        Hess_ptr_ = new MatrixSymIdent;
      if( current_k_ != s->k() ) {
        H_ptr = const_cast<MatrixSymIdent*>(dynamic_cast<const MatrixSymIdent*>(Hess_ptr_.get()));
        TEUCHOS_TEST_FOR_EXCEPT( !( H_ptr ) ); // Should not be null!
        H_ptr->init_setup(n-r,1.0);
        grad = 0.0;
      }
      break;
    }
      case OBJ_RSQP: // grad = qp_grad, Hess = rHL_k
    {
      grad.bind( s->qp_grad().get_k(0)() );
      Hess_ptr_ = Hess_ptr_t( &s->rHL().get_k(0), false ); // don't delete memory!
      break;
    }
      defaut:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Not a valid option
  }

  //
  // Solve the null space subproblem
  //

  Workspace<value_type>  wz_ws(wss,n-r),Zwz_ws(wss,n);
  DVectorSlice                 wz(&wz_ws[0],wz_ws.size());
  DVectorSlice                 Zwz(&Zwz_ws[0],Zwz_ws.size());
  value_type                  qp_eta      = 0;

  bool throw_qp_failure = false;

  if( bl.nz() == 0 && bu.nz() == 0 && m-r == 0 ) {
    //
    // Just solve the unconstrainted QP
    //
    // wz = - inv(Hess)*grad
#ifdef _WINDOWS
    const MatrixFactorized &Hess = dynamic_cast<const MatrixFactorized&>(*Hess_ptr_);
#else
    const MatrixFactorized &Hess = dyn_cast<const MatrixFactorized>(*Hess_ptr_);
#endif
    AbstractLinAlgPack::V_InvMtV( &wz, Hess, BLAS_Cpp::no_trans, grad );
    DenseLinAlgPack::Vt_S(&wz,-1.0);
    // Zwz = Z*wz
    LinAlgOpPack::V_MtV( &Zwz, Z_k, BLAS_Cpp::no_trans, wz );
  }
  else {

    //
    // Set the arguments to the QP subproblem
    //

    DVectorSlice 			qp_g		= grad;
    const MatrixOp& 	qp_G 		= *Hess_ptr_;
    const value_type		qp_etaL 	= 0.0;
    SpVectorSlice			qp_dL(NULL,0,0,n-r);	// If nz() == 0 then no simple bounds
    SpVectorSlice			qp_dU(NULL,0,0,n-r);
    const MatrixOp		*qp_E		= NULL;
    BLAS_Cpp::Transp		qp_trans_E	= BLAS_Cpp::no_trans;
    DVectorSlice             qp_b;
    SpVectorSlice			qp_eL(NULL,0,0,n);
    SpVectorSlice			qp_eU(NULL,0,0,n);
    const MatrixOp		*qp_F		= NULL;
    BLAS_Cpp::Transp		qp_trans_F	= BLAS_Cpp::no_trans;
    DVectorSlice				qp_f;
    DVectorSlice				qp_d		= wz;
    SpVector				*qp_nu		= NULL;
    SpVector				*qp_mu		= NULL;
    DVectorSlice				qp_Ed;
    DVectorSlice				qp_lambda;

    SpVector _nu_wz, _nu_Dwz, 	// Possible storage for multiplers for separate inequality
      _nu;				// constriants for wz.
    DVector _Dwz;				// Possible storage for D*wz computed by QP solver?

    //
    // Determine if we can use simple bounds on wz.
    // 
    // If we have a variable reduction null space matrix
    // (with any choice for Y) then:
    // 
    // w = Z*wz + (1-eta) * Y*wy
    // 
    // [ w(dep)   ]  = [ D ] * wz  + (1-eta) * [ Ywy(dep)   ]
    // [ w(indep) ]    [ I ]                   [ Ywy(indep) ]
    // 
    // For a cooridinate decomposition (Y = [ I ; 0 ]) then Ywy(indep) = 0 and
    // in this case the bounds on d(indep) become simple bounds on pz even
    // with the relaxation.
    // 
    // Otherwise, we can not use simple variable bounds and implement the
    // relaxation properly.
    // 

    const ZVarReductMatrix
      *Zvr = dynamic_cast<const ZVarReductMatrix*>( &Z_k );
    Range1D
      indep 	= Zvr ? Zvr->indep() : Range1D(),
      dep		= Zvr ? Zvr->dep()   : Range1D();

    const bool
      use_simple_wz_bounds = ( Zvr!=NULL && DenseLinAlgPack::norm_inf(Ywy(indep))==0.0 );

    if( use_simple_wz_bounds ) {

      // Set simple bound constraints on pz
      qp_dL.bind( bl(indep) );
      qp_dU.bind( bu(indep) );
      qp_nu = &( _nu_wz = s->nu().get_k(0)(indep) );	// warm start?
    
      // Set general inequality constraints for D*pz
      qp_E = &Zvr->D();
      qp_b.bind( Ywy(dep) );
      qp_eL.bind( bl(dep) );
      qp_eU.bind( bu(dep) );
      qp_mu = &( _nu_Dwz = s->nu().get_k(0)(dep) );	// warm start?
      _Dwz.resize(r);
      qp_Ed.bind(_Dwz());	// Part of Zwz will be updated directly!

    }
    else {

      // Set general inequality constraints for Z*pz
      qp_E = &Z_k;
      qp_b.bind( Ywy() );
      qp_eL.bind( bl() );
      qp_eU.bind( bu() );
      qp_mu = &(_nu = s->nu().get_k(0));	// warm start??
      qp_Ed.bind(Zwz);	// Zwz
    }

    // Set the general equality constriants (if they exist)
    DVector q(m-r);
    Range1D undecomp = s->con_undecomp();
    if( m > r ) {
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement when needed!
    }

    // Setup the rest of the arguments

    QPSolverRelaxed::EOutputLevel
      qp_olevel;
    switch( olevel ) {
        case PRINT_NOTHING:
        qp_olevel = QPSolverRelaxed::PRINT_NONE;
        break;
        case PRINT_BASIC_ALGORITHM_INFO:
        qp_olevel = QPSolverRelaxed::PRINT_BASIC_INFO;
        break;
        case PRINT_ALGORITHM_STEPS:
        qp_olevel = QPSolverRelaxed::PRINT_BASIC_INFO;
        break;
        case PRINT_ACTIVE_SET:
        qp_olevel = QPSolverRelaxed::PRINT_ITER_SUMMARY;
        break;
        case PRINT_VECTORS:
        qp_olevel = QPSolverRelaxed::PRINT_ITER_VECTORS;
        break;
        case PRINT_ITERATION_QUANTITIES:
        qp_olevel = QPSolverRelaxed::PRINT_EVERY_THING;
        break;
        default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }

    //
    // Solve the QP
    // 
    const QPSolverStats::ESolutionType
      solution_type =
      qp_solver().solve_qp(
        int(olevel) == int(PRINT_NOTHING) ? NULL : &out
        , qp_olevel
        , algo->algo_cntr().check_results()
        ? QPSolverRelaxed::RUN_TESTS :  QPSolverRelaxed::NO_TESTS
        , qp_g, qp_G, qp_etaL, qp_dL, qp_dU
        , qp_E, qp_trans_E, qp_E ? &qp_b : NULL
        , qp_E ? &qp_eL : NULL, qp_E ? &qp_eU : NULL
        , qp_F, qp_trans_F, qp_F ? &qp_f : NULL
        , NULL
        , &qp_eta, &qp_d
        , qp_nu
        , qp_mu, qp_E ? &qp_Ed : NULL
        , qp_F ? &qp_lambda : NULL, NULL
        );

    //
    // Check the optimality conditions for the QP
    //
    std::ostringstream omsg;
    if(		qp_testing() == QP_TEST
        || ( qp_testing() == QP_TEST_DEFAULT && algo->algo_cntr().check_results() )  )
    {
      if( int(olevel) >= int(PRINT_ALGORITHM_STEPS) ) {
        out	<< "\nChecking the optimality conditions of the reduced QP subproblem ...\n";
      }
      if(!qp_tester().check_optimality_conditions(
        solution_type
        , int(olevel) == int(PRINT_NOTHING) ? NULL : &out
        , int(olevel) >= int(PRINT_VECTORS) ? true : false
        , int(olevel) >= int(PRINT_ITERATION_QUANTITIES) ? true : false
        , qp_g, qp_G, qp_etaL, qp_dL, qp_dU
        , qp_E, qp_trans_E, qp_E ? &qp_b : NULL
        , qp_E ? &qp_eL : NULL, qp_E ? &qp_eU : NULL
        , qp_F, qp_trans_F, qp_F ? &qp_f : NULL
        , NULL
        , &qp_eta, &qp_d
        , qp_nu
        , qp_mu, qp_E ? &qp_Ed : NULL
        , qp_F ? &qp_lambda : NULL, NULL
        ))
      {
        omsg << "\n*** Alert! at least one of the QP optimality conditions did not check out.\n";
        if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
          out << omsg.str();
        }
        throw_qp_failure = true;
      }
    }

    if( solution_type !=  QPSolverStats::OPTIMAL_SOLUTION ) {
      if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
        out << "\nCould not solve the QP!\n";
      }
      return false;
    }

    //
    // Set the solution
    //
    if( use_simple_wz_bounds ) {
      // Set Zwz
      Zwz(dep)   = _Dwz;
      Zwz(indep) = wz;
    }
    else {
      // Everything should already be set!
    }

    // Cut back Ywy = (1-eta) * Ywy
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    if( fabs(qp_eta - 0.0) > eps ) {
      if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
        out
          << "\n*** Alert! the QP was infeasible (eta = "<<qp_eta<<").  Cutting back Ywy_k = (1.0 - eta)*Ywy  ...\n";
      }
      DenseLinAlgPack::Vt_S( &Ywy , 1.0 - qp_eta );
    }
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\n||wz||inf    = " << DenseLinAlgPack::norm_inf(wz);
    out	<< "\n||Zwz||2     = " << DenseLinAlgPack::norm_2(Zwz);
    if(qp_eta > 0.0) out << "\n||Ypy||2 = " << DenseLinAlgPack::norm_2(Ywy);
    out << std::endl;
  }
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nwz = \n" << wz;
    out << "\nZwz = \n" << Zwz;
    if(qp_eta > 0.0) out << "\nYwy = \n" << Ywy;
  }
  if( qp_eta == 1.0 ) {
    std::ostringstream omsg;
    omsg
      << "FeasibilityStepReducedStd_Strategy::compute_feasibility_step(...) : "
      << "Error, a QP relaxation parameter\n"
      << "of eta = " << qp_eta << " was calculated and therefore it must be assumed\n"
      << "that the NLP's constraints are infeasible\n"
      << "Throwing an InfeasibleConstraints exception!\n";
    if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << omsg.str();
    }
    throw_qp_failure = true;
//		throw InfeasibleConstraints(omsg.str());
  }

  //
  // Set the full step
  //
  // w = Ywy + Zwz
  //
  DenseLinAlgPack::V_VpV( w, Ywy, Zwz );

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\n||w||inf    = " << DenseLinAlgPack::norm_inf(*w);
    out << std::endl;
  }
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nw = \n" << *w;
  }

  current_k_ = s->k();

  if( throw_qp_failure )
    return false;

*/
  TEUCHOS_TEST_FOR_EXCEPT(true);

  return true;
}

void FeasibilityStepReducedStd_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
  out << L << "*** Computes feasibility steps by solving a constrained QP using the range and null\n"
    << L << "*** space decomposition\n"
    << L << "begin quais-range space step: \"" << typeName(quasi_range_space_step()) << "\"\n";

   quasi_range_space_step().print_step( out, L + "  " );

  out << L << "end quasi-range space step\n"
    << L << "if quasi-range space step failed then\n"
    << L << "  this feasibility step fails also!\n"
    << L << "end\n"
    << L << "Ywy = v\n"
    << L << "*** Set the bounds for bl <= Z*wz <= bu\n"
    << L << "bl = d_bounds_k.l - (xo - x_k) - Ywy\n"
    << L << "bu = d_bounds_k.u - (xo - x_k) - Ywy\n"
    << L << "Set bl(i) = bu(i) for those nu_k(i) != 0.0\n"
    << L << "if (qp_objective == OBJ_MIN_FULL_STEP) and (current_k != k) then\n"
    << L << "  grad = Z_k'*Ywy\n"
    << L << "  Hess = Z_k'*Z_k\n"
    << L << "elseif (qp_objective == OBJ_MIN_NULL_SPACE_STEP) and (current_k != k) then\n"
    << L << "  grad = 0\n"
    << L << "  Hess = I\n"
    << L << "elseif (qp_objective == OBJ_RSQP) and (current_k != k) then\n"
    << L << "  grad = qp_grad_k\n"
    << L << "  Hess = rHL_k\n"
    << L << "end\n"
    << L << "if check_results == true then\n"
    << L << "  assert that bl and bu are valid and sorted\n"
    << L << "end\n"
    << L << "etaL = 0.0\n"
    << L << "*** Determine if we can use simple bounds on pz or not\n"
    << L << "if Z_k is a variable reduction null space matrix and norm(Ypy_k(indep),0) == 0 then\n"
    << L << "  use_simple_wz_bounds = true\n"
    << L << "else\n"
    << L << "  use_simple_wz_bounds = false\n"
    << L << "end\n"
    << L << "*** Setup QP arguments\n"
    << L << "qp_g = qp_grad_k\n"
    << L << "qp_G = rHL_k\n"
    << L << "if use_simple_wz_bounds == true then\n"
    << L << "  qp_dL = bl(indep),  qp_dU = bu(indep)\n"
    << L << "  qp_E  = Z_k.D,      qp_b  = Ywy(dep)\n"
    << L << "  qp_eL = bl(dep),    qp_eU = bu(dep)\n"
    << L << "else\n"
    << L << "  qp_dL = -inf,       qp_dU = +inf\n"
    << L << "  qp_E  = Z_k,        qp_b  = Ywy\n"
    << L << "  qp_eL = bl,         qp_eU = bu\n"
    << L << "end\n"
    << L << "if m > r then\n"
    << L << "  qp_F  = V_k,        qp_f  = c_k(undecomp) + Gc_k(undecomp)'*Ywy\n"
    << L << "else\n"
    << L << "  qp_F  = empty,      qp_f  = empty\n"
    << L << "end\n"
    << L << "Use a warm start given the active-set in nu_k\n"
    << L << "Solve the following QP to compute qp_d, qp_eta, qp_Ed = qp_E * qp_d\n"
    << L << ",qp_nu, qp_mu and qp_lambda (" << typeName(qp_solver()) << "):\n"
    << L << "  min    qp_g' * qp_d + 1/2 * qp_d' * qp_G * qp_d + M(eta)\n"
    << L << "  qp_d <: R^(n-r)\n"
    << L << "  s.t.\n"
    << L << "         etaL  <=  qp_eta\n"
    << L << "         qp_dL <= qp_d <= qp_dU                          [qp_nu]\n"
    << L << "         qp_eL <= qp_E * qp_d + (1-eta)*qp_b  <= qp_eU   [qp_mu]\n"
    << L << "         qp_F * d_qp + (1-eta) * qp_f = 0                [qp_lambda]\n"
    << L << "if (qp_teing==QP_TEST) or (fd_deriv_testing==QP_TEST_DEFAULT\n"
    << L << "and check_results==true) then\n"
    << L << "  Check the optimality conditions of the above QP\n"
    << L << "  if the optimality conditions do not check out then\n"
    << L << "    set throw_qp_failure = true\n"
    << L << "  end\n"
    << L << "end\n"
    << L << "*** Set the solution to the QP subproblem\n"
    << L << "wz  = qp_d\n"
    << L << "eta = qp_eta\n"
    << L << "if use_simple_wz_bounds == true then\n"
    << L << "  Zwz(dep)   = qp_Ed,  Zwz(indep) = wz\n"
    << L << "else\n"
    << L << "  Zwz = qp_Ed\n"
    << L << "end\n"
    << L << "if eta > 0 then set Ywy = (1-eta) * Ywy\n"
    << L << "if QP solution is suboptimal then\n"
    << L << "  throw_qp_failure = true\n"
    << L << "elseif QP solution is primal feasible (not optimal) then\n"
    << L << "  throw_qp_failure = primal_feasible_point_error\n"
    << L << "elseif QP solution is dual feasible (not optimal) then\n"
    << L << "  find max u s.t.\n"
    << L << "    d_bounds_k.l <= (xo - x) + u*(Ywy+Zwz) <= d_bounds_k.u\n"
    << L << "  alpha_k = u\n"
    << L << "  throw_qp_failure = true\n"
    << L << "end\n"
    << L << "if eta == 1.0 then\n"
    << L << "  The constraints are infeasible!\n"
    << L << "  throw_qp_failure = true\n"
    << L << "end\n"
    << L << "current_k = k\n"
    << L << "w = Zwz + Ywy\n"
    << L << "if (throw_qp_failure == true) then\n"
    << L << "  The feasibility step computation has failed!\n"
    << L << "end\n"
    ;
}

} // end namespace MoochoPack
