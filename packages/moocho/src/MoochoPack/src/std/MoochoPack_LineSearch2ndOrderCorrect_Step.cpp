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

#include <ostream>
#include <typeinfo>

#include "MoochoPack_LineSearch2ndOrderCorrect_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack_print_vector_change_stats.hpp"
#include "ConstrainedOptPack_MeritFuncCalc1DQuadratic.hpp"
#include "ConstrainedOptPack_MeritFuncCalcNLP.hpp"
#include "ConstrainedOptPack_MeritFuncCalcNLE.hpp"
#include "ConstrainedOptPack_MeritFuncNLESqrResid.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack/src/max_near_feas_step.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

namespace MoochoPack {

LineSearch2ndOrderCorrect_Step::LineSearch2ndOrderCorrect_Step(
  const direct_ls_sqp_ptr_t&			direct_ls_sqp
  ,const merit_func_ptr_t&			merit_func
  ,const feasibility_step_ptr_t&		feasibility_step
  ,const direct_ls_newton_ptr_t&		direct_ls_newton
  ,value_type							eta
  ,ENewtonOutputLevel					newton_olevel
  ,value_type							constr_norm_threshold
  ,value_type							constr_incr_ratio
  ,int								after_k_iter
  ,EForcedConstrReduction				forced_constr_reduction
  ,value_type                         forced_reduct_ratio
  ,value_type							max_step_ratio
  ,int								max_newton_iter
  )
  :direct_ls_sqp_(direct_ls_sqp)
  ,merit_func_(merit_func)
  ,feasibility_step_(feasibility_step)
  ,direct_ls_newton_(direct_ls_newton)
  ,eta_(eta)
  ,newton_olevel_(newton_olevel)
  ,constr_norm_threshold_(constr_norm_threshold)
  ,constr_incr_ratio_(constr_incr_ratio)
  ,after_k_iter_(after_k_iter)
  ,forced_constr_reduction_(forced_constr_reduction)
  ,forced_reduct_ratio_(forced_reduct_ratio)
  ,max_step_ratio_(max_step_ratio)
  ,max_newton_iter_(max_newton_iter)
{}

bool LineSearch2ndOrderCorrect_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{

/* ToDo: Upate the below code!

  using std::setw;

  using DenseLinAlgPack::dot;
  using DenseLinAlgPack::norm_inf;
  using DenseLinAlgPack::V_VpV;
  using DenseLinAlgPack::V_VmV;
  using DenseLinAlgPack::Vp_StV;
  using DenseLinAlgPack::Vt_S;

  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_MtV;

  using AbstractLinAlgPack::max_near_feas_step;

  using ConstrainedOptPack::print_vector_change_stats;

  typedef LineSearch2ndOrderCorrect_Step	this_t;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  NLP			&nlp	= algo.nlp();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();
  out << std::boolalpha;

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  // /////////////////////////////////////////
  // Set references to iteration quantities
  //
  // Set k+1 first then go back to get k to ensure
  // we have backward storage.
  
  DVector
    &x_kp1 = s.x().set_k(+1).v();
  value_type
    &f_kp1 = s.f().set_k(+1);
  DVector
    &c_kp1 = s.c().set_k(+1).v();

  const value_type
    &f_k = s.f().get_k(0);
  const DVector
    &c_k = s.c().get_k(0).v();
  const DVector
    &x_k = s.x().get_k(0).v();
  const DVector
    &d_k = s.d().get_k(0).v();
  value_type
    &alpha_k = s.alpha().get_k(0);

  // /////////////////////////////////////
  // Compute Dphi_k, phi_kp1 and phi_k

  // Dphi_k
  const value_type
    Dphi_k = merit_func().deriv();
  if( Dphi_k >= 0 ) {
    throw LineSearchFailure( "LineSearch2ndOrderCorrect_Step::do_step(...) : " 
      "Error, d_k is not a descent direction for the merit function " );
  }

  // ph_kp1
  value_type
    &phi_kp1 = s.phi().set_k(+1) = merit_func().value( f_kp1, c_kp1 );

  // Must compute phi(x) at the base point x_k since the penalty parameter may have changed.
  const value_type
    &phi_k = s.phi().set_k(0) = merit_func().value( f_k, c_k );

  // //////////////////////////////////////
  // Setup the calculation merit function

  // Here f_kp1, and c_kp1 are updated at the same time the
  // line search is being performed.
  nlp.set_f( &f_kp1 );
  nlp.set_c( &c_kp1 );
  MeritFuncCalcNLP
    phi_calc( &merit_func(), &nlp );

  // //////////////////////////////////////////////////
  // Concider 2nd order correction if near solution?
  
  bool considering_correction = false;
  {
    const value_type
      small_num = std::numeric_limits<value_type>::min(),
      nrm_c_k   = s.c().get_k(0).norm_inf(),
      nrm_c_kp1 = s.c().get_k(+1).norm_inf();
    const bool
      test_1 = nrm_c_k <= constr_norm_threshold(),
      test_2 = (nrm_c_kp1/(1.0 + nrm_c_k)) < constr_incr_ratio(),
      test_3 = s.k() >= after_k_iter();
    considering_correction = test_1 && test_2 && test_3;
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "\n||c_k||inf = " << nrm_c_k << (test_1 ? " <= " : " > " )
        << "constr_norm_threshold = " << constr_norm_threshold()
        << "\n||c_kp1||inf/(1.0+||c_k||inf) = "
        << "(" << nrm_c_kp1 << ")/(" << 1.0 << " + " << nrm_c_k << ") = "
        << ( nrm_c_kp1 / (1.0 + nrm_c_k ) ) << (test_2 ? " <= " : " > " )
        << "constr_incr_ratio = " << constr_incr_ratio()
        << "\nk = " << s.k() << (test_3 ? " >= " : " < ")
        << "after_k_iter = " << after_k_iter()
        << (considering_correction
          ? ("\nThe computation of a 2nd order correction for x_kp1 = x_k + alpha_k*d_k + alpha_k^2*w"
             " will be considered ...\n")
          : "\nThe critera for considering a 2nd order correction has not been met ...\n" );
    }
  }

  // //////////////////////////////
  // See if we can take a full step

  bool chose_point = false;

  const value_type frac_phi = phi_k + eta() * Dphi_k;
  const bool armijo_test = phi_kp1 <= frac_phi;
  if( armijo_test ) {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "\nAccepting full step x_kp1 = x_k + d_k\n";
    }
    chose_point = true;	// The point meets the Armijo test.
  }

  // This is storage for function and gradient evaluations for
  // the trial newton points and must be remembered for latter
  value_type f_xdww;
  DVector     c_xdww;
  DVector w(x_kp1.size()),		// Full correction after completed computation.
       xdww(x_kp1.size());	// Will be set to xdw + sum( w(newton_i), newton_i = 1... )
                // where w(itr) is the local corrections for the current
                // newton iteration.
  bool use_correction = false;
  bool newton_failed  = true;
  
  bool considered_correction = ( considering_correction && !chose_point );
  if( considered_correction ) {

    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "\nConsidering whether to compute a 2nd order correction for\n"
          << "x_kp1 = x_k + alpha_k * d_k + alpha_k^2 * w ...\n";
    }

    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      const value_type obj_descent = dot( s.Gf().get_k(0)(), d_k() );
      out	<< "\nGf_k' * d_k = " << obj_descent << std::endl;
      if( obj_descent >= 0.0 ) {
        out	<< "\nWarning, this may not work well with Gf_k'*d_k >= 0.0\n";
      }
    }

    // Merit function for newton line searches
    ConstrainedOptPack::MeritFuncNLESqrResid
      phi_c;

    DVector
      xdw = x_kp1;	// Will be set to x + d + sum(w(i),i=1..itr-1)
              //     where w(i) are previous local corrections
    value_type
      phi_c_x    = phi_c.value( c_k() ),
      phi_c_xd   = phi_c.value( c_kp1() ),
      phi_c_xdw  = phi_c_xd,		// No correction is computed yet so w = 0
      phi_c_xdww = phi_c_xdw,
      nrm_d	   = norm_inf( d_k() );

    // Merit function for newton line searches
    nlp.set_f( &(f_xdww = f_kp1) );
    nlp.set_c( &(c_xdww = c_kp1) );
    ConstrainedOptPack::MeritFuncCalcNLE
      phi_c_calc( &phi_c, &nlp );

    DVector wy(s.con_decomp().size());	// Range space wy (see latter).

    const bool sufficient_reduction =
      phi_c_xd < forced_reduct_ratio() * phi_c_x;
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "\nphi_c(c(x_k+d_k)) = " << phi_c_xd << (sufficient_reduction ? " <= " : " > ")
        << "forced_reduct_ratio* phi_c(c(x_k)) = " << forced_reduct_ratio() << " * " << phi_c_x
        << " = " << (forced_reduct_ratio()*phi_c_x)
        << (sufficient_reduction
          ? "\nNo need for a 2nd order correciton, perform regular line search ... \n"
          : "\nWe need to try to compute a correction w ...\n" );
    }
    if(sufficient_reduction) {
      use_correction = false;
    }
    else {
      // Try to compute a second order correction term.

      // Set print level newton iterations
      ENewtonOutputLevel  use_newton_olevel;
      if( newton_olevel() == PRINT_USE_DEFAULT ) {
        switch(olevel) {
            case PRINT_NOTHING:
            case PRINT_BASIC_ALGORITHM_INFO:
            use_newton_olevel = PRINT_NEWTON_NOTHING;
            break;
            case PRINT_ALGORITHM_STEPS:
            case PRINT_ACTIVE_SET:
            use_newton_olevel = PRINT_NEWTON_SUMMARY_INFO;
            break;
            case PRINT_VECTORS:
            use_newton_olevel = PRINT_NEWTON_STEPS;
            break;
            case PRINT_ITERATION_QUANTITIES:
            use_newton_olevel = PRINT_NEWTON_VECTORS;
            break;
            default:
            TEUCHOS_TEST_FOR_EXCEPT(true);
        }
      }
      else {
        use_newton_olevel = newton_olevel();
      }
      EJournalOutputLevel inner_olevel;
      switch(use_newton_olevel) {
          case PRINT_NEWTON_NOTHING:
          case PRINT_NEWTON_SUMMARY_INFO:
          inner_olevel = PRINT_NOTHING;
          break;
          case PRINT_NEWTON_STEPS:
          inner_olevel = PRINT_ALGORITHM_STEPS;
          break;
          case PRINT_NEWTON_VECTORS:
          if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES )
            inner_olevel = PRINT_ITERATION_QUANTITIES;
          else if( (int)olevel >= (int)PRINT_ACTIVE_SET )
            inner_olevel = PRINT_ACTIVE_SET;
          else
            inner_olevel = PRINT_VECTORS;
          break;
          default:
          TEUCHOS_TEST_FOR_EXCEPT(true);
      }
      
      // Print header for summary information
      const int dbl_min_w = 21;
      const int dbl_w = std::_MAX(dbl_min_w,int(out.precision()+8));
      if( use_newton_olevel == PRINT_NEWTON_SUMMARY_INFO ) {
        out	<< "\nStarting Newton iterations ...\n"
          << "\nphi_c_x   = "	<< phi_c_x 
          << "\nphi_c_xd  = "	<< phi_c_xd
          << "\n||d_k||nf = "	<< nrm_d << "\n\n"
          << setw(5)			<< "it"
          << setw(dbl_w)		<< "||w||inf"
          << setw(dbl_w)		<< "u"
          << setw(dbl_w)		<< "step_ratio"
          << setw(5)			<< "lsit"
          << setw(dbl_w)		<< "a"
          << setw(dbl_w)		<< "phi_c_xdww"
          << setw(dbl_w)		<< "phi_c_xdww-phi_c_x"
          << setw(dbl_w)		<< "phi_c_xdww-phi_c_xd\n"
          << setw(5)			<< "----"
          << setw(dbl_w)		<< "--------------------"
          << setw(dbl_w)		<< "-------------------"
          << setw(dbl_w)		<< "-------------------"
          << setw(5)			<< "----"
          << setw(dbl_w)		<< "-------------------"
          << setw(dbl_w)		<< "-------------------"
          << setw(dbl_w)		<< "-------------------"
          << setw(dbl_w)		<< "-------------------\n";
      }

      // Perform newton feasibility iterations
      int newton_i;
      for( newton_i = 1; newton_i <= max_newton_iter(); ++newton_i ) {
        
        if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_STEPS ) {
          out << "\n**** newton_i = " << newton_i << std::endl;
        }

        // Compute a feasibility step
        if(!feasibility_step().compute_feasibility_step(
          out,inner_olevel,&algo,&s,xdw,nlp.c()(),&w() ))
        {
          if( (int)use_newton_olevel == (int)PRINT_NEWTON_SUMMARY_INFO ) {
            out << "\nCould not compute feasible direction!\n";
          }
          break; // exit the newton iterations
        }

        value_type
          nrm_w = norm_inf(w());				

        if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_STEPS ) {
          out << "\n||w||inf = " << nrm_w << std::endl;
        }

        if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_VECTORS ) {
          out << "\nw = " << w();
        }

        // ////////////////////////////////
        // Cutting back w

        value_type a = 1.0;	// This is the alpha for your linesearch

        // Cut back w to be in the relaxed bounds.
        std::pair<value_type,value_type>
          u_steps = max_near_feas_step( s.x().get_k(0)(), w()
            , algo.nlp().xl(), algo.nlp().xu()
            , algo.algo_cntr().max_var_bounds_viol() );
        const value_type u = u_steps.first;

        if( u < a ) {
          if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_STEPS ) {
            out << "\nCutting back w = (a=u) * w to be within relaxed bounds:\n"
              << "u = " << u << std::endl;
          }
          a = u;
        }

        // Cut back step so x+d+sum(w(i)) is not too far from x+d
        value_type
          step_ratio = nrm_w / nrm_d;
        if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_STEPS ) {
          out << "\nstep_ratio = ||w||inf/||d||inf = " << step_ratio
              << std::endl;
        }
        if( a * step_ratio > max_step_ratio() ) {
          const value_type aa = a*(max_step_ratio()/step_ratio);
          if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_STEPS ) {
            out << "\na*step_ratio = " << a*step_ratio
                << " > max_step_ratio = " << max_step_ratio() << std::endl
              << "Cutting back a = a*max_step_ratio/step_ratio = "
                << aa << std::endl;
          }
          a = aa;
        }

        // /////////////////////////////////////////////////
        // Perform a line search along xdww = xdw + a * w
        
        if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_STEPS ) {
          out << "\nPerform linesearch along xdww = xdw + a*w\n"
            << "starting from a = " << a << " ...\n";
        }

        xdww = xdw();									// xdww = xdw + a*w
        Vp_StV( &xdww(), a, w() );
        phi_c.calc_deriv(nlp.c());	// Set the directional derivative at c(xdw)
        phi_c_xdww = phi_c_calc( xdww() );	// phi_c_xdww = phi(xdww)
        const DVectorSlice xdw_w[2] = { xdw(), w() };
        MeritFuncCalc1DQuadratic
          phi_c_calc_1d( phi_c_calc, 1 , xdw_w, &xdww() );
        bool ls_okay = false;
        try {
          ls_okay = direct_ls_newton().do_line_search(
            phi_c_calc_1d,phi_c_xdw
            ,&a,&phi_c_xdww
            , (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_STEPS 
            ? &out : 0												
            );
        }
        catch(const DirectLineSearch_Strategy::NotDescentDirection& excpt ) {
          if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO ) {
            out << "\nThe line search object throw the exception:" << typeName(excpt) << ":\n"
              << excpt.what() << std::endl;
          }
          ls_okay = false;
        }
        // Note that the last value c(x) computed by the line search is for
        // xdw + a*w.

        // Print line for summary output
        if( use_newton_olevel == PRINT_NEWTON_SUMMARY_INFO ) {
          out	<< setw(5)			<< newton_i
            << setw(dbl_w)		<< nrm_w
            << setw(dbl_w)		<< u
            << setw(dbl_w)		<< step_ratio
            << setw(5)			<< direct_ls_newton().num_iterations()
            << setw(dbl_w)		<< a
            << setw(dbl_w)		<< phi_c_xdww
            << setw(dbl_w)		<< (phi_c_xdww-phi_c_x)
            << setw(dbl_w)		<< (phi_c_xdww-phi_c_xd) << std::endl;
        }

        if(!ls_okay) {
          if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO ) {
            out << "\nThe line search failed so forget about computing a correction ...\n";
          }
          use_correction = false;
          newton_failed = true;
          break;
        }

        // See if this point is okay
        bool good_correction = false;
        switch( forced_constr_reduction() ) {
          case CONSTR_LESS_X_D: {
            good_correction = ( phi_c_xdww < forced_reduct_ratio()*phi_c_xd );
            if( good_correction
              && (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO )
            {
              out << "\nphi_c(c(x_k+d_k+w)) = " << phi_c_xdww
                << " < forced_reduct_ratio * phi_c(c(x_k+d_k)) = "
                << forced_reduct_ratio() << " * " << phi_c_xd
                << " = " << (forced_reduct_ratio()*phi_c_xd) << std::endl;
            }
            break;
          }
          case CONSTR_LESS_X: {
            good_correction = ( phi_c_xdww < forced_reduct_ratio()*phi_c_x );
            if( good_correction
              && (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO )
            {
              out << "\nphi_c(c(x_k+d_k+w)) = " << phi_c_xdww
                << " < forced_reduct_ratio * phi_c(c(x_k)) = "
                << forced_reduct_ratio() << " * " << phi_c_x
                << " = " << (forced_reduct_ratio()*phi_c_x) << std::endl;
            }
            break;
          }
          default:
            TEUCHOS_TEST_FOR_EXCEPT(true);
        }

        if(good_correction) {
          if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO ) {
            out << "\nAccept this point and compute our full correction w ... \n";
          }
          // Compute the full correction and do a curved linesearch
          // w = xdww - x_kp1
          V_VmV( &w(), xdww(), x_kp1() );
          if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO ) {
            out << "\n||w||inf = " << norm_inf(w()) << std::endl;
          }
          if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_VECTORS ) {
            out << "\nw = " << w();
          }
          use_correction = true;
          newton_failed  = false;
          break;
        }

        // Else perform another newton iteration.
        xdw       = xdww;
        phi_c_xdw = phi_c_xdww;

      }	// end for
      if( !use_correction ) {
        newton_failed  = true;
        if( forced_constr_reduction() == CONSTR_LESS_X_D ) {
          if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
            out	<< "\nDam! This is really bad!\n"
              << "We where only looking for point phi_c(c(x_k+d_k+w)"
              << " < phi_c(c(x_k+k_k) and we could not find it\n"
              << " in the aloted number of newton iterations!\n"
              << "Perhaps the Gc_k did not give us a descent direction?\n"
              << "Just perform a standard line search from here ...\n";
          }
          use_correction = false;
        }
        else {
          if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
            out	<< "\nWe where looking for point phi_c(c(x_k+d_k+w))"
              << " < phi_c(c(x_k)) and we could not find it.\n";
          }
          if( phi_c_xdww < phi_c_xd ) {
            if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
              out	<< "However, we did find a point less than phi_c(c(x_k+d_k))\n"
                << "so lets use the correction anyway.\n";
            }
            // Compute the full correction and do a curved linesearch
            // w = xdww - x_kp1
            V_VmV( &w(), xdww(), x_kp1() );
            if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO ) {
              out << "\n||w||inf = " << norm_inf(w()) << std::endl;
            }
            if( (int)use_newton_olevel >= (int)this_t::PRINT_NEWTON_VECTORS ) {
              out << "\nw = " << w();
            }
            use_correction = true;
          }
          else {
            if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
              out	<< "Dam! We did not even find a point less than phi_c(c(x_k+d_k))\n"
                << "just perform a standard line search along x_k + alpha_k*d_k.\n";
            }
            use_correction = false;
          }
        }
      }
    }	// end else from if phi_c_xdw > phi_c_x
  } // end considered_correction

  // //////////////////////////
  // Set up for the line search

  if( considered_correction ) {
    if( use_correction ) {
      // We are using the correction so setup the full step for the
      // NLP linesearch to come.
      Vp_V( &x_kp1(), w() );      // Set point to x_kp1 = x_k + d_k + w
      nlp.calc_f(x_kp1(),false);  // same as last call to calc_c(x)
      f_kp1 = nlp.f();            // Here f and c where computed at x_k+d_k+w
      c_kp1 = nlp.c()();
      phi_kp1 = merit_func().value( f_kp1, c_kp1 );
    }
    else {
      // Just pretend the second order correction never happened
      // and we don't need to do anything.
    }
    // Set back the base point
    nlp.set_f( &f_kp1 );
    nlp.set_c( &c_kp1 );
  }

  // //////////////////////
  // Do the line search

  if( !chose_point ) {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      if( use_correction ) {
        out	<< "\nPerform a curved linesearch along:\n"
          << "x_kp1 = x_k + alpha_k * d_k + alpha_k^2 * w ...\n";
      }
      else {
        out	<< "\nPerform a standard linesearch along:\n"
          << "x_kp1 = x_k + alpha_k * d_k ...\n";
      }
    }
    const DVectorSlice xdw[3] = { x_k(), d_k(), w() };
    MeritFuncCalc1DQuadratic
      phi_calc_1d( phi_calc, (use_correction?2:1) , xdw, &x_kp1() );
    if( !direct_ls_sqp().do_line_search( phi_calc_1d, phi_k, &alpha_k, &phi_kp1
      , static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ?
        &out : static_cast<std::ostream*>(0)	)		)
    {
      // If the line search failed but the value of the merit function is less than
      // the base point value then just accept it and move on.  This many be a
      // bad thing to do.

      const value_type
        scaled_ared		= (s.phi().get_k(0) - s.phi().get_k(+1))/s.phi().get_k(0),
        keep_on_frac	= 1.0e-10;	// Make adjustable?
      bool keep_on = scaled_ared < keep_on_frac;

      if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO )
      {
        out
          << "\nThe maximum number of linesearch iterations has been exceeded "
          << "(k = " << algo.state().k() << ")\n"
          << "(phi_k - phi_kp1)/phi_k = " << scaled_ared;
//				if(keep_on) {
//					out
//						<< " < " << keep_on_frac
//						<< "\nso we will accept to step and move on.\n";
//				}
//				else {
//					out
//						<< " > " << keep_on_frac
//						<< "\nso we will reject the step and declare a line search failure.\n";
//				}
      }
//
//			if( keep_on ) return true;
      
      throw LineSearchFailure( "LineSearch2ndOrderCorrect_Step::do_step(): "
                   "Error, Line search failure" );
    }
  }

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    out	<< "\nalpha_k      = "	<< alpha_k				<< std::endl
      << "\n||x_kp1||inf = "	<< norm_inf( x_kp1 )	<< std::endl
      << "\nf_kp1        = "	<< f_kp1				<< std::endl
      << "\n||c_kp1||inf = "	<< norm_inf(c_kp1)		<< std::endl
      << "\nphi_kp1      = "	<< phi_kp1				<< std::endl;
  }

  if( (int)olevel >= (int)PRINT_VECTORS ) {
    out << "\nx_kp1 =\n"	<< x_kp1
      << "\nc_kp1 =\n"	<< c_kp1;
  }

*/
  TEUCHOS_TEST_FOR_EXCEPT(true);

  return true;
}

void LineSearch2ndOrderCorrect_Step::print_step(
  const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out	<< L << "*** Calculate a second order correction when near solution.\n"
    << L << "*** If we can compute a correction w then perform a curved\n"
    << L << "*** line search along x_kp1 = x_k + alpha_k * d_k + alpha_k^2 * w.\n"
    << L << "default: eta                     = " << eta() << std::endl
    << L << "         constr_norm_threshold   = " << constr_norm_threshold() << std::endl
    << L << "         constr_incr_ratio       = " << constr_incr_ratio() << std::endl
    << L << "         after_k_iter            = " << after_k_iter() << std::endl
    << L << "         forced_constr_reduction = " << (forced_constr_reduction()== CONSTR_LESS_X_D
                              ? "CONSTR_LESS_X_D\n"
                              : "CONSTR_LESS_X\n" )
    << L << "         forced_reduct_ratio     = " << forced_reduct_ratio() << std::endl
    << L << "         max_step_ratio          = " << max_step_ratio() << std::endl
    << L << "         max_newton_iter         = " << max_newton_iter() << std::endl
    << L << "begin definition of NLP merit function phi.value(f(x),c(x)):\n";
  
  merit_func().print_merit_func( out, L + "  " );
  
  out	<< L << "end definition\n"
    << L << "Dphi_k = phi.deriv()\n"
    << L << "if Dphi_k >= 0 then\n"
    << L << "  throw line_search_failure\n"
    << L << "end\n"
    << L << "phi_kp1 = phi_k.value(f_kp1,c_kp1)\n"
    << L << "phi_k = phi.value(f_k,c_k)\n"
    << L << "if (norm(c_k,inf) <= constr_norm_threshold)\n"
    << L << "and (norm(c_kp1,inf)/(norm(c_k,inf)+small_num) <= constr_incr_ratio)\n"
    << L << "and (k >= after_k_iter) then\n"
    << L << "considering_correction = ( (norm(c_k,inf) <= constr_norm_threshold)\n"
    << L << "  and (norm(c_kp1,inf)/(1.0 + norm(c_k,inf)) <= constr_incr_ratio)\n"
    << L << "  and (k >= after_k_iter) )\n"
    << L << "chose_point = false\n"
    << L << "if phi_kp1 < phi_k + eta * Dphi_k then\n"
    << L << "  chose_point = true\n"
    << L << "else\n"
    << L << "if (considering_correction == true) and (chose_point == false) then\n"
    << L << "  considered_correction = true\n"
    << L << "  begin definition of c(x) merit function phi_c.value(c(x)):\n";

  ConstrainedOptPack::MeritFuncNLESqrResid().print_merit_func(
    out, L + "    " );

  out	<< L << "  end definition\n"
    << L << "  xdw = x_kp1;\n"
    << L << "  phi_c_x = phi_c.value(c_k);\n"
    << L << "  phi_c_xd = phi_c.value(c_kp1);\n"
    << L << "  phi_c_xdw = phi_c_xd;\n"
    << L << "  phi_c_xdww = phi_c_xdw;\n"
    << L << "  if phi_c_xd < forced_reduct_ratio * phi_c_x then\n"
    << L << "    *** There is no need to perform a correction.\n"
    << L << "    use_correction = false;\n"
    << L << "  else\n"
    << L << "    *** Lets try to compute a correction by performing\n"
    << L << "    *** a series of newton steps to compute local steps w\n"
    << L << "    for newton_i = 1...max_newton_itr\n"
    << L << "      begin feasibility step calculation: \"" << typeName(feasibility_step()) << "\"\n";

  feasibility_step().print_step( out, L + "        " );

  out << L << "      end feasibility step calculation\n"
    << L << "      Find the largest positive step u where the slightly\n"
    << L << "      relaxed variable bounds:\n"
    << L << "        xl - delta <= xdw + u * w <= xu + delta\n"
    << L << "      are strictly satisfied (where delta = max_var_bounds_viol).\n"
    << L << "      a = min(1,u);\n"
    << L << "      step_ratio = norm(w,inf)/norm(d,inf);\n"
    << L << "      a = min(a,max_step_ratio/step_ratio);\n"
    << L << "      Perform line search for phi_c.value(c(xdww = xdw+a*w))\n"
    << L << "      starting from a and compute:\n"
    << L << "        a,\n"
    << L << "        xdww = xdw + a * w,\n"
    << L << "        phi_c_xdww = phi_c.value(c(xdww))\n"
    << L << "      print summary information;\n"
    << L << "      if line search failed then\n"
    << L << "        use_correction = false;\n"
    << L << "        exit for loop;\n"
    << L << "      end\n"
    << L << "      *** Determine if this is sufficent reduction in c(x) error\n"
    << L << "      if forced_constr_reduction == CONSTR_LESS_X_D then\n"
    << L << "        good_correction = (phi_c.value(c(xdww)) < forced_reduct_ratio*phi_c_xd);\n"
    << L << "      else if forced_constr_reduction == CONSTR_LESS_X then\n"
    << L << "        good_correction = (phi_c.value(c(xdww)) < forced_reduct_ratio*phi_c_x);\n"
    << L << "      end\n"
    << L << "      if good_correction == true then\n"
    << L << "        w = xdww - (x_k+d_k);\n"
    << L << "        use_correction = true;\n"
    << L << "        exit for loop;\n"
    << L << "      end\n"
    << L << "      *** This is not sufficient reduction is c(x) error yet.\n"
    << L << "      xdw = xdww;\n"
    << L << "      phi_c_xdw = phi_c_xdww;\n"
    << L << "    end\n"
    << L << "    if use_correction == false then\n"
    << L << "      if forced_constr_reduction == CONSTR_LESS_X_D then\n"
    << L << "        *** Dam! We could not find a point phi_c_xdww < phi_c_xd.\n"
    << L << "        *** Perhaps Gc_k does not give a descent direction for phi_c!\n"
    << L << "      else if forced_constr_reduction == CONSTR_LESS_X then\n"
    << L << "        if phi_c_dww < phi_c_xd then\n"
    << L << "           *** Accept this correction anyway.\n"
    << L << "           use_correction = true\n"
    << L << "        else\n"
    << L << "          *** Dam! we could not find any reduction in phi_c so\n"
    << L << "          *** Perhaps Gc_k does not give a descent direction for phi_c!\n"
    << L << "      end\n"
    << L << "    end\n"
    << L << "  end\n"
    << L << "end\n"
    << L << "if chose_point == false then\n"
    << L << "  if use_correction == true then\n"
    << L << "    Perform line search along x_kp1 = x_k + alpha_k * d_k + alpha_k^2 * w\n"
    << L << "  else\n"
    << L << "    Perform line search along x_kp1 = x_k + alpha_k * d_k\n"
    << L << "  end\n"
    << L << "  begin direct line search : \"" << typeName(direct_ls_sqp()) << "\"\n";

  direct_ls_sqp().print_algorithm( out, L + "    " );

  out
    << L << "  end direct line search\n"
    << L << "  if maximum number of linesearch iterations are exceeded then\n"
    << L << "    throw line_search_failure\n"
    << L << "  end\n"
    << L << "end\n";
}

}	// end namespace MoochoPack

#endif // 0
