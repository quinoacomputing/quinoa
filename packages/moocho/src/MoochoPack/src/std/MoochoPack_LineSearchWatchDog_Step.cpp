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

#include "MoochoPack_LineSearchWatchDog_Step.hpp"
#include "MoochoPack_MoochoAlgorithmStepNames.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack_MeritFuncCalc1DQuadratic.hpp"
#include "ConstrainedOptPack_MeritFuncCalcNLP.hpp"
#include "ConstrainedOptPack_print_vector_change_stats.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"

namespace {
  const int NORMAL_LINE_SEARCH = -1;
}

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

MoochoPack::LineSearchWatchDog_Step::LineSearchWatchDog_Step(
      const direct_line_search_ptr_t&	direct_line_search
    , const merit_func_ptr_t&			merit_func
    , value_type						eta
    , value_type						opt_kkt_err_threshold
    , value_type						feas_kkt_err_threshold
    )
  :
      direct_line_search_(direct_line_search)
    , merit_func_(merit_func)
    , eta_(eta)
    , opt_kkt_err_threshold_(opt_kkt_err_threshold)
    , feas_kkt_err_threshold_(feas_kkt_err_threshold)
    , watch_k_(NORMAL_LINE_SEARCH)
{}

bool MoochoPack::LineSearchWatchDog_Step::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{
  using DenseLinAlgPack::norm_inf;
  using DenseLinAlgPack::V_VpV;
  using DenseLinAlgPack::Vp_StV;
  using DenseLinAlgPack::Vt_S;

  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_MtV;

  using ConstrainedOptPack::print_vector_change_stats;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  NLP			&nlp	= algo.nlp();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();
  out << std::boolalpha;

  // print step header.
  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
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

  // ////////////////////////////////
  // Use Watchdog near the solution

  if( watch_k_ == NORMAL_LINE_SEARCH ) {
    const value_type
      opt_kkt_err_k	= s.opt_kkt_err().get_k(0),
      feas_kkt_err_k	= s.feas_kkt_err().get_k(0);
    if( opt_kkt_err_k <= opt_kkt_err_threshold() && feas_kkt_err_k <= feas_kkt_err_threshold() ) {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out	<< "\nopt_kkt_err_k = " << opt_kkt_err_k << " <= opt_kkt_err_threshold = "
            << opt_kkt_err_threshold() << std::endl
          << "\nfeas_kkt_err_k = " << feas_kkt_err_k << " <= feas_kkt_err_threshold = "
            << feas_kkt_err_threshold() << std::endl
          << "\nSwitching to watchdog linesearch ...\n";
      }
      watch_k_ = 0;
    }
  }

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    out	<< "\nTrial point:\n"
      << "phi_k   = " << phi_k << std::endl
      << "Dphi_k  = " << Dphi_k << std::endl
      << "phi_kp1 = " << phi_kp1 << std::endl;
  }

  bool	ls_success = true,
      step_return = true;

  switch( watch_k_ ) {
    case 0:
    {
      // Take  a full step
      const value_type phi_cord = phi_k + eta() * Dphi_k;
      const bool accept_step = phi_kp1 <= phi_cord;

      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out	<< "\n*** Zeroth watchdog iteration:\n"
          << "\nphi_kp1 = " << phi_kp1 << ( accept_step ? " <= " : " > " )
            << "phi_k + eta * Dphi_k = " << phi_cord << std::endl;
      }

      if( phi_kp1 > phi_cord ) {
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out	<< "\nAccept this increase for now but watch out next iteration!\n";
        }
        // Save this initial point
        xo_		= x_k;
        fo_		= f_k;
        nrm_co_	= norm_inf( c_k );
        do_		= d_k;
        phio_	= phi_k;
        Dphio_	= Dphi_k;
        phiop1_	= phi_kp1;
        // Slip the update of the penalty parameter
        const value_type mu_k = s.mu().get_k(0);
        s.mu().set_k(+1) = mu_k;
        // Move on to the next step in the watchdog procedure
        watch_k_ = 1;
      }
      else {
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out	<< "\nAll is good!\n";
        }
        // watch_k_ stays 0
      }
      step_return = true;
      break;
    }
    case 1:
    {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out	<< "\n*** First watchdog iteration:\n"
          << "\nDo a line search to determine x_kp1 = x_k + alpha_k * d_k ...\n";
      }
      // Now do a line search but and we require some type of reduction
      const DVectorSlice xd[2] = { x_k(), d_k() };
      MeritFuncCalc1DQuadratic phi_calc_1d( phi_calc, 1, xd, &x_kp1() );
      ls_success = direct_line_search().do_line_search( phi_calc_1d, phi_k
        , &alpha_k, &phi_kp1
        , (int)olevel >= (int)PRINT_ALGORITHM_STEPS ?
          &out : static_cast<std::ostream*>(0)	);

      // If the linesearch failed then the rest of the tests will catch this.

      value_type phi_cord = 0;
      bool test1, test2;

      if(		( test1 = ( phi_k <= phio_ ) )
        || ( test2 = phi_kp1 <= ( phi_cord = phio_ + eta() * Dphio_ ) )		)
      {
        // We will accept this step and and move on.
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out
            << "\nphi_k = " << phi_k << ( test1 ? " <= " : " > " )
              << "phi_km1 = " << phio_ << std::endl
            << "phi_kp1 = " << phi_kp1 << ( test2 ? " <= " : " > " )
              << "phi_km1 + eta * Dphi_km1 = " << phi_cord << std::endl
            << "This is a sufficent reduction so reset watchdog.\n";
        }
        watch_k_ = 0;
        step_return = true;
      }
      else if ( ! ( test1 = ( phi_kp1 <= phio_ ) ) ) {
        // Even this reduction is no good!
        // Go back to original point and do a linesearch from there.
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out
            << "\nphi_kp1 = " << phi_kp1 << " > phi_km1 = " << phio_ << std::endl
            << "This is not good reduction in phi so do linesearch from x_km1\n"
            << "\n* Go back to x_km1: x_kp1 = x_k - alpha_k * d_k ...\n";
        }

        // Go back from x_k to x_km1 for iteration k:
        //
        // x_kp1 = x_km1
        // x_kp1 = x_k - alpha_km1 * d_km1
        //
        // A negative sign for alpha is an indication that we are backtracking.
        //
        s.alpha().set_k(0)			= -1.0;
        s.d().set_k(0).v()			= do_;
        s.f().set_k(+1)				= fo_;

        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out << "Output iteration k ...\n"
            << "k = k+1\n";
        }

        // Output these iteration quantities
        algo.track().output_iteration( algo );	// k
        // Transition to iteration k+1
        s.next_iteration();

        // Take the step from x_k = x_km2 to x_kp1 for iteration k (k+1):
        //
        // x_kp1 = x_km2 + alpha_n * d_km2
        // x_kp1 = x_k   + alpha_n * d_km1
        // x_kp1 = x_n
        //				
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out << "\n* Take the step from x_k = x_km2 to x_kp1 for iteration k (k+1)\n"
            << "Find: x_kp1 = x_k + alpha_k * d_k = x_km2 + alpha_k * d_km2\n ...\n";
        }

        // alpha_k = 1.0
        value_type &alpha_k = s.alpha().set_k(0) = 1.0;
        
        // /////////////////////////////////////
        // Compute Dphi_k and phi_k

        // x_k
        const DVector &x_k								= xo_;

        // d_k
        const DVector &d_k = s.d().set_k(0).v()			= do_;

        // Dphi_k
        const value_type &Dphi_k						= Dphio_;

        // phi_k
        const value_type &phi_k = s.phi().set_k(0)		= phio_;

        // Here f_kp1, and c_kp1 are updated at the same time the
        // line search is being performed.
        algo.nlp().set_f( &s.f().set_k(+1) );
        algo.nlp().set_c( &s.c().set_k(+1).v() );
        phi_calc.set_nlp( algo.get_nlp() );

        // ////////////////////////////////////////
        // Compute x_xp1 and ph_kp1 for full step

        // x_kp1 = x_k + alpha_k * d_k
        DVector &x_kp1 = s.x().set_k(+1).v();
        V_VpV( &x_kp1, x_k, d_k );

        // phi_kp1
        value_type &phi_kp1 = s.phi().set_k(+1)			= phiop1_;

        const DVectorSlice xd[2] = { x_k(), d_k() };
        MeritFuncCalc1DQuadratic phi_calc_1d( phi_calc, 1, xd, &x_kp1() );
        ls_success = direct_line_search().do_line_search(
            phi_calc_1d, phi_k
          , &alpha_k, &phi_kp1
          , (int)olevel >= (int)PRINT_ALGORITHM_STEPS ?
            &out : static_cast<std::ostream*>(0)	);

        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out << "\nOutput iteration k (k+1) ...\n"
            << "k = k+1 (k+2)\n"
            << "Reinitialize watchdog algorithm\n";
        }

        // Output these iteration quantities
        algo.track().output_iteration( algo );	// (k+1)
        // Transition to iteration k+1 (k+2)
        s.next_iteration();

        watch_k_ = 0; // Reinitialize the watchdog

        // Any update for k (k+2) should use the last updated value
        // which was for k-2 (k) since there is not much info for k-1 (k+1).
        // Be careful here and make sure this is square with other steps.

        algo.do_step_next( EvalNewPoint_name );
        step_return = false;	// Redirect control
      }
      else {
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out
            << "phi_kp1 = " << phi_kp1 << " <= phi_km1 = " << phio_ << std::endl
            << "\nAccept this step but do a linesearch next iteration!\n";
        }
        // Slip the update of the penalty parameter
        const value_type mu_k = s.mu().get_k(0);
        s.mu().set_k(+1) = mu_k;
        // Do the last stage of the watchdog procedure next iteration.
        watch_k_ = 2;
        step_return = true;
      }
      break;
    }
    case NORMAL_LINE_SEARCH:
    case 2:
    {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        if( watch_k_ == 2 ) {
          out	<< "\n*** Second watchdog iteration:\n"
            << "Do a line search to determine x_kp1 = x_k + alpha_k * d_k ...\n";
        }
        else {
          out	<< "\n*** Normal linesearch:\n"
            << "Do a line search to determine x_kp1 = x_k + alpha_k * d_k ...\n";
        }
      }

      const DVectorSlice xd[2] = { x_k(), d_k() };
      MeritFuncCalc1DQuadratic phi_calc_1d( phi_calc, 1, xd, &x_kp1() );
      ls_success = direct_line_search().do_line_search( phi_calc_1d, phi_k
        , &alpha_k, &phi_kp1
        , (int)olevel >= (int)PRINT_ALGORITHM_STEPS ?
          &out : static_cast<std::ostream*>(0)	);

      if( watch_k_ == 2 )
        watch_k_ = 0;

      step_return = true;
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Only local programming error
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\nalpha    = " << s.alpha().get_k(0) << "\n";
    out << "\nphi_kp1 = " << s.phi().get_k(+1) << "\n";
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nd_k = \n" << s.d().get_k(0)();
    out << "\nx_kp1 = \n" << s.x().get_k(+1)();
  }

  if( !ls_success )
    throw LineSearchFailure("LineSearchWatchDog_Step::do_step(): Line search failure");

  return step_return;

}

void MoochoPack::LineSearchWatchDog_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out	<< L << "*** Use the Watchdog linesearch when near solution.\n"
    << L << "default: opt_kkt_err_threshold = 0.0\n"
    << L << "         feas_kkt_err_threshold = 0.0\n"
    << L << "         eta = 1.0e-4\n"
    << L << "         watch_k = NORMAL_LINE_SEARCH\n"
    << L << "begin definition of NLP merit function phi.value(f(x),c(x)):\n";

  merit_func().print_merit_func( out, L + "    " );
  
  out	<< L << "end definition\n"
    << L << "Dphi_k = phi.deriv()\n"
    << L << "if Dphi_k >= 0 then\n"
    << L << "    throw line_search_failure\n"
    << L << "end\n"
    << L << "phi_kp1 = phi_k.value(f_kp1,c_kp1)\n"
    << L << "phi_k = phi.value(f_k,c_k)\n"
    << L << "if watch_k == NORMAL_LINE_SEARCH then\n"
    << L << "    if opt_kkt_err <= opt_kkt_err_threshold\n"
    << L << "      and feas_kkt_err <= feas_kkt_err_threshold then\n"
    << L << "        *** Start using watchdog from now on!\n"
    << L << "        watch_k = 0\n"
    << L << "    end\n"
    << L << "end\n"
    << L << "if watch_k == 0 then\n"
    << L << "    *** Zeroth watchdog iteration\n"
    << L << "    if phi_kp1 >= phi_k + eta * Dphi_k then\n"
    << L << "        *** Accept this increase for now but watch out next iteration!\n"
    << L << "        *** Save the first point\n"
    << L << "        xo     = x_k\n"
    << L << "        fo     = f_k\n"
    << L << "        nrm_co = norm_inf_c_k\n"
    << L << "        do     = d_k\n"
    << L << "        phio   = phi_k\n"
    << L << "        Dphio  = Dphi_k\n"
    << L << "        phiop1 = phi_kp1\n"
    << L << "        *** Skip the update of the penalty parameter next iteration.\n"
    << L << "        mu_kp1 = mu_k\n"
    << L << "        *** Continue with next step in watchdog\n"
    << L << "        watch_k = 1\n"
    << L << "    else\n"
    << L << "        *** This is a good step so take it!\n"
    << L << "    end\n"
    << L << "else if watch_k == 1 then\n"
    << L << "    *** First watchdog iteration\n"
    << L << "    Do line search for: x_kp1 = x_k + alpha_k + d_k\n"
    << L << "        -> alpha_k, x_kp1, f_kp1, c_kp1, phi_kp1\n"
    << L << "    if ( phi_k <= phio ) or ( phi_kp1 <= phio + eta * Dphio ) then\n"
    << L << "        *** We will accept this step and reinitialize the watchdog\n"
    << L << "        watch_k = 0\n"
    << L << "    else if ( phi_kp1 > phio ) then\n"
    << L << "        *** This reduction is no good!\n"
    << L << "        *** Go back from x_k to x_km1 for this iteration (k)\n"
    << L << "        alpha_k        = -1.0\n"
    << L << "        d_k            = do\n"
    << L << "        f_kp1          = fo\n"
    << L << "        Output this iteration (k)\n"
    << L << "        k = k+1\n"
    << L << "        *** Go from x_k = x_km2 to x_kp1 for this iteration (k+1)\n"
    << L << "        alpha_k        = 1\n"
    << L << "        x_k            = xo\n"
    << L << "        d_k            = do\n"
    << L << "        Dphi_k         = Dphio\n"
    << L << "        phi_k          = phio\n"
    << L << "        x_kp1          = x_k + d_k\n"
    << L << "        phi_kp1        = phiop1\n"
    << L << "        Do line search for: x_kp1 = x_k + alpha_k + d_k\n"
    << L << "            -> alpha_k, x_kp1, f_kp1, c_kp1, phi_kp1\n"
    << L << "        Output this iteration (k+1)\n"
    << L << "        k = k+1\n"
    << L << "        *** Any updates for k (k+2) should use the last updated value\n"
    << L << "        *** which was for k-2 (k) since there is not much info for k-1 (k+1).\n"
    << L << "        *** Be careful here and make sure this works with other steps.\n"
    << L << "        goto EvalNewPoint\n"
    << L << "    else\n"
    << L << "        *** Accept this reduction but do a linesearch next iteration!\n"
    << L << "        *** Skip the update of the penalty parameter next iteration.\n"
    << L << "        mu_kp1 = mu_k\n"
    << L << "        *** Continue with next step in watchdog\n"
    << L << "        watch_k = 2\n"
    << L << "    end\n"
    << L << "else if ( watch_k == 2 ) then\n"
    << L << "    *** Second watchdog iteration\n"
    << L << "    Do line search for: x_kp1 = x_k + alpha_k + d_k\n"
    << L << "        -> alpha_k, x_kp1, f_kp1, c_kp1, phi_kp1\n"
    << L << "    *** Reset the watchdog algorithm\n"
    << L << "    watch_k = 1\n"
    << L << "else if ( watch_k == NORMAL_LINE_SEARCH ) then\n"
    << L << "    Do line search for: x_kp1 = x_k + alpha_k + d_k\n"
    << L << "        -> alpha_k, x_kp1, f_kp1, c_kp1, phi_kp1\n"
    << L << "    begin direct line search : \""
        << typeName(direct_line_search()) << "\"\n";

  direct_line_search().print_algorithm( out, L + "    " );

  out
    << L << "    end direct line search\n"
    << L << "end\n"
    << L << "if maximum number of linesearch iterations are exceeded then\n"
    << L << "    throw line_search_failure\n"
    << L << "end\n";
}

#endif // 0
