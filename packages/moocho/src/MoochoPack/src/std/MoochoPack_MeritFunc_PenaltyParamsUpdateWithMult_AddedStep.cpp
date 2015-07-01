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

#include "MoochoPack_MeritFunc_PenaltyParamsUpdateWithMult_AddedStep.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack_MeritFuncPenaltyParams.hpp"
#include "ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"

namespace {

typedef MoochoPack::value_type value_type;
inline value_type max(value_type v1, value_type v2)
{	return (v1 > v2) ? v1 : v2; }

}

namespace MoochoPack {

MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::MeritFunc_PenaltyParamsUpdateWithMult_AddedStep(
      const merit_func_ptr_t& merit_func, value_type small_mu, value_type min_mu_ratio
    , value_type mult_factor, value_type kkt_near_sol )
  : merit_func_(merit_func), near_solution_(false)
    , small_mu_(small_mu), min_mu_ratio_(min_mu_ratio), mult_factor_(mult_factor)
    , kkt_near_sol_(kkt_near_sol), norm_inf_mu_last_(0.0)
{}

bool MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type
  , poss_type assoc_step_poss)
{
  using DenseLinAlgPack::norm_inf;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  MeritFuncPenaltyParams
    *params = dynamic_cast<MeritFuncPenaltyParams*>(&merit_func());
  if( !params ) {
    std::ostringstream omsg;
    omsg
      << "MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::do_step(...), Error "
      << "The class " << typeName(&merit_func()) << " does not support the "
      << "MeritFuncPenaltyParams iterface\n";
    out << omsg.str();
    throw std::logic_error( omsg.str() );
  }

  MeritFuncNLPDirecDeriv
    *direc_deriv = dynamic_cast<MeritFuncNLPDirecDeriv*>(&merit_func());
  if( !direc_deriv ) {
    std::ostringstream omsg;
    omsg
      << "MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::do_step(...), Error "
      << "The class " << typeName(&merit_func()) << " does not support the "
      << "MeritFuncNLPDirecDeriv iterface\n";
    out << omsg.str();
    throw std::logic_error( omsg.str() );
  }

  bool perform_update = true;

  if( s.mu().updated_k(0) ) {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\nmu_k is already updated by someone else?\n";
    }
    const value_type mu_k = s.mu().get_k(0);
    if( mu_k == norm_inf_mu_last_ ) {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nmu_k " << mu_k << " == norm_inf_mu_last = " << norm_inf_mu_last_
          << "\nso we will take this as a signal to skip the update.\n";
      }
      perform_update = false;
    }
    else {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nmu_k " << mu_k << " != norm_inf_mu_last = " << norm_inf_mu_last_
          << "\nso we will ignore this and perform the update anyway.\n";
      }
    }		
  }
  if(perform_update) {

    if ( s.lambda().updated_k(0) ) {

      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nUpdate the penalty parameter...\n";
      }

      const DVector
        &lambda_k = s.lambda().get_k(0).cv();

      if( params->mu().size() != lambda_k.size() )
        params->resize( lambda_k.size() );
      DVectorSlice
        mu = params->mu();

      const value_type
        max_lambda	= norm_inf( lambda_k() ),
        mult_fact	= (1.0 + mult_factor_);

      if(near_solution_) {
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out << "\nNear solution, forcing mu(j) >= mu_old(j)...\n";
        }
        DVector::const_iterator	lb_itr = lambda_k.begin();
        DVectorSlice::iterator	mu_itr = mu.begin();
        for( ; lb_itr != lambda_k.end(); ++mu_itr, ++ lb_itr )
          *mu_itr = max( max( *mu_itr, mult_fact * ::fabs(*lb_itr) ), small_mu_ );
      }
      else {
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out << "\nNot near solution, allowing reduction in mu(j) ...\n";
        }
        DVector::const_iterator	lb_itr = lambda_k.begin();
        DVectorSlice::iterator	mu_itr = mu.begin();
        for( ; lb_itr != lambda_k.end(); ++mu_itr, ++ lb_itr ) {
          const value_type lb_j = ::fabs(*lb_itr);
          *mu_itr = max(
                  (3.0 * (*mu_itr) + lb_j) / 4.0	
                , max( mult_fact * lb_j, small_mu_ )
                );
        }
        value_type kkt_error = s.opt_kkt_err().get_k(0) + s.feas_kkt_err().get_k(0);
        if(kkt_error <= kkt_near_sol_) {
          if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
            out << "\nkkt_error = " << kkt_error << " <= kkt_near_sol = "
                << kkt_near_sol_ << std::endl
              << "Switching to forcing mu_k >= mu_km1 in the future\n";
          }
          near_solution_ = true;
        }
      }

      // Force the ratio
      const value_type
          max_mu	= norm_inf( mu() ),
          min_mu	= min_mu_ratio_ * max_mu;
      for(DVectorSlice::iterator mu_itr = mu.begin(); mu_itr != mu.end(); ++mu_itr)
        *mu_itr = max( (*mu_itr), min_mu );	

      s.mu().set_k(0) = norm_inf_mu_last_ = max_mu;

      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nmax(|mu(j)|) = " << (*std::max_element( mu.begin(), mu.end() ))
          << "\nmin(|mu(j)|) = " << (*std::min_element( mu.begin(), mu.end() ))
            << std::endl;
      }

      if( (int)olevel >= (int)PRINT_VECTORS ) {
        out << "\nmu = \n" << mu;
      }
    }
    else {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nDon't have the info to update penalty parameter so just use the last updated...\n";
      }
    }
  }

  // In addition also compute the directional derivative
  direc_deriv->calc_deriv( s.Gf().get_k(0)(), s.c().get_k(0)(), s.d().get_k(0)() );

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    out << "\nmu_k = " << s.mu().get_k(0) << "\n";
  }

  return true;
}

void MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Update the penalty parameter for the merit function to ensure\n"
    << L << "*** a descent direction a directional derivatieve.\n"
    << L << "*** phi is a merit function object that uses the penalty parameter mu.\n"
    << L << "default: near_solution = false\n"
    << L << "         small_mu = " << small_mu_ << std::endl
    << L << "         min_mu_ratio = " << min_mu_ratio_ << std::endl
    << L << "         mult_factor = " << mult_factor_ << std::endl
    << L << "         kkt_near_sol = " << kkt_near_sol_ << std::endl
    << L << "perform_update = true\n"
    << L << "if mu_k is already updated then\n"
    << L << "    if mu_k == norm_inf_mu_last then\n"
    << L << "        *** We will use this as a signal to skip the update\n"
    << L << "        perform_update = false\n"
    << L << "    else\n"
    << L << "        *** We will perform the update anyway\n"
    << L << "    end\n"
    << L << "if perform_update == true then\n"
    << L << "    if lambda_k is updated then\n"
    << L << "        max_lambda = norm(lambda_k,inf)\n"
    << L << "        mult_fact = (1+mult_factor)\n"
    << L << "        mu = phi.mu()\n"
    << L << "        if near_solution == true\n"
    << L << "            for j = 1...m\n"
    << L << "                mu(j) = max(max(mu(j),mult_fact*abs(lambda_k(j))),small_mu)\n"
    << L << "            end\n"
    << L << "        else\n"
    << L << "            for j = 1...m\n"
    << L << "                mu(j) = max( ( 3.0 * mu(j) + abs(lambda_k(j)) ) / 4.0\n"
    << L << "                            , max( 1.001 * abs(lambda_k(j)) , small_mu ) )\n"
    << L << "            end\n"
    << L << "            kkt_error = opt_kkt_err_k + feas_kkt_err_k\n"
    << L << "            if kkt_error <= kkt_near_sol then\n"
    << L << "                near_solution = true\n"
    << L << "            end\n"
    << L << "        end\n"
    << L << "        min_mu = min_mu_ratio * norm(mu,inf)\n"
    << L << "        for j = 1...m\n"
    << L << "            mu(j) = max( mu(j), min_mu )\n"
    << L << "        end\n"
    << L << "    else\n"
    << L << "        *** Don't have the information to perform the update.\n"
    << L << "    end\n"
    << L << "end\n"
    << L << "phi.calc_deriv(Gf_k,c_k,d_k)\n";
}

// Overridden from MeritFunc_PenaltyParamUpdate_AddedStep

void MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::small_mu( value_type small_mu )
{
  small_mu_ = small_mu;
}

value_type MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::small_mu() const
{
  return small_mu_;
}

void MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::min_mu_ratio( value_type min_mu_ratio )
{
  min_mu_ratio_ = min_mu_ratio;
}

value_type MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::min_mu_ratio() const
{
  return min_mu_ratio_;
}

void MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::mult_factor( value_type mult_factor )
{
  mult_factor_ = mult_factor;
}

value_type MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::mult_factor() const
{
  return mult_factor_;
}

void MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::kkt_near_sol( value_type kkt_near_sol )
{
  kkt_near_sol_ = kkt_near_sol;
}

value_type MeritFunc_PenaltyParamsUpdateWithMult_AddedStep::kkt_near_sol() const
{
  return kkt_near_sol_;
}


}	// end namespace MoochoPack

#endif // 0
