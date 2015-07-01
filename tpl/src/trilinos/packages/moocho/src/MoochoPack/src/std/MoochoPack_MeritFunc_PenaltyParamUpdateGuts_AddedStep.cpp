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

#include "MoochoPack_MeritFunc_PenaltyParamUpdateGuts_AddedStep.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_MeritFuncPenaltyParam.hpp"
#include "ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"

namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
} // end namespace

namespace MoochoPack {

MeritFunc_PenaltyParamUpdateGuts_AddedStep::MeritFunc_PenaltyParamUpdateGuts_AddedStep(
  value_type                     small_mu
  ,value_type                    mult_factor
  ,value_type                    kkt_near_sol
  )
  :near_solution_(false)
  ,small_mu_(small_mu)
  ,mult_factor_(mult_factor)
  ,kkt_near_sol_(kkt_near_sol)
{}

bool MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  NLP			&nlp	= algo.nlp();
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  const size_type
    //n  = nlp.n(),
    m  = nlp.m();
  IterQuantityAccess<MeritFuncNLP>
    &merit_func_nlp_iq = s.merit_func_nlp();

  if( !merit_func_nlp_iq.updated_k(0) ) {
    const int merit_func_k_last_updated = merit_func_nlp_iq.last_updated();
    if( merit_func_k_last_updated != IterQuantity::NONE_UPDATED ) {
      MeritFuncNLP
        &merit_func_nlp_k_last = merit_func_nlp_iq.get_k(merit_func_k_last_updated);
      merit_func_nlp_iq.set_k(0) = merit_func_nlp_k_last;
    }
    else {
      merit_func_nlp_iq.set_k(0); // Just use default constructor
    }
    MeritFuncNLP
      &merit_func_nlp_k = merit_func_nlp_iq.get_k(0);
    MeritFuncPenaltyParam
      *param = dynamic_cast<MeritFuncPenaltyParam*>(&merit_func_nlp_k);
    TEUCHOS_TEST_FOR_EXCEPTION(
      !param, std::logic_error
      ,"MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(...), Error "
      << "The class " << typeName(merit_func_nlp_k) << " does not support the "
      << "MeritFuncPenaltyParam iterface" );
    MeritFuncNLPDirecDeriv
      *direc_deriv = dynamic_cast<MeritFuncNLPDirecDeriv*>(&merit_func_nlp_k);
    TEUCHOS_TEST_FOR_EXCEPTION(
      !direc_deriv, std::logic_error
      ,"MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(...), Error "
      << "The class " << typeName(merit_func_nlp_k) << " does not support the "
      << "MeritFuncNLPDirecDeriv iterface" );
    value_type  new_mu = 0.0;
    value_type  min_mu = 0.0;
    if ( this->min_mu(s,&min_mu) ) {
      // Update the penalty parameter as defined in the fortran rSQP code (EXACT2())
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nUpdate the penalty parameter...\n";
      }
      value_type
        mu_km1 = param->mu(),
        mult_fact = (1.0 + mult_factor_);
      if(near_solution_) {
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out << "\nNear solution, forcing mu_k >= mu_km1...\n";
        }
        new_mu = my_max( my_max( mu_km1, mult_fact * min_mu ), small_mu_ );
      }
      else {
        if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
          out << "\nNot near solution, allowing reduction in mu ...\n";
        }
        new_mu =	my_max(
          (3.0 * mu_km1 + min_mu) / 4.0	
          , my_max( mult_fact * min_mu, small_mu_ )
          ); 
        value_type
          kkt_error = s.opt_kkt_err().get_k(0) + s.feas_kkt_err().get_k(0);
        if(kkt_error <= kkt_near_sol_) {
          if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
            out << "\nkkt_error = " << kkt_error << " <= kkt_near_sol = "
              << kkt_near_sol_ << std::endl
              << "Switching to forcing mu_k >= mu_km1 in the future\n";
          }
          near_solution_ = true;
        }
      }
    }
    else {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nDon't have the info to update penalty parameter so just use the last updated...\n";
      }
      new_mu = param->mu();
    }
    // Set the penalty parameter
    param->mu( new_mu );
    // In addition also compute the directional derivative
    direc_deriv->calc_deriv(
      s.Gf().get_k(0)
      ,m  ? &s.c().get_k(0) : NULL
      ,NULL   // h
      ,NULL   // hl
      ,NULL   // hu
      ,s.d().get_k(0)
      );
    
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\nmu = " << new_mu << "\n";
    }
  }
  return true;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Update the penalty parameter for the merit function to ensure\n"
    << L << "*** a descent direction a directional derivatieve.\n"
    << L << "*** phi is a merit function object that uses the penalty parameter mu.\n"
    << L << "default: near_solution = false\n"
    << L << "         small_mu = " << small_mu_ << std::endl
    << L << "         mult_factor = " << mult_factor_ << std::endl
    << L << "         kkt_near_sol = " << kkt_near_sol_ << std::endl
    << L << "if merit_func_nlp_k is not already updated then\n"
    << L << "    if some merit_func_nlp_k(?) has been udpated then\n"
    << L << "        merit_func_nlp_k = merit_func_nlp_k(last_udpated)\n"
    << L << "    else\n"
    << L << "        merit_func_nlp_k = default construction\n"
    << L << "    end\n"
    << L << "    if merit_func_nlp_k does not support MeritFuncPenaltyParam throw excpetion\n"
    << L << "    if merit_func_nlp_k does not support MeritFuncNLPDirecDeriv throw excpetion\n"
    ;
              print_min_mu_step( out, L + "    " ); 
  out
    << L << "   mu_new = merit_func_nlp_k.mu()\n"
    << L << "   if update_mu == true then\n"
    << L << "       mu_last = merit_func_nlp_k.mu()\n"
    << L << "       mult_fact = 1.0 + mult_factor\n"
    << L << "       if near_solution == true\n"
    << L << "           mu_new = max( max( mu_last, mult_fact*min_mu ), small_mu )\n"
    << L << "       else\n"
    << L << "           mu_new = max(   ( 3.0 * mu_last + min_mu ) / 4.0\n"
    << L << "                           , max( mult_fact * min_mu , small_mu ) )\n"
    << L << "           kkt_error = opt_kkt_err_k + feas_kkt_err_k\n"
    << L << "           if kkt_error <= kkt_near_sol then\n"
    << L << "               near_solution = true\n"
    << L << "           end\n"
    << L << "       end\n"
    << L << "   else\n"
    << L << "       mu_new = merit_func_nlp_k.mu()\n"
    << L << "   end\n"
    << L << "   merit_func_nlp_k..mu(mu_new)\n"
    << L << "   merit_func_nlp_k.calc_deriv(Gf_k,c_k,h_k,hl,hu,d_k)\n"
    << L << "end\n"
    ;
}

// Overridden from MeritFunc_PenaltyParamUpdate_AddedStep

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::small_mu( value_type small_mu )
{
  small_mu_ = small_mu;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::small_mu() const
{
  return small_mu_;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::mult_factor( value_type mult_factor )
{
  mult_factor_ = mult_factor;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::mult_factor() const
{
  return mult_factor_;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::kkt_near_sol( value_type kkt_near_sol )
{
  kkt_near_sol_ = kkt_near_sol;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::kkt_near_sol() const
{
  return kkt_near_sol_;
}

}	// end namespace MoochoPack
