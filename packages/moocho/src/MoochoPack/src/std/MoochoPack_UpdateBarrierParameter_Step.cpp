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
#include <iostream>
#include <math.h>

#include "MoochoPack_UpdateBarrierParameter_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "MoochoPack_IpState.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "Teuchos_Assert.hpp"

#include "OptionsFromStreamPack_StringToBool.hpp"

#define min(a,b) ( (a < b) ? a : b )
#define max(a,b) ( (a > b) ? a : b )

namespace MoochoPack {

UpdateBarrierParameter_Step::UpdateBarrierParameter_Step(
  const value_type init_barrier_parameter,
  const value_type tau_mu,
  const value_type theta_mu,
  const value_type tau_epsilon,
  const value_type theta_epsilon,
  const value_type e_tol_max
  )
  :
  init_barrier_parameter_(init_barrier_parameter),
  tau_mu_(tau_mu),
  theta_mu_(theta_mu),
  tau_epsilon_(tau_epsilon),
  theta_epsilon_(theta_epsilon),
  e_tol_max_(e_tol_max)
  {}
  

bool UpdateBarrierParameter_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
  {
  using Teuchos::dyn_cast;
  using IterationPack::print_algorithm_step;

  NLPAlgo            &algo   = dyn_cast<NLPAlgo>(_algo);
  IpState             &s      = dyn_cast<IpState>(_algo.state());
  NLP                 &nlp    = algo.nlp();

  if (!nlp.is_initialized())
    {
    nlp.initialize(false);
    }
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();
  
  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
    }

  
  ///***********************************************************
  // Get iteration quantities
  ///***********************************************************
  IterQuantityAccess<value_type>& e_tol_iq = s.e_tol();
  IterQuantityAccess<value_type>& mu_iq = s.barrier_parameter();

  ///***********************************************************
  // Check values and initialize, if necessary
  ///***********************************************************
  /*	if (mu_iq.last_updated() == IterQuantity::NONE_UPDATED)
    {
    mu_iq.set_k(0) = init_barrier_parameter_;
    e_tol_iq.set_k(0) = Calculate_e_tol(mu_iq.get_k(0));

    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
      {
      out << "\nInitializing barrier parameter (mu) and sub problem tolerance (e_tol) ...\n";
      }
    }
    else*/
    {
    ///***********************************************************
    // if tolerance is satisfied, calculate new barrier_parameter 
    //  and e_tol, otherwise update mu and e_tol from last iter
    ///***********************************************************
    const value_type opt_err = s.opt_kkt_err().get_k(0);
    const value_type feas_err = s.feas_kkt_err().get_k(0);
    const value_type comp_err_mu = s.comp_err_mu().get_k(0);

    const value_type mu_km1 = mu_iq.get_k(-1);
    if (e_tol_iq.last_updated() == IterQuantity::NONE_UPDATED)
      {
      // First time through, don't let mu change
      mu_iq.set_k(0,-1);
      e_tol_iq.set_k(0) = Calculate_e_tol(mu_iq.get_k(0));
          }
    else
      {
      const value_type e_tol_km1 = e_tol_iq.get_k(-1);
      bool sub_prob_converged = (opt_err < e_tol_km1 && feas_err < e_tol_km1 && comp_err_mu < e_tol_km1); 
      if (sub_prob_converged)
        {
        // Calculate new mu and e_tol
        value_type& mu_k = mu_iq.set_k(0);
        mu_k = min(tau_mu_*mu_km1, pow(mu_km1, theta_mu_));
        e_tol_iq.set_k(0) = Calculate_e_tol(mu_k);
        
        if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
          {
          out << "\nSub-problem converged!\n"
            << " Updating barrier parameter (mu) and sub problem tolerance (e_tol) ...\n";
          }
        }
      else
        {
        if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
          {
          out << "\nSub-problem not-converged!\n"
            << " Keeping existing barrier parameter (mu) and sub problem tolerance (e_tol) ...\n";
          }
        mu_iq.set_k(0,-1);
        e_tol_iq.set_k(0,-1);
        }
      }

    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
      {
      out << "\nbarrier_parameter (mu) = " << mu_iq.get_k(0)
        << "\nsub problem tolerance (e_tol) = " << e_tol_iq.get_k(0)  << std::endl;
      }
    }

  return true;
  }


void UpdateBarrierParameter_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
  {
  //const NLPAlgo   &algo = rsqp_algo(_algo);
  //const NLPAlgoState  &s    = algo.rsqp_state();
  out << L << "# Update the interior point barrier parameter (mu)\n"
    << L << "if (KKTerror < e_tol) then\n"
    << L << "   mu_kp1 = min(tau_mu*mu_k, mu_k^theta_mu)\n"
    << L << "   e_tol_kp1 = min(e_tol_max, tau_epsilon*min(mu_k, mu_k^theta_epsilon))\n"
    << L << "else\n"
    << L << "   mu_kp1 = mu_k\n"
    << L << "   e_tol_kp1 = e_tol_k\n"
    << L << "end;\n";
  }

value_type UpdateBarrierParameter_Step::Calculate_e_tol(value_type mu)
  {	
  value_type e_tol = tau_epsilon_*min(mu, pow(mu, theta_epsilon_));
  e_tol = min(e_tol_max_, e_tol);

  return e_tol;
  }


namespace {

const int local_num_options = 5;

enum local_EOptions 
  {
    TAU_MU,
    THETA_MU,
    TAU_EPSILON,
    THETA_EPSILON,
    E_TOL_MAX
  };

const char* local_SOptions[local_num_options] = 
  {
    "tau_mu",
    "theta_mu",
    "tau_epsilon",
    "theta_epsilon",
    "e_tol_max"
  };

}

 
UpdateBarrierParameter_StepSetOptions::UpdateBarrierParameter_StepSetOptions(
  UpdateBarrierParameter_Step* target
  , const char opt_grp_name[] )
  :
  OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions ),
  OptionsFromStreamPack::SetOptionsToTargetBase< UpdateBarrierParameter_Step >( target )
  {
  }

void UpdateBarrierParameter_StepSetOptions::setOption( 
  int option_num, const std::string& option_value )
  {
  using OptionsFromStreamPack::StringToBool;
  
  typedef UpdateBarrierParameter_Step target_t;
  switch( (local_EOptions)option_num ) 
    {
    case TAU_MU:
      target().tau_mu(std::atof(option_value.c_str()));
      break;
    case THETA_MU:
      target().theta_mu(std::atof(option_value.c_str()));
      break;
    case TAU_EPSILON:
      target().tau_epsilon(std::atof(option_value.c_str()));
      break;
    case THETA_EPSILON:
      target().theta_epsilon(std::atof(option_value.c_str()));
      break;
    case E_TOL_MAX:
      target().e_tol_max(std::atof(option_value.c_str()));
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
    }
  }

} // end namespace MoochoPack 
