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

#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "NLPInterfacePack_NLPBarrier.hpp"
#include "MoochoPack_PreProcessBarrierLineSearch_Step.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

#define min(a,b) ( (a < b) ? a : b )
#define max(a,b) ( (a > b) ? a : b )

namespace MoochoPack {

PreProcessBarrierLineSearch_Step::PreProcessBarrierLineSearch_Step(
  Teuchos::RCP<NLPInterfacePack::NLPBarrier> barrier_nlp,
  const value_type tau_boundary_frac
  )
  :
  tau_boundary_frac_(tau_boundary_frac),
  barrier_nlp_(barrier_nlp),
  filter_(FILTER_IQ_STRING)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !barrier_nlp_.get(),
    std::logic_error,
    "PreProcessBarrierLineSearch_Step given NULL NLPBarrier."
    );
}
  

bool PreProcessBarrierLineSearch_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  using Teuchos::dyn_cast;
  using IterationPack::print_algorithm_step;
    using AbstractLinAlgPack::assert_print_nan_inf;
  using AbstractLinAlgPack::fraction_to_boundary;
  using AbstractLinAlgPack::fraction_to_zero_boundary;
  using LinAlgOpPack::Vp_StV;

  NLPAlgo            &algo   = dyn_cast<NLPAlgo>(_algo);
  IpState             &s      = dyn_cast<IpState>(_algo.state());
  NLP                 &nlp    = algo.nlp();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();
  
  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
  {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
  }

  const value_type& mu_k = s.barrier_parameter().get_k(0);

  // if using filter and u changed, clear filter
  if (filter_.exists_in(s))
  {
    if ( s.barrier_parameter().updated_k(-1) )
    {
      const value_type mu_km1 = s.barrier_parameter().get_k(-1);
      if (mu_k != mu_km1)
      {
        if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
        {
          out << "\nBarrier Parameter changed - resetting the filter ...\n";
        }
        // reset the filter
        MoochoPack::Filter_T &filter_k = filter_(s).set_k(0);
        filter_k.clear();
      }
    }
  }

  // Update the barrier parameter in the NLP
  barrier_nlp_->mu(s.barrier_parameter().get_k(0));
    
  // Calculate the barrier k terms
  barrier_nlp_->set_Gf( &(s.grad_barrier_obj().set_k(0)) );
  barrier_nlp_->calc_Gf(s.x().get_k(0), true);

  barrier_nlp_->set_f( &(s.barrier_obj().set_k(0)) );
  barrier_nlp_->calc_f(s.x().get_k(0), true);


  // Calculate the k+1 terms
  // Get iteration quantities...
  value_type& alpha_k = s.alpha().set_k(0);
  value_type& alpha_vl_k = s.alpha_vl().set_k(0);
  value_type& alpha_vu_k = s.alpha_vu().set_k(0);

  const Vector& x_k = s.x().get_k(0);
  VectorMutable& x_kp1 = s.x().set_k(+1);

  const Vector& d_k = s.d().get_k(0);
  const Vector& dvl_k = s.dvl().get_k(0);
  const Vector& dvu_k = s.dvu().get_k(0);

  const Vector& vl_k = s.Vl().get_k(0).diag();
  VectorMutable& vl_kp1 = s.Vl().set_k(+1).diag();

  const Vector& vu_k = s.Vu().get_k(0).diag();
  VectorMutable& vu_kp1 = s.Vu().set_k(+1).diag();

  alpha_k = fraction_to_boundary(
    tau_boundary_frac_, 
    x_k, 
    d_k,
    nlp.xl(),
    nlp.xu()
    );

  alpha_vl_k = fraction_to_zero_boundary(
    tau_boundary_frac_,
    vl_k,
    dvl_k
    );

  alpha_vu_k = fraction_to_zero_boundary(
    tau_boundary_frac_,
    vu_k,
    dvu_k
    );

  TEUCHOS_TEST_FOR_EXCEPT( !( alpha_k <= 1.0 && alpha_vl_k <= 1.0 && alpha_vu_k <= 1.0 ) );
  TEUCHOS_TEST_FOR_EXCEPT( !( alpha_k >= 0.0 && alpha_vl_k >= 0.0 && alpha_vu_k >= 0.0 ) );

  x_kp1 = x_k;
  Vp_StV( &x_kp1, alpha_k, d_k);

  alpha_vl_k = alpha_vu_k = min(alpha_vl_k, alpha_vu_k);

  vl_kp1 = vl_k;
  Vp_StV( &vl_kp1, alpha_vl_k, dvl_k);

  vu_kp1 = vu_k;
  Vp_StV( &vu_kp1, alpha_vu_k, dvu_k);
 

    IterQuantityAccess<VectorMutable>
    *c_iq   = nlp.m() > 0 ? &s.c() : NULL;

    if (assert_print_nan_inf(x_kp1, "x", true, NULL))
  {
    // Calcuate f and c at the new point.
    barrier_nlp_->unset_quantities();
    barrier_nlp_->set_f( &s.barrier_obj().set_k(+1) );
    if (c_iq) {
      barrier_nlp_->set_c( &c_iq->set_k(+1) );
      barrier_nlp_->calc_c( x_kp1, true );
    }
    barrier_nlp_->calc_f( x_kp1, false ); 
    barrier_nlp_->unset_quantities();
  }
  
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) )
  {
    out << "\nalpha_vl_k = " << alpha_vl_k
      << "\nalpha_vu_k = " << alpha_vu_k
      << "\nalpha_k    = " << alpha_k
      << std::endl;
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
  {
    out << "\nvl_kp1 = \n" << vl_kp1
      << "\nvu_kp1 = \n" << vu_kp1
      << "\nx_kp1 = \n" << x_kp1;
  }

  return true;
}

void PreProcessBarrierLineSearch_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
{
  //const NLPAlgo   &algo = rsqp_algo(_algo);
  //const NLPAlgoState  &s    = algo.rsqp_state();
  out << L << "*** calculate alpha max by the fraction to boundary rule\n"
    << L << "ToDo: Complete documentation\n";
}

namespace {

const int local_num_options = 1;

enum local_EOptions 
{
  TAU_BOUNDARY_FRAC
};

const char* local_SOptions[local_num_options] = 
{
  "tau_boundary_frac"
};

}

 
PreProcessBarrierLineSearch_StepSetOptions::PreProcessBarrierLineSearch_StepSetOptions(
  PreProcessBarrierLineSearch_Step* target
  , const char opt_grp_name[] )
  :
  OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions ),
  OptionsFromStreamPack::SetOptionsToTargetBase< PreProcessBarrierLineSearch_Step >( target )
{}

void PreProcessBarrierLineSearch_StepSetOptions::setOption( 
  int option_num, const std::string& option_value )
{
  typedef PreProcessBarrierLineSearch_Step target_t;
  switch( (local_EOptions)option_num ) 
  {
    case TAU_BOUNDARY_FRAC:
      target().tau_boundary_frac(std::atof(option_value.c_str()));
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

} // end namespace MoochoPack 
