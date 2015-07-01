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

#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_PreEvalNewPointBarrier_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"

#include "OptionsFromStreamPack_StringToBool.hpp"

#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

PreEvalNewPointBarrier_Step::PreEvalNewPointBarrier_Step(
  const value_type relative_bound_push,
  const value_type absolute_bound_push
  )
  :
  relative_bound_push_(relative_bound_push),
  absolute_bound_push_(absolute_bound_push)
  {}

bool PreEvalNewPointBarrier_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
  {
  using Teuchos::dyn_cast;
  using IterationPack::print_algorithm_step;
  using AbstractLinAlgPack::force_in_bounds_buffer;

  NLPAlgo            &algo   = dyn_cast<NLPAlgo>(_algo);
  IpState             &s      = dyn_cast<IpState>(_algo.state());
  NLP                 &nlp    = algo.nlp();
  NLPFirstOrder   *nlp_foi = dynamic_cast<NLPFirstOrder*>(&nlp);
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();
  
  if(!nlp.is_initialized())
    nlp.initialize(algo.algo_cntr().check_results());

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
    }

  IterQuantityAccess<value_type>     &barrier_parameter_iq = s.barrier_parameter();
  IterQuantityAccess<VectorMutable>  &x_iq  = s.x();

  if( x_iq.last_updated() == IterQuantity::NONE_UPDATED ) 
    {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
      {
      out << "\nInitialize x with x_k = nlp.xinit() ...\n"
        << " and push x_k within bounds.\n";
      }
    VectorMutable& x_k = x_iq.set_k(0) = nlp.xinit();
  
    // apply transformation operator to push x sufficiently within bounds
    force_in_bounds_buffer(relative_bound_push_, 
                 absolute_bound_push_,
                 nlp.xl(),
                 nlp.xu(),
                 &x_k);

    // evaluate the func and constraints
    IterQuantityAccess<value_type>
      &f_iq    = s.f();
    IterQuantityAccess<VectorMutable>
      &Gf_iq   = s.Gf(),
      *c_iq    = nlp.m() > 0 ? &s.c() : NULL;
    IterQuantityAccess<MatrixOp>
      *Gc_iq   = nlp_foi ? &s.Gc() : NULL;

    using AbstractLinAlgPack::assert_print_nan_inf;
    assert_print_nan_inf(x_k, "x", true, NULL); // With throw exception if Inf or NaN!

    // Wipe out storage for computed iteration quantities (just in case?) : RAB: 7/29/2002
    if(f_iq.updated_k(0))
      f_iq.set_not_updated_k(0);
    if(Gf_iq.updated_k(0))
      Gf_iq.set_not_updated_k(0);
    if (c_iq)
      {
      if(c_iq->updated_k(0))
        c_iq->set_not_updated_k(0);
      }
    if (nlp_foi)
      {
      if(Gc_iq->updated_k(0))
        Gc_iq->set_not_updated_k(0);
      }
    }

  if (barrier_parameter_iq.last_updated() == IterQuantity::NONE_UPDATED)
    {
    barrier_parameter_iq.set_k(-1) = 0.1; // RAB: 7/29/2002: This should be parameterized! (allow user to set this!)
    }

  // Print vector information
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
    {
    out << "x_k =\n" << x_iq.get_k(0);
    }

  return true;
  }


void PreEvalNewPointBarrier_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
  {
  //const NLPAlgo   &algo = rsqp_algo(_algo);
  //const NLPAlgoState  &s    = algo.rsqp_state();
  out << L << "# Evaluate information specific to primal / dual barrier algorithms\n"
    << L << "if (x never updated) then\n"
    << L << "  x_k = nlp.xinit()\n"
    << L << "  force_in_bounds(x_k)\n"
    << L << "  set f_k not updated\n"
    << L << "  set Gf_k not updated\n"
    << L << "  if (m > 0) then\n"
    << L << "    set c_k not updated\n"
    << L << "    set Gc_k not updated\n"
    << L << "end\n";
  }


namespace {

const int local_num_options = 2;

enum local_EOptions 
  {
  RELATIVE_BOUND_PUSH,
  ABSOLUTE_BOUND_PUSH
  };

const char* local_SOptions[local_num_options] = 
  {
  "relative_bound_push",
  "absolute_bound_push"
  };

}

 
PreEvalNewPointBarrier_StepSetOptions::PreEvalNewPointBarrier_StepSetOptions(
  PreEvalNewPointBarrier_Step* target
  , const char opt_grp_name[] )
  :
  OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions ),
  OptionsFromStreamPack::SetOptionsToTargetBase< PreEvalNewPointBarrier_Step >( target )
  {
  }

void PreEvalNewPointBarrier_StepSetOptions::setOption( 
  int option_num, const std::string& option_value )
  {
  using OptionsFromStreamPack::StringToBool;

  typedef PreEvalNewPointBarrier_Step target_t;
  switch( (local_EOptions)option_num ) 
    {
    case RELATIVE_BOUND_PUSH:
      target().relative_bound_push(std::atof(option_value.c_str()));
      break;
    case ABSOLUTE_BOUND_PUSH:
      target().absolute_bound_push(std::atof(option_value.c_str()));
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
    }
  }

}  // end namespace MoochoPack
