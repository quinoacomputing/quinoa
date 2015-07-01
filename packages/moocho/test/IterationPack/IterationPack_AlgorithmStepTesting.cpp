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

#include <iomanip>

#include "IterationPack_AlgorithmStepTesting.hpp"
#include "IterationPack_Algorithm.hpp"
#include "IterationPack_print_algorithm_step.hpp"

namespace {
char step_type_name[3][15] = { "DO_MAIN_STEP", "DO_PRE_STEP" , "DO_POST_STEP" };
} // namespace

namespace IterationPack {

bool AlgorithmStepTesting::do_step(Algorithm& algo, poss_type step_poss, EDoStepType type
    , poss_type assoc_step_poss)
{
  print_algorithm_step( algo, step_poss, type, assoc_step_poss
    , algo.track().journal_out() );
  return true;
}

void AlgorithmStepTesting::initialize_step(
  Algorithm& algo, poss_type step_poss, EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  print_step_poss( algo, step_poss, type, assoc_step_poss );
  algo.track().journal_out() << "\" : initialize_step(...) called\n";
}

void AlgorithmStepTesting::inform_updated(
  Algorithm& algo, poss_type step_poss, EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  print_step_poss( algo, step_poss, type, assoc_step_poss );
  algo.track().journal_out() << "\" : inform_step(...) called\n";
}

void AlgorithmStepTesting::finalize_step(
  Algorithm& algo, poss_type step_poss, EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  print_step_poss( algo, step_poss, type, assoc_step_poss );
  algo.track().journal_out() << "\" : finalize_step(...) called\n";
}

void AlgorithmStepTesting::print_step( const Algorithm& algo, poss_type step_poss, EDoStepType type
  , poss_type assoc_step_poss ,std::ostream& out, const std::string& leading_str ) const
{
  algo.track().journal_out()
    << std::endl << leading_str << step_poss << ", " << step_type_name[type];
  if(type == DO_MAIN_STEP) {
    algo.track().journal_out()
      << ", \"" << algo.get_step_name(step_poss) << "\"";
  }
  else {
    EAssocStepType _type = (type == DO_PRE_STEP ? PRE_STEP : POST_STEP );
    algo.track().journal_out()
      << ", " << assoc_step_poss
      << ", \"" << algo.get_assoc_step_name(step_poss,_type,assoc_step_poss) << "\"";
  }
  algo.track().journal_out()
    << "\" : print_step(algo,step_poss,type,assoc_step_poss,out) called\n";
}

// private

void AlgorithmStepTesting::print_step_poss(
  const Algorithm& algo, poss_type step_poss, EDoStepType type
  ,poss_type assoc_step_poss
  ) const
{
  algo.track().journal_out()
    << std::endl << step_poss << ", " << step_type_name[type];
  if(type == DO_MAIN_STEP) {
    algo.track().journal_out()
      << ", \"" << algo.get_step_name(step_poss) << "\"";
  }
  else {
    EAssocStepType _type = (type == DO_PRE_STEP ? PRE_STEP : POST_STEP );
    algo.track().journal_out()
      << ", " << assoc_step_poss
      << ", \"" << algo.get_assoc_step_name(step_poss,_type,assoc_step_poss) << "\"";
  }
}

}	// end namespace IterationPack 
