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

// #define RELEASE_TRACE

#include <iostream> // for debugging Release version.
#include <typeinfo>

#include "MoochoPack_NLPAlgo.hpp"

namespace MoochoPack {

NLPAlgo::NLPAlgo()
  : algo_cntr_(NULL), nlp_(NULL), first_step_poss_(1)
{}

// Overridden form rSQPAlgoInteface

const NLPAlgoState& NLPAlgo::retrieve_state() const
{
  return dynamic_cast<const NLPAlgoState&>(state());
}

NLPSolverClientInterface::EFindMinReturn
NLPAlgo::dispatch() {
  switch( do_algorithm(first_step_poss_) ) {
    case IterationPack::TERMINATE_TRUE:
      return NLPSolverClientInterface::SOLUTION_FOUND;
    case IterationPack::TERMINATE_FALSE:
      return NLPSolverClientInterface::ALGORITHMIC_ERROR;
    case IterationPack::MAX_ITER_EXCEEDED:
      return NLPSolverClientInterface::MAX_ITER_EXCEEDED;
    case IterationPack::MAX_RUN_TIME_EXCEEDED:
      return NLPSolverClientInterface::MAX_RUN_TIME_EXCEEDED;
    case IterationPack::INTERRUPTED_TERMINATE_TRUE:
      return NLPSolverClientInterface::SOLUTION_FOUND;
    case IterationPack::INTERRUPTED_TERMINATE_FALSE:
      return NLPSolverClientInterface::ALGORITHMIC_ERROR;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return NLPSolverClientInterface::SOLUTION_FOUND;	// will never be called.
}

void NLPAlgo::interface_print_algorithm(std::ostream& out) const {
  print_steps(out);
  print_algorithm(out);
}

void NLPAlgo::interface_set_algo_timing( bool algo_timing ) {
  set_algo_timing(algo_timing);
}

bool NLPAlgo::interface_algo_timing() const {
  return algo_timing();
}

void NLPAlgo::interface_print_algorithm_times( std::ostream& out ) const {
  print_algorithm_times(out);
}

// Overridden from Algorithm.

void NLPAlgo::print_algorithm(std::ostream& out) const {
  out
    << "\n*** NLP ***\n"
    << typeName(*get_nlp()) << "\n";

  Algorithm::print_algorithm(out);
}

}	// end namespace MoochoPack
