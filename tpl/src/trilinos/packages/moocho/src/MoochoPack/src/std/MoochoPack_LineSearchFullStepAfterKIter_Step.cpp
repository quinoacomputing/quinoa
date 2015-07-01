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

#include "MoochoPack_LineSearchFullStepAfterKIter_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"

bool MoochoPack::LineSearchFullStepAfterKIter_Step::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{
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

  const bool take_full_step = s.k() > full_steps_after_k();

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "\nk = " << s.k() << ( take_full_step ? " > " : " < ")
          << "full_steps_after_k = " << full_steps_after_k() << std::endl;
  }

  if( !take_full_step ) {
    return line_search().do_step(_algo,step_poss,type,assoc_step_poss);
  }
  else {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "\nKeep the full step...\n";
    }
  }

  return true;
}

void MoochoPack::LineSearchFullStepAfterKIter_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out	<< L << "*** Start using full steps after full_steps_after_k iterations.\n"
    << L << "default: full_steps_after_k = very big\n";
  out	<< L << "if k < full_steps_after_k then\n";
  line_search().print_step(algo,step_poss,type,assoc_step_poss,out,L + "    " );
  out	<< L << "end\n";
}

#endif // 0
