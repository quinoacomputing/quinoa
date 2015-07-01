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

#ifndef LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H
#define LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H

#include <limits>

#include "rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "MiStandardAggregationMacros.h"

namespace MoochoPack {

/** \brief Changes from a line search step to just taking full steps after
  * full_steps_after_k iterations.
  */
class LineSearchFullStepAfterKIter_Step : public LineSearch_Step {
public:

  /// <<std comp>> members for the line search step
  STANDARD_COMPOSITION_MEMBERS(LineSearch_Step,line_search);

  /** \brief . */
  LineSearchFullStepAfterKIter_Step(
        const line_search_ptr_t&	line_search			= 0
      , int						full_steps_after_k
                      = std::numeric_limits<int>::max()	)
    : line_search_(line_search)
      , full_steps_after_k_(full_steps_after_k)
  {}

  /// 
  void full_steps_after_k( int full_steps_after_k )
  {	full_steps_after_k_ = full_steps_after_k; }
  /** \brief . */
  value_type full_steps_after_k() const
  {	return full_steps_after_k_; }

  // ////////////////////
  // Overridden

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
  int		full_steps_after_k_;
  
};	// end class LineSearchFullStepAfterKIter_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H
