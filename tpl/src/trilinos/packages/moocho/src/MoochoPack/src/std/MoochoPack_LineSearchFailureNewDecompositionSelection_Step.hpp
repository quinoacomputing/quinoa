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

#ifndef LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H
#define LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H

#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Directs the selection of a new decomposition if the line search fails.
  *
  * If the delegated line search Step object throws a \c LineSearchFailure
  * exception, then this object directs the selection of a new
  * decomposition.  If the very next iteration also results in a linesearch
  * failure then we must quit.
  */
class LineSearchFailureNewDecompositionSelection_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /// <<std comp>> members for LineSearch object.
  STANDARD_COMPOSITION_MEMBERS( IterationPack::AlgorithmStep, line_search_step );

  /// <<std comp>> members for Decomposition Select Strategy object.
  STANDARD_COMPOSITION_MEMBERS( NewDecompositionSelection_Strategy, new_decomp_strategy );

  /** \brief . */
  LineSearchFailureNewDecompositionSelection_Step(
    const line_search_step_ptr_t        &line_search_step
    ,const new_decomp_strategy_ptr_t    &new_decomp_strategy
    );

  /** @name Overridden from AlgorithmStep */
  //@{
  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
  //@}

private:

  int last_ls_failure_k_;

  // not defined and not to be called
  LineSearchFailureNewDecompositionSelection_Step();
  LineSearchFailureNewDecompositionSelection_Step(
    const LineSearchFailureNewDecompositionSelection_Step&);
  LineSearchFailureNewDecompositionSelection_Step& operator=(
    const LineSearchFailureNewDecompositionSelection_Step&);

};	// end class LineSearchFailureNewDecompositionSelection_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H
