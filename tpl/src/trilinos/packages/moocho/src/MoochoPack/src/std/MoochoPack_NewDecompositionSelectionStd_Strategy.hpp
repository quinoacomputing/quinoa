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

#ifndef NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H
#define NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H

#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"
#include "MoochoPack_DecompositionSystemHandlerSelectNew_Strategy.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Just force the decomposition system object to select a new
 * decomposition and let everyone else fend for themselves.
 */
class NewDecompositionSelectionStd_Strategy
  : public NewDecompositionSelection_Strategy
{
public:

  /// «std comp» members for range/null decomposition handler
  STANDARD_COMPOSITION_MEMBERS( DecompositionSystemHandlerSelectNew_Strategy, decomp_sys_handler );

  /** \brief . */
  NewDecompositionSelectionStd_Strategy(
    const decomp_sys_handler_ptr_t   &decomp_sys_handler
    );

  /** @name Overridden from NewDecompositionSelection_Strategy */
  //@{

  bool new_decomposition(
    NLPAlgo& algo, Algorithm::poss_type step_poss
    ,IterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
    );
  /** \brief . */
  void print_new_decomposition(
    const NLPAlgo& algo, Algorithm::poss_type step_poss
    ,IterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
    ,std::ostream& out, const std::string& leading_str
    ) const;

  //@}

private:

  // Not defined and not to be called
  NewDecompositionSelectionStd_Strategy();

};	// end class NewDecompositionSelectionStd_Strategy

}	// end namespace MoochoPack 

#endif	// NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H
