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

#ifndef CHECK_DECOMPOSITION_FROM_PY_STEP_H
#define CHECK_DECOMPOSITION_FROM_PY_STEP_H

#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Check if the decomposition is going singular and if it is
 * select a new decomposition.
 *
 * This steps checks if the decomposition is going singular by looking
 * to see if the computation for the range space step looks bad.
 *
 * In particular we want to know how \c cond(R) is changing.  We know that:
 \verbatim
 py = -inv(R) * c
 --> ||py|| <= ||inv(R)|| * ||c||
 --> ||py|| / ||c|| <= ||inv(R)||
 --> ( ||py|| / ||c|| ) * ||R|| <= ||inv(R)|| * ||R|| = cond(R)
 \endverbatim
 *
 * If we assume <tt>||R|| > 1</tt> we know that <tt>cond(R) > beta = ||py||/||c||</tt>
 * so if \c beta is very large then cond(R) is even larger.  Since the best decomposition
 * we can find may be fairly illconditioned we may not want to use and absolute measure
 * of \c beta to determine if the decomposition is illconditioned.  Instead we will look
 * at the change in \c beta between iterations and if \c beta gets very much
 * larger than its minimum value (i.e. <tt>beta / beta_min > max_cond_change</tt> 
 * ( default = 10000 ) ) then we will select a new decomposition.  Also if
 * <tt>beta > max_cond</tt> ( default = <tt>0.01 / mach_eps</tt> ) then we know we will be
 * computing inacurate solutions so we must select a new decomposition.
 */
class CheckDecompositionFromPy_Step
  : public IterationPack::AlgorithmStep // doxygen needs full name
{
public:

  /** \brief . */
  STANDARD_COMPOSITION_MEMBERS( NewDecompositionSelection_Strategy, new_decomp_strategy );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_decomposition_cond_change_frac );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_cond );

  /** \brief . */
  CheckDecompositionFromPy_Step(
    const new_decomp_strategy_ptr_t   &new_decomp_strategy
    ,value_type                       max_decomposition_cond_change_frac = 100.0
    );

  /// Call the reset initialization of all defaults.
  void reset();

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

  value_type  beta_min_;

  // Not defined and not to be called
  CheckDecompositionFromPy_Step();

}; // end class CheckDecompositionFromPy_Step

} // end namespace MoochoPack 

#endif // CHECK_DECOMPOSITION_FROM_PY_STEP_H
