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

#ifndef QP_FAILURE_REINIT_REDUCED_HESSIAN_STEP_H
#define QP_FAILURE_REINIT_REDUCED_HESSIAN_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Directs the algorithm to reinitalize the reduced Hessian on the event
 * of a QP failure.
 *
 * If the delegated Step object throws a \c QPFailure exception
 * then this Step object wipes out all reduced Hessian info rHL
 * for the current and previous iterations and then directs the algorithm
 * back to the ReducedHessian step (see the printed step description).
 */
class QPFailureReinitReducedHessian_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** \brief . */
  STANDARD_COMPOSITION_MEMBERS( IterationPack::AlgorithmStep, null_space_step );

  /** \brief . */
  QPFailureReinitReducedHessian_Step( const null_space_step_ptr_t& null_space_step );

  // ////////////////////
  // Overridden

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
  int last_qp_failure_k_;

  // not defined and not to be called
  QPFailureReinitReducedHessian_Step();
  QPFailureReinitReducedHessian_Step(const QPFailureReinitReducedHessian_Step&);
  QPFailureReinitReducedHessian_Step& operator=(const QPFailureReinitReducedHessian_Step&);
  
};	// end class QPFailureReinitReducedHessian_Step

}	// end namespace MoochoPack 

#endif	// QP_FAILURE_REINIT_REDUCED_HESSIAN_STEP_H
