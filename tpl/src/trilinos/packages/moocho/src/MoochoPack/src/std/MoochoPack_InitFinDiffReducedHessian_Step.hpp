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

#ifndef INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_H
#define INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_H

#include "IterationPack_AlgorithmStep.hpp"
#include "MoochoPack_quasi_newton_stats.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Initializes the reduced hessian using a single finite difference
 * along the null space of the constraints.
 *
 * A single finite difference correction is computed along:\\
 *
 * x_fd = x_k + u * Z * e
 *
 * The step length is set to u = step_scale / ||Z*e||inf.  The
 * step length is cut back if the point x_fd is outside the
 * relaxed variable bounds.
 *
 * The finite difference is then computed as:
 *
 * rGf_fd = ( Z_k' * g(x_k + u * Z*e) - rGf_k ) / u
 *
 * The diagonal terms of the reduced hessian are then set
 * as:
 *
 * diag(i) = max( ||rGf_fd||inf , smallest_ele )  if initialization_method == SCALE_IDENTITY\\
 * diag(i) = max( rGf_fd(i)     , smallest_ele )  if initialization_method == SCALE_DIAGONAL\\
 * diag(i) = max( abs(rGf_fd(i)), smallest_ele )  if initialization_method == SCALE_DIAGONAL_ABS\\
 *
 * Where:
 *
 * smallest_ele = max( ||rGf_fd||inf / max_cond , min_diag )
 *
 * Since the matrix is diagonal the diagonal is equal to the eigenvalues of
 * the matrix.  Therefore you can show that the condition number measured in
 * any norm is max(diag(i))/min(diag(i)).  therefore we just need
 * to limit the smallest diagonal as diag(i) > max(diag(i)) / max_cond.
 */
class InitFinDiffReducedHessian_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** @name Initializers/constructors */
  //@{

  /** \brief . */
  enum EInitializationMethod { SCALE_IDENTITY, SCALE_DIAGONAL, SCALE_DIAGONAL_ABS };

  /** \brief . */
  InitFinDiffReducedHessian_Step(
    EInitializationMethod   initialization_method  = SCALE_IDENTITY
    ,value_type             max_cond               = 1e+1
    ,value_type             min_diag               = 1e-8
    ,value_type             step_scale             = 1e-1
    );

  /// The initialization method for setting the diagonal
  STANDARD_MEMBER_COMPOSITION_MEMBERS(EInitializationMethod,initialization_method);

  /// Maximum condition (l2 norm) for the intial matrix = (max(diag(i))/min(diag(i)).
  STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,max_cond);

  /// The absolute minimum value of a diagonal element
  STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,min_diag);

  /// The scaling of the step length u = step_scale / ||Z*e||inf
  STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,step_scale);

  //@}

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
  quasi_newton_stats_iq_member	quasi_newton_stats_;

};	// end class ReducedHessianBFGS_Step

}	// end namespace MoochoPack 

#endif	// INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_H
