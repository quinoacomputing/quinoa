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

#ifndef NULL_SPACE_STEP_WITH_INEQU_STD_STEP_H
#define NULL_SPACE_STEP_WITH_INEQU_STD_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "IterationPack_CastIQMember.hpp"
#include "MoochoPack_d_bounds_iter_quant.hpp"
#include "MoochoPack_qp_solver_stats.hpp"
#include "MoochoPack_act_set_stats.hpp"
#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "ConstrainedOptPack_QPSolverRelaxedTester.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Solves the reduced QP subproblem with bounds and/or general
 * inequalities.
 *
 * ToDo: Finish documentation.
 */
class TangentialStepWithInequStd_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /// QP solver
  STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxed, qp_solver );

  /// QP solver tester
  STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxedTester, qp_tester );

  /** \brief Set the ratio of the number of inequality constraints in the
   * active-set of the last two calls before a warm start is attempted.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warm_start_frac );

  /** \brief . */
  enum EQPTesting { QP_TEST_DEFAULT, QP_TEST, QP_NO_TEST };

  /** \brief Set how and if the QP solution is tested.
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPTesting, qp_testing );

  /** \brief Determine if a QPFailure exception is thrown if the QP solver
   * returns PRIMAL_FEASIBLE_POINT.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, primal_feasible_point_error );

  /** \brief Determine if a \c QPFailure exception is thrown if the QP solver
   * returns \c DUAl_FEASIBLE_POINT.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dual_feasible_point_error );

  /// Construct and initialize
  TangentialStepWithInequStd_Step(
    const qp_solver_ptr_t       &qp_solver
    ,const qp_tester_ptr_t      &qp_tester
    ,value_type                 warm_start_frac             = 0.8
    ,EQPTesting                 qp_testing                  = QP_TEST_DEFAULT
    ,bool                       primal_feasible_point_error = true
    ,bool                       dual_feasible_point_error   = true
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
  IterationPack::CastIQMember<VectorMutable>  dl_iq_;
  IterationPack::CastIQMember<VectorMutable>  du_iq_;
  qp_solver_stats_iq_member                                qp_solver_stats_;
  act_set_stats_iq_member                                  act_set_stats_;

  // not defined and not to be called
  TangentialStepWithInequStd_Step();

};	// end class TangentialStepWithInequStd_Step

}	// end namespace MoochoPack 

#endif	// NULL_SPACE_STEP_WITH_INEQU_STD_STEP_H
