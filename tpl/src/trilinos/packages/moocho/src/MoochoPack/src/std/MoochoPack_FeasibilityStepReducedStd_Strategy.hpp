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

#ifndef FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H
#define FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H

#include "MoochoPack_FeasibilityStep_Strategy.hpp"
#include "MoochoPack_QuasiRangeSpaceStep_Strategy.hpp"
#include "MoochoPack_d_bounds_iter_quant.hpp"
#include "IterationPack_CastIQMember.hpp"
#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "ConstrainedOptPack_QPSolverRelaxedTester.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Implements the feasibility step computation for reduced space SQP.
 */
class FeasibilityStepReducedStd_Strategy : public FeasibilityStep_Strategy
{
public:

  /// <<std comp>> members for the qp solver
  STANDARD_COMPOSITION_MEMBERS( QuasiRangeSpaceStep_Strategy, quasi_range_space_step );

  typedef ConstrainedOptPack::QPSolverRelaxedTester
    QPSolverRelaxedTester;

  /// QP solver
  STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxed, qp_solver );

  /// Comparision object compatible with Gc
  STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxedTester, qp_tester );
    
  /** \brief . */
  enum EQPObjective {
    OBJ_MIN_FULL_STEP           ///< min 1/2 * (Y*wy + Z*wz)'*(Y*wy + Z*wz)
    ,OBJ_MIN_NULL_SPACE_STEP    ///< min 1/2 * wz'*wz
    ,OBJ_RSQP                   ///< min qp_grad_k'*wz + 1/2 * wz'*rHL_k*wz
  };

  /** \brief Set what is used for the QP objective.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPObjective, qp_objective );

  /** \brief . */
  enum EQPTesting {
    QP_TEST_DEFAULT     ///< Decide based on olevel input to <tt>compute_feasibility_step(...)</tt>
    ,QP_TEST            ///< Perform the tests
    ,QP_NO_TEST         ///< Don't perform the tests
  };

  /** \brief Set how and if the QP solution is tested.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPTesting, qp_testing );

  /// Construct and initialize
  FeasibilityStepReducedStd_Strategy(
    const quasi_range_space_step_ptr_t   &quasi_range_space_step
    ,const qp_solver_ptr_t               &qp_solver
    ,const qp_tester_ptr_t               &qp_tester
    ,EQPObjective                        qp_objective     = OBJ_MIN_NULL_SPACE_STEP
    ,EQPTesting                          qp_testing       = QP_TEST_DEFAULT
    );

  // ////////////////////////////////////////////
  // Overridden from FeasibilityStep_Strategy

  /** \brief Computes a feasibility step by computing simple quasi-range and null space components.
   *
   * ToDo: Finish documentation!
   *
   */
   bool compute_feasibility_step(
    std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,const Vector& xo, const Vector& c_xo, VectorMutable* w
      );

  /** \brief . */
  void print_step( std::ostream& out, const std::string& leading_str ) const;

private:

  IterationPack::CastIQMember<VectorMutable>  dl_iq_;
  IterationPack::CastIQMember<VectorMutable>  du_iq_;
  int                                                      current_k_;
  Teuchos::RCP<const MatrixOp>            Hess_ptr_;
  VectorSpace::vec_mut_ptr_t                               grad_store_;
  DMatrix                                                Hess_store_;

}; // end class FeasibilityStepReducedStd_Strategy

} // end namespace MoochoPack

#endif // FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H
