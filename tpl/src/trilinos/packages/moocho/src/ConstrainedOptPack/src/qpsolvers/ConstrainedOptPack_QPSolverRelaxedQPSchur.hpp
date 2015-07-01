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

#ifndef QP_SOLVER_RELAXED_QP_SCHUR_H
#define QP_SOLVER_RELAXED_QP_SCHUR_H

#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "ConstrainedOptPack_QPSchur.hpp"
#include "ConstrainedOptPack_QPInitFixedFreeStd.hpp"
#include "ConstrainedOptPack_MatrixSymHessianRelaxNonSing.hpp"
#include "ConstrainedOptPack_ConstraintsRelaxedStd.hpp"
#include "ConstrainedOptPack_MatrixSymAddDelBunchKaufman.hpp"
#include "AbstractLinAlgPack_VectorMutableDense.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Solves Quadratic Programming (QP) problems using QPSchur.
 *
 * This is the only subclass needed for QPSchur.  All of the specifics of how the
 * initial KKT system is formed is delegated to a strategy object of type
 * \c InitKKTSystem (see below).
 */
class QPSolverRelaxedQPSchur : public QPSolverRelaxed {
public:

  /** \brief Interface for the object that forms the initial KKT system {abstract}.
   *
   * Note that this interface is set up such that the relaxation variable
   * must always be initially fixed (and rightly so to avoid illconditioning).
   */
  class InitKKTSystem {
  public:
    /** \brief . */
    typedef std::vector<size_type> i_x_free_t;
    /** \brief . */
    typedef std::vector<size_type> i_x_fixed_t;
    /** \brief . */
    typedef std::vector<EBounds>   bnd_fixed_t;
    /** \brief . */
    typedef std::vector<size_type> j_f_decomp_t;
    /** \brief . */
    typedef Teuchos::RCP<const MatrixSymOpNonsing>
      Ko_ptr_t;
    /** \brief . */
    virtual ~InitKKTSystem() {}
    /** \brief Initializes the KKT system.
     *
     * Let the following permutation matrices define the selection of the
     * initial KKT system:
     *
     * <tt>Q = [ Q_R, Q_X ]</tt> : Initially fixed <tt>Q_R</tt> and free <tt>Q_X</tt> variables
     *
     * <tt>P = [ P_d, P_u ]</tt> : Decomposed <tt>P_d</tt> and undecomposed <tt>P_u</tt> constraints
     *
     * Given the definitions of <tt>Q</tt> and <tt>P</tt> above, this function will return
     * the initial KKT system:
     \verbatim

     Ko = [ Q_R'*G*Q_R         Q_R'*op(F')*P_d ]
          [ P_d'*op(F)*Q_R     0               ]

     fo = [ -Q_R'*g - Q_R'*G*Q_X*b_X    ]
          [ -P_d'f - P_d'*op(F)*Q_X*b_X ]

     b_X = ??? (see below)
     \endverbatim
     *
     * @param  g    [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  G    [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  dL   [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  dU   [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  F    [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  trans_f
     *              [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  f    [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  d    [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  nu   [in] See <tt>QPSolverRelaxed::solve_qp()</tt>
     * @param  n_R  [out] Number of initially free variables.
     * @param  i_x_free
     *              [out] array (size <tt>n_R</tt> or <tt>0</tt>):
     *              If <tt>i_x_free.size() > 0</tt> then <tt>i_x_free[l-1], l = 1...n_R</tt>
     *              defines the matrix <tt>Q_R</tt> as:<br>
     *              <tt>Q_R(:,l) = e(i_x_free[l-1]), l = 1...n_R</tt><br>
     *              If <tt>i_x_free.size() == 0</tt> then <tt>i_x_free</tt> is implicitly
     *              identity and <tt>Q_R</tt> is defiend as:<br>
     *              <tt>Q_R(:,l) = e(l), l = 1...n_R</tt><br>
     *              The ordering of these indices is significant.
     * @param  i_x_fixed
     *              [out] array (size <tt>n_X</tt>):
     *              <tt>i_x_fixed[l-1], l = 1...n_X</tt> defines the matrix <tt>Q_X</tt> as:<br>
     *              <tt>Q_X(:,l) = e(i_x_fixed[l-1]), l = 1...n_X</tt><br>
     *              The ordering of these indices is significant.
     * @param  bnd_fixed
     *             [out] array (size <tt>n_X</tt>):
     *             <tt>bnd_fixed[l-1], l = 1...n_X</tt> defines the initial active set as:<br>
     *\verbatim
                           / LOWER : b_X(l) = dL(i_x_fixed[l-1])
     bnd_fixed[l-1] = |  UPPER : b_X(l) = dU(i_x_fixed[l-1])
                       \ EQUALITY : b_X(l) = dL(i) = dU(i) (i = i_x_fixed[l-1])
     \endverbatim
     * @param  j_f_decomp
     *             [out] array (size <tt>m</tt>):
     *             <tt>j_f_decomp[p-1], p = 1...m</tt> defines the decomposed equalities included
     *             in <tt>Ko</tt> as:<br>
     *             <tt>P_d(:,p) = e(j_f_decomp[p-1]), p = 1...m</tt><br>
     *             The ordering of these indices is significant and are not necessarily
     *             sorted in assending or decending order.
     * @param  b_X [out] vector (size <tt>n_X</tt>):
     *             Initial varaible bounds (see <tt>bnd_fixed</tt> above).  Note that
     *             the relaxation variable is always one of the initially fixed
     *             variables.
     * @param  Ko  [in/out] Initial KKT matrix (size <tt>(n_R+m) x (n_R+m)</tt>).
     *             On output, Ko will contain a possibly dynamically allocated nonsingular
     *             matrix object that represents Ko.  In input, if <tt>Ko->get() != NULL</tt>,
     *             and no other objects have a reference to this object (based on
     *             <tt>Ko->count()</tt>, and it is of the  proper type, then this matrix may be reused.
     * @param  fo  [out] vector (size <tt>n_R + m</tt>) of the rhs for the initial KKT system.
     */
    virtual void initialize_kkt_system(
      const Vector    &g
      ,const MatrixOp   &G
      ,value_type           etaL
      ,const Vector   *dL
      ,const Vector   *dU
      ,const MatrixOp   *F
      ,BLAS_Cpp::Transp     trans_F
      ,const Vector   *f
      ,const Vector   *d
      ,const Vector   *nu
      ,size_type            *n_R
      ,i_x_free_t           *i_x_free
      ,i_x_fixed_t          *i_x_fixed
      ,bnd_fixed_t          *bnd_fixed
      ,j_f_decomp_t         *j_f_decomp
      ,DVector               *b_X
      ,Ko_ptr_t             *Ko
      ,DVector               *fo
      ) const = 0;

  }; // end class InitKKTSystem

  /** \brief Interface for the object that can reform an initial KKT system
   * dynamically {abstract}.
   *
   * This interface allows the definition of the initial KKT system to
   * be changed on the fly.  This may not be possible for may different
   * QPs so this is an optional interface.  Allowing a redefinition of the
   * initial KKT system may allow QPs with more degrees of freedom and lots
   * of changes to the initial active set to be efficiently solved.
   */
  class ReinitKKTSystem : public InitKKTSystem {
  public:
    // ToDo: Create method reinitailze_kkt_system(...)
  }; // end class ReinitKKTSystem

  /** \brief Strategy object that sets up the initial KKT system.
   */
  STANDARD_COMPOSITION_MEMBERS( InitKKTSystem, init_kkt_sys );

  /** \brief Constraints object.
   */
  STANDARD_COMPOSITION_MEMBERS( QPSchurPack::ConstraintsRelaxedStd, constraints );

  /** \brief Set the maximum number of QP iterations as <tt>max_itr = max_qp_iter_frac * n</tt>.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac );

  /** \brief Set the maximum real run-time in minutes.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_real_runtime );

  /** \brief Policy used to select a violated constraint.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(
    QPSchurPack::ConstraintsRelaxedStd::EInequalityPickPolicy
    ,inequality_pick_policy
    );

  /// Output level
  enum ELocalOutputLevel {
     USE_INPUT_ARG				= -1	// Use the value input to solve_qp(...)
    ,NO_OUTPUT					= 0		//
    ,OUTPUT_BASIC_INFO			= 1		// values sent to QPSchur::solve_qp(...)
    ,OUTPUT_ITER_SUMMARY		= 2		// ...
    ,OUTPUT_ITER_STEPS			= 3
    ,OUTPUT_ACT_SET				= 4
    ,OUTPUT_ITER_QUANTITIES		= 5
  };

  /** \brief Set the output level for QPSchur.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ELocalOutputLevel, print_level );

  /** \brief Set the feasibility tolerance for the bound constriants.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, bounds_tol );

  /** \brief Set the feasibility tolerance for the general inequality constraints.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, inequality_tol );

  /** \brief Set the feasibility tolerance for the general equality constriants.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, equality_tol );

  /** \brief Set a looser feasibility tolerance ( > feas_tol )
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, loose_feas_tol );

  /** \brief Set the tolerence where a scaled Langrange multiplier is considered
   * degenerate.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, dual_infeas_tol );

  /** \brief Set the tolerence for the size of the step in the primal space that is considered
   * to be a near infinite step.  This is used to determine if the KKT
   * system is near singular.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, huge_primal_step );

  /** \brief Set the tolerence for the size of the step in the dual space that is considered
   * to be a near infinite step.  This is used to determine if the constriants
   * are infeasible.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, huge_dual_step );

  /** \brief <<std member comp>> members for the Big M parameter used in the objective.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, bigM );

  /** \brief <<std member comp>> members for the warning tolerance for tests.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol );

  /** \brief <<std member comp>> members for the error tolerance for tests.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol );

  /** \brief Set the minimum number of refinement iterations to perform
   * when using iterative refinement.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, iter_refine_min_iter );
    
  /** \brief Set the maximum number of refinement iterations to perform
   * when using iterative refinement.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, iter_refine_max_iter );

  /** \brief Set the maxinum scaled tolerance the residual of the optimality conditions
   * must be before terminating iterative refinement.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, iter_refine_opt_tol );

  /** \brief Set the maxinum scaled tolerance the residual of the feasibility conditions
   * must be before terminating iterative refinement.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, iter_refine_feas_tol );

  /** \brief Set whether iterative refinement is automatically used once the solution
   * is found.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, iter_refine_at_solution );

  /** \brief Set the relative tolerance for pivots in the schur complement under
   * which a waning will be printed (see MatrixSymAddDelUpdateable) for
   * near singular updates.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, pivot_warning_tol );

  /** \brief Set the relative tolerance for pivots in the schur complement under
   * which a singularity exception will be thrown (see MatrixSymAddDelUpdateable)
   * for singular updates.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, pivot_singular_tol );

  /** \brief Set the relative tolerance for pivots in the schur complement over
   * which a wrong inertia exception will be throw (see MatrixSymAddDelUpdateable)
   * for updates with the wrong inertia.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, pivot_wrong_inertia_tol );

  /** \brief Set whether equality constriants are to be added to the active set
   * initialy to the schur complement or not.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, add_equalities_initially );

  /** \brief . */
  QPSolverRelaxedQPSchur(
    const init_kkt_sys_ptr_t&    init_kkt_sys       = Teuchos::null
    ,const constraints_ptr_t&    constraints        = Teuchos::rcp(new QPSchurPack::ConstraintsRelaxedStd)
    ,value_type                  max_qp_iter_frac   = 10.0
    ,value_type                  max_real_runtime   = 1e+20
    ,QPSchurPack::ConstraintsRelaxedStd::EInequalityPickPolicy
                                 inequality_pick_policy
                                     = QPSchurPack::ConstraintsRelaxedStd::ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY
    ,ELocalOutputLevel           print_level             = USE_INPUT_ARG  // Deduce from input arguments
    ,value_type                  bounds_tol              = -1.0           // use default
    ,value_type                  inequality_tol          = -1.0	          // use default
    ,value_type                  equality_tol            = -1.0	          // use default
    ,value_type                  loose_feas_tol          = -1.0	          // use default
    ,value_type                  dual_infeas_tol         = -1.0	          // use default
    ,value_type                  huge_primal_step        = -1.0	          // use defalut
    ,value_type                  huge_dual_step          = -1.0	          // use default
    ,value_type                  bigM                    = 1e+10
    ,value_type                  warning_tol             = 1e-10
    ,value_type                  error_tol               = 1e-5
    ,size_type                   iter_refine_min_iter    = 1
    ,size_type                   iter_refine_max_iter    = 3
    ,value_type                  iter_refine_opt_tol     = 1e-12
    ,value_type                  iter_refine_feas_tol    = 1e-12
    ,bool                        iter_refine_at_solution = true
    ,value_type                  pivot_warning_tol       = 1e-8
    ,value_type                  pivot_singular_tol      = 1e-11
    ,value_type                  pivot_wrong_inertia_tol = 1e-11
    ,bool                        add_equalities_initially= true
    );

  /** \brief . */
  ~QPSolverRelaxedQPSchur();

  /** @name Overridden from QPSolverRelaxed */
  //@{

  /** \brief . */
  QPSolverStats get_qp_stats() const;
  /** \brief . */
  void release_memory();

  //@}

protected:

  /** @name Overridden from QPSolverRelaxed */
  //@{

  /** \brief . */
  QPSolverStats::ESolutionType imp_solve_qp(
    std::ostream* out, EOutputLevel olevel, ERunTests test_what
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector* dL, const Vector* dU
    ,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
    ,const Vector* eL, const Vector* eU
    ,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
    ,value_type* obj_d
    ,value_type* eta, VectorMutable* d
    ,VectorMutable* nu
    ,VectorMutable* mu, VectorMutable* Ed
    ,VectorMutable* lambda, VectorMutable* Fd
    );

  //@}

private:

  // ////////////////////////////
  // Private data members

  QPSolverStats                    qp_stats_;
  QPSchur                          qp_solver_;
  QPSchurPack::QPInitFixedFreeStd  qp_;
  MatrixSymHessianRelaxNonSing     G_relaxed_;
  VectorMutableDense         bigM_vec_;
  MatrixSymAddDelBunchKaufman      schur_comp_;
  DVector                           g_relaxed_;
  DVector                           b_X_;
  InitKKTSystem::Ko_ptr_t          Ko_;
  DVector                           fo_;

}; // end class QPSolverRelaxedQPSchur

} // end namespace ConstrainedOptPack

#endif // QP_SOLVER_RELAXED_QP_SCHUR_H
