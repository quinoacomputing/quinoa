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

#ifndef QPSCHUR_H
#define QPSCHUR_H

#include <ostream>
#include <map>
#include <vector>

#include "ConstrainedOptPack_Types.hpp"
#include "ConstrainedOptPack_MatrixSymAddDelUpdateableWithOpNonsingular.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_MatrixSymAddDelUpdateable.hpp"
#include "AbstractLinAlgPack_MatrixOpSerial.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "StopWatchPack_stopwatch.hpp"

namespace ConstrainedOptPack {

namespace QPSchurPack {

/// Utility class for a ranged check vector
template < class T >
class vector_one_based_checked : public std::vector<T>
{
  typedef vector_one_based_checked this_t;
public:
  /// one based indexing
  T& operator()( typename this_t::size_type i )
  {
#ifdef LINALGPACK_CHECK_RANGE
      return this->at(i-1);
#else
      return this->operator[](i-1);
#endif
  }
  /// one based indexing
  T operator()( typename this_t::size_type i ) const
  {
#ifdef LINALGPACK_CHECK_RANGE
      return this->at(i-1);
#else
      return this->operator[](i-1);
#endif
  }
}; // end class vector_one_based

class Constraints;
class QP;

/** \brief Represents the QP to be solved by QPSchur {abstract}.
 *
 * In order to solve a QP, clients must define subclasses
 * for this interface and the \c Constraints interface
 * defined later.  This is where the specialized properties
 * of the QP are exploited.  For certain types of QPs, standard
 * implementation classes can be defined.
 *
 * Here the QP is of the form:
 \verbatim

  (1.a)	min		g'*x + 1/2*x'*G*x
  (1.b)	s.t.	A'*x = c
  (1.c)			cl_bar <= A_bar'*x <= cu_bar

  where:
    x <: R^n
    g <: R^n
    G <: R^(n x n)
    A <: R^(n x m)
    c <: R^m
    A_bar <: R^(n x m_bar)
    cl_bar, cu_bar <: R^m_bar
 \endverbatim
 * Above, <tt>cl_bar <= A_bar'*x <= cu_bar</tt> may represent variable bounds, general
 * inequalities and equality constraints and these constraints are represented
 * by the class \c Constraints.
 *
 * When solving the above QP using the schur complement QP solver we start out with
 * a KKT system with a subset of variables initially fixed at a bound:
 \verbatim

  [ G_RR    G_RX     A_F     0  ] [ x_R    ]   [ - g_R ]
  [ G_RX'   G_XX     A_X     I  ] [ x_X    ]   [ - g_X ]
  [ A_R'    A_X'      0      0  ]*[ lambda ] = [   c   ]
  [           I       0      0  ] [ mu_X   ]   [   b_X ]
 \endverbatim
 * We can simplify the above system by solving for the initially
 * fixed variables and removing them from the initial KKT system to give:
 \verbatim

  x_X = b_X

  [ G_RR     A_R  ]   [ x_R    ]   [ - g_R - G_RX*b_X ]
  [ A_R'      0   ] * [ lambda ] = [    c - A_X'*b_X  ]
  \_______________/   \________/   \__________________/
          Ko              vo                fo

 mu_X = - g_X - G_RX'*x_R - G_X*b_X - A_X*lambda

  where:
      n_X = n - n_R
    x_R  = Q_R'*x        <: R^n_R
    g_R  = Q_R'*g        <: R^n_R
    G_RR = Q_R'*G*Q_R    <: R^(n_R x n_R)
    G_RX = Q_R'*G*Q_X    <: R^(n_R x n_X)
    G_XX = Q_X'*G*Q_X    <: R^(n_X x n_X)
    A_R  = Q_R'*A        <: R^(n_R x m)
    A_X  = Q_X'*A        <: R^(n_X x m)
    Q_R                  <: R^(n x n_R)
    Q_X                  <: R^(n x n_X)
 \endverbatim
 * This class is an interface for encapsulating the QP.  Operations are available
 * for accessing <tt>g</tt>, <tt>G</tt>, <tt>A</tt>, <tt>Ko</tt>, <tt>vo</tt>, and
 * <tt>fo</tt> as well as the mapping
 * matrices <tt>Q_R</tt> and <tt>Q_X</tt> (both ordered by row). Also, operations are available
 * for accessing data structures that describe the set of initially fixed and free
 * variables.  See the \c Constraints interface for how to access the constraints
 * in (1.c) and the matrix <tt>A_bar</tt>.
 */
class QP {
public:

  // /////////////////
  // Public Types

  /** \brief . */
  typedef vector_one_based_checked<EBounds>		x_init_t;
  /** \brief . */
  typedef vector_one_based_checked<size_type>		l_x_X_map_t;
  /** \brief . */
  typedef vector_one_based_checked<size_type>		i_x_X_map_t;

  /** \brief . */
  typedef QPSchurPack::Constraints				Constraints;

  // /////////////////
  // Public Interface

  /** \brief . */
  virtual ~QP()
  {}

  // ///////////////////////////////////////
  // Initial active set independent members 

  /** \brief . */
  virtual size_type n() const = 0;
  /** \brief . */
  virtual size_type m() const = 0;
  /** \brief . */
  virtual const DVectorSlice g() const = 0;
  /** \brief . */
  virtual const MatrixSymOp& G() const = 0;
  /// If m == 0 then don't call this, it may throw an exception or worse.
  virtual const MatrixOp& A() const = 0;

  // /////////////////////////////////////
  // Initial active set specific members

  /** \brief . */
  virtual size_type n_R() const = 0;

  /** \brief Return the status of a variable initially.
   *
   * For 1 <= i <= n:
   \verbatim
               / FREE      : x(i) is initially free
               | LOWER     : x(i) is initially fixed at xl(i)
   x_init(i) = | UPPER     : x(i) is initially fixed at xu(i)
               \ EQUALITY  : x(i) fixed at xl(i) = xu(i) always
   \endverbatim
   */
  virtual const x_init_t& x_init() const = 0;

  /** \brief Map from full x(i) to initially fixed x_X(l).
   *
   * For 1 <= i <= n:
   * 
   \verbatim
                   / l : x(i) = x_X(l) = b_X(l) initially (1 <= l <= n_X)
   l_x_X_map(i) =  |
                   \ 0 : otherwise
   \endverbatim
   *
   */
  virtual const l_x_X_map_t& l_x_X_map() const = 0;

  /** \brief Map from initially fixed x_X(l) to full x(i).
   *
   * For 1 <= l <= n_X:
   * 
   \verbatim
   i_x_X_map(l) = i : x(i) = x_X(l) = b_X(l) initially (1 <= i <= n)
   \endverbatim
   *
   */
  virtual const i_x_X_map_t& i_x_X_map() const = 0;

  /** \brief . */
  /* The bounds of the initially fixed variables.
   *
   * For 1 <= l <= n_X:
   *
   \verbatim
             / xl(i_x_X_map(l))                     : if x_init(i_x_X_map(l)) == LOWER
   b_X(l) =  | xu(i_x_X_map(l))                     : if x_init(i_x_X_map(l)) == UPPER
             \ xl(i_x_X_map(l)) = xu(i_x_X_map(l))  : if x_init(i_x_X_map(l)) == EQUALITY
   \endverbatim
   *
   */
  virtual const DVectorSlice b_X() const = 0;

  /// (Q_R().ordered_by() == BY_ROW)
  virtual const GenPermMatrixSlice& Q_R() const = 0;

  /// (Q_X().ordered_by() == BY_ROW)
  virtual const GenPermMatrixSlice& Q_X() const = 0;

  /** \brief . */
  virtual const MatrixSymOpNonsing& Ko() const = 0;

  /** \brief . */
  virtual const DVectorSlice fo() const = 0;

  // //////////////////////////////////////////////////////////
  // Additional constaints for cl_bar <= A_bar'*x <= cu_bar

  /** \brief . */
  virtual Constraints& constraints() = 0;

  /** \brief . */
  virtual const Constraints& constraints() const = 0;

  /** \brief Dump the definition of the QP to a stream.
   *
   * This function is only to be used for debugging small problems.
   */
  virtual void dump_qp( std::ostream& out );

};	// end class QP

/** \brief Represents the extra constraints in the QP to be satisfied
 * by the schur complement QP solver QPSchur {abstract}.
 *
 * This class is only ment to be used in conjunction with the class \c QP
 * and \c QPSchur.  Its interface is designed to be minimal with respect to
 * the needs of the <tt>QPSchur</tt> solver.  However, this interface may be useful
 * for any primal-dual QP solver.
 *
 * These constraints are:
 \verbatim

  (1.c)	cl_bar <= A_bar'*x <= cu_bar

  where:
    A_bar <: R^(n x m_bar)
    cl_bar, cu_bar <: R^m_bar
 \endverbatim
 *
 * These constraints are also partitioned as:
 \verbatim

  s.t.
    [     xl     ]    [   I       ]       [     xu     ]
    [  cl_breve  ] <= [  A_breve' ]*x  <= [  cu_breve  ]

  where:
    I <: R^(n x n)
    xl, xu <: R^n, are the variable bounds for variables that have bounds (sparse)
    A_breve <: R^(n x m_breve), is the Jacobian for the general constraints
    cl_breve, cu_breve <: R^m_breve, are bounds for general constraints
 \endverbatim
 *
 * Here <tt>m_bar = n + m_breve</tt>
 *
 * Above, some of the bounds in <tt>xl</tt>, <tt>xu</tt>, <tt>cl_breve</tt>,
 * and <tt>cu_breve</tt> may be <tt>-inf</tt> or <tt>+inf</tt> and will
 * therefore never be violated and never be added to the active set.
 * Also, some of the lower and upper bounds may be equal which turns those
 * inequality constraints into equality constraints (or fixed variables).
 */
class Constraints {
public:

  /** \brief . */
  enum EPickPolicy { ANY_VIOLATED, MOST_VIOLATED };

  /** \brief . */
  virtual ~Constraints() {}

  /** \brief . */
  virtual size_type n() const = 0;
  
  /** \brief . */
  virtual size_type m_breve() const = 0;

  /** \brief . */
  virtual const MatrixOp& A_bar() const = 0;
  
  /// Set the policy used to pick a violated constraint.
  virtual void pick_violated_policy( EPickPolicy pick_policy ) = 0;
  /** \brief . */
  virtual EPickPolicy pick_violated_policy() const = 0;

  /** \brief Pick a violated constraint.
   *
   * @param	x			 [in] Trial point to pick a violated constraint at.
   * @param	j_viol		 [out] Indice of violated constraint.  j_viol = 0 if
   *							 no constraint is violated by more that some tolerance.
   * @param	constr_val	 [out] The value if the violated constraint a_bar(j)'*x.
   * @param	viol_bnd_val [out] The value if the violated bound.
   * @param	norm_2_constr[out] The 2 norm of the violated constraint ||a_bar(j)||2
   * @param	bnd			 [out] Classification of the bound being violated.
   * @param	can_ignore	 [out] True if the constraint can be ignored if it is linearly
   *							 dependent.
   */
  virtual void pick_violated(
     const DVectorSlice& x, size_type* j_viol, value_type* constr_val
    ,value_type* viol_bnd_val, value_type* norm_2_constr, EBounds* bnd, bool* can_ignore
    ) const = 0;

  /** \brief Inform to ignore the jth constraint the next time pick_violated(...) is called.
   */
  virtual void ignore( size_type j ) = 0;

  /** \brief Return the bound for a constraint.
   *
   * @param	j	[in] Indice of the constraint of the bound to obtain.
   * @param	bnd	[in] Which bound to obtain (UPPER or LOWER).
   * @return
   *		xl(j) [ 0 < j < n, bnd == LOWER ]<br>
   *		xu(j) [ 0 < j < n, bnd == UPPER ]<br>
   *		cl_breve(j - n) [ n + 1 < j < n + m_breve, bnd == LOWER ]<br>
   *		cu_breve(j - n) [ n + 1 < j < n + m_breve, bnd == UPPER ]<br>
   */
  virtual value_type get_bnd( size_type j, EBounds bnd ) const = 0;

};	// end class Constraints

}	// end namespace QPSchurPack 

/** \brief Solves a Quadratic Program with a dual QP method using a schur complement
 * factorization.
 *
 * See the paper "QPSchur: A Primal-Dual Active-Set Quadratic Programming
 * Algorithm Using a Schur Complement Factorization Method" for a description
 * of what this class does.
 */
class QPSchur {
public:

  /** @name Public Types */
  //@{

  /** \brief . */
  typedef QPSchurPack::QP               QP;
  /** \brief . */
  typedef MatrixSymAddDelUpdateable     MSADU;
  /// Thrown if a test failed
  class TestFailed : public std::logic_error
  {public: TestFailed(const std::string& what_arg) : std::logic_error(what_arg) {}};
  /// Thrown if constraints are inconsistant (no feasible region)
  class InconsistantConstraintsException : public std::logic_error
  {public: InconsistantConstraintsException(const std::string& what_arg) : std::logic_error(what_arg) {}};
  /// Thrown if there is some numerical instability
  class NumericalInstabilityException : public std::runtime_error
  {public: NumericalInstabilityException(const std::string& what_arg) : std::runtime_error(what_arg) {}};
  /// Thrown if during the course of the primal-dual iteration a non-dual feasible point if found.
  class DualInfeasibleException : public NumericalInstabilityException
  {public: DualInfeasibleException(const std::string& what_arg)
    : NumericalInstabilityException(what_arg) {}};
  /// Enumeration for if to run internal tests or not.
  enum ERunTests { RUN_TESTS, NO_TESTS };
  /// solve_qp return values
  enum ESolveReturn {
    OPTIMAL_SOLUTION
    ,MAX_ITER_EXCEEDED
    ,MAX_RUNTIME_EXEEDED_FAIL
    ,MAX_RUNTIME_EXEEDED_DUAL_FEAS
    ,MAX_ALLOWED_STORAGE_EXCEEDED
    ,INFEASIBLE_CONSTRAINTS
    ,NONCONVEX_QP
    ,DUAL_INFEASIBILITY
    ,SUBOPTIMAL_POINT
  };
  /// Output level
  enum EOutputLevel {
     NO_OUTPUT					= 0
    ,OUTPUT_BASIC_INFO			= 1
    ,OUTPUT_ITER_SUMMARY		= 2
    ,OUTPUT_ITER_STEPS			= 3
    ,OUTPUT_ACT_SET				= 4
    ,OUTPUT_ITER_QUANTITIES		= 5
  };
  /// Value for near degenerate lagrange multipliers
  static value_type DEGENERATE_MULT;

  //@}

  /** @name Public Member functions */
  //@{

  /// Schur complement matrix object S_hat
  STANDARD_COMPOSITION_MEMBERS( MatrixSymAddDelUpdateableWithOpNonsingular, schur_comp );

  /** \brief Set the maximum number of primal-dual QP iterations to take.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, max_iter );

  /** \brief Set the maximum wall clock runtime (in minutes).
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_real_runtime );

  /** \brief Set the feasibility tolerance for the constriants.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_tol );

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
    
  /** \brief Set whether a singular initial schur complement will attempted to be
   * salvaged by adding as many nonsingular rows/cols as possible.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, salvage_init_schur_comp );

  /** \brief Set the tolerances to use when updating the schur complement.
   */
  void pivot_tols( MSADU::PivotTolerances pivot_tols );
  /** \brief . */
  MSADU::PivotTolerances pivot_tols() const;

  /** \brief . */
  virtual ~QPSchur() {}

  /** \brief . */
  QPSchur(
    const schur_comp_ptr_t&   schur_comp           = Teuchos::null
    ,size_type                max_iter             = 100
    ,value_type               max_real_runtime     = 1e+20
    ,value_type               feas_tol             = 1e-8
    ,value_type               loose_feas_tol       = 1e-6
    ,value_type               dual_infeas_tol      = 1e-12
    ,value_type               huge_primal_step     = 1e+20
    ,value_type               huge_dual_step       = 1e+20
    ,value_type               warning_tol          = 1e-10
    ,value_type               error_tol            = 1e-5
    ,size_type                iter_refine_min_iter = 1
    ,size_type                iter_refine_max_iter = 3
    ,value_type               iter_refine_opt_tol  = 1e-12
    ,value_type               iter_refine_feas_tol = 1e-12
    ,bool                     iter_refine_at_solution = true
    ,bool                     salvage_init_schur_comp = true
    ,MSADU::PivotTolerances   pivot_tols = MSADU::PivotTolerances( 1e-8,1e-11,1e-11 )
    );

  /** \brief Solve a QP.
   *
   * If the initial schur complement turns out to have the wrong inertia then
   * the QP is nonconvex, and the exception \c WrongInteriaUpdateExecption will be thrown.
   * Otherwise, unless some other strange exception is thrown, this function
   * will return normally (see return).
   *
   * @param	qp	[in] The abstraction for the QP being solved
   * @param	num_act_change
   * 			[in] The number of changes to the
   *					active set before the primal-dual QP algorithm
   *					starts.
   * @param	ij_act_change
   * 			[in] Array (size num_act_change): specifying
   *					how to initialize the active set.  If i = -ij_act_change(s)
   *					> 0 then the initially fixed variable x(i) is to be
   *					freed.  If j = ij_act_change(s) > 0 then the constraint
   *					a_bar(j)'*x is to be added to the active set to the
   *					bound bnd(s).  The order of these changes can significantly
   *					effect the performance of the algorithm if these changes
   *					are not part of the optimal active set.  Put changes that
   *					you are sure about earlier in the list and those that you
   *					are not a sure about later.
   * @param	bnd [in] Array (size num_act_change):  bnd(s) gives which bound to
   *					make active.  If ij_act_change(s) < 0 then this is ignored.
   * @param	out [out] output stream.  Iteration information is printed according
   *					to output_level.  If <tt>output_level == NO_OUTPUT</tt> then <tt>out</tt> may
   *					be <tt>NULL</tt>.  If <tt>out==NULL</tt>, then output_level is forced to <tt>NO_OUTPUT</tt>
   * @param	output_level
   * 			[in] Specifies the level of output (see \c EOutputLevel).
   *	@param	test_what
   *				[in] Determines if internal validation tests are performed.
   *					The optimality conditions for the QP are not checked
   *					internally, since this is something that client can
   *					(and should) do independently.
   *					RUN_TESTS : As many validation/consistency tests
   *						are performed internally as possible.  If a test
   *						fails then a TestFailed execption will be thrown.
   *					NO_TEST : No tests are performed internally.  This is
   *						to allow the fastest possible execution.
   * @param	x	[out] vector (size qp.n()): Solution or current iteration value
   * @param	mu 	[out] sparse vector (size qp.n()): Optimal lagrange multipliers for
   * 				bound constraints.  On output mu->is_sorted() == true.
   * @param	lambda
   * 			[out] vector (size q.m()): Optimal lagrange multipliers for
   *					equality constraints.
   * @param	lambda_breve
   * 			[out] sparse vector (size qp.constraints().m_breve()) for the active
   *					constraints in A_breve.  On output lambda_breve->is_sorted() == true.
   *	@param	iter [out] The number of warm start drops and primal-dual iterations.
   *	@param	num_adds [out] The number of updates to the active set where a constraint
   *					was added.  These do not include initially fixed variables.
   *	@param	num_drops [out] The number of updates to the active set where a constraint
   *					was dropped.  These include constraints dropped during a warm start
   *					as well as during the primal-dual QP iterations.
   *
   * @return  Returns the type of point that x, mu , lambda and lambda_breve represents.
   */
  virtual ESolveReturn solve_qp(
    QP& qp
    ,size_type num_act_change, const int ij_act_change[], const EBounds bnds[]
    ,std::ostream *out, EOutputLevel output_level, ERunTests test_what
    ,DVectorSlice* x, SpVector* mu, DVectorSlice* lambda, SpVector* lambda_breve
    ,size_type* iter, size_type* num_adds, size_type* num_drops
    );

  //@}

  /** \brief Represents the matrix U_hat. 
   *
   * This matrix is only ment to be an aggregate of an <tt>ActiveSet</tt>
   * object and is only managed by the <tt>ActiveSet</tt> object.  It is made
   * public so that clients can developed specialized implementations
   * if needed.
   */
  class U_hat_t : public MatrixOpSerial {
  public:
    /// Construct uninitialized
    U_hat_t();
    /// Initialize.
    void initialize( 
       const MatrixSymOp		*G
      ,const MatrixOp			*A
      ,const MatrixOp			*A_bar
      ,const GenPermMatrixSlice	*Q_R
      ,const GenPermMatrixSlice	*P_XF_hat
      ,const GenPermMatrixSlice	*P_plus_hat
      );
    /** \brief . */
    const MatrixSymOp& G() const
    {	return *G_;	}
    /** \brief . */
    const MatrixOp* A() const
    {	return A_;	}
    /** \brief . */
    const MatrixOp& A_bar() const
    {	return *A_bar_;	}
    /** \brief . */
    const GenPermMatrixSlice& Q_R() const
    {	return *Q_R_; }
    /** \brief . */
    const GenPermMatrixSlice& P_XF_hat() const
    {	return *P_XF_hat_;	}
    /** \brief . */
    const GenPermMatrixSlice& P_plus_hat() const
    {	return *P_plus_hat_;	}
    
    /** @name Overridden from MatrixBase */
    //@{{

    /** \brief . */
    size_type rows() const;
    /** \brief . */
    size_type cols() const;

    //@}

    /** @name Overridden from MatrixOpSerial */
    //@{

    /** \brief . */
    void Vp_StMtV(
      DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
      ,const DVectorSlice& vs_rhs2, value_type beta
      ) const;
    /** \brief . */
    void Vp_StMtV(
      DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
      ,const SpVectorSlice& sv_rhs2, value_type beta
      ) const;

    //@}
    
  private:
    const MatrixSymOp	    *G_;
    const MatrixOp	        *A_;
    const MatrixOp	        *A_bar_;
    const GenPermMatrixSlice	*Q_R_;
    const GenPermMatrixSlice	*P_XF_hat_;
    const GenPermMatrixSlice	*P_plus_hat_;

  };	// end class U_hat_t

  /** \brief Represents and manages the active set for the QPSchur algorithm.
   *
   * This is a concrete type that encapsulates the maintaince of the active set and
   * abstracts quantities associated with it.
   *
   * At each iteration the dual active-set QP algorithm must solve the system:
   \verbatim

    [ Ko       U_hat ]   [    v  ]   [  fo   ]
    [ U_hat'   V_hat ] * [ z_hat ] = [ d_hat ]
   \endverbatim
   *
   * Above, \c U_hat contains the updates to the KKT system for adding constraints
   * to the active set and freeing variables that where initially fixed
   * and therefore left out of <tt>Ko</tt>.
   * 
   * This object maintains references to objects that represent the current
   * augmented KKT system:
   \verbatim

   MatrixOp                      : U_hat        <: R^((n_R+m) x q_hat)
   MatrixSymOp                   : V_hat        <: R^(q_hat x q_hat)
   MatrixSymOpNonsing        : S_hat        <: R^(q_hat x q_hat)
   GenPermMatrixSlice                : P_XF_hat     <: R^(n x q_hat)             (q_F_hat nonzeros)
   GenPermMatrixSlice                : P_FC_hat     <: R^(q_hat x q_hat)         (q_C_hat nonzeros)
   GenPermMatrixSlice                : P_plus_hat   <: R^((n+m_breve) x q_hat)   (q_plus_hat nonzeros)
   GenPermMatrixSlice                : Q_XD_hat     <: R^(n x q_D_hat)           (q_D_hat nonzeros)
   DVector                            : d_hat        <: R^(q_hat)
   DVector                            : z_hat        <: R^(q_hat)
   \endverbatim
   */
  class ActiveSet {
  public:

    // /////////////////////
    // Public types

    /** \brief . */
    typedef QPSchurPack::QP            QP;
    /** \brief . */
    typedef MatrixSymAddDelUpdateable  MSADU;

    // /////////////////////
    // Public interface

    /** \brief «std comp» members for schur complement matrix S_hat.
     *
     * Warning: Resetting schur_comp will cause a reinitialization to
     * an empty active set.
     */
    STANDARD_COMPOSITION_MEMBERS( MatrixSymAddDelUpdateableWithOpNonsingular, schur_comp );

    /** \brief Set the tolerances to use when updating the schur complement.
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( MSADU::PivotTolerances, pivot_tols );

    /** \brief . */
    ActiveSet(
      const schur_comp_ptr_t   &schur_comp
      ,MSADU::PivotTolerances  pivot_tols = MSADU::PivotTolerances( 1e-6,1e-8,1e-8 )
      );

    /** @name Update the active set. */
    //@{

    /** \brief Initialize with an additional active set.
     *
     * If the initial schur complement is not full rank
     * then an <tt>LDConstraintException</tt> exception will be thrown.
     * The active set will contain all of the constraints it
     * can such that the schur complement is nonsingular.
     */
    void initialize( 
      QP& qp, size_type num_act_change, const int ij_act_change[]
      ,const EBounds bnds[], bool test, bool salvage_init_schur_comp
      ,std::ostream *out, EOutputLevel output_level );

    /** \brief Reinitialize the schur complement factorization for the current active set
     *
     * ToDo: Finish documentation
     */
    void refactorize_schur_comp();

    /** \brief Add a constraint to the active set then refactorize the schur complemnt
     * (if forced).
     *
     * ToDo: Finish documentation
     *
     * If the new KKT system is singular then the exeption
     * \c MatrixSymAddDelUpdateable::SingularUpdateException will be thrown
     * but the old KKT system will be kept intact.
     *
     * If the reduced Hessian for the new KKT system does not have the
     * correct inertia then the exception
     * MatrixSymAddDelUpdateable::WrongInertiaUpdateException
     * will be thrown but the old KKT system will be kept intact.
     *
     * @return Returns true if any output was sent to *out.
     */
    bool add_constraint(
      size_type ja, EBounds bnd_ja, bool update_steps
      ,std::ostream *out, EOutputLevel output_level
      ,bool force_refactorization = true
      ,bool allow_any_cond = false );

    /** \brief Drop a constraint from the active set then refactorize the schur
     * complement (if forced).
     *
     * ToDo: Finish documentation
     *
     * Returns true if any output was sent to *out.
     */
    bool drop_constraint(
      int jd, std::ostream *out, EOutputLevel output_level
      ,bool force_refactorization = true, bool allow_any_cond = false );

    /** \brief Drop a constraint from, then add a constraint to the active set
     * and refactorize the schur complement.
     *
     * ToDo: Finish documentation
     *
     * Returns true if any output was sent to *out.
     */
    bool drop_add_constraints(
      int jd, size_type ja, EBounds bnd_ja, bool update_steps
      ,std::ostream *out, EOutputLevel output_level );

    //@}

    /** @name access the QP */
    //@{

    /** \brief . */
    QP& qp();
    /** \brief . */
    const QP& qp() const;

    //@}

    /** @name Access the active sets quantities. */
    //@{

    /** \brief Return the total size of the schur complement.
     *
     * q_hat = q_plus_hat + q_F_hat + q_C_hat.
     */
    size_type q_hat() const;

    /** \brief Return the number of constraints from A_bar added
     * to the active set.
     */
    size_type q_plus_hat() const;

    /** \brief Return the number of variables that where
     * initially fixed but are currently free or
     * fixed to another bound.
     */
    size_type q_F_hat() const;

    /** \brief Return the number of variables that where
     * initially fixed but are currently
     * fixed to another bound.
     */
    size_type q_C_hat() const;

    /** \brief Return the number of variables that where
     * initially fixed and are still currently
     * fixed to their intial bounds.
     */
    size_type q_D_hat() const;

    /** \brief Returns -i for row & column of S_bar for an initially
     * fixed variable left out of Ko that became free and returns
     * j for the constraint a(j)'*x that was added to the active
     * set.
     *
     * <tt>1 <= s <= q_hat</tt>
     */
    int ij_map( size_type s ) const;

    /** \brief Map from a constraint or initially fixed variable
     * to a row and column in the schur complement S_bar.
     *
     * To determine if an initially fixed variable x(i) is now
     * free call s_map(-i).  If s_map(-i) returns zero then
     * x(i) is still fixed.  Otherwise s_map(-i) returns the
     * row and column in S_bar for this change in the
     * active set.
     *
     * To determine if a constraint a(j)'*x is part of the
     * active set call s_map(j).  If s_map(j) returns zero
     * then a(j)'*x is not part of the active set.
     * Otherwise s_map(j) returns the row and column
     * in S_bar for this change in the active set.
     */
    size_type s_map( int ij ) const;

    /** \brief Returns ||a(j)||2 where j = ij_map(s).
     * 
     * If ij_map(s) < 0, the this function returns zero.
     * 
     * 1 <= s <= q_hat
     */
    value_type constr_norm( size_type s ) const;

    /** \brief Return which bound is active for the active constraint.
     */
    EBounds bnd( size_type s ) const;

    /** \brief Returns the indice of x_X(l) of the initially fixed variables
     * that are still fixed at their original bounds.
     *
     * i <= k <= q_D_hat
     */
    size_type l_fxfx( size_type k ) const;

    /** \brief . */
    const U_hat_t& U_hat() const;
    /** \brief . */
    const MatrixSymOpNonsing& S_hat() const;
    /** \brief . */
    const GenPermMatrixSlice& P_XF_hat() const;
    /** \brief . */
    const GenPermMatrixSlice& P_FC_hat() const;
    /** \brief . */
    const GenPermMatrixSlice& P_plus_hat() const;
    /** \brief . */
    const GenPermMatrixSlice& Q_XD_hat() const;
    /** \brief . */
    const DVectorSlice d_hat() const;
    /** \brief . */
    DVectorSlice z_hat();
    /** \brief . */
    const DVectorSlice z_hat() const;
    /** \brief . */
    DVectorSlice p_z_hat();
    /** \brief . */
    const DVectorSlice p_z_hat() const;
    /** \brief . */
    DVectorSlice mu_D_hat();
    /** \brief . */
    const DVectorSlice mu_D_hat() const;
    /** \brief . */
    DVectorSlice p_mu_D_hat();
    /** \brief . */
    const DVectorSlice p_mu_D_hat() const;

    /** \brief Determine if a constriant was an initially fixed variable.
     *
     * This function will return true if:
     * 
     * j <= n && x_init(j) != FREE
     * 
     * This is just a function of convienience
     * 
     */
    bool is_init_fixed( size_type j ) const;

    /// Returns true if all the degrees of freedom of the QP are used up
    bool all_dof_used_up() const;

    //@}

  private:

    // ///////////////////////////
    // Private types

    /** \brief . */
    typedef std::vector<int>			ij_map_t;
    /** \brief . */
    typedef std::map<int,size_type>		s_map_t;
    /** \brief . */
    typedef std::vector<EBounds>		bnds_t;
    /** \brief . */
    typedef std::vector<int>			l_fxfx_t;
    /** \brief . */
    typedef std::vector<size_type>		P_row_t;
    /** \brief . */
    typedef std::vector<size_type>		P_col_t;

    // ///////////////////////////
    // Private data members

    bool				initialized_;
    bool				test_;
    QP*					qp_;	// QP being solved.
    const QP::x_init_t	*x_init_;
    size_type			n_;
    size_type			n_R_;
    size_type			m_;
    size_type			m_breve_;
    size_type			q_plus_hat_;
    size_type			q_F_hat_;
    size_type			q_C_hat_;
    ij_map_t			ij_map_;
//		s_map_t				s_map_;
    DVector				constr_norm_;
    bnds_t				bnds_;
    l_fxfx_t            l_fxfx_;
    U_hat_t				U_hat_;
    //
    // for s = 1...q_hat
    //
    //                     /  e(i)    if i > 0 (where: i = -ij_map(s))
    // [P_XF_hat](:,s)   = |
    //                     \  0       otherwise
    //                     
    GenPermMatrixSlice	P_XF_hat_;		// \hat{P}^{XF} \in \Re^{n \times \hat{q}}
    P_row_t				P_XF_hat_row_;	// i
    P_row_t				P_XF_hat_col_;	// s
    //
    // for s = 1...q_hat
    //
    //                   /  e(sd)   if 0 < j <= n && is_init_fixed(j)
    //                   |          (where: j = ij_map(s), sd = s_map(-j))
    // [P_FC_hat](:,s) = |
    //                   \  0       otherwise
    //
    GenPermMatrixSlice	P_FC_hat_;		// {\tilde{P}^{F}}^{T} \hat{P}^{C} \in \Re^{\hat{q} \times \hat{q}}
    P_row_t				P_FC_hat_row_;	// sd
    P_row_t				P_FC_hat_col_;	// s
    //
    // for s = 1...q_hat
    //
    //                     /  e(j)    if j > 0 && !is_init_fixed(j) (where: j = ij_map(s))
    // [P_plus_hat](:,s) = |
    //                     \  0       otherwise
    //
    GenPermMatrixSlice	P_plus_hat_;	// \hat{P}^{(+)} \in \Re^{n+\breve{m} \times \hat{q}^{D}}
    P_row_t				P_plus_hat_row_;	// j
    P_row_t				P_plus_hat_col_;	// s
    //
    // for k = 1...q_D_hat
    //
    // [Q_XD_hat](:,k) = e(i)  (where is_init_fixed(i) && s_map(-i) == 0)
    //
    GenPermMatrixSlice	Q_XD_hat_;		// \hat{Q}^{XD} \in \Re^{n_X \times \hat{q}^{D}}
    P_row_t				Q_XD_hat_row_;	// i
    P_row_t				Q_XD_hat_col_;	// k
    //
    DVector				d_hat_;			// \hat{d}
    DVector				z_hat_;			// \hat{z}
    DVector				p_z_hat_;
    DVector				mu_D_hat_;		// \hat{\mu}^{D}
    DVector				p_mu_D_hat_;	// p^{\hat{\mu}^{D}}

    // ///////////////////////////
    // Private member functions

    //
    void assert_initialized() const;

    // Assert in range.
    void assert_s( size_type s) const;

    // Reinitialize P_XF_hat, P_plus_hat, Q_XD_hat, and U_hat
    void reinitialize_matrices(bool test);
    
    // Remove an element from the augmented KKT system.
    // This does not update P_plus_hat, P_XF_hat or any
    // of the dimensions.  Returns true if *out was
    // written to.
    bool remove_augmented_element(
      size_type sd, bool force_refactorization
      ,MatrixSymAddDelUpdateable::EEigenValType eigen_val_drop
      ,std::ostream *out, EOutputLevel output_level
      ,bool allow_any_cond );

    // not defined and not to be called.
    ActiveSet();

  };	// end class ActiveSet

  /// Return a reference to the active set object
  const ActiveSet& act_set() const;

  /// Dump all the active set quantities for debugging
  static void dump_act_set_quantities( const ActiveSet& act_set, std::ostream& out
    , bool print_S_hat = true );

protected:

  // /////////////////////////
  // Protected types

  /** \brief . */
  enum EPDSteps { PICK_VIOLATED_CONSTRAINT, UPDATE_ACTIVE_SET, COMPUTE_SEARCH_DIRECTION
    , COMPUTE_STEP_LENGTHS, TAKE_STEP };

  // ///////////////////////////
  // Protected Member functions

  /** \brief Run the algorithm from a dual feasible point.
   *
   * By default, the algorithm should start with
   * first_step = PICK_VIOLATED_CONSTRAINT if we are starting
   * with a dual feasible point.
   */
  virtual
  ESolveReturn qp_algo(
    EPDSteps first_step
    ,std::ostream *out, EOutputLevel output_level, ERunTests test_what
    ,const DVectorSlice& vo, ActiveSet* act_set, DVectorSlice* v
    ,DVectorSlice* x, size_type* iter, size_type* num_adds, size_type* num_drops
    ,size_type* iter_refine_num_resid, size_type* iter_refine_num_solves
    ,StopWatchPack::stopwatch* timer
    );

  /** \brief Set the values in x for all the variables.
   */
  virtual void set_x( const ActiveSet& act_set, const DVectorSlice& v, DVectorSlice* x );

  /// Map from the active set to the sparse multipliers for the inequality constraints
  virtual void set_multipliers(
    const ActiveSet& act_set, const DVectorSlice& v
    ,SpVector* mu, DVectorSlice* lambda, SpVector* lambda_breve );

  /// Determine if time has run out and if we should return.
  bool timeout_return( StopWatchPack::stopwatch*timer, std::ostream *out, EOutputLevel output_level ) const;

  /** \brief . */
  enum EIterRefineReturn {
    ITER_REFINE_NOT_PERFORMED    // Did not even perform it (iter_refine_max_iter == 0)
    ,ITER_REFINE_ONE_STEP        // Only performed one step and the status is not known.
    ,ITER_REFINE_NOT_NEEDED      // Convergence tolerance was already satisfied
    ,ITER_REFINE_IMPROVED        // Did not converge but it was improved
    ,ITER_REFINE_NOT_IMPROVED    // Tried iterative refinement but no improvement
    ,ITER_REFINE_CONVERGED       // Performed iterative refinement and converged!
  };
  /** \brief Perform iterative refinement on the augmented KKT system for the current active set.
   \verbatim

   [   Ko     U_hat ] [ v ] + [ ao * bo ]
   [ U_hat'   V_hat ] [ z ]   [ aa * ba ]
   \endverbatim
   * Returns \c true if iterative refinement satisfied the convergence criteria.
   */
    EIterRefineReturn iter_refine(
    const ActiveSet      &act_set
    ,std::ostream        *out
    ,EOutputLevel        output_level
    ,const value_type    ao  // Only used if bo != NULL
    ,const DVectorSlice   *bo // If NULL then assumed to be zero!
    ,const value_type    aa  // Only used if q_hat > 0
    ,const DVectorSlice   *ba // If NULL then assumed to be zero!  Not accessed if q_hat > 0
    ,DVectorSlice         *v
    ,DVectorSlice         *z  // Can be NULL if q_hat > 0
    ,size_type           *iter_refine_num_resid
    ,size_type           *iter_refine_num_solves
    );

private:

  // /////////////////////////
  // Private data members

  ActiveSet  act_set_;  // The active set.

};	// end class QPSchur

}	// end namespace ConstrainedOptPack 

#endif	// QPSCHUR_H
