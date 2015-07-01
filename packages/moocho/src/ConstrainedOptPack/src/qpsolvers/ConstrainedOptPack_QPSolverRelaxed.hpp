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

#ifndef QP_SOLVER_RELAXED_H
#define QP_SOLVER_RELAXED_H

#include "ConstrainedOptPack_QPSolverStats.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Solves Quadratic Programs (QPs) of several different forms while
 * allowing a relaxation of the constraints.
 *
 * The formulation for the QP being solved is:
 \verbatim
  
  (1.1)  min          g'*d + 1/2*d'*G*d + M(eta)
         d <: R^n
         
         s.t.
  (1.2)               etaL <=  eta
  (1.3)               dL   <=  d                       <= dU
  (1.4)               eL   <=  op(E)*d - b*eta         <= eU
  (1.5)                        op(F)*d + (1 - eta) * f  = 0
 \endverbatim
 * The relaxation is used to ensure that the QP will have a solution
 * (<tt>eta = 1</tt>, <tt>d = 0</tt> guaranteed if <tt>dL <= 0 <= dU</tt>
 * and <tt>eL <= b <= eU</tt>).  If the function <tt>M(eta)</tt> in the
 * objective is large enough, then the constraint <tt>etaL <= eta</tt> will be active
 * if a feasible region exists.  The form of the function <tt>M(eta)</tt> is not
 * specified by this interface but is defined as appropriate for each individual
 * QP solver method and implementation. 
 *
 * The Lagrangian for the QP in (1) is:
 \verbatim
  
  L = g' * d + 1/2 * d' * G * d + M(eta)
       + kappa * (etaL - eta)
       + nuL' * (dL - d)
       + nuU' * (d - dU)
       + muL' * (eL - op(E)*d - b*eta)
       + muU' * (op(E)*d - b*eta - eU)
     + lambda' * (op(F)*d + (1 - eta) * f)
 \endverbatim
 * The optimality conditions for this QP are:
 \verbatim

  Linear dependence of gradients:
    
  (2)  d(L)/d(d) = g + G*d - nuL + nuU + op(E)'*(- muL + muU) + op(F)'*lambda
                 = g + G*d + nu + op(E)'*mu + op(F)'*lambda = 0
    
       where: nu = nuU - nuL, mu = muU - muL

  (3)  d(L)/d(eta) = d(M)/d(eta) - kappa - b'*mu - f'*lambda = 0
    
  Feasibility:
    
  (4.1)  etaL <=  eta
  (4.2)  dL   <=  d                       <= dU
  (4.3)  eL   <=  op(E)*d - b*eta         <= eU
  (4.4)  op(F)*d + (1 - eta) * f  = 0

  Complementarity:

  (5.1)  nu(i) * (dL - d)(i) = 0, if nu(i) <= 0, i = 1...n
  (5.2)  nu(i) * (d - dU)(i) = 0, if nu(i) >= 0, i = 1...n
  (5.3)  mu(j) * (eL - op(E)*d + b*eta)(j) = 0, if mu(j) <= 0, j = 1...m_in
  (5.4)  mu(j) * (op(E)*d - b*eta - eU)(j) = 0, if mu(j) >= 0, j = 1...m_in

  Nonnegativity of Lagrange Multipliers for Inequality Constraints:

  (6.1)  nu(i) <= 0 if (dL - d)(i) = 0, i = 1...n
  (6.2)  nu(i) >= 0 if (d - dU)(i) = 0, i = 1...n
  (6.3)  mu(j) <= 0 if (eL - op(E)*d - b*eta)(j) = 0, j = 1...m_in
  (6.4)  mu(j) >= 0 if (op(E)*d - b*eta - eU)(j) = 0, j = 1...m_in
 \endverbatim
 * The optimal <tt>d</tt> and <tt>eta</tt> are determined as well as the lagrange multipliers
 * for the constriants:
 \verbatim

  nu :       dL <= d <= dU
 
  mu :       eL <= op(E)*d - b*eta <= eU

  lambda :   op(F)*d + (1 - eta) * f  = 0
 \endverbatim
 * The lagrange multiper <tt>kappa</tt> for the constraint <tt>etaL <= eta</tt> is not returned
 * since if this constraint is not active, then <tt>kappa == 0</tt> and all of the multiplier
 * estimates will be off because of the arbitrarily large value of <tt>d(M)/d(eta)</tt> in
 * the optimality condition (3).
 */
class QPSolverRelaxed {
public:

  /** @name Public Types */
  //@{

  /// Thrown if the QP is unbounded.
  class Unbounded : public std::logic_error
  {public: Unbounded(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if the QP is infeasible.
  class Infeasible : public std::logic_error
  {public: Infeasible(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if there is invalid input
  class InvalidInput : public std::logic_error
  {public: InvalidInput(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if a test failed
  class TestFailed : public std::logic_error
  {public: TestFailed(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Enumeration for the amount of output to create from \c solve_qp().
  enum EOutputLevel {
    PRINT_NONE			= 0,
    PRINT_BASIC_INFO	= 1,
    PRINT_ITER_SUMMARY	= 2,
    PRINT_ITER_STEPS	= 3,
    PRINT_ITER_ACT_SET	= 4,
    PRINT_ITER_VECTORS	= 5,
    PRINT_EVERY_THING	= 6
    };

  /// Enumeration for if to run internal tests or not.
  enum ERunTests { RUN_TESTS, NO_TESTS };

  //@}

  /** @name Initializers */
  //@{

  /// Set the scalar that will be used to identigy infinite bounds
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, infinite_bound );

  /** \brief . */
  QPSolverRelaxed();

  /** \brief . */
  virtual ~QPSolverRelaxed() {}

  //@}

  /** @name Interface methods with default implementations */
  //@{

  /** \brief Solve the QP.
   *
   *	@param	out			[out] If out != NULL then output is printed to this stream
   *							depending on the value of <tt>olevel</tt>.
   *	@param	olevel		[in] Determines the amount of output to print to *out.
   *							The exact type of output is determined by the implementing
   *							subclass but here is the sugguested behavior:
   *							PRINT_NONE : Don't print anything (same as out == NULL).
   *							PRINT_BASIC_INFO : Only print basic information about
   *								the solution of the QP.  Amount of output = O(1).
   *							PRINT_ITER_SUMMARY : Prints a summary of each iteration
   *								in the algorithm.  Amount of output = O(num_iter).
   *							PRINT_ITER_STEPS : Prints output about each iteration
   *								in more detail than PRINT_ITER_SUMMARY but is still
   *								O(num_iter).
   *							PRINT_ITER_VECTORS : Prints out important vectors computed
   *								in each QP iteration.  Mainly useful for debugging.
   *								Amount of output = O((num_iter)(n)).
   *							PRINT_EVERY_THING : Print out nearly every important
   *								quantity within each QP iteration including
   *								vectors and matrices.  Mainly useful for debugging.
   *								Amount of output = O((num_iter)(n-m)^2)+O((num_iter)(n))
   *	@param	test_what	[in] Determines if internal validation tests are performed.
   *							The optimality conditions for the QP are not checked
   *							internally, this is something that client can (and should)
   *							do independently (see QPSolverRelaxedTester).
   *							RUN_TESTS : As many validation/consistency tests
   *								are performed internally as possible.  If a test
   *								fails then a TestFailed execption will be thrown.
   *								The subclasses determine what the tests are and
   *								what failing a test means.
   *							NO_TEST : No tests are performed internally.  This is
   *								to allow the fastest possible execution.
   *	@param	g			[in] vector (size n): objective 1st order
   *	@param	G			[in] matrix (size n x n): objective, second order Hessian
   *	@param	etaL		[in] scalar: Lower bound for relaxation variable (0 usually)
   *	@param	dL			[in] vector (size n): lower variable bounds
   *	@param	dU			[in] vector (size n): upper variable bounds
   *	@param	E			[in] matrix (op(E)) (size m_in x n): inequality constraint Jacobian matrix 
   *	@param	trans_E		[in] E is transposed? 
   *	@param	b			[in] vector (size m_in): relaxation vector for inequalities
   *	@param	eL			[in] vector (size m_in): lower inequality bounds
   *	@param	eU			[in] vector (size m_in): upper inequality bounds
   *	@param	F			[in] matrix (op(F) size m_eq x n): equality constraint Jacobian matrix
   *	@param	trans_F		[in] F is transposed? 
   *	@param	f			[in] vector (size m_eq): equality constraint right hand side
   *	@param	obj_d		[out] If obj_d != NULL on input, then obj_d will be set with the
   *							value of obj_d = g'*d + 1/2*d'*G*d for the value of d
   *							returned.
   *	@param	eta			[out] scalar:  Relaxation variable
   *	@param	d			[in/out] vector (size n):  On input, it contains an intial estimate
   *                         of the solution.  On output it is the estimate of the solution
   *                         (see the return value).
   *	@param	nu			[in/out] vector (size n):  Lagrange multipilers
   *							for variable bounds.  On input it contains the
   *							estimate of the active set and multiplier values.  If
   *							nu->nz() == 0 on input then there is no estimate for the
   *							active-set.  Note that having nu->nz() > 0 on input is not
   *							a commandment to perform a warm start.  Ultimatly this
   *							decision is up to the subclass and lower level
   *							subclass/client interactions.
   *							On output nu contains the active-set for the
   *							returned solution.
   *	@param	mu			[in/out] vector (size m_in):  Lagrange multipilers
   *							for the general inequality constriants.
   *							On input it contains the
   *							estimate of the active set and multiplier values.  If
   *							mu->nz() == 0 on input then there is no estimate for the
   *							active-set.  Note that having mu->nz() > 0 on input is not
   *							a commandment to perform a warm start.  Ultimatly this
   *							decision is up to the subclass and lower level
   *							subclass/client interactions.
   *							On output mu contains the active-set for the
   *							returned solution.
   *	@param	Ed			[in/out] vector (size m_in) If Ed!=NULL on input, then on output
   *							Ed will contain the product opt(E)*d for the value of
   *							d on output.  This is included to save user from having
   *							to perform this computation again if it has already been
   *							done internally.
   *	@param	lambda		[out] vector (size m_eq):  Lagrange multipilers for equality
   *							constraints.
   *	@param	Fd			[in/out] vector (size m_eq) If Fd!=NULL on input, then on output
   *							Fd will contain the product opt(F)*d for the value of
   *							d on output.  This is included to save user from having
   *							to perform this computation again if it has already been
   *							done internally.
   *
   *	@return
   *		<tt>OPTIMAL_SOLUTION</tt> : Returned point satisfies the optimality conditions
   *			in (2)-(6) above.  This will generally be the case if the maximum
   *			number of QP iterations is not exceeded and none of the possible
   *			exeptions are thrown.
   *		<tt>PRIMAL_FEASIBLE_POINT</tt> : Returned point satisfies the feasibility
   *         and complementarity conditions in (4)-(5) above.
   *			For example, a primal, active-set QP algorithm may return this if
   *			the maxinum number of iterations has been exceeded durring the
   *			optimality phase (phase 2).
   *		<tt>DUAL_FEASIBLE_POINT</tt> : Returned point satisfies the optimality conditions
   *			in (2)-(4.1),(4.4),(5) and (6) but not the inequality constraints in (4.2),(4.3).
   *			For example, a dual, active-set QP algorithm might return this
   *			value if the maximum number of iterations has been exceeded.
   *		<tt>SUBOPTIMAL_POINT</tt> : Returned point does not accurately enough satisfy any of the
   *			optimality conditions in (2)-(6) above but the solution may still be
   *			of some use.  For example, an active-set (primal, or dual) QP algorithm
   *			might return this if there is some serious illconditioning in the QP
   *			and a solution satisfying the desired tolerance could not be found.
   *			Also, a primal-dual interior point QP solver might return this if the
   *			maxinum number of iterations is exceeded.  The returned solution may
   *			still be of some use to the client though.
   *
   * This function may throw many exceptions.  If there is some problem with
   * the QP definition the exceptions <tt>Unbounded</tt>, <tt>Infeasible</tt> or <tt>InvalidInput</tt>
   * may be thrown.  If the QP is illconditioned other exeptions may be thrown
   * by this function or perhaps warning messages may be printed to <tt>*out</tt> and
   * a value other than <tt>OPTIMAL_SOLUTION</tt> may be returned.
   *
   * After this function returns, <tt>this->get_qp_stats()</tt> can be called to return
   * some statistics for the QP just solved (or attempted to be solved).
   *
   * Note, the variable bounds can be removed by passing in <tt>dL.nz() == dU.nz() == 0</tt>.
   * 
   * By default this function calls the function <tt>this->solve_qp()</tt> which accepts
   * various sets of constraints.
   */
  virtual QPSolverStats::ESolutionType solve_qp(
    std::ostream* out, EOutputLevel olevel, ERunTests test_what
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector& dL, const Vector& dU
    ,const MatrixOp& E, BLAS_Cpp::Transp trans_E, const Vector& b
    ,const Vector& eL, const Vector& eU
    ,const MatrixOp& F, BLAS_Cpp::Transp trans_F, const Vector& f
    ,value_type* obj_d
    ,value_type* eta, VectorMutable* d
    ,VectorMutable* nu
    ,VectorMutable* mu, VectorMutable* Ed
    ,VectorMutable* lambda, VectorMutable* Fd
    );

  /** \brief Solve the QP without general equality constrants.
   *
   * By default this function calls <tt>solve_qp()</tt> which accepts
   * various sets of constraints.
   */
  virtual QPSolverStats::ESolutionType solve_qp(
    std::ostream* out, EOutputLevel olevel, ERunTests test_what
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector& dL, const Vector& dU
    ,const MatrixOp& E, BLAS_Cpp::Transp trans_E, const Vector& b
    ,const Vector& eL, const Vector& eU
    ,value_type* obj_d
    ,value_type* eta, VectorMutable* d
    ,VectorMutable* nu
    ,VectorMutable* mu, VectorMutable* Ed
    );

  /** \brief Solve the QP without general inequality constrants.
   *
   * By default this function calls <tt>solve_qp()</tt> which accepts
   * various sets of constraints.
   */
  virtual QPSolverStats::ESolutionType solve_qp(
    std::ostream* out, EOutputLevel olevel, ERunTests test_what
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector& dL, const Vector& dU
    ,const MatrixOp& F, BLAS_Cpp::Transp trans_F, const Vector& f
    ,value_type* obj_d
    ,value_type* eta, VectorMutable* d
    ,VectorMutable* nu
    ,VectorMutable* lambda, VectorMutable* Fd
    );


  /** \brief Solve the QP without general equality or inequality constrants (no relaxation
   * needed).
   *
   * By default this function calls <tt>solve_qp()</tt> which accepts
   * various sets of constraints.
   */
  virtual QPSolverStats::ESolutionType solve_qp(
    std::ostream* out, EOutputLevel olevel, ERunTests test_what
    ,const Vector& g, const MatrixSymOp& G
    ,const Vector& dL, const Vector& dU
    ,value_type* obj_d
    ,VectorMutable* d
    ,VectorMutable* nu
    );

  /** \brief This is a more flexible function where the client can
   * set different constraints to be included.
   *
   * The default implementation of this function validates the compatibily
     * of the input objects and that the proper sets of constraints are set
   * by calling \c validate_input() first.  Refere to the method
   * \c validate_input() to see how arguments are set.  After input validation,
   * \c print_qp_input() is called to print the QP input arguments.  Then
   * \c imp_solve_qp() is called which must be implemented by the subclass.
   * Finally, \c print_qp_output() is called to print the QP output.
   */
  virtual QPSolverStats::ESolutionType solve_qp(
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

  /** \brief Get the statistics of the last QP solved.
   *
   *	solution:	Returns the type of solution found.  See \c solve_qp()<br>
   *	num_iter:	Gives the number of QP iterations on output.<br>
   *	num_adds:	Gives the number of iterations where a variable
   *				was added to the active-set.  This does not include
   *				variables that where part of the initial estimate in nu
   *				for a warm start.<br>
   *	num_drops:	Gives the number QP iterations where a variable was
   *				dropped.<br>
   *	warm_start:	Returns if a warm start was performed.<br>
   *	infeas_qp:	Returns if eta > 0.0.
   */
  virtual QPSolverStats get_qp_stats() const = 0;

  /** \brief Release any memory that is being used.
   */
  virtual void release_memory() = 0;

  //@}

  /** @name Static utility functions */
  //@{

  /** \brief This is a (static) function that is used as a utility to
   * validate the input arguments to \c solve_qp().
   * 
   * The input arguments are validated as follows.
   * 
   * If the variable bounds are excluded then:
   * <tt>void(dL) == void(dU) == void(nu) == NULL</tt>
   *
   * If the general inequality constraints are excluded then:
   * <tt>void(E) == void(b) == void(eL) == void(eU) == void(mu) == NULL</tt>
   *
   * If the equality constraints are excluded then:
   * <tt>void(F) == void(f) == void(lambda) == NULL</tt>
   *
   * If the equality and inequality constraints are excluded then:
   * <tt>eta == NULL</tt>
   * <tt>void(E) == void(b) == void(eL) == void(eU) == void(mu) == NULL</tt>
   * <tt>void(F) == void(f) == void(lambda) == NULL</tt>
   * 
   * If any errors are found an std::invalid_argument exception
   * will be thrown.
   */
  static void validate_input(
    const value_type infinite_bound
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector* dL, const Vector* dU
    ,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
    ,const Vector* eL, const Vector* eU
    ,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
    ,const value_type* obj_d
    ,const value_type* eta, const Vector* d
    ,const Vector* nu
    ,const Vector* mu, const Vector* Ed
    ,const Vector* lambda, const Vector* Fd
    );

  /** \brief Utility (static) function for printing the input input/output arguments before
   * the QP solver is run.  The QP solver subclasses can call this function.
   *
   * @param out [out] stream printed to if <tt>out != NULL</tt>
   * @param olevel [in] Determines what is printed.
   *   \begin{description}
   *   \item[(int)olevel >= (int)PRINT_ITER_STEPS] Prints O(1) information about the arguments.
   *   \item[(int)olevel >= (int)PRINT_ITER_ACT_SET] Prints the contents of <tt>nu</tt>, <tt>mu</tt>, and <tt>lambda</tt>.
   *      Output is proportional to the number of active constraints O(nu->nz() + mu->nz() + lambda->dim()). 
   *   \item[(int)olevel >= (int)PRINT_ITER_VECTORS] Prints the contents of all the vectors.
   *      Output is proportional to O(d->dim()).
   *   \item[(int)olevel >= (int)PRINT_EVERY_THING] Prints the contents of all the vectors and matrices.
   *      Output could be as large as O(d->dim() * d->dim()) or larger.
   *   \end{description}
    */
  static void print_qp_input( 
    const value_type infinite_bound
    ,std::ostream* out, EOutputLevel olevel
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector* dL, const Vector* dU
    ,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
    ,const Vector* eL, const Vector* eU
    ,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
    ,value_type* eta, VectorMutable* d
    ,VectorMutable* nu
    ,VectorMutable* mu
    ,VectorMutable* lambda
    );

  /** \brief Utility (static) function for printing the output input/output arguments after
   * the QP solver is run.  The QP solver subclasses can call this function.
   *
   * @param out [out] stream printed to if <tt>out != NULL</tt>
   * @param olevel [in] Determines what is printed.
   *   \begin{description}
   *   \item[(int)olevel >= (int)PRINT_ITER_STEPS] Prints O(1) information about the arguments.
   *   \item[(int)olevel >= (int)PRINT_ITER_ACT_SET] Prints the contents of <tt>nu</tt>, <tt>mu</tt>, and <tt>lambda</tt>.
   *      Output is proportional to the number of active constraints O(nu->nz() + mu->nz() + lambda->dim()). 
   *   \item[(int)olevel >= (int)PRINT_ITER_VECTORS] Prints the contents of all the vectors.
   *      Output is proportional to O(d->dim()).
   *   \item[(int)olevel >= (int)PRINT_EVERY_THING] Prints the contents of all the vectors and matrices.
   *      Output could be as large as O(d->dim() * d->dim()) or larger.
   *   \end{description}
    */
  static void print_qp_output(
    const value_type infinite_bound
    ,std::ostream* out, EOutputLevel olevel
    ,const value_type* obj_d
    ,const value_type* eta, const Vector* d
    ,const Vector* nu
    ,const Vector* mu, const Vector* Ed
    ,const Vector* lambda, const Vector* Fd
    );

  //@}

protected:

  /** @name Pure virtual methods that must be overridden by subclass */
  //@{

  /** \brief Subclasses are to override this to implement the QP algorithm.
   *
   * Called by default implementations of \c solve_qp() methods.
   */
  virtual QPSolverStats::ESolutionType imp_solve_qp(
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
    ) = 0;

  //@}

};	// end class QPSovlerRelaxed

}	// end namespace ConstrainedOptPack

#endif	// QP_SOLVER_RELAXED_H
