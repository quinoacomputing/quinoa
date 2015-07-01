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

#ifndef QP_SOLVER_RELAXED_TESTER_H
#define QP_SOLVER_RELAXED_TESTER_H

#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Tests the optimality conditions of the output from a \c QPSolverRelaxed
 * object.
 *
 * For the given QP and its solution (if solved) this class tests
 * the optimality conditions.
 * 
 * The optimality conditions checked are:
 \verbatim

  Linear dependence of gradients:
    
  (2)  d(L)/d(d) = g + G*d - nuL + nuU + op(E)'*(- muL + muU) + op(F)'*lambda
                 = g + G*d + nu + op(E)'*mu + op(F)'*lambda = 0
    
       where: nu = nuU - nuL, mu = muU - muL

  Feasibility:
    
  (4.1)  etaL <=  eta
  (4.2)  dL   <=  d                       <= dU
  (4.3)  eL   <=  op(E)*d - b*eta         <= eU
  (4.4)  op(F)*d + (1 - eta) * f  = 0

  Complementarity:

  (5.1)  nu(i) * (dL - d)(i), if nu(i) <= 0, i = 1...n
  (5.2)  nu(i) * (d - dU)(i), if nu(i) >= 0, i = 1...n
  (5.3)  mu(j) * (eL - op(E)*d + b*eta)(j), if mu(j) <= 0, j = 1...m_in
  (5.4)  mu(j) * (op(E)*d - b*eta - eU)(j), if mu(j) >= 0, j = 1...m_in
 \endverbatim
 * The realtive error of each of these conditions is checked.  Specifically,
 * here is how the errors are computed which are compared to the error and warning
 * tolerances:
 \verbatim
  
    opt_err = || g + G*d + nu + op(E)'*mu + op(F)'*lambda ||inf / (1 + opt_scale)
                  
    feas_err = ||max(op(A)*x-b,0)||inf / ( 1 + ||op(A)*x||inf )

    comp_err(i) = gamma(i) * (op(A)*x - b)(i) / ( 1 + opt_scale + ||op(A).row(i)'*x||inf )
        ,for gamma(i) != 0

     where:
      op(A)*x <= b
      opt_scale = ||g||inf + ||G*d||inf + ||nu||inf + ||op(E)'*mu||inf + ||op(F)'*lambda||inf
                    + |(eta - etaL) * (b'*mu + f'*lambda)|
 \endverbatim
 * Above, <tt>op(A)*x <= b</tt> can represent any of the constraints in
 * (4.1)-(4.4).
 *
 * Any elements of <tt>opt_err(i) >= opt_warning_tol</tt> will result in an error
 * message printed to \c *out.  Any elements of <tt>opt_err(i) >= opt_error_tol</tt>
 * will cause the checks to stop and false to be returned from the function
 * \c check_optimiality_conditions().
 *
 * Any elements of <tt>feas_err(i) >= feas_warning_tol</tt> will result in an error
 * message printed to \c *out.  Any elements of <tt>feas_err(i) >= feas_error_tol</tt>
 * will cause the checks to stop and false to be returned from the function
 * \c check_optimiality_conditions().
 *
 * Any elements of <tt>comp_err(i) >= comp_warning_tol</tt> will result in an error
 * message printed to *out.  Any elements of <tt>comp_err(i) >= comp_error_tol</tt>
 * will cause the checks to stop and false to be returned from the function
 * \c check_optimiality_conditions().
 *
 * The goal of these tests is to first and foremost to catch gross programming
 * errors.  These tests can also be used to help flag and catch illconditioning
 * in the problem or instability in the QP solver.  The importance of such tests
 * can not be overstated.  The scalings above are done to try to adjust for
 * the scaling of the problem.  Note that we are accounting for very big numbers
 * but not for very small numbers very well and therefore tests may be conserative
 * in some cases.  At the very least we account for loss of precision due to 
 * catastrophic cancelation that occurs when subtracting large numbers and expecting
 * to get zero.  The purpose of including the term <tt>|b'*mu + f'*lambda|</tt> is to
 * account for the situation where the relaxation is needed and <tt>kappa != 0</tt>
 * and therefore, form the condition <tt>d(M)/d(eta) - kappa - b'*mu - f'*lambda = 0</tt>
 * (with <tt>d(M)/d(eta)</tt> very large), the multipliers \c mu and \c lambda will be
 * very large and will contribute to much roundoff errors.
 *
 * As shown above, the complementarity conditions (5.1)-(5.4) are specifically checked.
 * These should be satisfied for any solution type other than a \c SUBOPTIMAL_POINT
 * return value from \c QPSolverRelaxed::solve_qp().  Checking the complementarity
 * conditions for an active-set QP solver is just checking that the active constraints
 * are satisfied. Checking them for an iterior-point solver is critical to ensure that
 * the system was solved to satisfactory tolerance.
 * By scaling the active-constraint violation by the Langrange multiplier
 * we emphasis the feasibility of those constraints that have the greatest
 * impact on the objective function.  In other words, all things being equal, we are
 * more concerned with a tight feasibility tolerance for constraints with
 * larger lagrange multipliers than for those with smaller multipliers.
 * The complementarity error is also scaled by the inverse of the sums
 * of the optimality scaling opt_scale and the size of the constraint residual.
 * By scaling by the max term \c opt_scale in the linear dependence of gradients we are
 * trying to adjust the effect of the lagrange multiplier.  Therefore, if the gradient
 * of the objective <tt>g+G*d</tt> is large then \c opt_scale will account for this.
 */
class QPSolverRelaxedTester {
public:

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_warning_tol );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_error_tol );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_warning_tol );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_error_tol );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, comp_warning_tol );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, comp_error_tol );

  /** \brief . */
  QPSolverRelaxedTester(
    value_type    opt_warning_tol   = 1e-10
    ,value_type   opt_error_tol     = 1e-5
    ,value_type   feas_warning_tol  = 1e-10
    ,value_type   feas_error_tol    = 1e-5
    ,value_type   comp_warning_tol  = 1e-10
    ,value_type   comp_error_tol    = 1e-5
    );

  /** \brief . */
  virtual ~QPSolverRelaxedTester() {}

  /** \brief Check the optimality conditions for the solved (or partially solved) QP.
   *
   * The default implementation calls the function \c check_optimality_conditions()
   * which accepts various sets of constraints.
   *
   *	@param	solution_type
   *						[in] Value returned from \c QPSolverRelaxed::solve_qp().
   *							Even though all of the optimality conditions are
   *							checked the optimality conditions that are actually
   *							enforced is determined by this argument.
   *								OPTIMAL_SOLUTION : All of the optimality conditions
   *									are enforced.
   *								PRIMAL_FEASIBLE_POINT : Only the optimality conditions
   *								 	(4.1)-(4.4) are enforced.
   *								DUAL_FEASIBLE_POINT: Only the optimality condtions
   *								 	(2) and (6.1)-(6.4) are enforced.
   *								SUBOPTIMAL_POINT : None of the optimality conditions
   *								 	are enforced.
   *	@param	out			[out] If <tt>!=NULL</tt>, the output is sent to this stream.
   *	@param	print_all_warnings
   *						[in] If \c true, then all errors greater than \c warning_tol will
   *							be printed.
   *	@param	g			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	G			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	etaL		[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	dL			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	dU			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	E			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	trans_E		[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	b			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	eL			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	eU			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	F			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	trans_F		[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	f			[in] Input to \c QPSolverRelaxed::solve_qp().
   *	@param	obj_d		[in] Output from \c QPSolverRelaxed::solve_qp().
   *	@param	eta			[in] Output from \c QPSolverRelaxed::solve_qp().
   *	@param	d 			[in] Output from \c QPSolverRelaxed::solve_qp().
   *	@param	nu			[in] Output from \c QPSolverRelaxed::solve_qp().
   *	@param	mu 			[in] Output from \c QPSolverRelaxed::solve_qp().
   *	@param	Ed			[in] Output from \c QPSolverRelaxed::solve_qp().
   *	@param	Fd			[in] Output from \c QPSolverRelaxed::solve_qp().
   *
   * @return <tt>true</tt> if all of the errors are greater than the error tolerances
   * 	, otherwise it returns <tt>false</tt>
   */
  virtual bool check_optimality_conditions(
    QPSolverStats::ESolutionType solution_type
    ,const value_type infinite_bound
    ,std::ostream* out, bool print_all_warnings, bool print_vectors
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector& dL, const Vector& dU
    ,const MatrixOp& E, BLAS_Cpp::Transp trans_E, const Vector& b
    ,const Vector& eL, const Vector& eU
    ,const MatrixOp& F, BLAS_Cpp::Transp trans_F, const Vector& f
    ,const value_type* obj_d
    ,const value_type* eta, const Vector* d
    ,const Vector* nu
    ,const Vector* mu, const Vector* Ed
    ,const Vector* lambda, const Vector* Fd
    );

  /** \brief Check the optimality conditions without general equality constrants.
   */
  virtual bool check_optimality_conditions(
    QPSolverStats::ESolutionType solution_type
    ,const value_type infinite_bound
    ,std::ostream* out, bool print_all_warnings, bool print_vectors
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector& dL, const Vector& dU
    ,const MatrixOp& E, BLAS_Cpp::Transp trans_E, const Vector& b
    ,const Vector& eL, const Vector& eU
    ,const value_type* obj_d
    ,const value_type* eta, const Vector* d
    ,const Vector* nu
    ,const Vector* mu, const Vector* Ed
    );

  /** \brief Check the optimality conditions general inequality constrants.
   */
  virtual bool check_optimality_conditions(
    QPSolverStats::ESolutionType solution_type
    ,const value_type infinite_bound
    ,std::ostream* out, bool print_all_warnings, bool print_vectors
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector& dL, const Vector& dU
    ,const MatrixOp& F, BLAS_Cpp::Transp trans_F, const Vector& f
    ,const value_type* obj_d
    ,const value_type* eta, const Vector* d
    ,const Vector* nu
    ,const Vector* lambda, const Vector* Fd
    );


  /** \brief Check the optimality conditions without general equality or inequality
   * constrants (no relaxation needed).
   */
  virtual bool check_optimality_conditions(
    QPSolverStats::ESolutionType solution_type
    ,const value_type infinite_bound
    ,std::ostream* out, bool print_all_warnings, bool print_vectors
    ,const Vector& g, const MatrixSymOp& G
    ,const Vector& dL, const Vector& dU
    ,const value_type* obj_d
    ,const Vector* d
    ,const Vector* nu
    );

  /** \brief This is a more flexible function where the client can
   * set different constraints to be included.
   *
   */
  virtual bool check_optimality_conditions(
    QPSolverStats::ESolutionType solution_type
    ,const value_type infinite_bound
    ,std::ostream* out, bool print_all_warnings, bool print_vectors
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

protected:

  /** \brief Subclasses are to override this to implement the testing code.
   *
   * There is a default implementation that is very general and
   * should be considered good enough for most applications.
   */
  virtual bool imp_check_optimality_conditions(
    QPSolverStats::ESolutionType solution_type
    ,const value_type infinite_bound
    ,std::ostream* out, bool print_all_warnings, bool print_vectors
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

};	// end class QPSolverRelaxedTester

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// QP_SOLVER_RELAXED_TESTER_H
