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

#ifndef NLP_FIRST_ORDER_DIRECT_TESTER_H
#define NLP_FIRST_ORDER_DIRECT_TESTER_H

#include <iosfwd>

#include "NLPInterfacePack_Types.hpp"
#include "NLPInterfacePack_CalcFiniteDiffProd.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace NLPInterfacePack {

/** \brief Concrete class that tests the computed values of the
  * <tt>NLPDirect</tt> interface using finite differences.
  *
  * There are two options for testing the derivatives by finite differences.
  * Each option can be picked independently for the computations with
  * the objective \c f(x) and its gradient \c Gf and for the constraints
  * \a c(x) and its Jacobian <tt>Gc'</tt>.  The tests involving the objective
  * and the constraints will be discussed separatly.
  * 
  * For testing the gradient of the objective function, two options
  * are available.  The first option (<tt>Gf_testing_method==FD_COMPUTE_ALL</tt>)
  * is to compute \c FDGf by brute force which requries <tt>2*n</tt> evaluations of
  * \a f(x) using central differences.  Given <tt>FDGf</tt> the following comparison
  * is then made:
  \verbatim

    (1)    FDGf \approx Gf
  \endverbatim 
  * The other option (<tt>Gf_testing_method==FD_DIRECTIONAL</tt>) is to compute
  * random dot products of the form <tt>DFGf'*y</tt> where \c y is a randomly generated
  * vector.  Using central differences <tt>DFGf'*y</tt> can be computed using
  * two evaluations of \a f(x) per random \c y.  The number of random
  * <tt>y</tt>'s used is determined by the option <tt>num_fd_directions()</tt>.  So the
  * number of evaluations of \a f(x) for this option is <tt>2*num_fd_directions()</tt>.
  * 
  * The test for the quantity <tt>py = -inv(C)*c(con_decomp)</tt> is shown below:
  \verbatim

    (2)  - FDC * (-inv(C)*c) \approx c(c_decomp)
                 \_________/
                     py
  \endverbatim
  * Computing <tt>-FDC * py</tt> requires only two evaluations of
  * \a c(x) using central differences.  There is no other option needed
  * for this test.
  * 
  * Lastly, we have to test <tt>D = -inv(C)*N</tt>.  The first option
  * (<tt>Gc_testing_method==FD_COMPUTE_ALL</tt>) is to directly compute
  * \c N using central differences (<tt>2*(n-m)</tt> evaluations of \a c(x)) as
  * \c FDN and then perform the comparison:
  \verbatim 

    (3)  - FDC * (-inv(C)*N) \approx FDN
                 \_________/
                      D
  \endverbatim 
  * The matrix <tt>-FDC * D</tt> can be computed using <tt>2*(n-m)</tt> evaluations
  * with \a c(x) using central differences.
  * Therefore, the total number of evaluations with \a c(x) for
  * comparing (3) is <tt>4*(n-m)</tt>.  If <tt>n-m</tt> is not too large then this
  * is definitely the preferred method to use.
  * 
  * The other option for testing <tt>D = -inv(C)*N</tt> is to compute
  * directional derivatives using finite differences.  In this approach,
  * for the random vector \c y, we can compute:
  \verbatim

    (4)  - FDC * (-inv(C)*N) * y \approx FDN * y
                 \_________/
                      D
  \endverbatim 
  * Using central differences, (4) can be computed with 4 evaluations
  * of \a c(x).  The number of random <tt>y</tt>'s used is determined by the option
  * \c num_fd_directions().  So the number of evaluations of \a c(x) for this
  * option is <tt>4*num_fd_directions()</tt>.
  * 
  * The client can pick a set of tolerances to measure if the
  * values of the above comparisons are close enough to the finite difference
  * values.  Let's define the relative error between the computed value and the
  * finite difference value to be:
  \verbatim

  err(i) = | (h(i) - fdh(i)) | /  (||h||inf + ||fdh||inf + sqrt(epsilon))
  \endverbatim
  * The above error takes into account the relative sizes of the elements and also
  * allows one or both of the elements to be zero without ending up with <i>0/0</i>
  * or something like <tt>1e-16</tt> not comparing with zero.
  *
  * All errors <tt>err(i) >= warning_tol</tt> are reported to <tt>*out</tt> if
  * <tt>out != NULL</tt>.  The first error <tt>err(i) >= error_tol</tt> that is found
  * is reported to <tt>*out</tt> if <tt>out != NULL</tt> and immediatly
  * <tt>finite_diff_check()</tt> returns \c false.  If all errors <tt>err(i) < error_tol</tt>,
  * then <tt>finite_diff_check()</tt> will return \c true.
  *
  * Given these two tolerances the client can do many things:
  *
  * 1) Print out all the comparisons that are not equal by setting <tt>warning_tol
  *    <= epsilon</tt> and <tt>error_tol >> 1</tt>.
  *
  * 2) Print out all suspect comparisons by setting <tt>epsilon < warning_tol < 1</tt>
  *    and <tt>error_tol >> 1</tt>.
  *
  * 3) Just validate that the quantities are approximatly equal and report the first
  *    discrepency if not by setting <tt>epsilon < error_tol < 1</tt> and <tt>warning_tol
  *    >= error_tol</tt>.
  *
  * 4) Print out any suspect comparisons by setting <tt>epsilon < warning_tol < 1</tt>
  *    but also quit if the error is too large by setting <tt>1 > error_tol > warning_tol</tt>.
  *
  * The tolerances \c Gf_warning_tol and \c Gf_error_tol are applied to the tests for
  * \c Gf shown in (1) for instance.  The tolerances \c Gc_warning_tol and \c Gc_error_tol
  * are used for the comparisions (2), (3) and (4).
  * 
  * There is one minor hitch to this testing.  For many NLPs, there is a  * strict region of \a x where \a f(x) or \a c(x) are not defined.  In order to
  * help ensure that we stay out of these regions, variable bounds and a scalar
  * \c max_var_bounds_viol can be included so that the testing software
  * will never evaluate \a f(x) or \a c(x) outside the region:
  \verbatim 

   xl - max_var_bounds_viol <= x <= xu + max_var_bounds_viol
  \endverbatim
  * This is an important agreement made with the user.
  */
class NLPDirectTester {
public:

  /** \brief . */
  enum ETestingMethod {
    FD_COMPUTE_ALL
    ,FD_DIRECTIONAL
  };

  /** \brief . */
  STANDARD_COMPOSITION_MEMBERS( CalcFiniteDiffProd, calc_fd_prod );
  /// Members for option \c Gf_testing_method()
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ETestingMethod, Gf_testing_method );
  /// Members for option \c Gc_testing_method()
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ETestingMethod, Gc_testing_method );
  /// Members for option \c Gf_warning_tol()
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gf_warning_tol );
  /// Members for option \c Gf_error_tol()
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gf_error_tol );
  /// Members for option \c Gc_warning_tol()
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gc_warning_tol );
  /// Members for option \c Gc_error_tol()
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gc_error_tol );
  /// Members for option \c num_fd_directions()
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_fd_directions );
  /// Members for option \c dump_all()
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all );

  /// Constructor
  NLPDirectTester(
    const calc_fd_prod_ptr_t  &calc_fd_prod       = Teuchos::null
    ,ETestingMethod           Gf_testing_method   = FD_DIRECTIONAL
    ,ETestingMethod           Gc_testing_method   = FD_DIRECTIONAL
    ,value_type               Gf_warning_tol      = 1e-6
    ,value_type               Gf_error_tol        = 1e-1
    ,value_type               Gc_warning_tol      = 1e-6
    ,value_type               Gc_error_tol        = 1e-1
    ,size_type                num_fd_directions   = 1
    ,bool                     dump_all            = false
    );

  /** \brief This function takes an NLP object and its computed derivatives
   * and function values and validates
   * the functions and the derivatives by evaluating them
   * about the given point <tt>xo</tt>.
   * 
   * If all the checks as described in the
   * intro checkout then this function will return true, otherwise it
   * will return false.
   * 
   * If the finite difference steps are limited by relaxed variable
   * bounds then a warning message is printed and the derivatives
   * computed could be very inaccurate.
   *
   * @param  nlp     [in] %NLP object used to compute and test derivatives for.
   * @param  xo      [in] Point at which the derivatives are computed at.
   * @param  xl      [in] If != NULL then this is the lower variable bounds.
   * @param  xu      [in] If != NULL then this is the upper variable bounds.
   *	                If xl != NULL then xu != NULL must also be true
   *                 and visa-versa or a std::invalid_arguement exceptions
   *                 will be thrown.
   * @param  c       [in] Value of c(x) computed at xo.
   *                 If NULL, then none of the tests involving it will
   *                 be performed.
   * @param  Gf      [in] Gradient of f(x) computed at xo.
   *                 If NULL, then none of the tests involving it will
   *                 be performed.
   * @param  py      [in] Newton step <tt>py = -inv(C) * c(con_decomp)</tt>
   *                 If NULL, then none of the tests involving it will
   *                 be performed.
   * @param  rGf     [in] Reduced gradient of the objective function
   *                 </tt>rGf = Gf(var_indep) - D' * Gf(var_dep)</tt>.  If NULL,
   *                 then none of the tests involving it will be performed.
   * @param  GcU     [in]  Auxiliary jacobian matrix <tt>del(c(con_undecomp),x)</tt>.
   *                 If NULL, htne none of the tests involving it will be performed.
   * @param  D       [in] Direct sensitivity matrix <tt>D = -inv(C)*N</tt>.  If NULL,
   *                 none of the tests involving it will be performed.
   * @param  Uz      [in] <tt>Uz = F + E * D</tt>, which is the an auxiliary sensitivity matrix.
   *                 If NULL, then none of the tests involving it will be performed.
   * @param  print_all_warnings
   *                 [in] If true then all errors greater than warning_tol
   *                 will be printed if out!=NULL
   * @param  out     [in/out] If != null then some summary information is printed to it
   *                 and if a derivative does not match up then it prints which
   *                 derivative failed.  If <tt>out == 0</tt> then no output is printed.
   *
   * @return Returns <tt>true</tt> if all the derivatives comparisons are
   * within the error tolerances or returns false
   * otherwise.  This function will return false if any NaN or Inf values
   * where encountered.
   */
  bool finite_diff_check(
    NLPDirect         *nlp
    ,const Vector     &xo
    ,const Vector     *xl
    ,const Vector     *xu
    ,const Vector     *c
    ,const Vector     *Gf
    ,const Vector     *py
    ,const Vector     *rGf
    ,const MatrixOp   *GcU
    ,const MatrixOp   *D
    ,const MatrixOp   *Uz
    ,bool             print_all_warnings
    ,std::ostream     *out
    ) const;

};	// end class NLPDirectTester

}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_ORDER_DIRECT_TESTER_H
