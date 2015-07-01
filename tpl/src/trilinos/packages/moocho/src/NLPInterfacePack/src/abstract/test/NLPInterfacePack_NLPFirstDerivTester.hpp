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

#ifndef NLP_FIRST_DERIVATIVES_TESTER_H
#define NLP_FIRST_DERIVATIVES_TESTER_H

#include <iosfwd>

#include "NLPInterfacePack_Types.hpp"
#include "NLPInterfacePack_CalcFiniteDiffProd.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace NLPInterfacePack {

/** \brief Concrete class that tests the derivatives using finite differences.
 *
 * There are two options for testing the derivatives by finite differences.
 * 
 * The first option (<tt>fd_testing_method==FD_COMPUTE_ALL</tt>) is to compute all of
 * them as dense vectors and matrices. This option can be very expensive in runtime
 * and storage costs.  The amount of storage space needed is <tt>O(n*m)</tt> and
 * \c f(x) and \c c(x) will be computed <tt>O(n)</tt> times.
 * 
 * The other option (<tt>fd_testing_method==FD_DIRECTIONAL</tt>)
 * computes products of the form <tt>g'*v</tt> and compares them to
 * the finite difference computed value <tt>g_fd'*v</tt>.  This method
 * only costs <tt>O(n)</tt> storage and two function evaluations per
 * direction (assuming central differences are used.  The directions
 * <tt>v</tt> are computed randomly between <tt>[-1,+1]</tt> so that
 * they are well scaled and should give good results.  The option
 * <tt>num_fd_directions()</tt> determines how many random directions
 * are used.  A value of <tt>num_fd_directions() <= 0</tt> means that
 * a single finite difference direction of <tt>1.0</tt> will be used
 * for the test.
 *
 * This class computes the derivatives using a
 * <tt>CalcFiniteDiffProd</tt> object can can use up to fourth-order
 * (central) finite differences but can use as low as first-order
 * one-sided differences.
 *
 * The client can set the tolerances used to measure if the anylitical
 * values of \c Gf and \c Gc are close enough to the finite difference
 * values.  Let the function \a h(x) be <tt>f(x)</tt> or any
 * <tt>cj(x)</tt>, for <tt>j = 1...m</tt>.  Let <tt>gh(i) =
 * d(h(x))/d(x(i))</tt> and <tt>fdh(i) =
 * finite_diff(h(x))/d(x(i))</tt>.  Then let's define the relative
 * error between the anylitic value and the finite difference value to
 * be:

 \verbatim

     err(i) = |(gh(i) - fdh(i))| /  (||gh||inf + ||fdh||inf + (epsilon)^(1/4))
 \endverbatim

 * The above error takes into account the relative sizes of the
 * elements and also allows one or both of the elements to be zero
 * without ending up with <tt>0/0</tt> or something like
 * <tt>1e-16</tt> not comparing with zero.
 *
 * All errors <tt>err(i) >= warning_tol</tt> are reported to <tt>*out</tt> if
 * <tt>out != NULL</tt> and <tt>print_all_warnings==true</tt>.  Otherwise, if
 * <tt>out != NULL</tt>, only the number of elements and the maxinum violation of the
 * warning tolerance will be printed.  The first error <tt>err(i) >= error_tol</tt>
 * that is found is reported is reported to <tt>*out</tt> if <tt>out != NULL</tt> and
 * immediatly \c finite_diff_check() returns \c false.  If all errors
 * <tt>err(i) < error_tol</tt> then \c finite_diff_check() will return \c true.
 *
 * Given these two tolerances the client can do many things:
 * <ol>
 * <li> Print out all the comparisons that are not equal by setting warning_tol
 *    == 0.0 and error_tol = very_large_number.
 *
 * <li> Print out all suspect comparisons by setting epsilon < warning_tol < 1
 *    and error_tol = very_large_number.
 *
 * <li> Just validate that matrices are approximatly equal and report the first
 *    discrepency if not by setting epsilon < error_tol < 1 and warning_tol
 *    >= error_tol.
 *
 * <li> Print out any suspect comparisons by setting epsilon < warning_tol < 1
 *    but also quit if the error is too large by setting error_tol > warning_tol.
 * </ol>
 * There is one minor hitch to this testing.  For many NLPs, there is a
 * strict region of \a x where \a f(x) or \a c(x) are not defined.  In order to
 * help ensure that we stay out of these regions, variable bounds can be
 * included and a scalar \c max_var_bounds_viol so that the testing software
 * will never evaluate \a f(x) or \a c(x) outside the region:
 \verbatim
   
     xl - max_var_bounds_viol <= x <= xu + max_var_bounds_viol
 \endverbatim
 * This is an important agreement made with the user.
 */
class NLPFirstDerivTester {
public:

  /** \brief . */
  enum ETestingMethod {
    FD_COMPUTE_ALL
    ,FD_DIRECTIONAL
  };

  /** \brief . */
  STANDARD_COMPOSITION_MEMBERS( CalcFiniteDiffProd, calc_fd_prod );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ETestingMethod, fd_testing_method );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_fd_directions );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol );

  /// Constructor
  NLPFirstDerivTester(
    const calc_fd_prod_ptr_t  &calc_fd_prod      = Teuchos::rcp(new CalcFiniteDiffProd())
    ,ETestingMethod           fd_testing_method  = FD_DIRECTIONAL
    ,size_type                num_fd_directions  = 1
    ,value_type               warning_tol        = 1e-8
    ,value_type               error_tol          = 1e-3
    );

  /** \brief This function takes an NLP object and its computed derivatives
   * and function values and validates
   * the functions and the derivatives by evaluating them
   * about the given point <tt>x</tt>.  If all the checks as described in the
   * intro checkout then this function will return true, otherwise it
   * will return false.
   *
   * @param  nlp      [in] NLP object used to compute and test derivatives for.
   * @param  xo       [in] Point at which the derivatives are computed at.
   * @param  xl       [in] If != NULL then this is the lower variable bounds.
   * @param  xu       [in] If != NULL then this is the upper variable bounds.
   *                  If xl != NULL then xu != NULL must also be true
   *                  and visa-versa or a std::invalid_arguement exceptions
   *                  will be thrown.
   * @param  Gc       [in] A matrix object for the Gc computed at xo.
   *                  If Gc==NULL then this is not tested for.
   * @param  Gf       [in] Gradient of f(x) computed at xo.
   *                  If Gf==NULL then this is not tested for.
   * @param  print_all_warnings
   *                  [in] If true then all errors greater than warning_tol
   *                  will be printed if out!=NULL
   * @param  out      [in/out] If != null then some summary information is printed to it
   *                  and if a derivative does not match up then it prints which
   *                  derivative failed.  If <tt>out == 0</tt> then no output is printed.
   *
   * @return Returns <tt>true</tt> if all the derivatives check out, and false
   * otherwise.
   */
  bool finite_diff_check(
    NLP               *nlp
    ,const Vector     &xo
    ,const Vector     *xl
    ,const Vector     *xu
    ,const MatrixOp   *Gc
    ,const Vector     *Gf
    ,bool             print_all_warnings
    ,std::ostream     *out
    ) const;

private:

  /** \brief . */
  bool fd_check_all(
    NLP               *nlp
    ,const Vector     &xo
    ,const Vector     *xl
    ,const Vector     *xu
    ,const MatrixOp   *Gc
    ,const Vector     *Gf
    ,bool             print_all_warnings
    ,std::ostream     *out
    ) const;

  /** \brief . */
  bool fd_directional_check(
    NLP               *nlp
    ,const Vector     &xo
    ,const Vector     *xl
    ,const Vector     *xu
    ,const MatrixOp   *Gc
    ,const Vector     *Gf
    ,bool             print_all_warnings
    ,std::ostream     *out
    ) const;
};

}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_DERIVATIVES_TESTER_H
