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

#ifndef CALC_FINITE_DIFF_FIRST_DERIVATIVE_PRODUCT_H
#define CALC_FINITE_DIFF_FIRST_DERIVATIVE_PRODUCT_H

#include "NLPInterfacePack_Types.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace NLPInterfacePack {

/** \brief Strategy interface for computing the product of the derivatives of the functions
 * of an %NLP along given directions using finite differences.
  *
  * Specifically, this interface can be used to compute products with the finite
  * difference approximations of the gradient <tt>FDGf'*v</tt> and/or the Jacobian
  * <tt>FDGc'*v</tt> and/or the Jacobian <tt>FDGh'*v</tt>.  Doing so together may
  * take advantage of any shared calculations in the %NLP.  These products are
  * computed using finite differences.
  *
  * One of several different finite differencing schemes can be used.
  *
  * <ul>
  * <li> <b><tt>FD_ORDER_ONE</tt></b>: <tt>O(eps)</tt> one sided finite differences (two (one) evaluations of func).
         \verbatim

         grad(f(x),x)'*v \approx (1/eps) * ( -f(xo) + f(xo + eps * u) )\endverbatim
  * <li> <b><tt>FD_ORDER_TWO</tt></b>: <tt>O(eps^2)</tt> one sided finite differences (three (two) evaluations of func).
         \verbatim

         grad(f(x),x)'*v \approx (1/(2*eps)) * ( -3*f(xo) + 4*f(xo + eps*v) -f(xo + 2*eps*v))\endverbatim
  * <li> <b><tt>FD_ORDER_TWO_CENTRAL</tt></b>: <tt>O(eps^2)</tt> two sided central differences (two evaluations of func).
         \verbatim

         grad(f(x),x)'*v \approx (1/(2*eps)) * ( - f(xo - eps*v) + f(xo + eps*v) )\endverbatim
  * <li> <b><tt>FD_ORDER_TWO_AUTO</tt></b>: Use \c FD_ORDER_TWO_CENTRAL when not limited by bounds,
  *       otherwise use \c FD_ORDER_TWO.
  * <li> <b><tt>FD_ORDER_FOUR</tt></b>: <tt>O(eps^4)</tt> one sided finite differences (five (four) evaluations of func).
         \verbatim

         grad(f(x),x)'*v \approx (1/(12*eps)) * ( -25*f(xo) + 48*f(xo + eps*v) - 36*f(xo + 2*eps*v) + 16*f(xo + 3*eps*v)
                                                  - 3*f(xo + 4*eps*v))\endverbatim
  * <li> <b><tt>FD_ORDER_FOUR_CENTRAL</tt></b>: <tt>O(eps^4)</tt> two sided central differences (four evaluations of func).
         \verbatim

         grad(f(x),x)'*v \approx (1/(12*eps)) * ( f(xo - 2*eps*v) - 8*f(xo - eps*v) + 8*f(xo + eps*v) - f(xo + 2*eps*v))\endverbatim
  * <li> <b><tt>FD_ORDER_FOUR_AUTO</tt></b>: Use \c FD_ORDER_FOUR_CENTRAL when not limited by bounds,
  *       otherwise use \c FD_ORDER_FOUR.
  * </ul>
  *
  * The client can select the step sizes that are used to compute the finite differences.  First, the same step sizes \c eps
  * for all \a f(x), \c c(x) and \c h(x) can be selected using \c fd_step_size() with a positive value.  If \c fd_step_size()
  * returns a negative number, then a default value will be determined internally that is appropriate for the chosen
  * method (see \c fd_method_order()).  Whatever step size represented by \c fd_step_size() (or a default) will be scaled
  * by <tt>||xo||inf + 1.0</tt> if <tt>fd_step_select() == FD_STEP_ABSOLUTE</tt>.
  * Using the same step sizes for \a f(x), \a c(x) and \a h(x) an advantage
  * in that the %NLP implementation may be able exploit shared computations.  However, for some applications, it may be
  * advantagous to use different step lengths for \c f(x), \c c(x) and \c h(x).  These individual step lengths can be
  * set using \c fd_step_size_f(), \c fd_step_size_c() and \c fd_step_size_h() respectively.  These step lengths will also
  * be scaled the same as for \c fd_step_size() if <tt>fd_step_select() == FD_STEP_ABSOLUTE</tt>.  If any of these
  * individual step size functions returns a negative number, then the value for \c fd_step_size() will be used.
  *
  * The finite difference perturbations may be limited by the relaxed variable bounds:
  \verbatim

  xl - max_var_bounds_viol <= x <= xu + max_var_bounds_viol
  \endverbatim
  * if variable bounds are present.  If it is not expected that bounds will limit the step
  * length (or there are no bounds) then the central difference methods \c FD_ORDER_TWO_CENTRAL
  * and \c FD_ORDER_FOUR_CENTRAL should be preferred since they are more accurate.  The method
  * \c FD_ORDER_FOUR_AUTO will use fourth order central differences (\c FD_ORDER_FOUR_CENTRAL)
  * unless the variable bounds limit the minimum step size in which case fourth order
  * one sided differences (\c FD_ORDER_FOUR) will be used.  With one sided differences, the
  * implementation may be able to take a larger (in magnitude) negative step than a positive
  * step (and visa-versa) in which case the one sided methods would have an advantage.  Note
  * that the one situation where one sided differences is guaranteed to be able to take steps
  * away from the bounds is when <tt>xl - max_var_bounds_viol + fd_step_size <= xu</tt> and <tt>v = eta(j)</tt>
  * (i.e. <tt>eta(j)</tt> is the jth column of identity).  The situation <tt>v = eta(j)</tt> occurs
  * when the client is computing the full finite difference approximations to \c Gf, \c Gc and/or \c Gh
  * on variable at a time.
  */
class CalcFiniteDiffProd {
public:

  /** \brief . */
  enum EFDMethodOrder {
    FD_ORDER_ONE           ///< Use O(eps) one sided finite differences (cramped bounds)
    ,FD_ORDER_TWO          ///< Use O(eps^2) one sided finite differences (cramped bounds)
    ,FD_ORDER_TWO_CENTRAL  ///< Use O(eps^2) two sided central finite differences
    ,FD_ORDER_TWO_AUTO     ///< Use FD_ORDER_TWO_CENTRAL when not limited by bounds, otherwise use FD_ORDER_TWO
    ,FD_ORDER_FOUR         ///< Use O(eps^4) one sided finite differences (cramped bounds)
    ,FD_ORDER_FOUR_CENTRAL ///< Use O(eps^4) two sided central finite differences
    ,FD_ORDER_FOUR_AUTO    ///< Use FD_ORDER_FOUR_CENTRAL when not limited by bounds, otherwise use FD_ORDER_FOUR
  };
  /** \brief . */
  enum EFDStepSelect {
    FD_STEP_ABSOLUTE      ///< Use absolute step size <tt>fd_step_size</tt>
    ,FD_STEP_RELATIVE     ///< Use relative step size <tt>fd_step_size * ||xo||inf</tt>
  };
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDMethodOrder, fd_method_order );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDStepSelect, fd_step_select );
  /** \brief Pick the size of the finite difference step.
   *
   * If <tt>fd_step_size < 0</tt> then the implementation will try to
   * select it based on the order of method <tt>fd_method_order()</tt>
   * that is selected.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, fd_step_size );
  /** \brief Pick the minimum step size under which the finite difference product
   * will not be computed.
   *
   * If <tt>fd_step_size_min == 0</tt> then the finite difference computation
   * will always be performed.  If <tt>fd_step_size_min < 0</tt> then the
   * minimum step size will be determined internally.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, fd_step_size_min );
  /// Set the step size for \a f(x)
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, fd_step_size_f );
  /// Set the step size for \a c(x)
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, fd_step_size_c );

  /** \brief . */
  CalcFiniteDiffProd(
    EFDMethodOrder              fd_method_order  = FD_ORDER_FOUR_AUTO
    ,EFDStepSelect              fd_step_select   = FD_STEP_ABSOLUTE
    ,value_type                 fd_step_size     = -1.0
    ,value_type                 fd_step_size_min = -1.0
    ,value_type                 fd_step_size_f   = -1.0
    ,value_type                 fd_step_size_c   = -1.0
    );

  /** \brief . */
  virtual ~CalcFiniteDiffProd() {}

  /** \brief Compute the directional derivatives by finite differences.
   *
   * The computation may fail if \c NaN or \c Inf is encountered durring any
   * of the computations in which case a \c NaNInfException exception will be
   * thrown.  Otherwise the computation should be completed successfully.
   * 
   * The finite difference peturbations may be limited by the relaxed variable bounds
   \verbatim

   xl - max_var_bounds_viol <= x <= xu + max_var_bounds_viol
   \endverbatim
   * if variable bounds are present.  If these bounds do limit the finite difference step
   * size then a warning will be printed to *out (if <tt>out!=NULL</tt>) and the
   * derivatives may be very inaccurate.  If bounds do limit the steps taken then it is
   * advisable to use the one sided finite difference (FD_ORDER_ONE) and this implementation
   * can move away from the bounds.
   * 
   * @param  xo       [in] Base point for the unknown variables to compute derivatives at.
   * @param  xl       [in] If <tt>!= NULL</tt> then this is the lower variable bounds.
   * @param  xu       [in] If <tt>!= NULL</tt> then this is the upper variable bounds.
   *                  If <tt>xl != NULL</tt> then <tt>xu != NULL</tt> must also be true
   *                  and visa-versa or a <tt>std::invalid_arguement</tt> exception
   *                  will be thrown.
   * @param  v        [in] The vector for which to form the products with.
   * @param  fo       [in] If <tt>fo != NULL</tt> then <tt>*fo</tt> should be set to the
   *                  value of <tt>f(xo)</tt>.  Not useful for \c FD_ORDER_TWO_CENTRAL.
   * @param  co       [in] If <tt>co != NULL</tt> then <tt>*co</tt> should be set to the
   *                  value of <tt>c(xo)</tt>.  Not useful for \c FD_ORDER_TWO_CENTRAL.
   * @param  check_nan_inf
   *                  [in] If \c true, the the computed values will be checked for nan and inf.
   * @param  nlp      [in] Used to compute \a f(x), \a c(x) and \a h(x).  The current set
   *                  references to <tt>nlp->get_f()</tt>, <tt>nlp->get_c()</tt>
   *                  and <tt>nlp->get_h()</tt> will be preserved on output.  The
   *                  %NLP must be initialized before input.
   * @param  Gf_prod	[out] If <tt>!= NULL</tt>, the will contain the finite difference
   *                  computed product <tt>Gf'*v</tt> on output.  If <tt>== NULL</tt>,
   *                  then \a f(x) is not specifically computed.
   * @param  Gc_prod	[out] If <tt>!= NULL</tt>, the will contain the finite
   *                  difference computed product <tt>Gc'*v</tt>.  If <tt>== NULL</tt>,
   *                  then \a c(x) is not specifically computed.
   * @param  out      [in/out] If <tt>out != NULL</tt> then any waring or error messages are
   *                  output here.
   *
   * @return Returns <tt>true if the finite difference computations were performed.</tt>
   * Returns <tt>false</tt> if the finite difference computations were not performed
   * because the variables bounds limited the step length too much.
   *
   * ToDo: Discuss options!
   */
  virtual bool calc_deriv_product(
    const Vector       &xo
    ,const Vector      *xl
    ,const Vector      *xu
    ,const Vector      &v
    ,const value_type  *fo
    ,const Vector      *co
    ,bool              check_nan_inf
    ,NLP               *nlp
    ,value_type        *Gf_prod
    ,VectorMutable     *Gc_prod
    ,std::ostream      *out
    ,bool              trace    = false
    ,bool              dump_all = false
    ) const;

};	// end class CalcFiniteDiffProd

}	// end namespace NLPInterfacePack

#endif	// CALC_FINITE_DIFF_FIRST_DERIVATIVE_PRODUCT_H
