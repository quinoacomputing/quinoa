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

#ifndef NLP_OBJ_GRADIENT_H
#define NLP_OBJ_GRADIENT_H

#include "NLPInterfacePack_NLP.hpp"

namespace NLPInterfacePack {
/** \brief %NLP interface class that adds gradient information for the objective function {abstract}.
 *
 * <b>Overview:</b>
 *
 * This class adds the ability to compute the gradient of the objective function
 * \c Gf(x) to the basic information given in the \c NLP interface class.  Note that
 * \c Gf is in the vector space \c space_x().
 *
 * <b>Client Usage:</b>
 *
 * As with the <tt>NLP</tt> base interface, the <tt>initialize()</tt> method must be called before
 * the %NLP object can be used.   The method <tt>set_Gf()</tt> is used to set a pointer to a vector
 * to update when the gradient of the objective \c Gf is computed when <tt>calc_Gf()</tt> is called.
 *
 * The number of evaluations of \c Gf using <tt>calc_Gf()</tt> is returned by <tt>num_Gf_evals()</tt>.
 * 
 * <b>Subclass developer's notes:</b>
 *
 * <A NAME="must_override"></A>
 * In addition to the methods that must be overridden by the <tt>NLP</tt> interface
 * (<A HREF="classNLPInterfacePack_1_1NLP.html#must_override">see</A>) the following methods
 * must also be overridden: <tt>imp_calc_Gf()</tt>.
 *
 * <A NAME="should_override"></A>
 * In addition to the methods that should be overridden from <tt>%NLP</tt> by most subclasses
 * (<A HREF="classNLPInterfacePack_1_1NLP.html#should_override">see</A>), the following
 * additional methods should be overridden: \c initialize().
 *
 * The following methods should never have to be overridden by most subclasses except in some very
 * strange situations: \c set_Gf(), \c get_Gf(), \c Gf(), \c num_Gf_evals().
 */
class NLPObjGrad : virtual public NLP {
public:

  /** @name Constructors */
  //@{

  /// Initialize to no reference set to calculation quanities
  NLPObjGrad();

  //@}

  /** @name NLP initialization */
  //@{

  /** \brief Initialize the NLP for its first use.
    *
    * This function implementation should be called by subclass implementations
    * in order to reset counts for \c f(x), \c c(x), \c h(x) and \c Gf(x) evaluations.
    * This implementation calls <tt>this->NLP::initialize()</tt>
    *
    * Postconditions:<ul>
    * <li> See <tt>NLP::initialize()</tt>
    * <li> <tt>this->num_Gf_evals() == 0</tt>
    * </ul>
    */
  void initialize(bool test_setup);

  //@}

  /** @name Information */
  //@{

  /** \brief Determine if the objective gradient is supported or not.
   *
   * The default implementation returns <tt>true</tt>.
   */
  virtual bool supports_Gf() const;

  /** \brief Determine if the objective gradient product is supported or not.
   *
   * The default implementation returns <tt>true</tt>.
   */
  virtual bool supports_Gf_prod() const;

  //@}

  /** @name <<std aggr>> members for the gradient of the objective function Gf(x) */
  //@{

  /** \brief Set a pointer to a vector to be updated when <tt>this->calc_Gf()</tt> is called.
   *
   * @param  Gf  [in] Pointer to gradient vector.  May be \c NULL.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->supports_Gf()</tt>
   * <li> [<tt>Gf != NULL</tt>] <tt>Gf->space().is_compatible(*this->space_x()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get_Gf() == Gf</tt>
   * </ul>
   */
  virtual void set_Gf(VectorMutable* Gf);
  /** \brief Return pointer passed to <tt>this->set_Gf()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->supports_Gf()</tt>
   * </ul>
   */
  virtual VectorMutable* get_Gf();
  /** \brief Returns non-<tt>const</tt> <tt>*this->get_Gf()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->supports_Gf()</tt>
   * <li> <tt>this->get_Gf() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual VectorMutable& Gf();
  /** \brief Returns <tt>const</tt> <tt>*this->get_Gf()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->supports_Gf()</tt>
   * <li> <tt>this->get_Gf() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual const Vector& Gf() const;

  //@}

  /** @name Unset calculation quantities */
  //@{
  
  /** \brief Call to unset all storage quantities (both in this class and all subclasses).
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> See <tt>NLP::unset_quantities()</tt>
   * <li> <tt>this->get_Gf() == NULL</tt>
   * </ul>
   *
   * This method must be called by all subclasses that override it.
   */
  void unset_quantities();

  //@}

  /** @name Calculation Members */
  //@{

  /** \brief Update the vector for \c Gf at the point \c x and put it in the stored reference.
   *
   * @param  x     [in] Point at which to calculate the gradient of the objective <tt>Gf(x)</tt>.
   * @param  newx  [in] (default \c true) If \c false, the values in \c x are assumed to be the same as
   *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
   *               If \c true, the values in \c x are assumed to not be the same as the last call to a
   *               <tt>this->calc_*(x,newx)</tt> member.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->supports_Gf()</tt>
   * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>this->get_Gf() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->Gf()</tt> is updated to \c Gf(x)
   * </ul>
   *
   * If <tt>set_multi_calc(true)</tt> was called then referenced storage for \c f and/or \c c
   * may also be updated but are not guaranteed to be.  But no other quanities from possible subclasses are allowed
   * to be updated as a side effect (i.e. no higher order derivatives).
   */ 
  virtual void calc_Gf(const Vector& x, bool newx = true) const;

  /** \brief Calculate the inner product <tt>Gf(x)'*d</tt> at the point <tt>x</tt> and put it in the stored reference.
   *
   * @param  x     [in] Base point
   * @param  d     [in] Direction to compute the product along.
   * @param  newx  [in] (default \c true) If \c false, the values in \c x are assumed to be the same as
   *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
   *               If \c true, the values in \c x are assumed to not be the same as the last call to a
   *               <tt>this->calc_*(x,newx)</tt> member.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->supports_Gf()</tt>
   * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>return</tt> gives the desired product.
   * </ul>
   *
   * If <tt>set_multi_calc(true)</tt> was called then referenced storage for \c f and/or \c c
   * may also be updated but are not guaranteed to be.  But no other quanities from possible subclasses are allowed
   * to be updated as a side effect (i.e. no higher order derivatives).
   */ 
  virtual value_type calc_Gf_prod(const Vector& x, const Vector& d, bool newx = true) const;

  //@}

  //@}

  /** @name Function evaluation counts. */
  //@{

  /** \brief Objective gradient evaluations count.
   *
   * This function can be called to find out how many evaluations
   * \c this->calc_Gf() the client requested since \c this->initialize() was called.
   */
  virtual size_type num_Gf_evals() const;
  
  //@}

  /** @name Protected types */
  //@{
  
  /** \brief Struct for gradient (objective), objective and constriants (pointers)
   */
  struct ObjGradInfo {
    /** \brief . */
    ObjGradInfo()
      : Gf(NULL), f(NULL), c(NULL)
    {}
    /** \brief . */
    ObjGradInfo( VectorMutable* Gf_in, const ZeroOrderInfo& first_order_info_in )
      : Gf(Gf_in), f(first_order_info_in.f), c(first_order_info_in.c)
    {}
    /// Pointer to gradient of objective function <tt>Gf</tt> (may be NULL if not set)
    VectorMutable*       Gf;
    /// Pointer to objective function <tt>f</tt> (may be NULL if not set)
    value_type*          f;
    /// Pointer to constraints residual <tt>c</tt> (may be NULL if not set)
    VectorMutable*       c;
  }; // end struct ObjGradInfo

  //@}

protected:

  /// Return objective gradient and zero order information.
  const ObjGradInfo obj_grad_info() const;

  /** @name Protected methods to be overridden by subclasses */
  //@{

  /** \brief Overridden to compute f(x) and perhaps c(x) (if multiple calculaiton = true).
   *
   * Preconditions:<ul>
   * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
   * <li> <tt>obj_grad_info.Gf != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*obj_grad_info.Gf</tt> is updated to \a Gf(x).
   * </ul>
   *
   * @param x       [in]  Unknown vector (size n).
   * @param  newx   [in] (default \c true) If \c false, the values in \c x are assumed to be the same as
   *                the last call to a <tt>this->imp_calc_*(x,newx)</tt> member.
   *                If \c true, the values in \c x are assumed to not be the same as the last call to a
   *                <tt>this->imp_calc_*(x,newx)</tt> member.
   * @param obj_grad_info
   *                [out] Pointers to \c f, \c c and \c Gf.
   *                On output <tt>*obj_grad_info.Gf</tt> is updated to \a Gf(x).
   *                Any of the other objects pointed to in
   *                \c obj_grad_info may be set if <tt>this->multi_calc() == true</tt> but are
   *                now guaranteed to be.
   */
  virtual void imp_calc_Gf(const Vector& x, bool newx, const ObjGradInfo& obj_grad_info) const = 0;

  //@}

private:

  mutable VectorMutable     *Gf_;
  mutable size_type         num_Gf_evals_;

};	// end class NLPObjGrad

// //////////////////
// Inline members

inline
const NLPObjGrad::ObjGradInfo NLPObjGrad::obj_grad_info() const
{
  return ObjGradInfo(Gf_,zero_order_info());
}

}	// end namespace NLPInterfacePack 

#endif // NLP_OBJ_GRADIENT_H
