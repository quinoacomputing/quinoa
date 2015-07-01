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

#ifndef NLP_SECOND_ORDER_INFO_H
#define NLP_SECOND_ORDER_INFO_H

#include "NLPInterfacePack_NLPFirstOrder.hpp"

namespace NLPInterfacePack {

/** \brief NLP second order information interface class {abstract}.
 *
 * <b>Overview:</b>
 *
 * This class adds second order inforamtion to the first order information
 * and basic information given in the <tt>NLPFirstOrder</tt> and base interfaces.
 *
 * Specifically the Hesssian of the Lagrangian is defined as:
 \verbatim

 HL = Hf + sum( Hc(j) * lambda(j), j = 1...m )
 \endverbatim
 * Where: <ul>
 * <li> \c Hf is the hessian of the objective function \a f(x)
 * <li> \c Hc(j) is the hessian of the \c jth equality constriant <i>c<sub>j</sub>(x)</i>
 * <li> \c lambda is the vector of lagrange multipliers for the equality
 *      constraints \a c(x) 
 * </ul>
 *
 * <b>Client Usage:</b>
 *
 * ToDo: Finish Documentation!
 * 
 * <b>Subclass developer's notes:</b>
 *
 * ToDo: Finish Documentation!
 *
 */
class NLPSecondOrder : virtual public NLPFirstOrder {
public:

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixSymOp> >    mat_sym_fcty_ptr_t;

  /** @name Constructors */
  //@{

  /// Initialize to no reference set to calculation quanities
  NLPSecondOrder();

  //@}

  /** @name NLP initialization */
  //@{
  
  /** \brief Initialize the NLP for its first use.
   *
   * This function implementation should be called by subclass implementations
   * in order to reset counts for \c f(x), \c c(x), \c h(x), \c Gf(x), \c Gc(x),
   * \c Gh(x) and \c HL(x) evaluations.  This implementation calls
   * <tt>this->NLPFirstOrder::initialize()</tt>
    *
   * Postconditions:<ul>
   * <li> See <tt>NLPFirstOrder::initialize()</tt>
   * <li> <tt>this->num_HL_evals() == 0</tt>
   * </ul>
   */
  void initialize(bool test_setup);

  //@}

  /** @name Matrix factory objects */
  //@{

  /** \brief Return a matrix factory object for creating <tt>HL</tt>.
   *
   * The returned matrix object may not support the creation of any
   * sub-matrix spaces (i.e. <tt>return->sub_space(rrng,crng).get() == NULL</tt>
   * for all <tt>rrng</tt> and <tt>crng</tt>).
   */
  virtual const mat_sym_fcty_ptr_t factory_HL() const = 0;

  //@}

  /** @name <<std aggr>> members for the Hessian of the Lagrangian HL */
  //@{

  /** \brief Set a pointer to a matrix object to be updated when <tt>this->calc_HL()</tt> is called.
   *
   * @param  HL  [in] Pointer to Hessian of the Lagrangian matrix.  May be \c NULL.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get_HL() == HL</tt>
   * </ul>
   */
  virtual void set_HL(MatrixSymOp* HL);
  /** \brief Return pointer passed to <tt>this->set_HL()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   */
  virtual MatrixSymOp* get_HL();
  /** \brief Returns non-<tt>const</tt> <tt>*this->get_HL()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_HL() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual MatrixSymOp& HL();
  /** \brief Returns <tt>const</tt> <tt>*this->get_HL()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_HL() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual const MatrixSymOp& HL() const;

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
   * <li> See <tt>NLPFirstOrder::unset_quantities()</tt>
   * <li> <tt>this->get_HL() == NULL</tt>
   * </ul>
   *
   * This method must be called by all subclasses that override it.
   */
  void unset_quantities();

  //@}

  /** @name Calculation Members */
  //@{

  /** \brief Update the matrix for <tt>HL</tt> at the point <tt>x</tt>, <tt>lambda</tt>,
   * <tt>lambdaI</tt> and put it in the stored reference.
   *
   * The referenced storage for <tt>f</tt>, <tt>c</tt>, <tt>Gf</tt> and <tt>Gc</tt>
   * may also be changed but are not guarentied to be.
   * But no other quanities from possible subclasses are allowed to be updated as a side effect.
   *
   * @param  x        [in] Unknown primal variables
   * @param  lambda   [in] Lagrange muitipliers for equality constriants.
   *                  If <tt>m() == 0</tt> then <tt>lambda</tt> must be <tt>NULL</tt>.  However, if
   *                  <tt>m() > 0</tt> then <tt>lambda == NULL</tt> is still allowed and is treated
   *                  as <tt>lambda = 0</tt>.
   * @param  newpoint [in] (default \c true) If \c false, the values in \c x, \c lambda and \c lambdaI
   *                  are the same as the last call to <tt>this->calc_HL()</tt>.
   *                  If \c true, then this is a new point.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>this->get_HL() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * <li> [<tt>this->m() == 0</tt>] <tt>lambda == NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>this->m() != 0 && lambda != 0</tt>] <tt>lambda->space().is_compatible(*this->space_c()) == true)</tt>
   *      (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->HL()</tt> is updated to \c HL(x)
   * </ul>
   */ 
  virtual void calc_HL(
    const Vector& x, const Vector* lambda, bool newpoint = true) const;

  //@}

  /** @name Number of function evaluations */
  //@{

  /** \brief Number of Hessian evaluations.
    *
    * This function can be called to find out how many evaluations
    * the client requested since \c initialize() was called.
    */
  virtual size_type num_HL_evals() const;

  //@}

protected:

  /** \brief Struct for zero, first and second order quantities (pointers)
   */
  struct SecondOrderInfo {
    /** \brief . */
    SecondOrderInfo()
      : HL(NULL), Gc(NULL), Gf(NULL), f(NULL), c(NULL)
      {}
    /** \brief . */
    SecondOrderInfo( MatrixSymOp* HL_in, const FirstOrderInfo& first_order_info )
      :HL(HL_in), Gc(first_order_info.Gc), Gf(first_order_info.Gf)
      ,f(first_order_info.f), c(first_order_info.c)
      {}
    /// Pointer to Hessiand of the Lagrangian <tt>HL</tt>) (may be NULL is not set)
    MatrixSymOp*        HL;
    /// Pointer to Hessian of the equality constraints <tt>Gc</tt> (may be NULL if not set)
    MatrixOp*           Gc;
    /// Pointer to gradient of objective function <tt>Gf</tt> (may be NULL if not set)
    VectorMutable*      Gf;
    /// Pointer to objective function <tt>f</tt> (may be NULL if not set)
    value_type*         f;
    /// Pointer to equality constraints residule <tt>c</tt> (may be NULL if not set)
    VectorMutable*      c;
  }; // end struct SecondOrderInfo

  /// Return objective gradient and zero order information.
  const SecondOrderInfo second_order_info() const;

  /** @name Protected methods to be overridden by subclasses */
  //@{

  /** \brief Overridden to compute <tt>Gc(x)</tt> and perhaps <tt>Gf(x)</tt>, <tt>f(x)</tt> and <tt>c(x)</tt>.
   *
   * @param x                     [in] Unknown vector (size n).
   * @param lambda                [in] Lagrange multipliers for equality constraints c(x).
   *                              Must be <tt>NULL</tt> if <tt>m() == 0</tt>.  If \c NULL, then
   *                              treated as <tt>lambda = 0</tt>.
   * @param  newpoint             [in] (default \c true) If \c false, the values in \c x, \c lambda and \c lambdaI
   *                              are the same as the last call to <tt>this->calc_HL()</tt>.
   *                              If \c true, then this is a new point.
   * @param second_order_info     [out] Pointers to \c HL, \c Gc, \c Gh, \c Gf, \c f, \c c and \c h
   *                              On output <tt>*second_order_info.HL</tt> is updated to \a HL(x).
   *                              Any of the other objects pointed to in \c second_order_info may
   *                              also be updated but are not guaranteed to be.
   *
   * Preconditions:<ul>
   * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
   * <li> <tt>second_order_info.HL != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*second_order_info.HL</tt> is updated to \a HL(x).
   * </ul>
   */
  virtual void imp_calc_HL(
    const Vector& x, const Vector* lambda, bool newpoint
    ,const SecondOrderInfo& second_order_info
    ) const = 0;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  Teuchos::AbstractFactory<AbstractLinAlgPack::MatrixSymOp>  *factory_HL;
#endif
  mutable MatrixSymOp   *HL_;
  mutable bool          num_HL_evals_;

};	// end class NLPSecondOrder

// //////////////////
// Inline members

inline
const NLPSecondOrder::SecondOrderInfo NLPSecondOrder::second_order_info() const
{
  return SecondOrderInfo(HL_,first_order_info());
}

}	// end namespace NLPInterfacePack 

#endif // NLP_SECOND_ORDER_INFO_H
