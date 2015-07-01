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

#ifndef NLP_FIRST_ORDER_DIRECT_H
#define NLP_FIRST_ORDER_DIRECT_H

#include "NLPInterfacePack_NLPObjGrad.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace NLPInterfacePack {

/** \brief Interface providing only direct first order sensitivity information.
 *
 * <b>Overview:</b>
 *
 * This interface defines a basis for the equality constriants and then only
 * certain linear systems with this basis are solved for.  This interface is
 * useful in reduced space SQP-type and other related optimization algorithms.
 *
 * Specifically, the variables are partitioned into dependent and independent
 * sets <tt>x = [ x_dep' x_indep' ]'</tt> and Jacobians of the constraints
 * <tt>c(x)</tt> at the point <tt>x</tt> are:

 \verbatim

  del(c,x) = Gc' = [ del(c(con_decomp))   ] = [ GcD' ] = [ GcDD'  GcDI' ] = [ C  N ]
                   [ del(c(con_undecomp)) ]   [ GcU' ]   [ GcUD'  GcUI' ]   [ E  F ]

    where:
      C <: R^(r x r) is nonsingular
      N <: R^(r x (n-r))
      E <: R^((m-r) x r)
      F <: R^((m-r) x (n-r))
 \endverbatim

 * This partitions the general equality constraints c(x) into two sets;
 * decomposed c(con_decomp) and undecomposed c(con_undecomp).  It is therefore
 * expected that sub-vectors and subspaces from
 * <tt>space_x().sub_space(var_dep)</tt>,
 * <tt>space_x().sub_space(var_indep)</tt>,
 * <tt>space_c().sub_space(con_decomp)</tt> and
 * <tt>space_c().sub_space(con_undecomp)</tt> can all be accessed.  Other
 * sub-vectors and sub-spaces may not be available (but the algorithm should
 * not need access to other sub-spaces).
 *
 * Free access to solves with the basis <tt>C</tt> is not given however and instead this interface
 * computes, for the current point \a x, the direct sensitivity matrice <tt>D = -inv(C)*N</tt>,
 * the auxiliary matrices <tt>Uz = F + E * D</tt> and <tt>GcU = [ GcUD; GcUI ] = [ E';  F' ]</tt>,
 * and the Newton step <tt>py = -inv(C)*c(con_decomp)</tt>.
 * In general, linear solves with the transpose with <tt>C</tt> are not possible and
 * therefore are not avalible.  A number of very specialized applications can only
 * provide this information but this is all that is needed by many numerical
 * optimization (and related) algorithms.
 *
 * <b>Client Usage:</b>
 *
 * The dimension of the basis matrix \c C is returned by \c r().  The ranges for the dependent and
 * independent varaibles are returned by \c var_dep() and \c var_indep().  The ranges for the
 * decomposed and undecomposed equality constraints are \c con_decomp() and \c con_undecomp().
 * Note that \c con_undecomp() will return an invalid range if there are no undecomposed equalities.
 *
 * Note that the matrix objects returned from \c factory_GcU(), \c factory_D()
 * and \c factory_Uz() can not be expected to be usable until they are
 * passed to the calculation routines or have been intialized in some other way.
 *
 * <b>Subclass Developer's Notes:</b>
 *
 * The default implementation of this interface assumes that there are no undecomposed
 * equality constraints (i.e. <tt>this->con_decomp().size() == this->m()).
 *
 * ToDo: Finish Documentation!
 */
class NLPDirect : virtual public NLPObjGrad
{
public:

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixOp> >               mat_fcty_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixSymOp> >            mat_sym_fcty_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixSymOpNonsing> > mat_sym_nonsing_fcty_ptr_t;

  /** \brief Initialize the factory objects for the special matrices for <tt>D'*D</tt> and <tt>S = I + D'*D</tt>.
   *
   * Postconditions:<ul>
   * <li>this->factory_transDtD().get() == factory_transDtD.get()</tt>
   * <li>this->factory_S().get() == factory_S.get()</tt>
   * </ul>
   */
  void set_factories(
    const mat_sym_fcty_ptr_t             &factory_transDtD
    ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
    );

  /** @name Dimensionality */
  //@{

  /** \brief Returns the number of decomposed equality constraints (<tt>r <= m</tt>).
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * The default implementation returns <tt>this->con_decomp().size()</tt>.
   * This implementation will work for all implementations.
   */
  virtual size_type r() const;

  //@}

  /** @name Ranges for dependent and independent variables and decomposed and undecomposed equalities
   */
  //@{

  /** \brief Return the range of dependent (i.e.\ basic) variables.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * The default implementation returns <tt>Range1D(1,this->m())</tt>.
   */
  virtual Range1D var_dep() const;
  /** \brief Return the range of independent (i.e.\ nonbasic) variables.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * The default implementation returns <tt>Range1D(this->m()+1,this->n())</tt>.
   */
  virtual Range1D var_indep() const;
  /** \brief Return the range of decomposed equality constraints.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * The default implementation returns <tt>Range1D(1,this->m())</tt>.
   */
  virtual Range1D con_decomp() const;
  /** \brief Return the range of undecomposed equality constraints.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * The default implementation returns <tt>Range1D::Invalid</tt>.
   */
  virtual Range1D con_undecomp() const;

  //@}

  /** @name Matrix factory objects */
  //@{
  
  /** \brief Return a matrix factory object for creating <tt>GcU</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * The default implementation is to return <tt>return.get() == NULL</tt>.
   * This is the proper implementation when <tt>m() == r()</tt>.
   * When <tt>m() > r()</tt> then the subclass must override this method to
   * return a valid matrix factory object.  Moreover, the returned
   * matrix object from <tt>this->factory_GcU()->create()->get_sub_view(rng,Range1D())</tt>
   * must be non-null for <tt>rng == this->var_dep()</tt> or <tt>rng == this->var_indep()</tt>.
   * This gives access to the matrices <tt>E'</tt> and <tt>F'</tt> as shown above.
   */
  virtual const mat_fcty_ptr_t factory_GcU() const;
  /** \brief Return a matrix factory object for <tt>D = -inv(C)*N</tt> {abstract}.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   */
  virtual const mat_fcty_ptr_t factory_D() const = 0;
  /** \brief Return a matrix factory object for <tt>Uz = F + E * D</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * The default implementation is to return <tt>return.get() == NULL</tt>.
   * This is the correct implementation when <tt>m() == r()</tt>.  However,
   * when <tt>m() > r()</tt> this method must be overridden to return a
   * non-null matrix factory object.
   */
  virtual const mat_fcty_ptr_t factory_Uz() const;
  /** \brief Return a matrix factory object for a mutable matrix compatible with <tt>GcU(var_dep)</tt>.
   *
   * This matrix factory object is designed to create mutable matrix objects compatible
   * with <tt>GcU(var_dep)</tt>.  For example, a matrix object <tt>Uy</tt> created by this matrix factory
   * can be used to compute <tt>Uy = Gc(var_dep,con_undecomp)' - Gc(var_indep,con_undecomp)'*D'</tt>
   * (this is needed by a orthogonal range/null decomposition.
   *
   * The default implementation is to return <tt>return.get() == NULL</tt>.
   * This is the correct implementation when <tt>m() == r()</tt>.  However,
   * when <tt>m() > r()</tt> this method must be overridden to return a
   * non-null matrix factory object.
   */
  virtual const mat_fcty_ptr_t factory_GcUD() const;

  /** \brief Returns a matrix factory for the result of <tt>J = D'*D</tt>
   * 
   * The resulting matrix is symmetric but is assumed to be singular.
   */
  virtual const mat_sym_fcty_ptr_t factory_transDtD() const;
  
  /** \brief Returns a matrix factory for the result of <tt>S = I + D'*D</tt>
   * 
   * The resulting matrix is symmetric and is guarrenteed to be nonsingular
   */
  virtual const mat_sym_nonsing_fcty_ptr_t factory_S() const;

  //@}

  /** @name Calculation members */
  //@{

  /** \brief Compute all of the needed quanities for direct sensitivities.
   *
   *	@param	x	[in] (dim == n()) Current value of unkowns.  This vector should
   *              have been created by <tt>this->space_x()->create_member()</tt>.
   *	@param	f 	[out] Value of <tt>f(x)</tt>.
   *				If f == NULL then this quantity is not computed.
   *	@param	c 	[in/out] (dim == m()) Value of the equality constraints \a c(x).
   *				If <tt>c == NULL</tt> then this quantity is not computed.
   *				If </tt>c != NULL and <tt>recalc_c == true</tt> then this quantity is recomputed.
   *				If </tt>c != NULL and <tt>recalc_c == false</tt> then this quantity is not
   *				recomputed and is used in the computation of \c py if requested (i.e. <tt>py != NULL</tt>).
   *				If <tt>c != NULL</tt> this this vector should have been created by
   *				<tt>this->space_c()->create_member()</tt>.
   *	@param  recalc_c
   *				[in] If \c true then \c c will be recomputed at \c x.
   *              If \c false then <tt>c</tt> will not be recomputed but will be used as stated above.
   *
   *	@param	Gf	[out] (dim == n()) Gradient of <tt>f(x)</tt>.
   *				If <tt>Gf == NULL</tt> then this quantity is not computed.  If <tt>Gf!=NULL</tt> this
   *				this vector should have been created by <tt>this->space_x()->create_member()</tt>.
   *	@param	py
   *				[out] (dim == r()) <tt>py = -inv(C)*c(con_decomp)</tt>.
   *				If <tt>py == NULL</tt> then this quantity is not computed.
   *				If <tt>recalc_c == false</tt> on input then the input <tt>c != NULL</tt> argument may
   *				be used	in the computation of \c py.  If <tt>py!=NULL</tt> this this vector should have
   *				been created by <tt>this->space_x()->sub_space(this->var_dep())->create_member()</tt>.
   *	@param	rGf
   *				[out] (dim == n()-r()) <tt>rGf = Gf(var_indep()) + D'*Gf(var_dep())</tt>,
   *				which is the reduced gradient of the objective function projected
   *				into the manifold of the decomposed equality constraints.  If <tt>rGf==NULL</tt>,
   *				this vector is not computed.  If <tt>rGf!=NULL</tt> then this vector
   *              should have been created by <tt>this->space_x(this->var_indep())->create_member()</tt>.
   *  @param  GcU [out] (dim = n x (m()-r())) Auxiliary jacobian matrix <tt>del(c(con_undecomp),x)</tt>.
   *              If m() == r() then <tt>GcU</tt> should be set to <tt>NULL</tt> on input.
   *              If GcU == NULL then this quantitiy is not computed.  If <tt>!=NULL</tt> this this matrix
   *              should have been created by <tt>this->factory_GcU()->create()</tt>.
   *	@param	D   [out] (dim = r() x (n()-r())) <tt>D = -inv(C)*N</tt>, which is the direct
   *              sensitivity of the constraints to the independent variables.
   *				If D == NULL then this quantity is not computed.  If <tt>!=NULL</tt> this this matrix
   *              should have been created by <tt>this->factory_D()->create()</tt>.
   *	@param	Uz   [out] (dim = (m()-r()) x (n()-r())) <tt>Uz = F + E * D</tt>, which is the an
   *              auxiliary sensitivity matrix.  If <tt>m() == r()</tt> then <tt>Uz</tt> should be set to
   *              <tt>NULL</tt> on input.  If <tt>Uz==NULL</tt> then this quantity is not computed.
   *              If <tt>!=NULL</tt> this this matrix should have been created by
   *              <tt>this->factory_Uz()->create()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> Lots more!
   * </ul>
   */
  virtual void calc_point(
    const Vector     &x
    ,value_type      *f
    ,VectorMutable   *c
    ,bool            recalc_c
    ,VectorMutable   *Gf
    ,VectorMutable   *py
    ,VectorMutable   *rGf
    ,MatrixOp        *GcU
    ,MatrixOp        *D
    ,MatrixOp        *Uz
    ) const = 0;

  /** \brief Calculate an approximate newton step given the Jacobian computed
   * for the last call to <tt>calc_point()</tt>.
   *
   * The idea behind this method is that with some applications it may be
   * much cheaper to compute an approximate Newton step for the constraints
   * given information computed during the last call to <tt>calc_point()</tt>.
   * It is assumed that this approximate solution <tt>py</tt> will still be a
   * descent direction for <tt>c(x)</tt>.  Some subclasses may have to perform an equal
   * amount of work as <tt>calc_point(...)</tt> to perform this calculation but those
   * are the breaks.
   *
   *	@param	x	[in] (dim == n()) current value of unkowns.
   *	@param	c 	[out] (dim == m()) Value of the constraints c(x)
   *				If c == NULL then this quantity is not computed.
   *				If c != NULL and recalc_c == true on input then this quantity is
   *				not recomputed and is used in the computation of
   *				py if requested (i.e. py!=NULL).
   *	@param  recalc_c
   *	@param	py
   *				[out] (size == r() on output) Approximate value of -inv(C)*c
   *				Note that py == NULL is not allowed here.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> Lots more.
   * </ul>
   */
  virtual void calc_semi_newton_step(
    const Vector    &x
    ,VectorMutable  *c
    ,bool           recalc_c
    ,VectorMutable  *py
    ) const = 0;

  //@}

  /** @name Overridden from NLP */
  //@{

  /** \brief Initialize the NLP for its first use.
    *
    * This function implementation should be called by subclass implementations
    * in order to reset counts for \c f(x), \c c(x), \c h(x) and \c Gf(x) evaluations.
    * This implementation calls <tt>this->NLPObjGrad::initialize()</tt>
    *
    * Postconditions:<ul>
    * <li> See <tt>NLPObjGrad::initialize()</tt>
    * </ul>
    */
  void initialize(bool test_setup);

  //@}

private:
  mat_sym_fcty_ptr_t             factory_transDtD_;
  mat_sym_nonsing_fcty_ptr_t     factory_S_;

};	// end class NLPDirect

}	// end namespace NLPInterfacePack

#endif   // NLP_FIRST_ORDER_DIRECT_H
