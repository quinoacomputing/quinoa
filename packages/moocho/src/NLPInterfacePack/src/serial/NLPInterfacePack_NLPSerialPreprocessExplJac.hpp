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

#ifndef NLP_SERIAL_PREPROCESS_EXPL_JAC_H
#define NLP_SERIAL_PREPROCESS_EXPL_JAC_H

#include <valarray>

#include "NLPInterfacePack_NLPSerialPreprocess.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_BasisSystemFactoryStd.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "Teuchos_AbstractFactory.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace NLPInterfacePack {

/** \brief NLP node subclass complementing \c NLPSerialPreprocess for explicit Jacobians.
 *
 * This subclass does a lot of work.  It has to consider several different
 * types of variability.  The matrices \c Gc and \c Gh that are computed must
 * take into consideration whether or not inequalities are converted
 * to equalities (<tt>convert_inequ_to_equ</tt>) and the permutation
 * of the entries according to the current basis selection.
 *
 \verbatim

    Gc = P_var * [  Gc_orig    Gh_orig   ] * P_equ'
                 [     0          -I     ]

 \endverbatim
 * This class also comes with a default implementation for the
 * <tt>BasisSystemPerm</tt> object which is created by a
 * <tt>BasisSystemPermFactory</tt> object that the client (or the
 * subclass) can specify.  The default implementation for this factory
 * object is from <tt>BasisSystemPermFactoryStd</tt> which uses the
 * <tt>AbstractLinAlgPack::BasisSystemPermDirectSparse</tt> subclass and
 * supports several different linear solvers by default.  The client
 * (or subclass) can augment the list of supported linear solvers
 * easily.
 *
 * ToDo: Finish documentation!
 *
 * <b>Subclass developers</b>
 *
 * Subclass developer's don't have to worry about slack variables or basis
 * permutations.  A concreate subclass just has to override the functions
 * that defined the original %NLP (see the tutorial example %NLP ???).
 *
 * In addition to the methods that must be overridden in \c NLPSerialPreprocess
 * (<A HREF="classNLPInterfacePack_1_1NLPSerialPreprocess.html#must_override">see</A>)
 * the following methods must be overridden as well: \c imp_Gc_nz_orig(), \c imp_Gh_nz_orig(),
 * \c imp_calc_Gc_orig(), \c imp_calc_Gh_orig().
 */
class NLPSerialPreprocessExplJac
  : virtual public NLPSerialPreprocess
  , virtual public NLPFirstOrder
{
public:

  /** @name Public types */
  //@{
  
  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixOp> >    factory_mat_ptr_t;

  //@}

  /** @name Constructors / initializers */
  //@{

  /// Set the <tt>BasisSystemFactory</tt> object used to create the basis system.
  STANDARD_COMPOSITION_MEMBERS( BasisSystemFactory, basis_sys_fcty );

  /** \brief Calls <tt>this->set_basis_sys_fcty()</tt> and <tt>this->set_mat_factories()</tt> methods.
   */
  NLPSerialPreprocessExplJac(
    const basis_sys_fcty_ptr_t  &basis_sys_fcty  = Teuchos::rcp(new BasisSystemFactoryStd())
    ,const factory_mat_ptr_t    &factory_Gc_full = Teuchos::null
    );

  /** \brief Initialize with matrix factory for original matrices \c Gc.
   *
   * This matrix type will be used for \c AbstractLinAlgPack::MatrixPermAggr::mat_orig()
   * returned by the initialized \c Gc.
   *
   * @param  factory_Gc_full
   *                [in] Smart pointer to matrix factory for \c Gc_full.  If
   *                <tt>factory_Gc_full.get() == NULL</tt> then the concrete matrix
   *                type ??? will be used as the default.
   */
  void set_factory_Gc_full( const factory_mat_ptr_t &factory_Gc_full );

  //@}

  /** @name Overridden public members from NLP */
  //@{

  /// Passes these options on to <tt>this->basis_sys_fcty().set_options(options)</tt>.
  void set_options( const options_ptr_t& options );
  /** \brief . */
  const options_ptr_t& get_options() const;
  /** \brief . */
  void initialize(bool test_setup);	
  /** \brief . */
  bool is_initialized() const;

  //@}

  /** @name Overridden public members from NLPFirstOrder */
  //@{
  
  /** \brief . */
  const mat_fcty_ptr_t factory_Gc() const;
  /// Calls <tt>basis_sys_fcty()->create()</tt>
  const basis_sys_ptr_t basis_sys() const;
  /// Validates the type of Gc is correct
  void set_Gc(MatrixOp* Gc);

  //@}

  /** @name Overridden public members from NLPVarReductPerm */
  //@{

  /** \brief . */
  bool get_next_basis(
    Permutation*  P_var,   Range1D* var_dep
    ,Permutation* P_equ,   Range1D* equ_decomp
    );
  /** \brief . */
  void set_basis(
    const Permutation   &P_var,   const Range1D  &var_dep
    ,const Permutation  *P_equ,   const Range1D  *equ_decomp
    );

  //@}


protected:

  /** @name Overridden protected members from NLPFirstOrder */
  //@{

  /** \brief . */
  void imp_calc_Gc(
    const Vector& x, bool newx
    ,const FirstOrderInfo& first_order_info
    ) const;
  
  //@}

  /** @name Protected types */
  //@{

  /** \brief Struct for zero and explicit first order quantities that subclass must fill in.
   *
   * When computing <tt>Gc</tt> and/or <tt>Gh</tt>, the subclass can
   * be instructed to set the row and columns index arrays by
   * setting <tt>Gc_ivect!=NULL</tt> and/or <tt>Gh_ivect!=NULL</tt>
   * respecitively.
   *
   * Objects of this type are passed on to subclasses and contain
   * pointers to quantities to be updated.  Note that %NLP
   * subclasses are not to resize the <tt>DVector</tt> or
   * <tt>std::valarray</tt> objects <tt>Gc_val</tt>,
   * <tt>Gc_ivect</tt>, <tt>Gc_jvect</tt>, <tt>Gh_val</tt>,
   * <tt>Gh_ivect</tt>, <tt>Gh_jvect</tt>, <tt>*Gf</tt>, <tt>*c</tt>
   * or </tt>h</tt> since the these will already be resized.
   *
   * The storage format for the gradient matrices <tt>Gc</tt> and
   * <tt>Gh</tt> use the coordinate data structure.  For <tt>Gc</tt>,
   * for instance, the elements are stored as:
   \verbatim

   for k = 0 ... Gc_nz
       Gc(Gc_ivect[k],Gc_jvect[k]) == Gc_val[k]
   \endverbatim
   * and all of the other matrix entries in <tt>Gc</tt> are
   * implicitly zero.
   *
   * In general, it is allowed for duplicate entries 
   * <tt>(Gc_ivect[k],Gc_jvect[k])</tt> to exist with the
   * convention that the corresponding <tt>Gc_val[k]</tt>
   * are to be added in matrix operations.  This is a relaxed
   * requirement that can make things much more complicated for
   * the code that accesses these matrix entries.
   */
  struct FirstOrderExplInfo {
    /** \brief . */
    typedef std::valarray<value_type>    val_t;
    /** \brief . */
    typedef std::valarray<index_type>    ivect_t;
    //
    typedef std::valarray<index_type>    jvect_t;
    /** \brief . */
    FirstOrderExplInfo()
      :Gc_val(NULL), Gc_ivect(NULL), Gc_jvect(NULL)
      ,Gh_val(NULL), Gh_ivect(NULL), Gh_jvect(NULL)
      ,f(NULL)
    {}
    /** \brief . */
    FirstOrderExplInfo(
      index_type* Gc_nz_in, val_t* Gc_val_in, ivect_t* Gc_ivect_in, jvect_t* Gc_jvect_in
      ,index_type* Gh_nz_in, val_t* Gh_val_in, ivect_t* Gh_ivect_in, jvect_t* Gh_jvect_in
      ,const ObjGradInfoSerial& obj_grad
      )
      :Gc_nz(Gc_nz_in), Gc_val(Gc_val_in), Gc_ivect(Gc_ivect_in), Gc_jvect(Gc_jvect_in)
      ,Gh_nz(Gh_nz_in), Gh_val(Gh_val_in), Gh_ivect(Gh_ivect_in), Gh_jvect(Gh_jvect_in)
      ,Gf(obj_grad.Gf), f(obj_grad.f), c(obj_grad.c), h(obj_grad.h)
    {}
    /** \brief . */
    size_type*    Gc_nz;
    /** \brief . */
    val_t*        Gc_val;
    /** \brief . */
    ivect_t*      Gc_ivect;
    /** \brief . */
    jvect_t*      Gc_jvect;
    /** \brief . */
    size_type*    Gh_nz;
    /** \brief . */
    val_t*        Gh_val;
    /** \brief . */
    ivect_t*      Gh_ivect;
    /** \brief . */
    jvect_t*      Gh_jvect;
    /** \brief . */
    DVector*       Gf;
    /** \brief . */
    value_type*   f;
    /** \brief . */
    DVector*       c;
    /** \brief . */
    DVector*       h;
  }; // end struct FirstOrderExplInfo

  //@}

  /** @name Pure virtual template methods to be defined by subclasses */
  //@{

  /** \brief Return the number of nonzero elements in \c Gc before elements are removed for fixed variables.
    *
    * The value returned from this method before the first time \c imp_calc_Gc() is called
    * is an upper estimate of the number of nonzeros.  To get the actual number
    * of nonzeros, call this function again after \c imp_calc_Gc() has been called.
    */
  virtual size_type imp_Gc_nz_orig() const = 0;

  /** \brief Return the number of nonzero elements in \c Gh before elements are removed for fixed variables.
    *
    * The value returned from this method before the first time \c imp_calc_Gh() is called
    * is an upper estimate of the number of nonzeros.  To get the actual number
    * of nonzeros, call this function again after \c imp_calc_Gh() has been called.
    */
  virtual size_type imp_Gh_nz_orig() const = 0;

  /** \brief Calculate the COOR matrix for the gradient for all of the
   * <tt>c(x)</tt> constaints in the original %NLP.
   *
   * @param x_full  [in] Unknown vector (size n_full).
   * @param newx    [in] True if is a new point.
   * @param first_order_expl_info
   *                [out] Pointers to zero and first order quantities .
   *                On output, <tt>*first_order_expl_info.Gc_nz</tt> must be set to the actual
   *                number of nonzero elements in <tt>Gc</tt> and the array of nonzero entry
   *                values <tt>*first_order_expl_info.Gc_val</tt> must also be set.
   *                The nonzero structure must also be set in the arrays
   *                <tt>*first_order_expl_info.Gc_ivect</tt> and
   *                <tt>*first_order_expl_info.Gc_jvect</tt> if
   *                <tt>first_order_expl_info.Gc_ivect != NULL</tt>.
   *                In addition, any of the other quantities pointed to in
   *                <tt>first_order_expl_info</tt> may be set on
   *                output, but are not guaranteed to be.
   *
   * Preconditions:<ul>
   * <li> <tt>first_order_expl_info.Gc_nz != NULL</tt>
   * <li> <tt>first_order_expl_info.Gc_val != NULL</tt>
   * <li> <tt>(first_order_expl_info.Gc_ivect != NULL) == (first_order_expl_info.Gc_jvect != NULL)</tt> 
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*first_order_expl_info.Gc_nz</tt> is updated to number of nonzero elements set in
   *      <tt>*first_order_expl_info.Gc_val</tt>.
   * <li> <tt>(*first_order_expl_info.Gc_val)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gc_nz</tt>
   *      is set to the nonzero entry values in \c Gc.
   * <li> [<tt>first_order_expl_info.Gc_ivect != NULL</tt>]
   *      <tt>(*first_order_expl_info.Gc_ivect)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gc_nz</tt>
   *      is set to the row indexes for the nonzero entires in \c Gc.
   * <li> [<tt>first_order_expl_info.Gc_jvect != NULL</tt>]
   *      <tt>(*first_order_expl_info.Gc_jvect)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gc_nz</tt>
   *      is set to the column indexes for the nonzero entires in \c Gc.
   * </ul>
   *
   * Note that duplicate entires with the same row and column indexes are allowed.  In this case, the
   * matrix entries are considered to be summed.
   */
  virtual void imp_calc_Gc_orig(
    const DVectorSlice& x_full, bool newx
    , const FirstOrderExplInfo& first_order_expl_info
    ) const = 0;

  /** \brief Calculate the COOR matrix for the gradient for all of the
   * <tt>h(x)</tt> constaints in the original %NLP.
   *
   * @param x_full  [in] Unknown vector (size n_full).
   * @param newx    [in] True if is a new point.
   * @param first_order_expl_info
   *                [out] Pointers to zero and first order quantities .
   *                On output, <tt>*first_order_expl_info.Gh_nz</tt> must be set to the actual
   *                number of nonzero elements in <tt>Gh</tt> and the array of nonzero entry
   *                values <tt>*first_order_expl_info.Gh_val</tt> must also be set.
   *                The nonzero structure must also be set in the arrays
   *                <tt>*first_order_expl_info.Gh_ivect</tt> and
   *                <tt>*first_order_expl_info.Gh_jvect</tt> if
   *                <tt>first_order_expl_info.Gh_ivect != NULL</tt>.
   *                In addition, any of the other quantities pointed to in
   *                <tt>first_order_expl_info</tt> may be set on
   *                output, but are not guaranteed to be.
   *
   * Preconditions:<ul>
   * <li> <tt>first_order_expl_info.Gh_nz != NULL</tt>
   * <li> <tt>first_order_expl_info.Gh_val != NULL</tt>
   * <li> <tt>(first_order_expl_info.Gh_ivect != NULL) == (first_order_expl_info.Gh_jvect != NULL)</tt> 
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*first_order_expl_info.Gh_nz</tt> is updated to number of nonzero elements set in
   *      <tt>*first_order_expl_info.Gh_val</tt>.
   * <li> <tt>(*first_order_expl_info.Gh_val)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gh_nz</tt>
   *      is set to the nonzero entry values in \c Gh.
   * <li> [<tt>first_order_expl_info.Gh_ivect != NULL</tt>]
   *      <tt>(*first_order_expl_info.Gh_ivect)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gh_nz</tt>
   *      is set to the row indexes for the nonzero entires in \c Gh.
   * <li> [<tt>first_order_expl_info.Gh_jvect != NULL</tt>]
   *      <tt>(*first_order_expl_info.Gh_jvect)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gh_nz</tt>
   *      is set to the column indexes for the nonzero entires in \c Gh.
   * </ul>
   *
   * Note that duplicate entires with the same row and column indexes are allowed.  In this case, the
   * matrix entries are considered to be summed.
   */
  virtual void imp_calc_Gh_orig(
    const DVectorSlice& x_full, bool newx
    , const FirstOrderExplInfo& first_order_expl_info
    ) const = 0;

  //@}

  /** @name Protected member functions for subclasses to use */
  //@{

  /// Assert if we have been initizlized (throws UnInitialized)
  void assert_initialized() const;

  /** \brief . */
  const FirstOrderExplInfo first_order_expl_info() const;

  //@}

private:

  // ////////////////////////////////////////
  // Private data members
  
  bool initialized_;              // Flag for if the NLP has has been properly initialized
  bool test_setup_;               // Flag for if to test the setup of things or not
  options_ptr_t options_;         // The options being used

  factory_mat_ptr_t   factory_Gc_full_;
  mat_fcty_ptr_t      factory_Gc_;

  mutable size_type   Gc_nz_orig_;    // Number of nonzeros in the original NLP Gc
  mutable size_type   Gh_nz_orig_;    // Number of nonzeros in the original NLP Gh
  mutable size_type   Gc_nz_full_;    // Number of nonzeros in the full NLP Gc
  mutable size_type   Gh_nz_full_;    // Number of nonzeros in the full NLP Gh
  mutable FirstOrderExplInfo::val_t    Gc_val_orig_;   // Storage for explicit nonzeros of full Gc
  mutable FirstOrderExplInfo::ivect_t  Gc_ivect_orig_;
  mutable FirstOrderExplInfo::jvect_t  Gc_jvect_orig_;
  mutable FirstOrderExplInfo::val_t    Gh_val_orig_;   // Storage for explicit nonzeros of orig Gh
  mutable FirstOrderExplInfo::ivect_t  Gh_ivect_orig_;
  mutable FirstOrderExplInfo::jvect_t  Gh_jvect_orig_;

  mutable bool                         Gc_perm_new_basis_updated_;  // Flag for if a new basis was set!

  // ////////////////////////////
  // Private member functions

  //
  void imp_calc_Gc_or_Gh(
    bool calc_Gc
    ,const Vector& x, bool newx
    ,const FirstOrderInfo& first_order_info
    ) const;

  //
  void imp_fill_jacobian_entries(
    size_type           n             // [in]
    ,size_type          n_full        // [in]
    ,bool               load_struct   // [in] If true, then the structure is loaded also
    ,const index_type   col_offset    // [in] Offset for filled column indexes
    ,const value_type   *val_full     // [in] Values (!=NULL)
    ,const value_type   *val_full_end // [in] Values end (!=NULL)
    ,const index_type   *ivect_full   // [in] Row indexes (!=NULL)
    ,const index_type   *jvect_full   // [in] Column indexes (!=NULL)
    ,index_type         *nz           // [in/out] Number of nonzeros added (!=NULL)            
    ,value_type         *val_itr      // [out] Values to fill (!=NULL)
    ,index_type         *ivect_itr    // [out] Row indexes (can be NULL if load_struct == false)
    ,index_type         *jvect_itr    // [out] Column indexes  (can be NULL if load_struct == false)
    ) const;

};	// end class NLPSerialPreprocessExplJac

// ///////////////////////////
// inline members

inline
const NLPSerialPreprocessExplJac::FirstOrderExplInfo
NLPSerialPreprocessExplJac::first_order_expl_info() const
{
  return FirstOrderExplInfo(
    &Gc_nz_orig_
    ,&Gc_val_orig_
    ,&Gc_ivect_orig_
    ,&Gc_jvect_orig_
    ,&Gh_nz_orig_
    ,&Gh_val_orig_
    ,&Gh_ivect_orig_
    ,&Gh_jvect_orig_
    ,obj_grad_orig_info()
    );
}

}	// end namespace NLPInterfacePack 

#endif // NLP_SERIAL_PREPROCESS_EXPL_JAC_H
