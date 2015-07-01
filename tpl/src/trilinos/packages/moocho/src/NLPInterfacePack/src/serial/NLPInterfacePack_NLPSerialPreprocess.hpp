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

#ifndef NLP_FULL_TO_REDUCED_H
#define NLP_FULL_TO_REDUCED_H

#include <valarray>

#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "NLPInterfacePack_NLPVarReductPerm.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_VectorMutableDense.hpp"
#include "AbstractLinAlgPack_PermutationSerial.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_IVector.hpp"

namespace NLPInterfacePack {

/** \brief %NLP node implementation subclass for preprocessing and basis manipulation.
 *
 * This is an implementation node class that takes an original NLP and transforms
 * it by:
 * <ul>
 * <li> Converting general inequalities to equalities with slack variables
 * <li> Removing variables fixed by bounds
 * <li> Converting general inequalities with cramped bounds into general equalities
 * <li> Reordering the quantities according to the current basis selection
 *      (by implementing the \c NLPVarReductPerm interface).
 * </ul>
 *
 * <b>Original %NLP formulation</b>
 *
 * The original %NLP (as specified by the subclass) takes the form:
 *
 \verbatim

     min    f_orig(x_orig)
     s.t.   c_orig(x_orig) = 0
            hl_orig <= h(x_orig) <= hu_orig
            xl_orig <= x_orig <= xu_orig
  where:
          x_orig         <: REAL^n_orig
          f_orig(x_orig) <: REAL^n_orig -> REAL
          c_orig(x_orig) <: REAL^n_orig -> REAL^m_orig
          h_orig(x_orig) <: REAL^n_orig -> REAL^mI_orig
 \endverbatim
 *
 * <b>Conversion of general inequalities to equalities using slack variables</b>
 *
 * The original %NLP formulation above
 * is transformed by adding slack variables <tt>s_orig <: REAL^mI_orig</tt>,
 * defining a new <tt>x_full = [ x_orig; s_orig ]</tt> and forming the new %NLP:
 \verbatim

     min    f_full(x_full)
     s.t.   c_full(x_full) = 0
            xl_full <= x_full <= xu_full

  where:

          x_full           = [ x_orig ]  n_orig
                             [ s_orig ]  mI_orig

          f_full(x_full)   = f_orig(x_orig)

          c_full(x_full)   = [ c_orig(x_orig)          ]  m_orig
                             [ h_orig(x_orig) - s_orig ]  mI_orig

          xl_full          = [ xl_orig ]  n_orig
                             [ hl_orig ]  mI_orig

          xu_full          = [ xu_orig ]  n_orig
                             [ hu_orig ]  mI_orig
 \endverbatim
 * Note that in this case, the Jacobian of the new equality constraints
 * and the gradient of the new objective become:
 \verbatim

    Gc_full = [  Gc_orig    Gh_orig   ]  n_orig
              [     0          -I     ]  mI_orig

                  m_orig     mI_orig

    Gf_full = [ Gf_orig ]  n_orig
              [    0    ]  mI_orig
 \endverbatim
 * It is up to the subclass to implement \c imp_calc_Gc()
 * and \c imp_calc_Gh() in a way that is consistent with the above
 * transformation while also considering basis permutations (see 
 * \c NLPSerialPreprocessExplJac).  As for the gradient
 * \c Gc_full, the subclass can actually include terms for the slack
 * variables in the objective function but the most common behavior
 * will be to just ignore slack variables in the subclass.
 *
 * <b>Preprocessing and basis manipulation</b>
 *
 * The initial basis selection is the original order (<tt>x_full = [ x_orig; s_orig ]</tt>)
 * with the variables fixed by bounds being removed,
 * and assumes there are no dependent equations (<tt>r == m</tt>).
 *
 * The implementations of the Jacobian matrices \c Gc and \c Gh are not determined here and
 * must be defined by an %NLP subclass (see \c NLPSerialPreprocessExplJac for example).
 *
 * This class stores the variable permutations and processing information in two parts.
 * In the first state, the fixed variables are removed as:
 \verbatim

    var_remove_fixed_to_full =  [ not fixed by bounds |  fixed by bounds  ]
                                [1 ..                n|n+1 ..        n_full]
 \endverbatim
 *
 * The mapping <tt>i_full = var_remove_fixed_to_full()(i_free_fixed)</tt> gives the index of the
 * original variable (\c i_full) for the sets of variables not fixed and fixed by bounds.
 *
 * The inverse mapping <tt>i_free_fixed = var_full_to_remove_fixed()(i_full)</tt> can be used
 * to determine if a variable is fixed by bounds or not..
 *
 * On top of this partitioning of free and fixed variables, there is a second stage which
 * is a permutation of the free variables into dependent and independent sets that is needed
 * by the client.
 \verbatim

    var_perm = [ dependent variables | independent variables ]
               [1..               n-r|n-r+1...              n]
 \endverbatim
 *
 * The mapping <tt>i_free_fixed = var_perm()(i_perm)</tt> is used to determine the index
 * of a free variable in \c var_remove_fixed_to_full() given its index (\c i_perm) 
 * for the current basis selection.
 *
 * For example, if \c x is the vector of variables for the current basis selection
 * and \c x_full is the vector of variables in the original order including
 * the fixed variables then the following is true:
 *
 * <tt>x(i) == x_full(var_remove_fixed_to_full()(var_perm()(i))), for i = 1...n</tt>
 *
 * The permutation <tt>equ_perm()</tt> gives the partitioning of the equality constraints
 * into decomposed and undecomposed equalities.  Decomposed inequality constraints are not
 * supported currently.
 *
 * <b>Subclass developers notes</b>
 *
 * Handling of multiple updates by subclasses: Here we discuss the protocol for the 
 * handling of multiple updates to quantities during the calculation of other quantities.
 * In order to simplify the implementation of subclasses as much as possible, storage
 * for all iteration quantities will be passed to the subclass in the methods
 * <tt>imp_calc_f_orig()</tt>, <tt>imp_calc_c_orig()</tt>, <tt>imp_calc_h_orig()</tt>
 * and <tt>imp_calc_Gf_orig()</tt> regardless of what quantities where set by the user
 * in the <tt>NLP</tt> interface.  The subclass can always find out what was set
 * by the client by calling <tt>get_f()</tt>, <tt>get_c()</tt>, <tt>get_Gf()</tt>
 * etc.  Therefore, in general, clients should just only compute what is required
 * in each call to <tt>imp_calc_xxx_orig()</tt> and only update other quantities
 * if it is absolutely free to do so (e.g. computing a function value when a gradient
 * is computed using AD) or is required to do so (e.g.an external interface that 
 * forces both <tt>f_orig(x_orig)</tt>, <tt>c_orig(x_orig)</tt> and <tt>h_orig(x_orig)</tt>
 * be computed at the same time).  It is up to the subclass to remember when a quantity
 * has already been computed so that it will not be computed again unnecessarily.  It is
 * always safe for the subclass to ignore these issues and just do what is easiest.
 * More careful implementations can be handled by the subclass by keeping track of
 * <tt>get_xxx()</tt> and <tt>newx</tt> and remembering when quantities are computed.
 *
 * <A NAME="must_override"></A>
 *
 * The following methods from the \c NLP interface must be overridden by the %NLP subclass:
 * \c max_var_bounds_viol(), \c set_multi_calc(), \c multi_calc().
 *
 * The following methods from the \c NLPVarReductPerm interface must be overridden by the %NLP subclass:
 * \c nlp_selects_basis().
 *
 * In addition, the methods from this interface that must be overridden are: \c imp_n_orig(),
 * \c imp_m_orig(), \c imp_mI_orig(), \c imp_xinit_orig(), \c imp_has_var_bounds(),
 * \c imp_xl_orig(), \c imp_xu_orig(), \c imp_hl_orig(), \c imp_hu_orig(), \c imp_calc_f_orig(),
 * \c imp_calc_c_orig(), \c imp_calc_h_orig() and \c imp_calc_Gf_orig().
 *
 * <A NAME="should_override"></A>
 *
 * The \c NLP method \c initialize() should also be overridden by all of the subclasses
 * (and call \c initialize() on its direct subclass).
 *
 * The following methods (with default implementations) may also be overridden by a subclass:
 * \c imp_get_next_basis() and \c imp_report_orig_final_solution().
 */
class NLPSerialPreprocess
  : virtual public NLPObjGrad
  , virtual public NLPVarReductPerm
{
public:

  /** @name Exceptions */
  //@{

  /// Thrown if xl(i) > xu(i)
  class InconsistantBounds : public std::logic_error
  {public: InconsistantBounds(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}

  /** \brief Default Constructor.
    *
    * This initalizes the basis to the first basis if the subclass specifies one and
    * if not picks to first \c r variables as the dependent variables and the last
    * <tt>n-r</tt> variables as the independent variables.  Also the default behavior
    * is to force the initial point in bounds.
    */
  NLPSerialPreprocess();

  /** \brief Gives the value of a Lagrange multipler for a fixed variable bound
   *.that has been preprocessed out of the problem.
   */
  static value_type fixed_var_mult(); 
  
  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void force_xinit_in_bounds(bool force_xinit_in_bounds);
  /** \brief . */
  bool force_xinit_in_bounds() const;
  /** \brief . */
  void initialize(bool test_setup);	
  /** \brief . */
  bool is_initialized() const;
  /** \brief . */
  size_type n() const;
  /** \brief . */
  size_type m() const;
  /** \brief . */
  vec_space_ptr_t space_x() const;
  /** \brief . */
  vec_space_ptr_t space_c() const;
  /** \brief . */
  size_type num_bounded_x() const;
  /** \brief . */
  const Vector& xl() const;
  /** \brief . */
  const Vector& xu() const;
  /** \brief . */
  const Vector& xinit() const;
  /** \brief . */
  void get_init_lagrange_mult(
    VectorMutable*   lambda
    ,VectorMutable*  nu
    ) const;
  /** \brief . */
  void scale_f( value_type scale_f );
  /** \brief . */
  value_type scale_f() const;
  /** \brief Overridden to permute the variables back into an order that is natural to the subclass.
    *
    * The default implementation of this function is to call the method
    * <tt>imp_report_full_final_solution(x_full,lambda_full,nu_full)</tt>.
    * This function translates from \c x, \c lambda and \c nu into the original
    * order with fixed variables added back to form \c x_full, \c lambda_full, \c lambdaI_full
    * and \c nu_full.
    */
  void report_final_solution(
    const Vector&    x
    ,const Vector*   lambda
    ,const Vector*   nu
    ,bool            is_optimal
    );
  /** \brief . */
  virtual size_type ns() const;
  /** \brief . */
  vec_space_ptr_t space_c_breve() const;
  /** \brief . */
  vec_space_ptr_t space_h_breve() const;
  /** \brief . */
  const Vector& hl_breve() const;
  /** \brief . */
  const Vector& hu_breve() const;
  /** \brief . */
  const Permutation& P_var() const;
  /** \brief . */
  const Permutation& P_equ() const;

  //@}

  /** @name Overridden public members from NLPVarReductPerm */
  //@{

  /** \brief . */
  const perm_fcty_ptr_t factory_P_var() const;
  /** \brief . */
  const perm_fcty_ptr_t factory_P_equ() const;
  /** \brief . */
  Range1D var_dep() const;
  /** \brief . */
  Range1D var_indep() const;
  /** \brief . */
  Range1D equ_decomp() const;
  /** \brief . */
  Range1D equ_undecomp() const;
  /** \brief . */
    bool nlp_selects_basis() const;
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
  /** \brief . */
  void get_basis(
    Permutation*  P_var,   Range1D* var_dep
    ,Permutation* P_equ,   Range1D* equ_decomp
    ) const;

  //@}

protected:

  /** @name Overridden protected members from NLP */
  //@{

  /** \brief . */
  void imp_calc_f(
    const Vector            &x
    ,bool                   newx
    ,const ZeroOrderInfo    &zero_order_info
    ) const;
  /** \brief . */
  void imp_calc_c(
    const Vector            &x
    ,bool                   newx
    ,const ZeroOrderInfo    &zero_order_info
    ) const;
  /** \brief . */
  void imp_calc_c_breve(
    const Vector            &x
    ,bool                   newx
    ,const ZeroOrderInfo    &zero_order_info_breve
    ) const;
  /** \brief . */
  void imp_calc_h_breve(
    const Vector            &x
    ,bool                   newx
    ,const ZeroOrderInfo    &zero_order_info_breve
    ) const;

  //@}

  /** @name Overridden protected members from NLPObjGrad */
  //@{

  /** \brief . */
  void imp_calc_Gf(
    const Vector            &x
    ,bool                   newx
    ,const ObjGradInfo      &obj_grad_info
    ) const;

  //@}

  /** @name Protected types */
  //@{

  /** \brief Struct for objective and constriants (pointer) as serial vectors.
   *
   * Objects of this type are passed on to subclasses and contain pointers to
   * quantities to be updated.  Note that %NLP subclasses are not to resize
   * the <tt>DVector</tt> objects <tt>*c</tt> or </tt>h</tt> since the
   * these will already be resized.
   */
  struct ZeroOrderInfoSerial {
  public:
    /** \brief . */
        ZeroOrderInfoSerial() : f(NULL)
    {}
    /** \brief . */
    ZeroOrderInfoSerial( value_type* f_in, DVector* c_in, DVector* h_in )
      : f(f_in), c(c_in), h(h_in)
    {}
    /// Pointer to objective function <tt>f</tt>   (may be NULL if not set)
    value_type*    f;
    /// Pointer to constraints residual <tt>c</tt> (may be NULL if not set)
    DVector*        c;
    /// Pointer to constraints residual <tt>h</tt> (may be NULL if not set)
    DVector*        h;
  }; // end struct ZeroOrderInfoSerial

  /** \brief Struct for serial gradient (objective), objective and constriants (pointers)
   *
   * Objects of this type are passed on to subclasses and contain
   * pointers to quantities to be updated.  Note that %NLP
   * subclasses are not to resize the <tt>DVector</tt> objects
   * <tt>*Gf</tt>, <tt>*c</tt> or </tt>h</tt> since the these will
   * already be resized.
   */
  struct ObjGradInfoSerial {
  public:
    /** \brief . */
    ObjGradInfoSerial()	: f(NULL)
    {}
    /** \brief . */
    ObjGradInfoSerial( DVector* Gf_in, const ZeroOrderInfoSerial& first_order_info_in )
      : Gf(Gf_in), f(first_order_info_in.f), c(first_order_info_in.c), h(first_order_info_in.h)
    {}
    /// Gradient of objective function <tt>Gf</tt> (may be NULL if not set)
    DVector*        Gf;
    /// Pointer to objective function <tt>f</tt>   (may be NULL if not set)
    value_type*    f;
    /// Pointer to constraints residual <tt>c</tt> (may be NULL if not set)
    DVector*        c;
    /// Pointer to constraints residual <tt>h</tt> (may be NULL if not set)
    DVector*        h;
  }; // end struct ObjGradInfoSerial

  //@}

  /** @name Pure virtual methods to be defined by subclasses */
  //@{

  /** \brief Return if the definition of the %NLP has changed since the last call to \c initialize()
   *
   * The default return is \c true.  This function is present in order to avoid
   * preprocessing when \c initialize() is called but nothing has changed.
   */
  virtual bool imp_nlp_has_changed() const { return true; }
  /// Return the number of variables in the original problem (including those fixed by bounds)
  virtual size_type imp_n_orig() const = 0;
  /// Return the number of general equality constraints in the original problem.
  virtual size_type imp_m_orig() const = 0;
  /// Return the number of general inequality constraints in the original problem.
  virtual size_type imp_mI_orig() const = 0;
  /// Return the original initial point (size \c imp_n_orig()).
  virtual const DVectorSlice imp_xinit_orig() const = 0;
  /// Return if the %NLP has bounds
  virtual bool imp_has_var_bounds() const = 0;
  /** \brief Return the original lower variable bounds (size \c imp_n_orig()).
   *
   * Only to be called if <tt>this->imp_has_var_bounds() == true</tt>.
   * A lower bound is considered free if it is less than or equal to:
   \verbatim

   <tt>-NLP::infinite_bound()</tt>
   \endverbatim
   */
  virtual const DVectorSlice imp_xl_orig() const = 0;
  /** \brief Return the original upper variable bounds (size \c imp_n_orig()).
   *
   * Only to be called if <tt>this->imp_has_var_bounds() == true</tt>.
   * An upper bound is considered free if it is greater than or equal to:
   \verbatim

   <tt>+NLP::infinite_bound()</tt>
   \endverbatim
   */
  virtual const DVectorSlice imp_xu_orig() const = 0;
  /** \brief Return the original lower general inequality bounds (size \c imp_mI_orig()).
   *
   * Only to be called if <tt>this->imp_mI_orig() == true</tt>.
   * A lower bound is considered free if it is equal to:
   * 
   * <tt>-NLP::infinite_bound()</tt>
   */
  virtual const DVectorSlice imp_hl_orig() const = 0;
  /** \brief Return the original upper general inequality bounds (size \c imp_mI_orig()).
   *
   * Only to be called if <tt>this->imp_mI_orig() == true</tt>.
   * An upper bound is considered free if it is equal to:
   * 
   * <tt>+NLP::infinite_bound()</tt>
   */
  virtual const DVectorSlice imp_hu_orig() const = 0;
  /** \brief Calculate the objective function for the original %NLP.
   */
  virtual void imp_calc_f_orig(
    const DVectorSlice           &x_full
    ,bool                        newx
    ,const ZeroOrderInfoSerial   &zero_order_info
    ) const = 0;
  /** \brief Calculate the vector for all of the general equality constaints in the original %NLP.
   */
  virtual void imp_calc_c_orig(
    const DVectorSlice           &x_full
    ,bool                        newx
    ,const ZeroOrderInfoSerial   &zero_order_info
    ) const = 0;
  /** \brief Calculate the vector for all of the general inequality constaints in the original %NLP.
   */
  virtual void imp_calc_h_orig(
    const DVectorSlice           &x_full
    ,bool                        newx
    ,const ZeroOrderInfoSerial   &zero_order_info
    ) const = 0;
  /** \brief Calculate the vector for the gradient of the objective in the original NLP.
   *
   * Note that the dimension of <tt>obj_grad_info.Gf->dim()</tt> is
   * <tt>n_orig + mI_orig</tt>.
   *
   * On input, if <tt>mI_orig > 0</tt>
   * then <tt><tt>(*obj_grad_info.Gf)(n_orig+1,n_orig+mI_orig)</tt> is
   * initialized to 0.0 (since slacks do not ordinarily do not appear
   * in the objective function).  However, the subclass can assign
   * (smooth) contributions for the slacks if desired.
   */
  virtual void imp_calc_Gf_orig(
    const DVectorSlice           &x_full
    ,bool                        newx
    ,const ObjGradInfoSerial     &obj_grad_info
    ) const = 0;
  /** \brief Return the next basis selection (default returns \c false).
   *
   * @param  var_perm_full
   *                   [out] (size = <tt>n_orig + mI_orig</tt>).
   *                   Contains the variable permutations (including slack variables possibly).
   * @param  equ_perm_full
   *                   [out] (size = <tt>m_orig + mI_orig)</tt>).
   *                   Contains the constriant permutations (including general inequalities
   *                   possibly).
   * @param  rank_full
   *                   [out] Returns the rank of the basis before fixed variables are taken out of
   *                   \c var_perm_full and \c equ_perm.
   * @param  rank      [out] Returns the rank of the basis after fixed variables are taken out of
   *                   \c var_perm_full and \c equ_perm.
   *                   If there are no fixed variables then <tt>rank</tt> should be equal to
   *                   <tt>rank_full</tt>.
   *
   * Postconditions:<ul>
   * <li> [<tt>return == true</tt>]
   *      <tt>var_perm_full(i) < var_perm_full(i+1)</tt>, for <tt>i = 1...rank_full-1</tt>
   * <li> [<tt>return == true</tt>]
   *      <tt>var_perm_full(i) < var_perm_full(i+1)</tt>, for <tt>i = rank_full...n_full-1</tt>
   * <li> [<tt>return == true</tt>]
   *      <tt>equ_perm_full(i) < equ_perm_full(i+1)</tt>, for <tt>i = 1...rank-1</tt>
   * <li> [<tt>return == true</tt>]
   *      <tt>equ_perm_full(i) < equ_perm_full(i+1)</tt>, for <tt>i = rank...m_full-1</tt>
   * </ul>
   *
   * This method will only be called if <tt>this->nlp_selects_basis() == true</tt>.
   *
   * The basis returned by the subclass must be sorted
   * <tt>var_perm_full = [ dep | indep ]</tt> and <tt>equ_perm_full
   * = [ equ_decomp | equ_undecomp ]</tt>.  The subclass should not
   * remove the variables fixed by bounds from \c var_perm_full as
   * they will be removed by this class as they are translated.  In
   * addition, the subclass can also include slack variables in the
   * basis (if mI_orig > 0>/tt>).  Therefore, a nonsingular basis before
   * fixed variables are removed may not be nonsingular once the fixed
   * variables are removed.  During the translation of <tt>var_perm_perm</tt>,
   * the variables fixed by bounds are removed by compacting
   * <tt>var_perm_full</tt> and adjusting the remaining indices.
   * For this to be correct with variables fixed by bounds, it is
   * assumed that the subclass knows which variables are fixed by
   * bounds and can construct <tt>var_perm_full</tt> so that after
   * the translated the basis will be nonsingular.  The first
   * <tt>rank</tt> entries in <tt>var_perm_full[1:rank_full]</tt>
   * left after the fixed variables have been removed give the
   * indices of the dependent (basic) variables and the remaining
   * variables in <tt>var_perm_full[rank_full+1,n_full]</tt> are the
   * indices for the independent (nonbasic) variables.  To simplify
   * things, it would be wise for the %NLP subclass not to put fixed
   * variables in the basis since this will greatly simplify
   * selecting a nonsingular basis.
   *
   * The first time this method is called, the subclass should
   * return the first suggested basis selection (even if it happens
   * to be identical to the original ordering).
   *
   * The default implementation returns <tt>false</tt> which implies
   * that the %NLP subclass has no idea what a good basis selection
   * should be..
   */
  virtual bool imp_get_next_basis(
    IVector      *var_perm_full
    ,IVector     *equ_perm_full
    ,size_type   *rank_full
    ,size_type   *rank
    );
  /** \brief To be overridden by subclasses to report the final solution in the
   * original ordering natural to the subclass.
   *
   * Note that the lagrange multipliers for fixed variables that have been
   * preprocessed out of the problem are not computed by the optimization
   * algorithm and are therefore not available.  These multipliers are
   * designated with the special value \c fixed_var_mult() but the numerical
   * value is not significant.
   *
   * The default implementation of this function is to do nothing.
   */
  virtual void imp_report_orig_final_solution(
    const DVectorSlice      &x_full
    ,const DVectorSlice     *lambda_orig
    ,const DVectorSlice     *lambdaI_orig
    ,const DVectorSlice     *nu_orig
    ,bool                   optimal
    )
  {}

  //@}

  /** @name Other protected implementation functions for subclasses to call */
  //@{

  /// Used by subclasses to set the state of the NLP to not initialized.
  void set_not_initialized();

  /// Assert if we have been initizlized (throws UnInitialized)
  void assert_initialized() const;

  /// Set the full x vector if <tt>newx == true</tt>
  void set_x_full(const DVectorSlice& x, bool newx, DVectorSlice* x_full) const;

  /// Give reference to current x_full
  DVectorSlice x_full() const;

  /** \brief . */
  const ZeroOrderInfoSerial zero_order_orig_info() const;

  /** \brief . */
  const ObjGradInfoSerial obj_grad_orig_info() const;
  
  /** \brief Permutation vector for partitioning free and fixed variables.
   *
   \verbatim

   var_remove_fixed_to_full =  [ not fixed by bounds |  fixed by bounds  ]
                               [1 ..                n|n + 1 ..     n_full]
   \endverbatim
   * The mapping <tt>i_full = var_remove_fixed_to_full()(i_free_fixed)</tt> gives the index of the
   * original variable (\c i_full) for the sets of variables not fixed and fixed (upper
   * and lower bounds where equal).
   */
  const IVector& var_remove_fixed_to_full() const;

  /** \brief Inverse permutation vector of \c var_remove_fixed_to_full().
   *
   * The inverse mapping <tt>i_free_fixed = var_full_to_remove_fixed()(i_full)</tt> can be used
   * to determine if a variable is free for fixed.
   */
  const IVector& var_full_to_remove_fixed() const;

  /** \brief Permutes from the compated variable vector (removing fixed variables) to the current
   * basis selection.
   *
   * On top of this partitioning of free and fixed variables, there is a permutation
   * of the free variables into dependent and independent variables that is needed
   * by the optimization algorithm.
   *
   \verbatim

   var_perm = [ dependent variables | independent variables ]
              [1..                 r|r+1..                 n]
   \endverbatim
   *
   * The mapping <tt>i_free_fixed = var_perm()(i_perm)</tt> is used to determine the index
   * of a free variable in \c var_remove_fixed_to_full() given its index (\c i_perm) being
   * used by the client.
   */
  const IVector& var_perm() const;

  /** \brief Permutes from the original constriant ordering to the current basis selection.
   *
   \verbatim

   equ_perm = [ decomposed equalities | undecomposed equalities ]
              [1..                   r|n-r+1...                n]
   \endverbatim
   *
   * The mapping <tt>j_full = equ_perm()(j_perm)</tt> is used to determine the index
   * of the constriant in c_full given its index \c i_perm being used by the NLP client.
   */
  const IVector& equ_perm() const;

  /** \brief Inverse of \c equ_perm()
   *
   * The mapping <tt>j_perm = inv_equ_perm()(j_full)</tt> is used to determine the index
   * \c j_perm of the constriant \c c being used by the client given the index in c_full.
   */
  const IVector& inv_equ_perm() const;

  // Perform the mapping from a full variable vector to the reduced permuted variable vector
  void var_from_full( DVectorSlice::const_iterator vec_full, DVectorSlice::iterator vec ) const;

  // Perform the mapping from a reduced permuted variable vector the full variable vector
  void var_to_full(DVectorSlice::const_iterator vec, DVectorSlice::iterator vec_full) const;

  // Perform the mapping from c_orig, h_orig, s_orig to the permuted constraint vector c
  void equ_from_full(
    const DVectorSlice   &c_orig
    ,const DVectorSlice  &h_orig
    ,const DVectorSlice  &s_orig
    ,DVectorSlice        *c_full
    ) const;

  //@}

private:

  // ///////////////////////////
  // Private data members

  mutable value_type					f_orig_;    // Computed by subclasses
  mutable DVector						c_orig_;    // ...
  mutable DVector						h_orig_;    // ...
  mutable DVector						Gf_full_;   // ...

  bool initialized_;
  // Flag for if the NLP has has been properly initialized

  bool force_xinit_in_bounds_;
  // Determine if the initial point will be adjusted between bounds

  value_type scale_f_;
  // Set the scaling of the objective function used.

  IVector		var_full_to_fixed_;
  // Holds the indices of variables that are fixed by bounds and those
  // that are not (Length = n_full_).  These partitions are not
  // necessarly sorted in assending order as var_perm and con_perm are.
  //
  //	var_full_to_fixed_ =  [	not fixed by bounds	| fixed by bounds	 ]
  //						  [1 ..				  n_|n_ + 1 ..	  n_full_]
  //

  IVector		inv_var_full_to_fixed_;
  // Inverse of var_full_to_fixed_.  If inv_var_full_to_fixed_(i) > n_ then this variable
  // is fixed between bounds, else inv_var_full_to_fixed_(i) is the indice of the 
  // variable in the unsorted x (not permuted to the current basis).

  IVector		var_perm_;
  // Variable permutations (length = n_) from the vector of unstorted variables not fixed
  // by bounds as defined by var_full_to_fixed_
  //
  // var_perm_	=	[ dependent variables | independent variables ]
  //					[1..                r_|r_+1...              n_]
  //

  IVector		equ_perm_;
  // Equality Constraint permutations (length = m_)
  //
  // equ_perm_	=	[ decomposed equalities | undecomposed equalities ]
  //					[1..                  r_|r_+1...                m_]
  //

  IVector		inv_equ_perm_;
  // Inverse of equ_perm
  //

  mutable DVector x_full_;
  DVector         xinit_full_;
  DVector         xl_full_;
  DVector         xu_full_;
  // The full vector (length = n_full_).  This vector may include
  // slack variables if mI_orig > 0:
  //
  //    [ x_orig; s_orig ]
  //

  perm_fcty_ptr_t            factory_P_var_;
  perm_fcty_ptr_t            factory_P_equ_;
  VectorSpaceSerial          space_x_;
  VectorSpaceSerial          space_c_;
  VectorSpaceSerial          space_c_breve_;
  VectorSpaceSerial          space_h_breve_;
  size_type                  num_bounded_x_;
  VectorMutableDense         xinit_; // Initial point of the shrunken NLP
  VectorMutableDense         xl_;    // Lower bounds of transformed NLP
  VectorMutableDense         xu_;    // Uppers bounds of transformed NLP
  VectorMutableDense         hl_breve_;// Lower bounds for general inequalities of transformed NLP
  VectorMutableDense         hu_breve_;// Uppers bounds for general inequalitiess of transformed NLP
  PermutationSerial          P_var_;
  PermutationSerial          P_equ_;
  size_type              n_orig_;  // Number of variables in the original NLP
  size_type              m_orig_;  // Number of general equality constraints in the original NLP
  size_type              mI_orig_; // Number of general inequality constraints in the original NLP
  size_type              n_full_;  // Number of variables in the transformed NLP (before fixed variables are removed)
  size_type              m_full_;  // Number of general equality constraints in the transformed NLP
  size_type              n_;       // Number of variables in the transformed NLP (with slacks and not fixed by bounds)
  size_type              r_;       // Number of independent equations in the transformed NLP

  int                    basis_selection_num_; // Number of the basis to select next

  // ///////////////////////////
  // Private member functions

  // Get the next basis (or first basis) from the NLP subclass and remove the
  // fixed variables.  Note that this function does not modify var_perm_, equ_perm_
  // or r_.  You must do that yourself by calling assert_and_set_basis.
  bool get_next_basis_remove_fixed(
    IVector* var_perm, IVector* equ_perm, size_type* rank );
  
  // Assert (throw std::length_error, NLPVarReductPerm::InvalidBasis) and set a basis selection
  // If &var_perm == &var_perm_ and/or &equ_perm == &equ_perm_ then the unneeded copy
  // is avoided.
  void assert_and_set_basis(
    const IVector& var_perm, const IVector& equ_perm, size_type rank );

  // Assert that there are bounds on the variables (throw NLP::NoBoundsOnVariables)
  void assert_bounds_on_variables() const;

  // Adjust initial point this->xinit_ to be within bound
  void do_force_xinit_in_bounds();

};	// end class NLPSerialPreprocess

// //////////////////////////////////////////////////
// Inline member functions

// protected

inline
void NLPSerialPreprocess::set_not_initialized()
{
  initialized_ = false;
}

inline
DVectorSlice NLPSerialPreprocess::x_full() const
{
  return x_full_();
}

inline
const NLPSerialPreprocess::ZeroOrderInfoSerial
NLPSerialPreprocess::zero_order_orig_info() const
{
  return ZeroOrderInfoSerial( &f_orig_, &c_orig_, &h_orig_ );
}

inline
const NLPSerialPreprocess::ObjGradInfoSerial
NLPSerialPreprocess::obj_grad_orig_info() const
{
  return ObjGradInfoSerial( &Gf_full_, zero_order_orig_info() );
}

inline
const IVector& NLPSerialPreprocess::var_remove_fixed_to_full() const
{
  return var_full_to_fixed_;
}

inline
const IVector& NLPSerialPreprocess::var_full_to_remove_fixed() const
{
  return inv_var_full_to_fixed_;
}

inline
const IVector& NLPSerialPreprocess::var_perm() const
{
  return var_perm_;
}

inline
const IVector& NLPSerialPreprocess::equ_perm() const
{
  return equ_perm_;
}

inline
const IVector& NLPSerialPreprocess::inv_equ_perm() const
{
  return inv_equ_perm_;
}

}	// end namespace NLPInterfacePack 

#endif // NLP_FULL_TO_REDUCED_H
