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

#ifndef NLP_VAR_REDUCT_PERM_H
#define NLP_VAR_REDUCT_PERM_H

#include "NLPInterfacePack_NLP.hpp"

namespace NLPInterfacePack {

/** \brief NLP interface class that adds variable and constriant permutations for variable
 * reduction basis selections.
 *
 * This class adds basis selection and manipulation.  This functionality is needed by many
 * optimization algorithms that categorize variables and constraints into specific sets
 * according to a basis selection.  To understand what sets these are, consider the following
 * equality constraints (from the \c NLP interface).
 \verbatim
  
            c(x)  = 0,    c(x) <: R^n -> R^m
 \endverbatim
 * This interface allows \a x and \a c(x) to be partitioned into different
 * sets.  The variables \a x are partitioned into a dependent set \c x(var_dep)
 * and an independent set \c x(var_dep) by the permutation \c P_var.  The equality
 * constraints \a c(x) are partitioned into decomposed \c c(equ_decomp) and undecomposed
 * \c c(equ_undecomp) sets by the permutation \c P_equ.  These permutations
 * permute from an original order to a new ordering.  For example:
 \verbatim

  Original Ordering      Permutation to new ordering                Partitioning
  -----------------    -------------------------------   -------------------------------------
       x_orig          P_var.permute(trans,x_orig,x)     -> x(var_dep),      x(var_indep)
       c_orig          P_equ.permute(trans,c_orig,c)     -> c(equ_decomp),   c(equ_undecomp)
 \endverbatim
 * Because of this partitioning, it is expected that the following vector sub-spaces will be
 * non-null: <tt>space_x()->sub_space(var_indep)</tt>, <tt>space_x()->sub_space(var_dep)</tt>,
 * <tt>space_c()->sub_space(equ_decomp)</tt>, <tt>space_c()->sub_space(equ_undecomp)</tt>.
 * Other subspaces may be non-null also but these are the only ones that are required
 * to be.
 *
 * After initialization, the %NLP subclass will be initialized to the first basis.
 * This basis may be the original ordering if \c P_var and \c P_equ
 * all return <tt>xxx_perm.is_identity()</tt>.
 * If the concrete %NLP is selecting the basis (<tt>nlp_selects_basis() == true</tt>) this
 * basis will be that first basis.  The first time that \c this->get_next_basis() is called
 * it will return this initial basis (which may not be the original ordering).
 *
 * The client can always see what this first basis is by calling <tt>this->get_basis()</tt>.
 * If a basis goes singular the client can request other basis selections from the %NLP by
 * calling <tt>this->get_next_basis()</tt> (which will return true if more basis selections
 * are available).
 * The client can also select a basis itself and then set that basis by calling
 * <tt>this->set_basis()</tt> to force the use of that basis selection.
 * In this way a valid basis is automatically selected after initialization so that
 * clients using another interface (\c NLP, \c NLPFirstOrder, or 
 * \c NLPSecondOrder) will be able to use the %NLP object without even knowing
 * about a basis selection.
 *
 * Below are some obviouls assertions about the basis selection: <ul>
 * <li> <tt>P_var.space().dim() == this->n()</tt>  (throw \c std::length_error)
 * <li> <tt>P_equ.space().dim() == this->m()</tt>  (throw \c std::length_error)
 * <li> <tt>var_dep.size() <= min( this->m() , this->n() )</tt> (throw \c InvalidBasis)
 * <li> <tt>var_dep.size() == equ_decomp.size()</tt> (throw \c InvalidBasis)
 * <li> Other obvious assertions on the basis selection?
 * </ul>
 */
class NLPVarReductPerm : virtual public NLP
{
public: 

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<Permutation> >         perm_fcty_ptr_t;

  /// Thrown if an invalid basis selection is made
  class InvalidBasis : public std::logic_error
  {public: InvalidBasis(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}	

  /** @name Abstract factories for Permutation objects */
  //@{

  /** \brief . */
  virtual const perm_fcty_ptr_t factory_P_var() const = 0;
  /** \brief . */
  virtual const perm_fcty_ptr_t factory_P_equ() const = 0;

  //@}

  /** @name Return ranges for the partitioning of variables and constraints */
  //@{

  /** \brief . */
  virtual Range1D var_dep() const = 0;
  /** \brief . */
  virtual Range1D var_indep() const = 0;
  /** \brief . */
  virtual Range1D equ_decomp() const = 0;
  /** \brief . */
  virtual Range1D equ_undecomp() const = 0;

  //@}

  /** @name Basis manipulation functions */
  //@{
  
  /// Returns true if the NLP can suggest one or more basis selections.
  virtual bool nlp_selects_basis() const = 0;

  /** \brief Returns the next basis the %NLP has and sets the NLP to the returned basis.
   *
   * @param  P_var        [out] Variable permutations defined as <tt>P_var'*x_old -> x_new = [ x(var_dep); x(var_indep) ]</tt>
   * @param  var_dep      [out] Range of dependent variables in <tt>x_new(var_dep)</tt>
   * @param  P_equ        [out] Equality constraint permutations defined as <tt>P_equ'*c_old -> c_new = [ c(equ_decomp); c(equ_undecomp) ]</tt>
   * @param  equ_decomp   [out] Range of decomposed equalities in <tt>c_new(equ_decomp)</tt>
   *
   * Postconditions:  The %NLP is set to the basis returned in the arguments and <tt>this->get_basis()</tt> will return the same basis.
   *
   * This member returns \c true if the %NLP has another basis to select, and is \c false if not.
   * If \c false is returned the client has the option of selecting another basis on its own
   * and passing it to the %NLP by calling <tt>this->set_basis()</tt>.
   */
  virtual bool get_next_basis(
    Permutation*  P_var,   Range1D* var_dep
    ,Permutation* P_equ,   Range1D* equ_decomp
    ) = 0;
  
  /** \brief Sets the basis the that the %NLP will use to permute the problem.
   *
   * @param  P_var        [in] Variable permutations defined as <tt>P_var'*x_old -> x_new = [ x(var_dep); x(var_indep) ]</tt>
   * @param  var_dep      [in] Range of dependent variables in <tt>x_new(var_dep)</tt>
   * @param  P_equ        [in] Equality constraint permutations defined as <tt>P_equ'*c_old -> c_new = [ c(equ_decomp); c(equ_undecomp) ]</tt>
   * @param  equ_decomp   [in] Range of decomposed equalities in <tt>c_new(equ_decomp)</tt>
   *
   * Preconditions: The input basis meets the <A ref=BasisAssertions>basis assertions</A> stated above or an \c InvalidBasis exceptin
   * is thrown.
   *
   * Postconditions:  The %NLP is set to the basis given in the arguments and <tt>this->get_basis()</tt> will return this same basis.
   */
  virtual void set_basis(
    const Permutation   &P_var,   const Range1D  &var_dep
    ,const Permutation  *P_equ,   const Range1D  *equ_decomp
    ) =  0;

  /** \brief Returns the basis selection currently being used by the NLP.
   *
   * @param  P_var        [out] Variable permutations defined as <tt>P_var'*x_old -> x_new = [ x(var_dep); x(var_indep) ]</tt>
   * @param  var_dep      [out] Range of dependent variables in <tt>x_new(var_dep)</tt>
   * @param  P_equ        [out] Equality constraint permutations defined as <tt>P_equ'*c_old -> c_new = [ c(equ_decomp); c(equ_undecomp) ]</tt>
   * @param  equ_decomp   [out] Range of decomposed equalities in <tt>c_new(equ_decomp)</tt>
   */
  virtual void get_basis(
    Permutation*  P_var,   Range1D* var_dep
    ,Permutation* P_equ,   Range1D* equ_decomp
    ) const = 0;
  
  //@}

private:

#ifdef DOXYGEN_COMPILE
  Teuchos::AbstractFactory<AbstractLinAlgPack::Permutation>    *factory_P_var;
  Teuchos::AbstractFactory<AbstractLinAlgPack::Permutation>    *factory_P_equ;
#endif	
  
}; // end class NLPVarReductPerm

} // end namespace NLPInterfacePack 

#endif // NLP_VAR_REDUCT_PERM_H
