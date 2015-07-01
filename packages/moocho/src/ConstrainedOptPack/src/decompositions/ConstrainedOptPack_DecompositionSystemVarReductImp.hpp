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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_IMP_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_IMP_H

#include "ConstrainedOptPack_DecompositionSystemVarReduct.hpp"
#include "AbstractLinAlgPack_BasisSystemTester.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Specialization node implementation subclass of \c DecompositionSystem for
 * variable reduction decompositions.
 *
 * This interface abstracts a variable reduction decomposition where:
 *
 \verbatim
  
  Gc' = [ C  N ] 
        [ E  F ]

  Z   = [ D ]
        [ I ]

  Uz  = F + E * D

      where:
           C = Gc(var_dep,equ_decomp)'     [nonsingular]
           N = Gc(var_indep,equ_decomp)'
           E = Gc(var_dep,equ_undecomp)'
           F = Gc(var_indep,equ_undecomp)'
           D = -inv(C) * N
 \endverbatim
 *
 * Above, \a C is a <tt>r x r</tt> nonsingular matrix.  Subclasses define
 * how \c Y is defined which in turn determines how \c R and \c Uy are
 * defined.
 *
 * The implementation of this subclass is completly determined by an aggregate
 * <tt>BasisSytem</tt> object.  Since the <tt>BasisSystem</tt> interface does
 * not allow any permutations of the basis to be performed ???.
 *
 * <b>Subclass implementors notes:</b>
 *
 * It is up to subclasses to override \c factory_R() and \c factory_Uy()
 * in order to define the types for \c R and \c Uy
 * respectively.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemVarReductImp : public DecompositionSystemVarReduct {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef DecompositionSystem	                               inherited;
  /** \brief . */
  typedef Teuchos::RCP<const BasisSystem>       basis_sys_ptr_t;

  //@}

  /** @name Constructors / initializers */
  //@{

  /// Set the BasisSystem tester object
  STANDARD_COMPOSITION_MEMBERS( BasisSystemTester, basis_sys_tester );

  /** \brief Construct a variable reduction decomposition.
   *
   * Calls <tt>this->initialize()</tt>.
   *
   * Preconditions:<ul>
   * <li> ???
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->factory_Z().get() != NULL</tt>
   * </ul>
   */
  DecompositionSystemVarReductImp(
    const VectorSpace::space_ptr_t     &space_x
    ,const VectorSpace::space_ptr_t    &space_c
    ,const basis_sys_ptr_t             &basis_sys
    ,const basis_sys_tester_ptr_t      &basis_sys_tester
    ,EExplicitImplicit                 D_imp
    ,EExplicitImplicit                 Uz_imp
    );

  //@}

  /** \brief Initialize.
   *
   * @param  space_x
   *             [in] DVector space for variables \c x.
   * @param  space_c
   *             [in] DVector space for general equalities \c c.
   * @param  basis_sys
   *             [in] The <tt>BasisSystem</tt> object that defines the
   *             variable reduction that the decomposition is based on.
   *             This object must be fully initialized before this
   *             method is called.  The object <tt>*basis_sys</tt> must
   *             not be altered while \c this is still using it.  It is
   *             allowed for <tt>basis_sys.get() == NULL</tt> in which
   *             case \c this will not be fully initialized.
   *
   * Preconditions:<ul>
   * <li> ToDo: Spell these out!
   * </ul>
   *
   * Postconditions:<ul>
   * <li> ToDo: Spell these out!
   * </ul>
   */
  void initialize(
    const VectorSpace::space_ptr_t     &space_x
    ,const VectorSpace::space_ptr_t    &space_c
    ,const basis_sys_ptr_t             &basis_sys
    );

  /** @name Access */
  //@{
  
  /** \brief . */
  const VectorSpace::space_ptr_t& space_x() const;
  /** \brief . */
  const VectorSpace::space_ptr_t& space_c() const;
  /** \brief . */
  const basis_sys_ptr_t& basis_sys() const;

  //@}

  /** @name Basis manipulation */
  //@{

  /** \brief Called by client to uninitialize decomposition matrices in prepairation
   * for selecting a different basis.
   *
   * @param  Z     [in/out] On output, \c Z will have all references to \c C and \c D removed.
   * @param  Y     [in/out] On output, \c Y will have all references to \c C and \c D removed.
   * @param  R     [in/out] On output, \c R will have all references to \c C and \c D removed.
   * @param  Uz    [in/out] If <tt>this->Uz_imp() == MAT_IMP_IMPLICIT</tt> then on output,
   *               \c Uz will have all references to \c C and \c D removed.  If
   *               <tt>this->Uz_imp() == MAT_IMP_EXPLICIT</tt> then \c Uz will be unaltered
   *               and it is expected that the client will initialize it properly before
   *               the next call to \c this->set_basis_matrices().
   * @param  Uy    [in/out] On output, \c Uy will have all references to \c C and \c D removed.
   * @param  C_ptr [out] On output, <tt>C_ptr->get() != NULL</tt> will point to the basis matrix to
   *               be updated by the client before the next call to \c this->set_basis_matrices().
   *               It is guarrenteed that \c *C_ptr->get() will not be referenced by any other entity
   *               so that changing the basis matrix object will not impact other objects in
   *               unexpected ways.
   * @param  D_ptr [out] If <tt>this->D_imp() == MAT_IMP_IMPLICIT</tt> then on output,
   *               <tt>D_ptr == NULL</tt> must be true.  If <tt>this->D_imp() == MAT_IMP_EXPLICIT</tt>
   *               then <tt>D_ptr->get() != NULL</tt> will point to the direct sensitivity matrix
   *               and it is expected that the client will initialize it properly before
   *               the next call to \c this->set_basis_matrices().
   *               It is guarrenteed that \c *D_ptr->get() will not be referenced by any other entity
   *               so that changing the basis matrix object will not impact other objects in
   *               unexpected ways.
   *
   * Preconditions:<ul>
   * <li> [<tt>this->D_imp() == MAT_IMP_IMPLICIT</tt>] <tt>D_ptr == NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>this->D_imp() == MAT_IMP_EXPLICIT</tt>] <tt>D_ptr != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>C_ptr->get() != NULL</tt>
   * <li> [<tt>this->D_imp() == MAT_IMP_EXPLICIT</tt>] <tt>D_ptr->get() != NULL</tt>
   * </ul>
   */
  void get_basis_matrices(
    std::ostream                                      *out
    ,EOutputLevel                                     olevel
    ,ERunTests                                        test_what
    ,MatrixOp                                         *Z
    ,MatrixOp                                         *Y
    ,MatrixOpNonsing                                  *R
    ,MatrixOp                                         *Uz
    ,MatrixOp                                         *Uy
    ,Teuchos::RCP<MatrixOpNonsing>       *C_ptr
    ,Teuchos::RCP<MatrixOp>              *D_ptr
    );

  /** \brief Set updated basis matrices along with a possibly updated basis system object.
   *
   * @param  C_ptr  [in] <tt>C_ptr.get()</tt> points to basis matrix object returned from
   *                \c this->uninitialize_matrices() which must be updated to current basis
   *                for current Jacobian matrices.
   * @param  D_ptr  [in] If <tt>this->D_imp() == MAT_IMP_EXPLICIT</tt>, then \c D_ptr->get()
   *                points to direct sensitivity matrix object returned from
   *                \c this->uninitialize_matrices() which must be updated to current basis
   *                for current Jacobian matrices.  If <tt>this->D_imp() == MAT_IMP_IMPLICIT</tt>
   *                then \c D_ptr must be \c NULL.
   * @param  Uz     [in] If <tt>this->D_imp() == MAT_IMP_EXPLICIT</tt>, then \c Uz points to the
   *                projected sensitivity matrix \c Uz which must be updated for the current
   *                basis for the current Jacobian matrices.  If <tt>this->Uz_imp() == MAT_IMP_IMPLICIT</tt>
   *                then \c Uz must be \c NULL.
   * @param  basis_sys
   *                [in] If the basis system has changed then set <tt>basis_sys.get() != NULL</tt> to
   *                pass in this new basis system object.
   *
   * Preconditions:<ul>
   * <li> <tt>C_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>this->D_imp() == MAT_IMP_IMPLICIT</tt>] <tt>D_ptr == NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>this->D_imp() == MAT_IMP_EXPLICIT</tt>] <tt>D_ptr != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>this->D_imp() == MAT_IMP_EXPLICIT</tt>] <tt>D_ptr->get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * It is expected that immediatly after this method is called that \c this->updated_decomp()
   * will be called to update the rest of the decomposition matrices for these basis matrices.
   */
  void set_basis_matrices(
    std::ostream                                       *out
    ,EOutputLevel                                      olevel
    ,ERunTests                                         test_what
    ,const Teuchos::RCP<MatrixOpNonsing>  &C_ptr
    ,const Teuchos::RCP<MatrixOp>         &D_ptr
    ,MatrixOp                                          *Uz
    ,const basis_sys_ptr_t                             &basis_sys   = Teuchos::null
    );

  /// Get the type of D matrix to be used or is being used (returns MAT_IMP_EXPLICIT or MAT_IMP_IMPLICIT only).
  EExplicitImplicit D_imp_used() const;	

  //@}

  /** @name Overridden from DecompositionSystem */
  //@{

  /// Returns <tt>this->space_x()->dim()</tt>.
  size_type n() const;
  /// Returns <tt>this->space_c()->dim()</tt>.
  size_type m() const;
  /// Returns <tt>this->basis_sys()->equ_decomp().size()</tt>.
  size_type r() const;
  /// Returns <tt>this->space_x()->sub_space(var_dep)</tt>
  const VectorSpace::space_ptr_t space_range() const;
  /// Returns <tt>this->space_x()->sub_space(var_indep)</tt>
  const VectorSpace::space_ptr_t space_null() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_Z() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_Uz() const;
  /** \brief
   *
   * Preconditions:<ul>
   * <li> <tt>this->space_x().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>this->space_c().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> [<tt>test_what == TEST</tt>] <tt>this->basis_sys_tester().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
   * </ul>
   */
  void update_decomp(
    std::ostream          *out
    ,EOutputLevel         olevel
    ,ERunTests            test_what
    ,const MatrixOp       &Gc
    ,MatrixOp             *Z
    ,MatrixOp             *Y
    ,MatrixOpNonsing      *R
    ,MatrixOp             *Uz
    ,MatrixOp             *Uy
    ,EMatRelations        mat_rel
    ) const;
  /** \brief . */
  void print_update_decomp(
    std::ostream& out, const std::string& leading_str ) const;

  //@}

  /** @name Overridden from DecompositionSystemVarReduct */
  //@{

  /** \brief . */
  Range1D var_indep() const;
  /** \brief . */
  Range1D var_dep() const;

  //@}

protected:

  /// Update D_imp_used
  virtual void update_D_imp_used(EExplicitImplicit *D_imp_used) const;

  /** \brief Overridden by subclasses to uninitialized Y, R and Uy then return C if referenced.
   *
   * Note that the returned smart pointer to \c C may have <tt>return.has_ownership() == false</tt>
   * in which case this will not be a shared resource with any other object (at least not
   * on this level).
   *
   * ToDo: Finish documentatation!
`	 */
  virtual	mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t	uninitialize_matrices(
    std::ostream                     *out
    ,EOutputLevel                    olevel
    ,MatrixOp                        *Y
    ,MatrixOpNonsing                 *R
    ,MatrixOp                        *Uy
    ) const = 0;

  /** \brief Overridden by subclasses to initialize Y, R and Uy given C and D.
   *
    * If C_ptr.has_ownership() == false, then the subclass implementation of this
   * method will use clone_mwons() to clone it so that the output matrices are
   * independent.
   *
   * ToDo: Finish documentatation!
`	 */
  virtual void initialize_matrices(
    std::ostream                                           *out
    ,EOutputLevel                                          olevel
    ,const mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t &C_ptr
    ,const mat_fcty_ptr_t::element_type::obj_ptr_t         &D_ptr
    ,MatrixOp                                              *Y
    ,MatrixOpNonsing                                       *R
    ,MatrixOp                                              *Uy
    ,EMatRelations                                         mat_rel
    ) const = 0;

  /** \brief Print the sub-algorithm by which the matrices Y, R, Uy and Uy are updated.
   *
   * ToDo: Finish documentatation!
   */
  virtual void print_update_matrices(
    std::ostream& out, const std::string& leading_str ) const = 0;

private:

  // //////////////////////////////////
  // Private data members
  
#ifdef DOXYGEN_COMPILE
  AbstractLinAlgPack::BasisSystem       *basis_sys;
  VectorSpace                           *space_x;
  VectorSpace                           *space_c;
  VectorSpace                           *space_range;
  VectorSpace                           *space_null;
#else
  basis_sys_ptr_t                       basis_sys_;
  VectorSpace::space_ptr_t              space_x_;
  VectorSpace::space_ptr_t              space_c_;
  VectorSpace::space_ptr_t              space_range_;
  VectorSpace::space_ptr_t              space_null_;
  mutable Teuchos::RCP<MatrixOpNonsing>  C_ptr_;
  mutable Teuchos::RCP<MatrixOp>             D_ptr_;
  mutable EExplicitImplicit                                   D_imp_used_;
#endif
  // //////////////////////////////////
  // Private member functions

  /// Allocate a new D_ptr matrix
  void alloc_new_D_matrix( 
    std::ostream                             *out
    ,EOutputLevel                            olevel
    ,Teuchos::RCP<MatrixOp> *D_ptr
    ) const;
  
  // not defined and not to be called!
  DecompositionSystemVarReductImp(const DecompositionSystemVarReductImp&);
  DecompositionSystemVarReductImp& operator=(const DecompositionSystemVarReductImp&);
  
};	// end class DecompositionSystemVarReductImp

// //////////////////////////////////////////
// Inline members

inline
const VectorSpace::space_ptr_t&
DecompositionSystemVarReductImp::space_x() const
{
  return space_x_;
}

inline
const VectorSpace::space_ptr_t&
DecompositionSystemVarReductImp::space_c() const
{
  return space_c_;
}

inline
const DecompositionSystemVarReductImp::basis_sys_ptr_t&
DecompositionSystemVarReductImp::basis_sys() const
{
  return basis_sys_;
}

inline
DecompositionSystemVarReductImp::EExplicitImplicit
DecompositionSystemVarReductImp::D_imp_used() const
{
  update_D_imp_used(&D_imp_used_);
  return D_imp_used_;
}

}	// end namespace ConstrainedOptPack

#endif	// DECOMPOSITION_SYSTEM_VAR_REDUCT_IMP_H
