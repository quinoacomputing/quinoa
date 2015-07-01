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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_STD_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_STD_H

#include <stdexcept>

#include "ConstrainedOptPack_DecompositionSystemVarReductPerm.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace ConstrainedOptPack {

/** \brief Concreate subclass of \c DecompositionSystemVarReductPerm that uses an
 * aggregate \c DecompostionSystemVarReductImp object.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemVarReductPermStd : public DecompositionSystemVarReductPerm {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<DecompositionSystemVarReductImp>    decomp_sys_imp_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<BasisSystemPerm>                    basis_sys_ptr_t;

  //@}

  /** @name Constructors / initializers */
  //@{

  /// Calls \c this->initialize().
  DecompositionSystemVarReductPermStd(
    const decomp_sys_imp_ptr_t&        decomp_sys_imp  = Teuchos::null
    ,const basis_sys_ptr_t&            basis_sys       = Teuchos::null
    ,bool                              basis_selected  = false
    ,EExplicitImplicit                 D_imp           = MAT_IMP_AUTO
    ,EExplicitImplicit                 Uz_imp          = MAT_IMP_AUTO
    );

  /// Initialize given decomposition system and basis system objects.
  void initialize(
    const decomp_sys_imp_ptr_t&        decomp_sys_imp
    ,const basis_sys_ptr_t&            basis_sys
    ,bool                              basis_selected  = false
    ,EExplicitImplicit                 D_imp           = MAT_IMP_AUTO
    ,EExplicitImplicit                 Uz_imp          = MAT_IMP_AUTO
    );

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  const decomp_sys_imp_ptr_t& decomp_sys_imp() const;
  /** \brief . */
  const basis_sys_ptr_t& basis_sys() const;

  //@}

  /** @name Overridden from DecompositionSystem */
  //@{

  /** \brief . */
  size_type n() const;
  /** \brief . */
  size_type m() const;
  /** \brief . */
  size_type r() const;
  /** \brief . */
  Range1D equ_decomp() const;
  /** \brief . */
  Range1D equ_undecomp() const;
  /** \brief . */
  const VectorSpace::space_ptr_t space_range() const;
  /** \brief . */
  const VectorSpace::space_ptr_t space_null() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_Z() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_Y() const;
  /** \brief . */
  const mat_nonsing_fcty_ptr_t factory_R() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_Uz() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_Uy() const;
  /** \brief . */
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

  /** @name Overridden from DecompositionSystemVarReductPerm */
  //@{

  /** \brief . */
  const perm_fcty_ptr_t   factory_P_var() const;
  /** \brief . */
  const perm_fcty_ptr_t   factory_P_equ() const;
  /** \brief . */
  bool has_basis() const;
  /** \brief . */
  void set_decomp(
    std::ostream          *out
    ,EOutputLevel         olevel
    ,ERunTests            test_what
    ,const Permutation    &P_var
    ,const Range1D        &var_dep
    ,const Permutation    *P_equ
    ,const Range1D        *equ_decomp
    ,const MatrixOp       &Gc
    ,MatrixOp             *Z
    ,MatrixOp             *Y
    ,MatrixOpNonsing      *R
    ,MatrixOp             *Uz
    ,MatrixOp             *Uy
    ,EMatRelations        mat_rel
    );
  /** \brief . */
  void select_decomp(
    std::ostream          *out
    ,EOutputLevel         olevel
    ,ERunTests            test_what
    ,const Vector         *nu
    ,MatrixOp             *Gc
    ,Permutation          *P_var
    ,Range1D              *var_dep
    ,Permutation          *P_equ
    ,Range1D              *equ_decomp
    ,MatrixOp             *Z
    ,MatrixOp             *Y
    ,MatrixOpNonsing      *R
    ,MatrixOp             *Uz
    ,MatrixOp             *Uy
    ,EMatRelations        mat_rel
    );

  //@}

private:

  // /////////////////////////
  // Private data members

  bool                        basis_selected_;  // True if a basis is currently selected
  decomp_sys_imp_ptr_t        decomp_sys_imp_;
  basis_sys_ptr_t             basis_sys_;

  // /////////////////////////
  // Private member functions

  /** \brief . */
  void assert_basis_selected() const;

  // Not defined and not to be called!
  DecompositionSystemVarReductPermStd();
  DecompositionSystemVarReductPermStd(const DecompositionSystemVarReductPermStd&);
  DecompositionSystemVarReductPermStd& operator=(const DecompositionSystemVarReductPermStd&);
  
};	// end class DecompositionSystemVarReductPermStd

// ///////////////////////////////////////
// Inline members

inline
const DecompositionSystemVarReductPermStd::decomp_sys_imp_ptr_t&
DecompositionSystemVarReductPermStd::decomp_sys_imp() const
{
  return decomp_sys_imp_;
}

inline
const DecompositionSystemVarReductPermStd::basis_sys_ptr_t&
DecompositionSystemVarReductPermStd::basis_sys() const
{
  return basis_sys_;
}

}	// end namespace ConstrainedOptPack

#endif // DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_STD_H
