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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H

#include <stdexcept>

#include "ConstrainedOptPack_DecompositionSystemVarReduct.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace ConstrainedOptPack {

/** \brief Specialization interface of \c DecompositonSystem that allows basis permutations.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemVarReductPerm : public DecompositionSystemVarReduct {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<Permutation> >         perm_fcty_ptr_t;

  //@}

  /** @name Constructors / initializers */
  //@{

  /** \brief . */
  DecompositionSystemVarReductPerm(
    EExplicitImplicit     D_imp    = MAT_IMP_AUTO
    ,EExplicitImplicit    Uz_imp   = MAT_IMP_AUTO
    )
    :DecompositionSystemVarReduct(D_imp,Uz_imp)
  {}

  //@}

  /** @name Permutation factories */
  //@{

  /** \brief . */
  virtual const perm_fcty_ptr_t   factory_P_var() const = 0;
  /** \brief . */
  virtual const perm_fcty_ptr_t   factory_P_equ() const = 0;

  //@}

  /** @name Setting or selecting a new decomposition */
  //@{

  /** \brief Query to see if a current basis is already selected.
   */
  virtual bool has_basis() const = 0;

  /** \brief Client selects the basis for <tt>Gc(:,con_decomp)'</tt>.
   *
   * ToDo: Finish documentation!
   */
  virtual void set_decomp(
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
    ,EMatRelations        mat_rel = MATRICES_INDEP_IMPS
    ) = 0;
  
  /** \brief Client asks decompostion system object to select the basis for <tt>Gc(:,con_decomp)'</tt>.
   *
   * ToDo: Finish documentation!
   */
  virtual void select_decomp(
    std::ostream              *out
    ,EOutputLevel             olevel
    ,ERunTests                test_what
    ,const Vector             *nu
    ,MatrixOp                 *Gc
    ,Permutation              *P_var
    ,Range1D                  *var_dep
    ,Permutation              *P_equ
    ,Range1D                  *equ_decomp
    ,MatrixOp                 *Z
    ,MatrixOp                 *Y
    ,MatrixOpNonsing          *R
    ,MatrixOp                 *Uz
    ,MatrixOp                 *Uy
    ,EMatRelations            mat_rel = MATRICES_INDEP_IMPS
    ) = 0;

  //@}
  
};	// end class DecompositionSystemVarReductPerm

}	// end namespace ConstrainedOptPack

#endif // DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H
