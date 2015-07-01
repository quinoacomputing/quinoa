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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H

#include "AbstractLinAlgPack_BasisSystem.hpp"

namespace AbstractLinAlgPack {

/** \brief Interface for setting and selecting a basis from the Jacobian
 * from a set of equations.
 *
 * ToDo: Finish documentation!
 */
class BasisSystemPerm : public BasisSystem {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<Permutation> >   perm_fcty_ptr_t;

  //@}


  /** \brief Required constructor (calls <tt>initialize()</tt>).
   */
  BasisSystemPerm(
    const mat_sym_fcty_ptr_t             &factory_transDtD
    ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
    )
    :BasisSystem(factory_transDtD,factory_S)
      {}

  /** @name Permutation factories */
  //@{

  /** \brief . */
  virtual const perm_fcty_ptr_t   factory_P_var() const = 0;
  /** \brief . */
  virtual const perm_fcty_ptr_t   factory_P_equ() const = 0;

  //@}

  /** @name Basis selection / manipulation */
  //@{

  /** \brief Factor a basis selected by the client.
   *
   * ToDo: Finish documentation!
   */
  virtual void set_basis(
    const Permutation          &P_var
    ,const Range1D             &var_dep
    ,const Permutation         *P_equ
    ,const Range1D             *equ_decomp
    ,const MatrixOp            &Gc
    ,MatrixOpNonsing           *C
    ,MatrixOp                  *D
    ,MatrixOp                  *GcUP
    ,EMatRelations             mat_rel = MATRICES_INDEP_IMPS
    ,std::ostream              *out    = NULL
    ) = 0;

  /** \brief Select a basis.
   *
   * ToDo: Finish documentation!
   */
  virtual void select_basis(
    const Vector               *nu
    ,MatrixOp                  *Gc
    ,Permutation               *P_var
    ,Range1D                   *var_dep
    ,Permutation               *P_equ
    ,Range1D                   *equ_decomp
    ,MatrixOpNonsing           *C
    ,MatrixOp                  *D
    ,MatrixOp                  *GcUP
    ,EMatRelations             mat_rel = MATRICES_INDEP_IMPS
    ,std::ostream              *out    = NULL
    ) = 0;
  
  //@}

private:
  // not defined and not to be called
  BasisSystemPerm();

}; // end class BasisSystemPerm

}  // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H
