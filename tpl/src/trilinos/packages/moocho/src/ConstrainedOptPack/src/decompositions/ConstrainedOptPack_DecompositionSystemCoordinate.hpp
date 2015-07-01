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

#ifndef DECOMPOSITION_SYSTEM_COORDINATE_H
#define DECOMPOSITION_SYSTEM_COORDINATE_H

#include "ConstrainedOptPack_DecompositionSystemVarReductImp.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Coordinate variable reduction subclass.
 *
 * This is the interface for the coordinate variable reduction decomposition
 * where:
 \verbatim

  Y = [ I ]   (class MatrixIdentConcatStd with MatrixZero)
      [ 0 ]

  R = Gc(:,con_decomp)'*Y = [ C   N ] * [ I ] = C
                                        [ 0 ]

  Uy = Gc(:,con_undecomp)'*Y = [ E  F ] * [ I ] = E
                                          [ 0 ]

 \endverbatim
 * The solution of the
 *
 * For now the copy constructor and the assignment operator are not defined.
 */
class DecompositionSystemCoordinate : public DecompositionSystemVarReductImp {
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief . */
  DecompositionSystemCoordinate(
    const VectorSpace::space_ptr_t     &space_x          = Teuchos::null
    ,const VectorSpace::space_ptr_t    &space_c          = Teuchos::null
    ,const basis_sys_ptr_t             &basis_sys        = Teuchos::null
    ,const basis_sys_tester_ptr_t      &basis_sys_tester = Teuchos::null
    ,EExplicitImplicit                 D_imp             = MAT_IMP_AUTO
    ,EExplicitImplicit                 Uz_imp            = MAT_IMP_AUTO
    );

  //@}

  /** @name Overridden from DecompositionSystem */
  //@{

  /** \brief . */
  const mat_fcty_ptr_t factory_Y() const;
  /** \brief . */
  const mat_nonsing_fcty_ptr_t factory_R() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_Uy() const;

  //@}

protected:

  /** @name Overridden from DecompositionSystemVarReductImp */
  //@{

  /** \brief . */
  mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t	uninitialize_matrices(
    std::ostream                                       *out
    ,EOutputLevel                                      olevel
    ,MatrixOp                                          *Y
    ,MatrixOpNonsing                                   *R
    ,MatrixOp                                          *Uy
    ) const;
  /** \brief . */
  void initialize_matrices(
    std::ostream                                           *out
    ,EOutputLevel                                          olevel
    ,const mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t &C
    ,const mat_fcty_ptr_t::element_type::obj_ptr_t         &D
    ,MatrixOp                                              *Y
    ,MatrixOpNonsing                                       *R
    ,MatrixOp                                              *Uy
    ,EMatRelations                                         mat_rel
    ) const;
  /** \brief . */
  void print_update_matrices(
    std::ostream& out, const std::string& leading_str ) const;

  //@}

};	// end class DecompositionSystemCoordinate

}	// end namespace ConstrainedOptPack

#endif	// DECOMPOSITION_SYSTEM_COORDINATE_H
