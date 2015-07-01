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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_H

#include "ConstrainedOptPack_DecompositionSystem.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"

namespace ConstrainedOptPack {

/** \brief Specialization of \c DecompositionSystem for variable reduction decompositions.
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
           C = Gc(var_dep,con_decomp)'     [nonsingular]
           N = Gc(var_indep,con_decomp)'
           E = Gc(var_dep,con_undecomp)'
           F = Gc(var_indep,con_undecomp)'
           D = -inv(C) * N
 \endverbatim
 *
 * This interface simply allows clients to determine how \c D and \c Uz
 * are implemented (implicitly or explicity).
 */
class DecompositionSystemVarReduct : public DecompositionSystem {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  enum EExplicitImplicit {
    MAT_IMP_EXPLICIT
    ,MAT_IMP_IMPLICIT
    ,MAT_IMP_AUTO
  };

  //@}

  /** @name Matrix representations */
  //@{

  /// Set whether to use explicit or implicit <tt>D = -inv(C)*N</tt> matrix.
  STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,D_imp);
  /// Set whether to use explicit or implicit <tt>Uz = F + E * D</tt> matrix.
  STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,Uz_imp);

    // ToDo: The above could be implemented as pure virtual funtions if needed later!

  //@}
    
  /** @name Constructors / initializers */
  //@{

  /** \brief . */
  DecompositionSystemVarReduct(
    EExplicitImplicit     D_imp    = MAT_IMP_AUTO
    ,EExplicitImplicit    Uz_imp   = MAT_IMP_AUTO
    )
    :D_imp_(D_imp), Uz_imp_(Uz_imp)
  {}

  //@}

  /** @name Variable partitions. */
  //@{

  /** \brief . */
  virtual Range1D var_indep() const = 0;
  /** \brief . */
  virtual Range1D var_dep() const = 0;

  //@}

private:

  // not defined and not to be called!
  DecompositionSystemVarReduct(const DecompositionSystemVarReduct&);
  DecompositionSystemVarReduct& operator=(const DecompositionSystemVarReduct&);

};	// end class DecompositionSystemVarReduct

}	// end namespace ConstrainedOptPack

#endif	// DECOMPOSITION_SYSTEM_VAR_REDUCT_H
