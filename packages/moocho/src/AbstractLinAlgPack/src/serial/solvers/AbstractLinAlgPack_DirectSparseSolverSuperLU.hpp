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

#ifdef SPARSE_SOLVER_PACK_USE_SUPERLU

#ifndef	DIRECT_SPARSE_SOLVER_SUPERLU_H
#define DIRECT_SPARSE_SOLVER_SUPERLU_H

#include <valarray>
#include <vector>
#include <string>

#include "AbstractLinAlgPack_DirectSparseSolverImp.hpp"
#include "AbstractLinAlgPack_SuperLUSolver.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Concreate sparse solver subclass that uses SuperLU.
 *
 * ToDo: Finish documentation!
 */
class DirectSparseSolverSuperLU : public DirectSparseSolverImp {
public:

  /** @name Control parameters */
  //@{

  // ToDo: Fill these in!

  //@}

  /** @name Constructors/initializers */
  //@{

  /** \brief Default constructor */
  DirectSparseSolverSuperLU();

  //@}

  /** @name Overridden from DirectSparseSolver */
  //@{

  /** \brief . */
  const basis_matrix_factory_ptr_t basis_matrix_factory() const;
  /** \brief . */
  void estimated_fillin_ratio( value_type estimated_fillin_ratio );

  //@}

protected:

  /** @name Protected types */
  //@{

  /** \brief Implements the BasisMatrix object for SuperLU.
   */
  class BasisMatrixSuperLU : public BasisMatrixImp {
  public:

    /** @name Overridden from BasisMatrixImp */
    //@{

    /** \brief . */
    Teuchos::RCP<BasisMatrixImp> create_matrix() const;
    /** \brief . */
    void V_InvMtV(
      VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
      ,const Vector& v_rhs2) const ;
    
    //@}

  }; // end class BasisMatrixSuperLU

  /** \brief Stores the factorization structure for SuperLU
   */
  class FactorizationStructureSuperLU : public FactorizationStructure {
  public:
    friend class DirectSparseSolverSuperLU;
    friend class BasisMatrixSuperLU;
  private:
    Teuchos::RCP<SuperLUPack::SuperLUSolver>
      superlu_solver_;
    Teuchos::RCP<SuperLUPack::SuperLUSolver::FactorizationStructure>
      fact_struct_;
    FactorizationStructureSuperLU();
  }; // end class FactorizationStructureSuperLU

  /** \brief Stores the factorization nonzeros for SuperLU
   */
  class FactorizationNonzerosSuperLU : public FactorizationNonzeros {
  public:
    friend class DirectSparseSolverSuperLU;
    friend class BasisMatrixSuperLU;
  private:
    Teuchos::RCP<SuperLUPack::SuperLUSolver::FactorizationNonzeros>
      fact_nonzeros_;
  }; // end class FactorizationNonzerosSuperLU

  //@}

  /** @name Overridden from DirectSparseSolverImp */
  //@{

  /** \brief . */
  const Teuchos::RCP<FactorizationStructure> create_fact_struc() const;
  /** \brief . */
  const Teuchos::RCP<FactorizationNonzeros> create_fact_nonzeros() const;
  /** \brief . */
  void imp_analyze_and_factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,FactorizationStructure                         *fact_struc
    ,FactorizationNonzeros                          *fact_nonzeros
    ,DenseLinAlgPack::IVector                            *row_perm
    ,DenseLinAlgPack::IVector                            *col_perm
    ,size_type                                      *rank
    ,std::ostream                                   *out
    );
  /** \brief . */
  void imp_factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,const FactorizationStructure                   &fact_struc
    ,FactorizationNonzeros                          *fact_nonzeros
    ,std::ostream                                   *out
    );

  //@}

private:

  // /////////////////////////////////
  // Private data members

  // ////////////////////////////////
  // Private member functions

};	// end class DirectSparseSolverSuperLU 

// ////////////////////////////////////////
// Inline members

}	// end namespace AbstractLinAlgPack 

#endif	// DIRECT_SPARSE_SOLVER_SUPERLU_H

#endif // SPARSE_SOLVER_PACK_USE_SUPERLU
