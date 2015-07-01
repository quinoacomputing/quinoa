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

#ifndef	ALAP_DIRECT_SPARSE_SOLVER_MA28_H
#define ALAP_DIRECT_SPARSE_SOLVER_MA28_H

#include <valarray>
#include <vector>
#include <string>

#include "AbstractLinAlgPack_DirectSparseSolverImp.hpp"
#include "AbstractLinAlgPack_MA28Solver.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Concreate sparse solver subclass that uses MA28.
 *
 * ToDo: Finish documentation!
 */
class DirectSparseSolverMA28 : public DirectSparseSolverImp {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  enum EScalingMethod { NO_SCALING, INITIAL_SCALING, SUCCESSIVE_SCALING };

  //@}

  /** @name Control parameters */
  //@{

  /// Pivot tolerance versus sparsity
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, u );

  /// If true, than an estimate of growth of the factors is computed
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, grow );

  /// Drop tolerance for an incomplete factorization
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, tol );

  /// Number of columns to search to reduce fill-in
  STANDARD_MEMBER_COMPOSITION_MEMBERS( index_type, nsrch );

  /// If true, then the largest entry encountered is returned
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, lbig );

  /// If true, then outputs from ma28 are printed to output stream
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_ma28_outputs );

  /// If output_file != "", then output from MA28 is sent to this file.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, output_file_name );

  //@}

  /** @name Constructors/initializers */
  //@{

  /** \brief Default constructor */
  DirectSparseSolverMA28(
    value_type          estimated_fillin_ratio  = 10.0
    ,value_type         u                       = 0.1
    ,bool               grow                    = false
    ,value_type         tol                     = 0.0
    ,index_type         nsrch                   = 4
    ,bool               lbig                    = false
    ,bool               print_ma28_outputs      = false
    ,const std::string& output_file_name        = ""
    );

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

  /** \brief Implements the BasisMatrix object for MA28.
   */
  class BasisMatrixMA28 : public BasisMatrixImp {
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

  }; // end class BasisMatrixMA28

  /** \brief Stores the factorization structure for MA28
   */
  class FactorizationStructureMA28 : public FactorizationStructure {
  public:
    /** \brief . */
    friend class DirectSparseSolverMA28;
    /** \brief . */
    friend class BasisMatrixMA28;
  private:
    // /////////////////////////////////////////
    // Private types
    typedef Teuchos::RCP<MatrixScaling_Strategy>    matrix_scaling_ptr_t;
    // /////////////////////////////////////////
    // Private data members
    mutable MA28_Cpp::MA28Solver ma28_; // Management of common block data
    // Keep a memory of the size of the system to check for consistent usage.
    index_type  m_;     // number of rows (keep for checks on consistancy)
    index_type  n_;     // number of columns ("")
    index_type  max_n_; // max(m_,n_)
    index_type  nz_;    // numner of non-zero elements in unfactorized matrix ("")
    index_type  licn_;  // size of icn_ and a_ (default = 3 * nz_)
    index_type  lirn_;  // size of irn_ (default = 3 * nz_)
    // Control parameters
    value_type	u_; // fill-in vs. stability ratio (default = 0.1)
    EScalingMethod	scaling_; // Scaling method
    // Matrix scaling
    matrix_scaling_ptr_t   matrix_scaling_;
    // Memory for factorization structure
    std::valarray<index_type>  ivect_;
    std::valarray<index_type>  jvect_;
    std::valarray<index_type>  icn_;
    std::valarray<index_type>  ikeep_;
    // Basis matrix selection
    IVector     row_perm_;
    IVector     col_perm_;
    index_type  rank_;
    // /////////////////////////////////////////
    // Private member functions
    /** \brief . */
    FactorizationStructureMA28();
  }; // end class FactorizationStructureMA28

  /** \brief Stores the factorization nonzeros for MA28
   */
  class FactorizationNonzerosMA28 : public FactorizationNonzeros {
  public:
    /** \brief . */
    friend class DirectSparseSolverMA28;
    /** \brief . */
    friend class BasisMatrixMA28;
  private:
    std::valarray<value_type>	a_; // holds the non-zeros of the factorized matrix 'a'
  }; // end class FactorizationNonzerosMA28

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
  // Private types

  /// Enumeration for iflag
  enum E_IFlag {
    SLOW_ITER_CONV                        = -17,
    MAXIT_REACHED                         = -16,
    MA28BD_CALLED_WITH_DROPPED            = -15,
    DUPLICATE_ELEMENTS                    = -14,
    NEW_NONZERO_ELEMENT                   = -13,
    N_OUT_OF_RANGE                        = -11,
    NZ_LE_ZERO                            = -10,
    LICN_LE_NZ                            = -9,
    LIRN_LE_NZ                            = -8,
    ERROR_DURRING_BLOCK_TRI               = -7,
    LICN_AND_LIRN_TOO_SMALL               = -6,
    LICN_TOO_SMALL                        = -5,
    LICN_FAR_TOO_SMALL                    = -4,
    LIRN_TOO_SMALL                        = -3,
    NUMERICALLY_SINGULAR                  = -2,
    STRUCTURALLY_SINGULAR                 = -1,
    SUCCESSFUL_DECOMP_ON_STRUCT_SINGULAR  =  1,
    SUCCESSFUL_DECOMP_ON_NUMER_SINGULAR   =  2
  };

  // /////////////////////////////////
  // Private data members

  value_type                          estimated_fillin_ratio_;
  Teuchos::RCP<std::ostream>  output_file_;
  int                                 file_output_num_;

  // ////////////////////////////////
  // Private member functions

  // Set MA28 control parameters
  void set_ma28_parameters( FactorizationStructureMA28* fs );

  // Print MA28 return parameters
  void print_ma28_outputs(
    bool                               ma28ad_bd
    ,index_type                        iflag
    ,const FactorizationStructureMA28  &fs
    ,const value_type                  w[]
    ,std::ostream                      *out
    );

  // Throw an exception for an iflag error
  void ThrowIFlagException(index_type iflag);

};	// end class DirectSparseSolverMA28 

// ////////////////////////////////////////
// Inline members

}	// end namespace AbstractLinAlgPack 

#endif	// ALAP_DIRECT_SPARSE_SOLVER_MA28_H
