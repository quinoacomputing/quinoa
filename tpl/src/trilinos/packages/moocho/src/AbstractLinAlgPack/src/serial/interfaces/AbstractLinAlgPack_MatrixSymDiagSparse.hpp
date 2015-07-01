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
#ifndef SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H
#define SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H

#include "AbstractLinAlgPack_MatrixSymOpSerial.hpp"
#include "AbstractLinAlgPack_MatrixConvertToSparse.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all serial symmetric diagonal matrices with
 * significant zeros along the diagonal.
 */
class MatrixSymDiagSparse
  : virtual public MatrixSymOpSerial
  , virtual public MatrixConvertToSparse
{
public:

  /** \brief <<std member comp>> members for how many updates to compute
    * at once in the operation M_MtMtM(....).
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_updates_at_once );

  /** \brief The default value of num_updates_at_once == 0 is set to allow
    * this class to determine the appropriate size internally.
    */
  MatrixSymDiagSparse();

  /** @name To be overridden by subclass */
  //@{

  /// Give access to the sparse diagonal
  virtual const SpVectorSlice diag() const = 0;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /** \brief . */
  size_type rows() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  std::ostream& output(std::ostream& out) const;

  //@}

  /** @name Overridden from MatrixOpSerial */
  //@{

  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;

  //@}

  /** @name Overridden from MatrixSymOpSerial */
  //@{

  /** \brief Computes the dense symmetric matrix B += a*op(A')*M*op(A).
   *
   * This matrix is computed using a set of rank-1 updates.
   *
   * Runtime ~ O( (m^2)*nz )
   *
   * Storage ~ O( num_updates_at_once * m )
   *
   * Where:<ul>
   * <li> <tt>n = A.rows() == this->rows()</tt>
   * <li> <tt>m = A.cols()</tt>
   * <li> <tt>nz = this->diag().nz()</tt>
   * </ul>
   *
   * Note that a necessary condition for \c B to be full rank is for
   * <tt>nz >= m</tt>.
   *
   * Also note that this default implementation is only for nonnegative
   * diagonal entries.
   */
  void Mp_StMtMtM( DMatrixSliceSym* sym_lhs, value_type alpha
    , EMatRhsPlaceHolder dummy_place_holder
    , const MatrixOpSerial& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
    , value_type beta ) const;

  //@}

  /** @name Overridden from MatrixConvertToSparse */
  //@{

  /** \brief . */
  index_type num_nonzeros(
    EExtractRegion        extract_region
    ,EElementUniqueness   element_uniqueness
    ) const;
  /** \brief . */
  void coor_extract_nonzeros(
    EExtractRegion                extract_region
    ,EElementUniqueness           element_uniqueness
    ,const index_type             len_Aval
    ,value_type                   Aval[]
    ,const index_type             len_Aij
    ,index_type                   Arow[]
    ,index_type                   Acol[]
    ,const index_type             row_offset
    ,const index_type             col_offset
    ) const;

  //@}

};	// end class MatrixSymDiagSparse

}	// end namespace AbstractLinAlgPack

#endif	// SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H
