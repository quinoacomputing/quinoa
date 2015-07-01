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

#ifndef MATRIX_CONVERT_TO_SPARSE_ENCAP_H
#define MATRIX_CONVERT_TO_SPARSE_ENCAP_H

#include "AbstractLinAlgPack_MatrixConvertToSparse.hpp"
#include "Teuchos_RCP.hpp"

namespace AbstractLinAlgPack {

/** \brief Sparse conversion subclass based on views of a \c MatrixExtractSparseElements object.
 *
 * ToDo:Finish documentation!
 */
class MatrixConvertToSparseEncap
  : virtual public MatrixConvertToSparse
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<const MatrixExtractSparseElements>  mese_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const IVector>                      i_vector_ptr_t;

  //@}

  /** @name Constructors/initializers */
  //@{

  /** \brief Construct to uninitialized.
   */
  MatrixConvertToSparseEncap();

  /** \brief Calls \c this->initialize().
   */
  MatrixConvertToSparseEncap(
    const mese_ptr_t           &mese
    ,const i_vector_ptr_t      &inv_row_perm
    ,const Range1D             &row_rng
    ,const i_vector_ptr_t      &inv_col_perm
    ,const Range1D             &col_rng
    ,const BLAS_Cpp::Transp    mese_trans
    ,const value_type          alpha = 1.0
    );

  /** \brief Initialize a permuted view of a sparse matrix.
   *
   * <tt>A = alpha * op( (P'*B*Q)(row_rng,col_rng) )</tt>
   *
   * ToDo: Finish documentation!
   */
  void initialize(
    const mese_ptr_t           &mese
    ,const i_vector_ptr_t      &inv_row_perm
    ,const Range1D             &row_rng
    ,const i_vector_ptr_t      &inv_col_perm
    ,const Range1D             &col_rng
    ,const BLAS_Cpp::Transp    mese_trans
    ,const value_type          alpha = 1.0
    );

  /** \brief Set uninitialized.
   *
   * ToDo: Finish documentation!
   */
  void set_uninitialized();

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  const mese_ptr_t& mese() const;
  /** \brief . */
  const i_vector_ptr_t& inv_row_perm() const;
  /** \brief . */
  const Range1D& row_rng() const;
  /** \brief . */
  const i_vector_ptr_t& inv_col_perm() const;
  /** \brief . */
  const Range1D& col_rng() const;
  /** \brief . */
  const BLAS_Cpp::Transp mese_trans() const;
  /** \brief . */
  const value_type alpha() const;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  const VectorSpace& space_rows() const;
  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;

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

private:

  typedef Teuchos::RCP<const VectorSpace> space_ptr_t;  

#ifdef DOXYGEN_COMPILE
  const MatrixExtractSparseElements    *mese;
  const DenseLinAlgPack::IVector            *inv_row_perm;
  Range1D                              row_rng;
  const DenseLinAlgPack::IVector            *inv_col_perm;
  Range1D                              col_rng;
#else
  mese_ptr_t          mese_;
  BLAS_Cpp::Transp    mese_trans_;
  value_type          alpha_;
  Range1D             row_rng_;
  Range1D             col_rng_;
  i_vector_ptr_t      inv_row_perm_;
  i_vector_ptr_t      inv_col_perm_;
  size_type           nz_full_;
  space_ptr_t         space_cols_;
  space_ptr_t         space_rows_;
#endif

};	// end class MatrixConvertToSparseEncap

// /////////////////////////////////////////////
// Inline members

// Access

inline
const MatrixConvertToSparseEncap::mese_ptr_t&
MatrixConvertToSparseEncap::mese() const
{
  return mese_;
}

inline
const MatrixConvertToSparseEncap::i_vector_ptr_t&
MatrixConvertToSparseEncap::inv_row_perm() const
{
  return inv_row_perm_;
}

inline
const Range1D& MatrixConvertToSparseEncap::row_rng() const
{
  return row_rng_;
}

inline
const MatrixConvertToSparseEncap::i_vector_ptr_t&
MatrixConvertToSparseEncap::inv_col_perm() const
{
  return inv_col_perm_;
}

inline
const Range1D& MatrixConvertToSparseEncap::col_rng() const
{
  return col_rng_;
}

inline
const BLAS_Cpp::Transp
MatrixConvertToSparseEncap::mese_trans() const
{
  return mese_trans_;
}

inline
const value_type MatrixConvertToSparseEncap::alpha() const
{
  return alpha_;
}
  
}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_CONVERT_TO_SPARSE_ENCAP_H
