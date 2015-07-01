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

#ifndef ALAP_MULTI_VECTOR_MUTABLE_H
#define ALAP_MULTI_VECTOR_MUTABLE_H

#include "AbstractLinAlgPack_MultiVector.hpp"

namespace AbstractLinAlgPack {

/** \brief Interface for a collection of mutable vectors (multi-vector, matrix).
 *
 * This interface extends the \c MutiVector interface an allows mutable access to
 * the constituent vectors.
 *
 * These vectors allow the modification of the matrix row by row, column by column,
 * and/or diagonal by diagonal.  Each of the views is transient and should be used
 * and discarded quickly.
 *
 * Note that the underlying matrix is only guaranteed to be modified after the smart reference
 * counted pointer returned from these methods is destoryed.  For example, consider the following code:
 \code

 void f( MultiVectorMutable* M, index_type i )
 {
  MultiVectorMutable::vec_mut_ptr_t
      row_i =M->row(i);
  *row_i = 0.0;
  // The underlying matrix may not be modified at this point.
  row_i = NULL;
  // Now the underlying matrix is guaranteed to be modified and
  // we can assume this in the following code.
  ...
 }
 \endcode
 * Default implementations of the const access methods \c row() \c col()
 * and \c diag() from \c MultiVector call the non-const methods defined
 * here and cast the pointers.
 *
 * Many of the default implementations of the linear algebra operations in
 * \c MatrixOp and the other matrix interfaces rely on the left hand side
 * matrix objects supporting the \c MultiVectorMutable interface.
 */
class MultiVectorMutable : virtual public MultiVector {
public:
  
  /** \brief . */
  using MultiVector::col;
  /** \brief . */
  using MultiVector::row;
  /** \brief . */
  using MultiVector::diag;

  /** \brief . */
  typedef Teuchos::RCP<VectorMutable>       vec_mut_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<MultiVectorMutable>  multi_vec_mut_ptr_t;

  /** @name Clone */
  //@{

  /** \brief Clone the non-const multi-vector object.
   *
   * The default implementation creates a new multi-vector
   * and then copies the values.
   */
  virtual multi_vec_mut_ptr_t mv_clone();

  //@}

  /** @name Provide mutable row, column and/or diagonal access */
  //@{

  /** \brief Get a mutable column vector.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->access_by() & COL_ACCESS</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>return.get() != NULL</tt>] <tt>space_cols().is_compatible(return->space()) == true</tt>
   * </ul>
   *
   * ToDo: Finish documentation!
   */
  virtual vec_mut_ptr_t col(index_type j) = 0;
  /** \brief Get a mutable row vector.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->access_by() & ROW_ACCESS</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>return.get() != NULL</tt>] <tt>space_rows().is_compatible(return->space()) == true</tt>
   * </ul>
   *
   * ToDo: Finish documentation!
   */
  virtual vec_mut_ptr_t row(index_type i) = 0;
  /** \brief Get a mutable diagonal vector.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->access_by() & DIAG_ACCESS</tt>] <tt>return.get() != NULL</tt>
   * </ul>
   *
   * ToDo: Finish documentation!
   */
  virtual vec_mut_ptr_t diag(int k) = 0;

  //@}

  /** @name Sub-view methods */
  //@{

  /** \brief Returns a mutable sub-view of the multi vector.
   *
   * ToDo: Finish documentation!
   *
   * The default implementation returns a \c MultiVectorMutableSubView object for
   * any valid arbitary sub-view.
   */
  virtual multi_vec_mut_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng);
  
  /** \brief Inlined implementation calls <tt>this->mv_sub_view(Range1D(rl,ru),Range1D(cl,cu))</tt>.
   */
  multi_vec_mut_ptr_t mv_sub_view(
    const index_type& rl, const index_type& ru
    ,const index_type& cl, const index_type& cu
    );

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  mat_mut_ptr_t clone();
  /** \brief . */
  void zero_out();
  /** \brief . */
  void Mt_S( value_type alpha );
  /** \brief . */
  MatrixOp& operator=(const MatrixOp& mwo_rhs);
  /** \brief . */
  bool Mp_StM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs
    ) const;
  /** \brief . */
  bool Mp_StM(
    value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs
    );

  //@}

  /** @name Overridden from MultiVector */
  //@{

  /** \brief . */
  multi_vec_ptr_t mv_clone() const;
  /** \brief . */
  vec_ptr_t col(index_type j) const;
  /** \brief . */
  vec_ptr_t row(index_type i) const;
  /** \brief . */
  vec_ptr_t diag(int k) const;
  /** \brief . */
  multi_vec_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const;

  //@}

}; // end class MultiVectorMutable

// //////////////////////////////////////////////////
// Inlined methods for MultiVector

inline
MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutable::mv_sub_view(
  const index_type& rl, const index_type& ru
  ,const index_type& cl, const index_type& cu
  )
{
  return this->mv_sub_view(Range1D(rl,ru),Range1D(cl,cu));
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MULTI_VECTOR_MUTABLE_H
