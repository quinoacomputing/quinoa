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

#ifndef ALAP_MULTI_VECTOR_H
#define ALAP_MULTI_VECTOR_H

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_RCP.hpp"

namespace AbstractLinAlgPack {

/** \brief . */
enum EApplyBy {
  APPLY_BY_ROW        ///<
  ,APPLY_BY_COL       ///<
};

/** \brief Apply a reduction/transformation operator column by column and
 * return an array of the reduction objects.
 *
 * ToDo: Finish documentation!
 */
void apply_op(
  EApplyBy                        apply_by
  ,const RTOpPack::RTOp           &primary_op
  ,const size_t                   num_multi_vecs
  ,const MultiVector*             multi_vecs[]
  ,const size_t                   num_targ_multi_vecs
  ,MultiVectorMutable*            targ_multi_vecs[]
  ,RTOpPack::ReductTarget*        reduct_objs[]         = NULL
  ,const index_type               primary_first_ele     = 1
  ,const index_type               primary_sub_dim       = 0
  ,const index_type               primary_global_offset = 0
  ,const index_type               secondary_first_ele   = 1
  ,const index_type               secondary_sub_dim     = 0
  );

/** \brief Apply a reduction/transformation operator column by column and reduce the intermediate
 * reduction objects into one reduction object.
 *
 * ToDo: Finish documentation!
 */
void apply_op(
  EApplyBy                        apply_by
  ,const RTOpPack::RTOp           &primary_op
  ,const RTOpPack::RTOp           &secondary_op
  ,const size_t                   num_multi_vecs
  ,const MultiVector*             multi_vecs[]
  ,const size_t                   num_targ_multi_vecs
  ,MultiVectorMutable*            targ_multi_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
  ,const index_type               primary_first_ele     = 1
  ,const index_type               primary_sub_dim       = 0
  ,const index_type               primary_global_offset = 0
  ,const index_type               secondary_first_ele   = 1
  ,const index_type               secondary_sub_dim     = 0
  );

/** \brief Interface for a collection of non-mutable vectors (multi-vector, matrix).
 *
 * This interface is quite restrictive in that it allows a client to access a
 * matrix by accessing rows, columns and/or diagonals.
 * The vector objects returned from the provided access methods
 * \c row(), \c col() and \c diag() are abstract vectors so there is
 * still good implementation flexibility but many <tt>%MatrixOp</tt>
 * implementations will not be able to support this interface.
 *
 * The primary purpose for this interface is to allow for convienent aggregations
 * of column vectors.  Such an orderly arrangement allows for better optimized
 * linear algebra operations such as matrix-matrix multiplication and the solution
 * of linear systems for multiple right hand sides.  Every application area (serial
 * parallel, out-of-core etc.) should be able to define at least one reasonbly
 * efficient implementation of a <tt>%MultiVector</tt> (or a <tt>%MultiVectorMutable</tt>)
 * subclass.
 *
 * The <tt>%MultiVector</tt> interface is derived from the \c MatrixOp 
 * interface and therefore a <tt>%MultiVector</tt> can be considered as a matrix
 * which has some interesting implications.  As an extended matrix interface, this
 * is somewhat of a "last resort" interface that allows many matrix operations to
 * have default implementations based on vector operations. None of the linear
 * algebra methods in <tt>%MatrixOp</tt> or any of the other matrix interfaces
 * have methods that directly accept <tt>%MultiVector</tt> objects.  However,
 * since <tt>%MultiVector</tt> is derived from <tt>%MatrixOp</tt>, a
 * <tt>%MultiVector</tt> object can be used anywere a <tt>%MatrixOp</tt>
 * object is accepted.  In fact, many of the default implementations for the
 * linear algebra methods in <tt>%MatrixOp</tt> test to see if the matrix
 * arguments support <tt>%MultiVector</tt> (or \c MultiVectorMutable</tt>)
 * and will fail if these interfaces are not supported.
 *
 * Note that only certain kinds of access may be preferred and it is allowed
 * for subclasses to return \c NULL vector objects for some types of access.  For
 * example, a matrix may be naturally oriented by column (for the primary role of
 * as a multi-vector) but row or diagonal access may be very inefficient.  For this
 * reason, the client should call the \c access_by() method which returns a bit
 * field that the client can compare to the constants \c ROW_ACCESS, \c COL_ACCESS
 * and \c DIAG_ACCESS.  The method \c access_by() only returns the types of access
 * that are guarrentted to be efficient, but
 * does not necessarily imply that a type of access is not supported.  For example,
 * <tt>(this->access_by() & ROW_ACCESS) == false</tt> this does not mean that 
 * <tt>this->row(i)</tt> will return \c NULL, it only means that row access will
 * be inefficient.  To determine if a certain type of access is even possible, check
 * the return for \c row(), \c col() and/or \c diag().  For example, if <tt>this->rows(1)</tt>
 * returns \c NULL, then this is a flag that row access for every row is not supported.
 * Diagonal access may be a different story.  For some matrix subclasses, only the
 * center diagonal my be easily accessable in which case \c diag(0) may return
 * <tt>!= NULL</tt> but \c diag(k) for all <tt>k != 0</tt> may return \c NULL.
 *
 * Note that since, this interface is derived from \c MatrixOp that it must
 * support the methods \c space_rows() and \c space_cols().  This does not imply
 * however that either of the access methods \c row() or \c col() must return
 * non-<tt>NULL</tt>.
 *
 * Examples of matrix implementations that can support this interface are a dense
 * BLAS compatible matrix (\c ROW_ACCESS, \c COL_ACCESS and \c DIAG_ACCESS), a
 * compressed column sparse matrix (\c COL_ACCESS only), a compressed row sparse matrix
 * (\c ROW_ACCESS only), a unordered sparse matrix with the diagonal explicity stored
 * (<tt>diag(0) != NULL</tt>) etc.
 *
 * Another very powerfull feature of this interface is the ability to
 * apply reduction/transformation operators over a sub-set of rows and
 * columns in a set of multi-vector objects.  The behavior is
 * identical as if the client extracted the rows or columns in a set
 * of multi-vectors and called <tt>Vector::apply_op()</tt> itself.
 * However, the advantage of using the multi-vector methods is that
 * there may be greater opertunity to exploit parallelism.  Also, the
 * intermediate reduction objects over a set of rows or columns can be
 * reduced by a secondary reduction object.
 *
 * ToDo: Finish documentation!
 */
class MultiVector : virtual public MatrixOp {
public:

  /** \brief . */
  using MatrixOp::clone;
  /** \brief . */
  using MatrixOp::Mp_StMtM;

  /** \brief . */
  typedef int  access_by_t;
  /** \brief . */
  enum {
    ROW_ACCESS    = 0x1 ///< 
    ,COL_ACCESS   = 0x2 ///<
    ,DIAG_ACCESS  = 0x4 ///<
  };
  /** \brief . */
  typedef Teuchos::RCP<const Vector>         vec_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const MultiVector>    multi_vec_ptr_t;

  /** @name Friends */
  //@{

  /** \brief . */
  friend void apply_op(
    EApplyBy                        apply_by
    ,const RTOpPack::RTOp           &primary_op
    ,const size_t                   num_multi_vecs
    ,const MultiVector*             multi_vecs[]
    ,const size_t                   num_targ_multi_vecs
    ,MultiVectorMutable*            targ_multi_vecs[]
    ,RTOpPack::ReductTarget*        reduct_objs[]
    ,const index_type               primary_first_ele
    ,const index_type               primary_sub_dim
    ,const index_type               primary_global_offset
    ,const index_type               secondary_first_ele
    ,const index_type               secondary_sub_dim
    );
  /** \brief . */
  friend void apply_op(
    EApplyBy                        apply_by
    ,const RTOpPack::RTOp           &primary_op
    ,const RTOpPack::RTOp           &secondary_op
    ,const size_t                   num_multi_vecs
    ,const MultiVector*             multi_vecs[]
    ,const size_t                   num_targ_multi_vecs
    ,MultiVectorMutable*            targ_multi_vecs[]
    ,RTOpPack::ReductTarget         *reduct_obj
    ,const index_type               primary_first_ele
    ,const index_type               primary_sub_dim
    ,const index_type               primary_global_offset
    ,const index_type               secondary_first_ele
    ,const index_type               secondary_sub_dim
    );

  //@}

  /** @name Clone */
  //@{

  /** \brief Clone the non-const multi-vector object.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>.
   */
  virtual multi_vec_ptr_t mv_clone() const;

  //@}

  /** @name Provide row, column and diagonal access as non-mutable vectors */
  //@{

  /** \brief Return a bit field for the types of access that are the most convenient.
   *
   * Postconditions:<ul>
   * <li> <tt>return & COL_ACCESS || return & ROW_ACCESS || return & DIAG_ACCESS</tt>
   * </ul>
   */
  virtual access_by_t access_by() const = 0;

  /** \brief Get a non-mutable column vector.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->access_by() & COL_ACCESS</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>return.get() != NULL</tt>] <tt>space_cols().is_compatible(return->space()) == true</tt>
   * </ul>
   */
  virtual vec_ptr_t col(index_type j) const = 0;
  /** \brief Get a non-mutable row vector.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->access_by() & ROW_ACCESS</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>return.get() != NULL</tt>] <tt>space_rows().is_compatible(return->space()) == true</tt>
   * </ul>
   */
  virtual vec_ptr_t row(index_type i) const = 0;
  /** \brief Get a non-mutable diagonal vector.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->access_by() & DIAG_ACCESS</tt>] <tt>return.get() != NULL</tt>
   * </ul>
   */
  virtual vec_ptr_t diag(int k) const = 0;

  //@}

  /** @name Sub-view methods */
  //@{

  /** \brief Returns a sub-view of the multi vector.
   *
   * ToDo: Finish documentation!
   *
   * The default implementation returns a \c MultiVectorSubView object for
   * any valid arbitary sub-view.
   */
  virtual multi_vec_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
  
  /** \brief Inlined implementation calls <tt>this->mv_sub_view(Range1D(rl,ru),Range1D(cl,cu))</tt>.
   */
  multi_vec_ptr_t mv_sub_view(
    const index_type& rl, const index_type& ru
    ,const index_type& cl, const index_type& cu
    ) const;

  //@}

protected:

  /** @name Collective apply_op() methods */
  //@{

  /** \brief Apply a reduction/transformation operator row by row, or column by column and return an array
   * of the reduction objects.
   *
   * Preconditions:<ul>
   * <li> [<tt>apply_by == APPLY_BY_COL</tt>] <tt>(this->access_by() & COL_ACCESS) == true)</tt> (throw <tt>???</tt>)
   * <li> [<tt>apply_by == APPLY_BY_ROW</tt>] <tt>(this->access_by() & ROW_ACCESS) == true)</tt> (throw <tt>???</tt>)
   * <li> ToDo: Finish!
   * </ul>
   *
   * The default implementation calls \c this->apply_op().
   *
   * ToDo: Finish documentation!
   */
  virtual void apply_op(
    EApplyBy apply_by, const RTOpPack::RTOp& primary_op
    ,const size_t num_multi_vecs,      const MultiVector*   multi_vecs[]
    ,const size_t num_targ_multi_vecs, MultiVectorMutable*  targ_multi_vecs[]
    ,RTOpPack::ReductTarget* reduct_objs[]
    ,const index_type primary_first_ele,   const index_type primary_sub_dim, const index_type primary_global_offset
    ,const index_type secondary_first_ele, const index_type secondary_sub_dim
    ) const;

  /** \brief Apply a reduction/transformation operator row by row, or column by column and reduce the intermediate
   * reduction objects into one reduction object.
   *
   * Preconditions:<ul>
   * <li> [<tt>apply_by == APPLY_BY_COL</tt>] <tt>(this->access_by() & COL_ACCESS) == true)</tt> (throw <tt>???</tt>)
   * <li> [<tt>apply_by == APPLY_BY_ROW</tt>] <tt>(this->access_by() & ROW_ACCESS) == true)</tt> (throw <tt>???</tt>)
   * <li> ToDo: Finish!
   * </ul>
   *
   * The default implementation calls \c this->apply_op().
   *
   * ToDo: Finish documentation!
   */
  virtual void apply_op(
    EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
    ,const size_t num_multi_vecs,      const MultiVector*   multi_vecs[]
    ,const size_t num_targ_multi_vecs, MultiVectorMutable*  targ_multi_vecs[]
    ,RTOpPack::ReductTarget* reduct_obj
    ,const index_type primary_first_ele, const index_type primary_sub_dim, const index_type primary_global_offset
    ,const index_type secondary_first_ele, const index_type secondary_sub_dim
    ) const;

  //@}

public:

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief Returns <tt>this->mv_clone()<tt>.
   */
  mat_ptr_t clone() const;

  /** \brief Returns <tt>this->mv_sub_view(row_rng,col_rng)</tt> casted to a MatrixOp.
   */
  mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;

  /** \brief Provides a specialized implementation for <tt>mwo_rhs1</tt> of type <tt>MatrixSymDiag</tt>.
   *
   * @return Returns <tt>true</tt> and implements the operation if
   * <tt>dynamic_cast<MatrixSymDiag>(&mwo_rhs1) != NULL
   * && op(*this).access_by() =& MultiVector::COL_ACCESS
   * && (mvm_lhs = dynamic_cast<MultiVectorMutable*>(&mwo_lhs)) != NULL
   * && mvm_lhs->access_by() & MultiVector::COL_ACCESS</tt>.
   * Otherwise, this function returns <tt>false</tt> and does not implement the operation.
   * or <tt>dynamic_cast<const MatrixSymDiag>(&mwo_rhs1) != NULL</tt>.
   *
   * The default implementation relies on column access of <tt>op(*this)</tt>
   * and <tt>mwo_lhs</tt> to implement this method.
   */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2
    ,value_type beta
    ) const;

  /** \brief Provides a specialized implementation for <tt>mwo_rhs2</tt> of type <tt>MatrixSymDiag</tt>.
   *
   * @return Returns <tt>true</tt> and implements the operation if
   * <tt>dynamic_cast<MatrixSymDiag>(&mwo_rhs1) != NULL
   * && op(*this).access_by() =& MultiVector::ROW_ACCESS
   * && (mvm_lhs = dynamic_cast<MultiVectorMutable*>(&mwo_lhs)) != NULL
   * && mvm_lhs->access_by() & MultiVector::ROW_ACCESS</tt>.
   * Otherwise, this function returns <tt>false</tt> and does not implement the operation.
   * or <tt>dynamic_cast<const MatrixSymDiag>(&mwo_rhs1) != NULL</tt>.
   *
   * The default implementation relies on row access of <tt>op(*this)</tt>
   * and <tt>mwo_lhs</tt> to implement this method.
   */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ,value_type beta
    ) const;

  //@}

private:
  
#ifdef DOXYGEN_COMPILE
  Vector         *rows;
  Vector         *columns;
  Vector         *diagonals;
#endif	

}; // end class MultiVector

// //////////////////////////////////////////////////
// Inlined functions

inline
MultiVector::multi_vec_ptr_t
MultiVector::mv_sub_view(
  const index_type& rl, const index_type& ru
  ,const index_type& cl, const index_type& cu
  ) const
{
  return this->mv_sub_view(Range1D(rl,ru),Range1D(cl,cu));
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MULTI_VECTOR_H
