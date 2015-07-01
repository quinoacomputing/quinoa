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

#ifndef MATRIX_COMPOSITE_STD_H
#define MATRIX_COMPOSITE_STD_H

#include <deque>

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Teuchos_RCP.hpp"
#include "ReleaseResource.hpp"

namespace AbstractLinAlgPack {

/** \brief Matrix class for matrices composed out of a set of other matrices and vectors.
 *
 * This matrix object represents:
 \verbatim

 M =   sum( a(j) * E(rM(j)) * op(P(j))*op(A(j))*op(Q(j))*F(cM(j)), j=1..num_mats )
     + sum( b(k) * E(rV(k)) * op(op(G(k))*v(k)) * F(cV(k)), k=1..num_vecs )
 \endverbatim
 * where: <br>
 * <tt>a(j)</tt> : Scalar for sub-matrix <tt>A(j)</tt><br>
 * <tt>P(j), Q(j)</tt> : GenPermMatrixSlice objects for sub-matrix <tt>A(j)</tt><br>
 * <tt>A(j)</tt> : <tt>MatrixOp</tt> object<br>
 * <tt>b(k)</tt> : Scalar for the sub-vector <tt>v(k)</tt> <br>
 * <tt>G(k)</tt> : GenPermMatrixSlice object for sub-vector <tt>v(k)</tt><br>
 * <tt>v(k)</tt> : <tt>Vector</tt> object <br>
 * <tt>E(r)</tt> : An appropriatly sized matrix object that moves the
 *     first row of aggregate sub-matrix or sub-vector to row \c r in \c M. <br>
 * <tt>F(c)</tt> : An appropriatly sized matrix object that moves the
 *     first column of aggregate sub-matrix or sub-vector to column \c c in \c M. <br>
 *
 * This formulation allows sub-matrices and sub-vectors to overlap.  However, bacause of the need
 * for the compatibility of the vector sub-spaces, the overlapped sub-matrices and sub-vectors
 * may not be compatible and therefore the matrix will be invalid and not work.  This can be
 * checked for by this matrix subclass.  As long as the sub-matrics and sub-vectors do not
 * overlap, the resulting matrix object will be perfectly compatible.
 *
 * In the above formulation is it possible for parts to be missing or replaced by simpler
 * means.  For example, a sub-matrix may forgo <tt>P(j)</tt> and <tt>Q(j)</tt> in which case
 * they are replaced (implicitly) by the identity matrix.  Or they can be replaced by <tt>Range1D</tt>
 * objects so that <tt>op(P(j))*op(A(j))*op(Q(j))*F(cM(j))</tt> becomes <tt>op(A(j))(rng_P(j),rng_Q(j))</tt>.
 * Also, a sub-matrix may only be composed of only a general permulation matrix <tt>P(j)</tt> and forgo
 * <tt>Q(j)</tt>.  Likewise a sub-vector entry can replace <tt>op(G(k)*v(k))</tt> using a range
 * as <tt>op(v(k)(rng_G(k)))</tt> or be a naked sub-vector <tt>op(v)</tt>.
 *
 * Before adding any sub-vectors or sub-matrices the method \c reinitialize() must be called
 * to set the shape of the matrix.  This will have the effect of wiping the matrix clean.
 * Then the methods \c add_vector() and \c add_matrix() can be called to add the constituent
 * sub-vectors and sub-matrices.  There are different forms of these methods for the special
 * cases described above.
 *
 * After all of the sub-vectors and sub-matrices have been added using the methods \c add_vector
 * and \c add_matrix() , the method \c finish_construction() must be called before \c this matrix
 * object will become fully initialized and ready for use.
 *
 * The cleanup for all of the matrix and vector objects is given to
 * <tt>\ref MemMngPack::ReleaseResource "ReleaseResouce"</tt> objects
 * so this is a very flexible aggregate matrix class from many perspactives.
 *
 * Access to the constituent sub-vectors and sub-matrices is given through STL iterators
 * with the element types <tt>SubVectorEntry</tt> and <tt>SubMatrixEntry</tt>.  The methods
 * \c vectors_begin() and \c vectors_end() return the iterators for the list of sub-vectors
 * while \c matrices_begin() and \c matrices_end() return the list of sub-matrices.  Note that
 * non-const iterators are returned (from non-const functions) as well as const iterators
 * (from const functions).  The client should not attempt to modify the elememnts returned
 * from these non-const iterators.  The reason that non-const iterators are returned is so
 * that specific sub-vectors and sub-matrices can be removed by calling \c remove_vector()
 * and \c remove_matrix().  If the client does modify any of the elements in any way or
 * calls \c add_vector(), \c add_matrix(), \c remove_vector() or \c remove_matrix() at
 * anytime, the current matrix is invalidated and \c finish_construction() must be called
 * to update the composite matrix object.
 *
 * Note that implementing the <tt>VectorSpace</tt> objects that are returned by
 * \c space_rows() and \c space_cols() may be non-trivial.  These methods return
 * <tt>VectorSpaceBlocked</tt> objects.  Therefore, in order for clients to
 * use vectors with this matrix object, it must create vectors compatible with these
 * vector space objects.  This can be a little tricky but should be doable.
 *
 * The method <tt>::sub_view(row_rng,col_rng)</tt> will return the simplest object it can given
 * the input \c row_rng, \c col_rng.  For example if the input ranges correspond to whole, unique
 * matrix, then that matrix will be returned and not some encapulation of that matrix.
 */
class MatrixComposite : public MatrixOp {
public:

  // ///////////////////////////////////
  // Public types

  /** \brief . */
  typedef Teuchos::RCP<
    MemMngPack::ReleaseResource>  release_resource_ptr_t;

  /** \brief Vector list entry for a sub-vector.
   *
   * ToDo: Finish Documentation!
   */
  struct SubVectorEntry {
    /** \brief . */
    typedef Teuchos::RCP<const GenPermMatrixSlice> GPMS_ptr_t;
    /** \brief . */
    SubVectorEntry(
      size_type r_l, size_type c_l, value_type beta
      ,const Range1D& rng_G
      ,const GPMS_ptr_t& G, const release_resource_ptr_t& G_release, BLAS_Cpp::Transp G_trans
      ,const Vector* v
      ,const release_resource_ptr_t& v_release, BLAS_Cpp::Transp v_trans
      )
      :r_l_(r_l),c_l_(c_l),beta_(beta),rng_G_(rng_G),G_(G),G_release_(G_release),G_trans_(G_trans)
      ,v_(v),v_release_(v_release),v_trans_(v_trans)
      {}
    /** \brief . */
    bool operator==(const SubVectorEntry v)
      {
        return 
          r_l_==v.r_l_ && c_l_==v.c_l_ && beta_==v.beta_
          && rng_G_==v.rng_G_ && G_.get()==v.G_.get() && G_release_.get()==v.G_release_.get() && G_trans_==v.G_trans_
          && v_==v.v_ && v_release_.get()==v.v_release_.get() && v_trans_==v.v_trans_;
      }
    /** \brief . */
    size_type                  r_l_,   ///< row of first element of vector in composite matrix.
                               c_l_;   ///< column of first element of vector in composite matrix.
    /** \brief . */
    value_type                 beta_;  ///< Scaling vector for vector elements
    /** \brief . */
    Range1D                    rng_G_; ///< rng_G_.size() > 0 => G_ is ignored, rng_G_.full_range() whole v!
    /** \brief . */
    GPMS_ptr_t                 G_;     ///< Will be non-identity if rng_G_.size() == 0.
    /** \brief . */
    release_resource_ptr_t     G_release_;
    /** \brief . */
    BLAS_Cpp::Transp           G_trans_;///< Determines op(G) == G (no_trans) or op(G) == G' (trans)
    /** \brief . */
    const Vector         *v_;     ///< Pointer to the vector (non-NULL)
    /** \brief . */
    release_resource_ptr_t     v_release_;
    /** \brief . */
    BLAS_Cpp::Transp           v_trans_;///< Determines op(v) = v (no_trans) or op(v) == v' (trans)
  }; // end struct SubVectorEntry

  /// Warning!  This could be changed to some other STL container!
  typedef std::deque<SubVectorEntry> vector_list_t;

  /** \brief Matrix list entry for a sub-matrix.
   *
   * ToDo: Finish Documentation!
   */
  struct SubMatrixEntry {
    /** \brief . */
    typedef Teuchos::RCP<const GenPermMatrixSlice> GPMS_ptr_t;
    /** \brief . */
    SubMatrixEntry(
      size_type r_l, size_type r_u, size_type c_l, size_type c_u, value_type alpha
      ,const Range1D& rng_P
      ,const GPMS_ptr_t& P, const release_resource_ptr_t& P_release, BLAS_Cpp::Transp P_trans
      ,const MatrixOp* A, const release_resource_ptr_t& A_release, BLAS_Cpp::Transp A_trans
      ,const Range1D& rng_Q
      ,const GPMS_ptr_t& Q, const release_resource_ptr_t& Q_release, BLAS_Cpp::Transp Q_trans
      )
      :r_l_(r_l),r_u_(r_u),c_l_(c_l),c_u_(c_u),alpha_(alpha),rng_P_(rng_P),P_(P),P_release_(P_release),P_trans_(P_trans)
      ,A_(A),A_release_(A_release),A_trans_(A_trans),rng_Q_(rng_Q),Q_(Q),Q_release_(Q_release),Q_trans_(Q_trans)
      {}
    /** \brief . */
    bool operator==(const SubMatrixEntry m)
      {
        return
          r_l_==m.r_l_ && r_u_==m.r_u_ && c_l_==m.c_l_ && c_u_==m.c_u_ && alpha_==m.alpha_
          && rng_P_==m.rng_P_ && P_.get()==m.P_.get() && P_release_.get()==m.P_release_.get() && P_trans_==m.P_trans_
          && A_==m.A_ && A_release_.get()==m.A_release_.get() && A_trans_==m.A_trans_
          && rng_Q_==m.rng_Q_ && Q_.get()==m.Q_.get() && Q_release_.get()==m.Q_release_.get() && Q_trans_==m.Q_trans_;
      }
    /** \brief . */
    size_type                  r_l_, r_u_, c_l_, c_u_;
    /** \brief . */
    value_type                 alpha_;
    /** \brief . */
    Range1D                    rng_P_;  // rng_P_.size() > 0 => P_ is ignored, rng_P_.full_range() => all rows op(A)
    /** \brief . */
    GPMS_ptr_t                 P_;
    /** \brief . */
    release_resource_ptr_t     P_release_;
    /** \brief . */
    BLAS_Cpp::Transp           P_trans_;
    /** \brief . */
    const MatrixOp         *A_;
    /** \brief . */
    release_resource_ptr_t     A_release_;
    /** \brief . */
    BLAS_Cpp::Transp           A_trans_;
    /** \brief . */
    Range1D                    rng_Q_; // rng_Q_.size() > 0 => Q_ is ignored, rng_Q_.full_range() => all columns op(A)
    /** \brief . */
    GPMS_ptr_t                 Q_;
    /** \brief . */
    release_resource_ptr_t     Q_release_;
    /** \brief . */
    BLAS_Cpp::Transp           Q_trans_;
  }; // end struct SubMatrixEntry

  /// Warning!  This could be changed to some other STL container!
  typedef std::deque<SubMatrixEntry> matrix_list_t;

  /** @name Constructors, initializers */
  //@{

  /** \brief Construct.
   *
   * Calls <tt>this->reinitalize()</tt>.
   */
  MatrixComposite( size_type rows = 0, size_type cols = 0 );

  /** \brief Initialize a sized (on unsized) zero matrix to start with.
   *
   * After calling this function the user can add the constituent matrices and
   * vectors using the \c add_matrix() and \c add_vector() methods.
   *
   * Postconditions:<ul>
   * <li> <tt>this->vectors_begin()  == this->vectors_end()</tt>
   * <li> <tt>this->matrices_begin() == this->matrices_end()</tt>
   * </ul>
   *
   */
  void reinitialize( size_type rows = 0, size_type cols = 0 );

  /** \brief Add a sub-vector beta*op(op(G)*v).
   *
   * ToDo : Finish Documentation!
   */
  void add_vector(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    beta
    ,const GenPermMatrixSlice      *G
    ,const release_resource_ptr_t  &G_release
    ,BLAS_Cpp::Transp              G_trans
    ,const Vector            *v
    ,const release_resource_ptr_t  &v_release
    ,BLAS_Cpp::Transp              v_trans
    );

  /** \brief Add a sub-vector beta*op(v(rng_G)).
   *
   * ToDo : Finish Documentation!
   */
  void add_vector(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    beta
    ,const Range1D                 &rng_G
    ,const Vector            *v
    ,const release_resource_ptr_t  &v_release
    ,BLAS_Cpp::Transp              v_trans
    );

  /** \brief Add a sub-vector beta*op(v)
   *
   * ToDo : Finish Documentation!
   */
  void add_vector(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    beta
    ,const Vector            *v
    ,const release_resource_ptr_t  &v_release
    ,BLAS_Cpp::Transp              v_trans
    );

  /** \brief Remove a sub-vector.
   *
   * Preconditions:<ul>
   * <li> <tt>this->vectors_begin() != this->vectors_end()</tt>
   * <li> <tt>this->vectors_begin() <= itr && itr < this->vectors_end()</tt>
   * </ul>
   */
  void remove_vector( vector_list_t::iterator itr );

  /** \brief Add a sub-matrix alpha*op(P)*op(A)*op(Q).
   *
   * ToDo : Finish Documentation!
   */
  void add_matrix(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    alpha
    ,const GenPermMatrixSlice      *P
    ,const release_resource_ptr_t  &P_release
    ,BLAS_Cpp::Transp              P_trans
    ,const MatrixOp            *A
    ,const release_resource_ptr_t  &A_release
    ,BLAS_Cpp::Transp              A_trans
    ,const GenPermMatrixSlice      *Q
    ,const release_resource_ptr_t  &Q_release
    ,BLAS_Cpp::Transp              Q_trans
    );

  /** \brief Add a sub-matrix alpha*op(A)(rng_P,rng_Q).
   *
   * ToDo : Finish Documentation!
   */
  void add_matrix(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    alpha
    ,const Range1D                 &rng_P
    ,const MatrixOp            *A
    ,const release_resource_ptr_t  &A_release
    ,BLAS_Cpp::Transp              A_trans
    ,const Range1D                 &rng_Q
    );

  /** \brief Add a sub-matrix alpha*op(A)(rng_P,:)*op(Q).
   *
   * ToDo : Finish Documentation!
   */
  void add_matrix(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    alpha
    ,const Range1D                 &rng_P
    ,const MatrixOp            *A
    ,const release_resource_ptr_t  &A_release
    ,BLAS_Cpp::Transp              A_trans
    ,const GenPermMatrixSlice      *Q
    ,const release_resource_ptr_t  &Q_release
    ,BLAS_Cpp::Transp              Q_trans
    );

  /** \brief Add a sub-matrix alpha*op(P)*op(A)(:,rng_Q)
   *
   * ToDo : Finish Documentation!
   */
  void add_matrix(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    alpha
    ,const GenPermMatrixSlice      *P
    ,const release_resource_ptr_t  &P_release
    ,BLAS_Cpp::Transp              P_trans
    ,const MatrixOp            *A
    ,const release_resource_ptr_t  &A_release
    ,BLAS_Cpp::Transp              A_trans
    ,const Range1D                 &rng_Q
    );

  /** \brief Add a sub-matrix alpha*op(A).
   *
   * ToDo : Finish Documentation!
   */
  void add_matrix(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    alpha
    ,const MatrixOp            *A
    ,const release_resource_ptr_t  &A_release
    ,BLAS_Cpp::Transp              A_trans
    );

  /** \brief Add a general permutation sub-matrix alpha*op(P).
   *
   * ToDo : Finish Documentation!
   */
  void add_matrix(
    size_type                      row_offset
    ,size_type                     col_offset
    ,value_type                    alpha
    ,const GenPermMatrixSlice      *P
    ,const release_resource_ptr_t  &P_release
    ,BLAS_Cpp::Transp              P_trans
    );

  /** \brief Remove a sub-matrix.
   *
   * Preconditions:<ul>
   * <li> <tt>this->matrices_begin() != this->matrices_end()</tt>
   * <li> <tt>this->matrices_begin() <= itr && itr < this->matrices_end()</tt>
   * </ul>
   */
  void remove_matrix( matrix_list_t::iterator itr );

  /** \brief Call to finish the construction process.
   *
   * This method must be called after all of the sub-vectors and sub-matrices have
   * been added and before <tt>this</tt> matrix object can be used.  This method will
   * validate that the constructed matrix is valid.  It is up to the client to pass
   * in vector space objects (presumably of type \c VectorSpaceBlocked) that
   * represent the rows and columns of this matrix.  It would be very complicated
   * for this matrix class to figure out how to construct these vector space
   * objects in general.
   *
   * Preconditions:<ul>
   * <li> <tt>space_rows.get() != NULL</tt> (throw \c std::invalid_argument)
   * <li> <tt>space_rows->dim() == cols</tt> where \c cols was passed to \c reinitialize()
   *      (throw \c std::invalid_argument).
   * <li> <tt>space_cols.get() != NULL</tt> (throw \c std::invalid_argument)
   * <li> <tt>space_cols->dim() == rows</tt> where \c rows was passed to \c reinitialize()
   *      (throw \c std::invalid_argument).
   * </ul>
   *
   * @param  space_cols  [in] Will represent the vector space returned by <tt>this->space_cols()</tt>.
   *                     The vector space object will be owned by \c this and must be not be modified
   *                     while \c this is in use.  Of course, the sub-spaces from this vector space
   *                     object must be compatible with the vector spaces for the constitient matrices
   *                     and vectors that were added.
   * @param  space_rows  [in] Will represent the vector space returned by <tt>this->space_rows()</tt>.
   *                     The vector space object will be owned by \c this and must be not be modified
   *                     while \c this is in use.  Of course, the sub-spaces from this vector space
   *                     object must be compatible with the vector spaces for the constitient matrices
   *                     and vectors that were added.
   */
  virtual void finish_construction(
    const VectorSpace::space_ptr_t&  space_cols
    ,const VectorSpace::space_ptr_t& space_rows
    );

  //@}

  /** @name Sub-vector, sub-matrix access (using iterators) */
  //@{

  /** \brief . */
  int                             num_vectors() const;
  /** \brief . */
  vector_list_t::iterator         vectors_begin();
  /** \brief . */
  vector_list_t::iterator         vectors_end();
  /** \brief . */
  vector_list_t::const_iterator   vectors_begin() const;
  /** \brief . */
  vector_list_t::const_iterator   vectors_end() const;
  /** \brief . */
  int                             num_matrices() const;
  /** \brief . */
  matrix_list_t::iterator         matrices_begin();
  /** \brief . */
  matrix_list_t::iterator         matrices_end();
  /** \brief . */
  matrix_list_t::const_iterator   matrices_begin() const;
  /** \brief . */
  matrix_list_t::const_iterator   matrices_end() const;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  const VectorSpace& space_rows() const;
  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
  /** \brief . */
  void Vp_StMtV(VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const Vector& v_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(VectorMutable* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const Vector& v_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(VectorMutable* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const SpVectorSlice& sv_rhs3, value_type beta) const;

  //@}

private:

  // ///////////////////////////////
  // private data members

  bool                       fully_constructed_;
  size_type                  rows_, cols_;
#ifdef DOXYGEN_COMPILE
  MatrixOp               *matrices;
    Vector               *vectors;
    VectorSpace                *space_cols;
    VectorSpace                *space_rows;
#else
  matrix_list_t              matrix_list_;
  vector_list_t              vector_list_;
  VectorSpace::space_ptr_t   space_cols_;
  VectorSpace::space_ptr_t   space_rows_;
#endif

  // ///////////////////////////////
  // private member functions

  void assert_fully_constructed() const;

}; // end class MatrixComposite

} // end namespace AbstractLinAlgPack

#endif // MATRIX_COMPOSITE_STD_H
