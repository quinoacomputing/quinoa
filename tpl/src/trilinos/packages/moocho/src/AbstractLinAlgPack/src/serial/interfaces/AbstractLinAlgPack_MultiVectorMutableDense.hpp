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

#ifndef MULTI_VECTOR_MUTABLE_DENSE_H
#define MULTI_VECTOR_MUTABLE_DENSE_H

#include "AbstractLinAlgPack_MatrixOpSerial.hpp"
#include "AbstractLinAlgPack_MatrixOpGetGMSMutable.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "ReleaseResource.hpp"

namespace AbstractLinAlgPack {

/** \brief <tt>MultiVectorMutable</tt> "Adapter" subclass for <tt>DenseLinAlgPack::DMatrixSlice</tt>
 * or <tt>DenseLinAlgPack::DMatrix</tt> object.
 *
 * This class can be used either as a view of a <tt>DenseLinAlgPack::DMatrixSlice</tt> object
 * or as a storage type for a <tt>DenseLinAlgPack::DMatrix</tt> object.
 *
 * To create a storage type with the dimensions of <tt>rows x cols</tt> just call the
 * constructor <tt>MatrixWithOpMutableDense(rows,cols)</tt> or after construction you
 * can call <tt>this->initialize(rows,cols)</tt>.
 *
 * To simply create a view of a matrix or its transpose, say \c M, without ownership
 * just call <tt>MatrixWithOpMutableDense(M(),M_trans,NULL)</tt> or after construction call
 * <tt>this->initialize(M(),M_trans,NULL)</tt>.
 *
 * Alternately, \c this can be given a matrix with the responsibility to delete any associated
 * memory by calling <tt>this->initialize()</tt> with a <tt>ReleaseResource</tt> object to
 * perform the deallocation.
 *
 * If \c this has been initialized by <tt>this->initialize(rows,cols)</tt> and
 * if the client really needs to get at the <tt>DenseLinAlgPack::DMatrix</tt> object
 * itself, then it can be obtained as:
 \code
 DMatrix* f( MultiVectorMutableDense* M )
 {
     return dynamic_cast<MemMngPack::ReleaseResource_ref_count_ptr<DMatrix&> >(*M.gms_release()).ptr.get();
 }
 \endcode
 * This is not pretty but it is not supposed to be.  Of course the above function will return NULL
 * the <tt>dynamic_cast<></tt> fails or if there is no allocated <tt>DMatrix</tt> object to get.
 */
class MultiVectorMutableDense
  : public AbstractLinAlgPack::MultiVectorMutable   // doxygen needs the full path
  , public MatrixOpSerial
  , public MatrixOpGetGMS
{
public:

  /** \brief . */
  typedef Teuchos::RCP<
    MemMngPack::ReleaseResource>  release_resource_ptr_t;

  /** @name Constructors / initializers */
  //@{

  /** \brief Calls <tt>this->initialize(rows,cols)</tt>.
   */
  MultiVectorMutableDense(
    const size_type                    rows = 0
    ,const size_type                   cols = 0
    );
  /** \brief Calls <tt>this->initialize(gms,gms_trans,gms_release)</tt>.
   */
  MultiVectorMutableDense(
    DMatrixSlice                       gms
    ,BLAS_Cpp::Transp                  gms_trans
    ,const release_resource_ptr_t&     gms_release
    );
  /** \brief Call <tt>this->initialize(v,v_release)</tt> with an allocated <tt>DenseLinAlgPack::DVector</tt>
   * object.
   */
  void initialize(
    const size_type                    rows
    ,const size_type                   cols
    );
  /** \brief Initialize with a dense matrix slice.
   *
   * Note that solve of the method overrides have to allocate a temporary memory
   * if <tt>gms_trans != no_trans</tt> (see \c get_gms_view() and \c output()).
   *
   * Postconditions:<ul>
   * <li> <tt>this->rows() = ( gms_trans = no_trans ? gms.rows() : gms.cols())</tt>
   * <li> <tt>this->cols() = ( gms_trans = no_trans ? gms.cols() : gms.rows())</tt>
   * </ul>
   */
  void initialize(
    DMatrixSlice                       gms
    ,BLAS_Cpp::Transp                  gms_trans
    ,const release_resource_ptr_t&     gms_release
    );

  //@}

  /** @name Access */
  //@{
  
  /** \brief Return the non-const dense matrix.
   *
   * Note that calling this method may result in the matrix implementation
   * being modified.  Therefore, no other methods on \c this object should be
   * called until the <tt>DMatrixSlice</tt> returned from this method is
   * discarded.
   */
  DMatrixSlice set_gms();
  /** \brief Return a const dense matrix.
   */
  const DMatrixSlice get_gms() const;
  /** \brief Return if underlying matrix is being viewed as the transpose or non-transposed.
   */
  BLAS_Cpp::Transp gms_trans() const;
  /** \brief Return a <tt>RCP<></tt> pointer to the object that will
   * release the associated resource.
   */
  const release_resource_ptr_t& gms_release() const;

  //@}

  /** @name Overridden from MatrixOpGetGMS */
  //@{

  /** \brief . */
  const DMatrixSlice get_gms_view() const;
  /** \brief . */
  void free_gms_view(const DMatrixSlice* gms_view) const;

  //@}

  /** @name Overridden from MatrixOpGetGMSMutable */
  //@{

  /** \brief . */
  DMatrixSlice get_gms_view();
  /** \brief . */
  void commit_gms_view(DMatrixSlice* gms_view);

  //@}
  
  /** @name Overridden from MatrixBase */
  //@{

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  void zero_out();
  /** \brief . */
  void Mt_S( value_type alpha );
  /** \brief . */
  MatrixOp& operator=(const MatrixOp& mwo_rhs);
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  bool Mp_StM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs
    ) const;
  /** \brief . */
  bool Mp_StM(
    value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs
    );
  /** \brief . */
  bool syrk(
     BLAS_Cpp::Transp M_trans, value_type alpha
    ,value_type beta, MatrixSymOp* sym_lhs
    ) const;
  /** \brief . */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2
    ,value_type beta ) const;
  /** \brief . */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ,value_type beta ) const;

  //@}

  /** @name Overridden from MultiVector */
  //@{

  /** \brief . */
  access_by_t access_by() const;

  //@}

  /** @name Overridden from MultiVectorMutable */
  //@{

  /** \brief . */
  vec_mut_ptr_t col(index_type j);
  /** \brief . */
  vec_mut_ptr_t row(index_type i);
  /** \brief . */
  vec_mut_ptr_t diag(int k);
  /** \brief . */
  multi_vec_mut_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng);

  //@}

  /** @name Overridden from MatrixOpSerial */
  //@{

  /** \brief . */
  void Vp_StMtV(
    DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(
    DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;

  // ToDo: Add more overrides as they are needed!

  //@}

protected:

  /** @name Overridden from MultiVector */
  //@{

  /** \brief . */
  void apply_op(
    EApplyBy apply_by, const RTOpPack::RTOp& primary_op
    ,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
    ,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
    ,RTOpPack::ReductTarget* reduct_objs[]
    ,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
    ,const index_type secondary_first_ele , const index_type secondary_sub_dim 
    ) const;
  /** \brief . */
  void apply_op(
    EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
    ,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
    ,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
    ,RTOpPack::ReductTarget *reduct_obj
    ,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
    ,const index_type secondary_first_ele , const index_type secondary_sub_dim 
    ) const;

  //@}

private:

  // ///////////////////////////////////////
  // Private data members
  
  DMatrixSlice              gms_;
  BLAS_Cpp::Transp          gms_trans_;
  release_resource_ptr_t    gms_release_;

}; // end class MultiVectorMutableDense			

// //////////////////////////////////////
// Inline members

inline
DMatrixSlice
MultiVectorMutableDense::set_gms()
{
  return gms_;
}

inline
const DMatrixSlice
MultiVectorMutableDense::get_gms() const
{
  return gms_;
}

inline
BLAS_Cpp::Transp
MultiVectorMutableDense::gms_trans() const
{
  return gms_trans_;
}

inline
const MultiVectorMutableDense::release_resource_ptr_t&
MultiVectorMutableDense::gms_release() const
{
  return gms_release_;
}

} // end namespace AbstractLinAlgPack

#endif // MULTI_VECTOR_MUTABLE_DENSE_H
