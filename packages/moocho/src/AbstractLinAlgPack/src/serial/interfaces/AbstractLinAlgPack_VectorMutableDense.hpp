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

#ifndef VECTOR_WITH_OP_MUTABLE_DENSE_H
#define VECTOR_WITH_OP_MUTABLE_DENSE_H

#include "AbstractLinAlgPack_VectorSpaceSerial.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorApplyOpSerialBase.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "Teuchos_RCP.hpp"
#include "ReleaseResource.hpp"

namespace AbstractLinAlgPack {

/** \brief DVector "Adaptor" subclass for <tt>DenseLinAlgPack::DVectorSlice</tt>
 * or <tt>DenseLinAlgPack::DVector</tt> objects.
 *
 * This class can be used either as a view of a <tt>DenseLinAlgPack::DVectorSlice</tt> object
 * or as a storage type for a <tt>DenseLinAlgPack::DVector</tt> object.
 *
 * To create a storage type with the dimension of \c dim just call the constructor
 * <tt>VectorMutableDense(dim)</tt> or after construction you can call
 * <tt>this->initialize(dim)</tt>.
 *
 * To simply create a view of a vector, say \c v, without ownership just call
 * <tt>VectorMutableDense(v(),NULL)</tt> or after construction call
 * <tt>this->initialize(v(),NULL)</tt>.
 *
 * Alternately, \c this can be given a vector with the responsibility to
 * delete any associated memory by calling <tt>this->initialize()</tt>
 * with a <tt>ReleaseResource</tt> object to perform the deallocation.
 *
 * If \c this has been initialized by <tt>this->initialize(dim)</tt> and if
 * the client really needs to get at the <tt>DenseLinAlgPack::DVector</tt> object
 * itself, then it can be obtained as:
 \code
 void f( VectorMutableDense* v )
     namespace rmp = MemMngPack;
     DVector &_v = *dynamic_cast<rmp::ReleaseResource_ref_count_ptr<DVector>&>(*v.vec_release()).ptr;

 \endcode
 * This is not pretty but it is not supposed to be.  Of course the above function will throw
 * an exception if the <tt>dynamic_cast<></tt> fails.
 */
class VectorMutableDense
  : virtual public VectorMutable
  , virtual private VectorApplyOpSerialBase
{
public:

  /** \brief . */
  typedef Teuchos::RCP<
    MemMngPack::ReleaseResource>  release_resource_ptr_t;

  /** @name Constructors/initializers */
  //@{

  /** \brief Calls <tt>this->initialize(dim)</tt>.
   */
  VectorMutableDense(
    const size_type                    dim = 0
    );
  /** \brief Calls <tt>this->initialize(v,v_release)</tt>.
   */
  VectorMutableDense(
    DVectorSlice                        v
    ,const release_resource_ptr_t&     v_release
    );
  /** \brief Call <tt>this->initialize(v,v_release)</tt> with an allocated <tt>DenseLinAlgPack::DVector</tt>
   * object.
   */
  void initialize(
    const size_type                    dim
    );
  /** \brief Initialize with a dense vector slice.
   */
  void initialize(
    DVectorSlice                        v
    ,const release_resource_ptr_t&     v_release
    );

  //@}

  /** @name Access */
  //@{
  
  /** \brief Return the non-const dense vector.
   *
   * Note that calling this method will result in the vector implementation
   * being modified.  Therefore, no other methods on \c this object should be
   * called until the <tt>DVectorSlice</tt> returned from this method is
   * discarded.
   *
   * Note that the underlying implementation calls <tt>this->has_changed()</tt>
   * before this method returns.
   */
  DVectorSlice set_vec();
  /** \brief Return a const dense vector.
   */
  const DVectorSlice get_vec() const;
  /** \brief Return a <tt>RCP<></tt> pointer to the object that will
   * release the associated resource.
   */
  const release_resource_ptr_t& vec_release() const;

  //@}

  /** @name Overriddenn from Vector */
  //@{

  /** \brief . */
  const VectorSpace& space() const;
  /** \brief . */
  void apply_op(
    const RTOpPack::RTOp& op
    ,const size_t num_vecs, const Vector* vecs[]
    ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
    ,RTOpPack::ReductTarget *reduct_obj
    ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
    ) const;
  /** \brief . */
  index_type dim() const;
  /** \brief . */
  value_type get_ele(index_type i) const;
  /** \brief . */
  void get_sub_vector( const Range1D& rng, RTOpPack::SubVector* sub_vec ) const;
  /** \brief . */
  void free_sub_vector( RTOpPack::SubVector* sub_vec ) const;

  //@}

  /** @name Overriddenn from VectorMutable */
  //@{

  /** \brief . */
  VectorMutable& operator=(value_type alpha);
  /** \brief . */
  VectorMutable& operator=(const Vector& v);
  /** \brief . */
  VectorMutable& operator=(const VectorMutable& v);
  /** \brief . */
  void set_ele( index_type i, value_type val );
  /** \brief . */
  vec_mut_ptr_t sub_view( const Range1D& rng );
  /** \brief . */
  void get_sub_vector( const Range1D& rng, RTOpPack::MutableSubVector* sub_vec );
  /** \brief . */
  void commit_sub_vector( RTOpPack::MutableSubVector* sub_vec );
  /** \brief . */
  void set_sub_vector( const RTOpPack::SparseSubVector& sub_vec );
  /** \brief . */
  void Vp_StMtV(
    value_type                       alpha
    ,const GenPermMatrixSlice        &P
    ,BLAS_Cpp::Transp                P_trans
    ,const Vector                    &x
    ,value_type                      beta
    );

  //@}

  /// Hack
  VectorMutableDense* operator&()
  {
    return this;
  }
private:

  // ///////////////////////////////////////
  // Private data members
  
  DVectorSlice              v_;
  release_resource_ptr_t    v_release_;
  VectorSpaceSerial         space_;

  // Not defined and not to be called
  //VectorMutableDense(const VectorMutableDense&);
  VectorMutableDense& operator=(const VectorMutableDense&);

}; // end class VectorMutableDense

// //////////////////////////////////////
// Inline members

inline
DVectorSlice
VectorMutableDense::set_vec()
{
  this->has_changed();
  return v_;
}

inline
const DVectorSlice
VectorMutableDense::get_vec() const
{
  return v_;
}

inline
const VectorMutableDense::release_resource_ptr_t&
VectorMutableDense::vec_release() const
{
  return v_release_;
}

} // end namespace AbstractLinAlgPack

#endif // VECTOR_WITH_OP_MUTABLE_DENSE_H
