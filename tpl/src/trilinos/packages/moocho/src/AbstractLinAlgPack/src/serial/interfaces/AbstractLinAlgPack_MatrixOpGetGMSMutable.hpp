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

#ifndef MATRIX_WITH_OP_GET_GMS_MUTABLE_H
#define MATRIX_WITH_OP_GET_GMS_MUTABLE_H

#include "AbstractLinAlgPack_MatrixOpGetGMS.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface that allows the extraction of a non-const <tt>DMatrixSlice</tt>
 * view of an abstract matrix.
 *
 * This interface is ment to be used by <tt>MatrixOp</tt> objects
 * that store all of their matrix elements in the local address space or can easily
 * access all of the elements from this process and can modify the elements in their
 * data structures.
 *
 * Subclasses that store a Fortran compatible dense dense matrix can implement
 * these methods without any dynamic memory allocations.  There is no default
 * implementation for these methods so subclasses that derive from this interface
 * must implement these methods.
 *
 * These methods should never be called directly.  Instead, use the helper
 * class type <tt>MatrixDenseMutableEncap</tt>.
 */
class MatrixOpGetGMSMutable : virtual public MatrixOpGetGMS {
public:

  /** \brief . */
  using MatrixOpGetGMS::get_gms_view;

  /** \brief Get a representation of the abstract matrixr in the form <tt>DenseLinAlgPack::DMatrixSlice</tt>.
   *
   * @return On ouput, \c return will be initialized to point to storage to the dense matrix elements.
   * The output from this function <tt>gms_view = this->get_gms_view()</tt> must be passed to
   * <tt>this->commit_gms_view(gms)</tt> to commit and free any memory that may have been allocated
   * and to ensure the that underlying abstract matrix object has been updated.
   * After <tt>this->commit_gms_view(gms_view)</tt> is called, \c gms_view must not be used any longer!
   *
   * Postconditions:<ul>
   * <li> <tt>return.rows() == this->rows()</tt>
   * <li> <tt>return.cols() == this->cols()</tt>
   * </ul>
   *
   * Warning!  If a subclass overrides this method, it must also override \c commit_gms_view().
   */
  virtual DMatrixSlice get_gms_view() = 0;

  /** \brief Commit changes to a view of a dense matrix initialized from <tt>this->get_gms_view()</tt>.
   *
   * @param  gms_view
   *              [in/out] On input, \c gms_view must have been initialized from \c this->get_gms_view().
   *              On output, \c gms_view will become invalid and must not be used.
   *
   * Preconditions:<ul>
   * <li> \c gms_view must have been initialized by \c this->get_gms_view)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> \c this is guaranteed to be updated.
   * <li> \c gms_view becomes invalid and must not be used any longer!
   * </ul>
   */
  virtual void commit_gms_view(DMatrixSlice* gms_view) = 0;

}; // end class MatrixOpGetGMSMutable

/** \brief Helper class type that simplifies the usage of the <tt>MatrixOpGetGMSMutable</tt> interface for clients.
 *
 * This takes care of worrying about if the <tt>MatrixOpGetGMSMutable</tt> interface is supported or not
 * and remembering to free the <tt>DMatrixSlice</tt> view properly.
 *
 * This class is only to be used on the stack as an automatic variable.  For example, to extract a
 * <tt>DMatrixSlice</tt> view of an abstract vector and use it to set the matrix to a scalar
 * one could write a function like:
 \code
 void assign( const value_type alpha, MatrixOpGetGMSMutable* mat_inout ) {
     MatrixDenseMutableEncap  gms_inout(*mat_inout);
   gms_inout() = alpha;
 }
 \endcode
 * In the above code, if the underlying <tt>MatrixOpGetGMSMutable</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method
 * <tt>MatrixOpGetGMSMutable::get_gms_view()</tt> then the above code will only have a constant
 * time overhead.
 */
class MatrixDenseMutableEncap {
public:

  /** \brief Construct a <tt>DMatrixSlice</tt> view from a <tt>MatrixOpGetGMSMutable</tt> object.
   */
  MatrixDenseMutableEncap( MatrixOpGetGMSMutable*  mat_get );
  /** \brief Construct a <tt>DMatrixSlice</tt> view from a <tt>MatrixOp</tt> object.
   *
   * If <tt>dynamic_cast<MatrixOpGetGMSMutable*>(mat) == NULL</tt> then a ???
   * exception is thrown.
   */
  MatrixDenseMutableEncap( MatrixOp* mat );
  /// Frees the <tt>DMatrixSlice</tt> view and commits the changes.
  ~MatrixDenseMutableEncap();
  /// Returns a non-const view of the <tt>DMatrixSlice</tt> view.
  DMatrixSlice operator()();
  /// Returns a const view of the <tt>DMatrixSlice</tt> view.
  const DMatrixSlice operator()() const;

private:

  MatrixOpGetGMSMutable     *mat_get_;
  DMatrixSlice                gms_view_;
  MatrixDenseMutableEncap();                                          // Not defined and not to be called!
  MatrixDenseMutableEncap(const MatrixDenseMutableEncap&);            // ""
  MatrixDenseMutableEncap& operator=(const MatrixDenseMutableEncap&); // ""

}; // end class MatrixDenseMutableEncap

// ///////////////////////////////////////////
// Inline members

// MatrixDenseMutableEncap

inline
MatrixDenseMutableEncap::MatrixDenseMutableEncap( MatrixOpGetGMSMutable*  mat_get )
  :mat_get_(mat_get)
  ,gms_view_(mat_get_->get_gms_view())
{}

inline
MatrixDenseMutableEncap::MatrixDenseMutableEncap( MatrixOp* mat )
  :mat_get_(&Teuchos::dyn_cast<MatrixOpGetGMSMutable>(*mat))
  ,gms_view_(mat_get_->get_gms_view())
{}

inline
MatrixDenseMutableEncap::~MatrixDenseMutableEncap()
{
  mat_get_->commit_gms_view(&gms_view_);
}

inline
DMatrixSlice MatrixDenseMutableEncap::operator()()
{
  return gms_view_;
}

inline
const DMatrixSlice MatrixDenseMutableEncap::operator()() const
{
  return gms_view_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_WITH_OP_GET_GMS_MUTABLE_H
