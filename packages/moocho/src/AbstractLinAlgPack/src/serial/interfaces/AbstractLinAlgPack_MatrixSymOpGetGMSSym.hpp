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

#ifndef MATRIX_SYM_WITH_OP_GET_GMS_SYM_H
#define MATRIX_SYM_WITH_OP_GET_GMS_SYM_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface that allows the extraction of a const <tt>DenseLinAlgPack::DMatrixSliceSym</tt>
 * view of an abstract matrix.
 *
 * This interface is ment to be used by <tt>MatrixSymOp</tt> objects
 * that store all of their matrix elements in the local address space or can easily
 * access all of the elements from this process.
 *
 * Subclasses that store a BLAS compatible symmetric dense matrix can implement
 * these methods without any dynamic memory allocations.  There is no default
 * implementation for these methods so subclasses that derive from this interface
 * must implement these methods.
 *
 * These methods should never be called directly.  Instead, use the helper
 * class type <tt>MatrixDenseSymEncap</tt>.
 */
class MatrixSymOpGetGMSSym
  : virtual public AbstractLinAlgPack::MatrixSymOp // doxygen needs full name
{
public:

  /** \brief Get a const view of the symmetric abstract matrix in the form <tt>DenseLinAlgPack::DMatrixSliceSym</tt>.
   *
   * @return On ouput, \c return will be initialized to point to storage to the symmetric dense
   *  matrix elements.
   * The output from this function <tt>sym_gms_view = this->get_sym_gms_view()</tt> must be passed to
   * <tt>this->free_sym_gms_view(gms)</tt> to free any memory that may have been allocated.
   * After <tt>this->free_gms_view(gms_view)</tt> is called, \c gms_view must not be used any longer!
   *
   * Postconditions:<ul>
   * <li> <tt>return.rows() == this->rows()</tt>
   * <li> <tt>return.cols() == this->cols()</tt>
   * </ul>
   *
   * Warning!  If a subclass overrides this method, it must also override \c free_sym_gms_view().
   */
  virtual const DenseLinAlgPack::DMatrixSliceSym get_sym_gms_view() const = 0;

  /** \brief Free a view of a symmetric dense matrix initialized from <tt>get_sym_gms_view()>/tt>.
   *
   * @param  sym_gms_view
   *              [in/out] On input, \c sym_gms_view must have been initialized from \c this->get_sym_gms_view().
   *              On output, \c sym_gms_view will become invalid and must not be used.
   *
   * Preconditions:<ul>
   * <li> \c sym_gms_view must have been initialized by \c this->get_sym_gms_view)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> \c sym_gms_view becomes invalid and must not be used any longer!
   * </ul>
   */
  virtual void free_sym_gms_view(const DenseLinAlgPack::DMatrixSliceSym* sym_gms_view) const = 0;

}; // end class MatrixSymOpGetGMSSym

/** \brief Helper class type that simplifies the usage of the <tt>MatrixSymOpGetGMSSym</tt> interface for clients.
 *
 * This takes care of worrying about if the <tt>MatrixSymOpGetGMSSym</tt> interface is supported or not
 * and remembering to free the <tt>DenseLinAlgPack::DMatrixSliceSym</tt> view properly.
 *
 * This class is only to be used on the stack as an automatic variable.  For example, to extract a
 * <tt>DenseLinAlgPack::DMatrixSliceSym</tt> view of an abstract vector and use it to call another function
 * one could write a function like:
 \code
 void call_func(const MatrixSymOpGetGMSSym& mat_in ) {
     func( MatrixDenseSymEncap(mat_in)() );
 }
 \endcode
 * In the above code, if the underlying <tt>MatrixSymOpGetGMSSym</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method
 * <tt>MatrixSymOpGetGMSSym::get_sym_gms_view()</tt> then the above code will only have a constant
 * time overhead.
 */
class MatrixDenseSymEncap {
public:

  /** \brief Construct a <tt>DenseLinAlgPack::DMatrixSliceSym</tt> view from a <tt>MatrixSymOpGetGMSSym</tt> object.
   */
  MatrixDenseSymEncap( const MatrixSymOpGetGMSSym&  mat_get );
  /** \brief Construct a <tt>DenseLinAlgPack::DMatrixSliceSym</tt> view from a <tt>MatrixSymOp</tt> object.
   *
   * If <tt>dynamic_cast<const MatrixSymOpGetGMSSym*>(&mat) == NULL</tt> then a ???
   * exception is thrown.
   */
  MatrixDenseSymEncap( const MatrixSymOp& mat );
  /// Frees the <tt>DenseLinAlgPack::DMatrixSliceSym</tt> view.
  ~MatrixDenseSymEncap();
  /// Returns a constant view of the <tt>DenseLinAlgPack::DMatrixSliceSym</tt> view.
  const DenseLinAlgPack::DMatrixSliceSym operator()() const;

private:

  const MatrixSymOpGetGMSSym     &mat_get_;
  const DenseLinAlgPack::DMatrixSliceSym          sym_gms_view_;
  MatrixDenseSymEncap();                                       // Not defined and not to be called!
  MatrixDenseSymEncap(const MatrixDenseSymEncap&);             // ""
  MatrixDenseSymEncap& operator=(const MatrixDenseSymEncap&);  // ""

}; // end class MatrixDenseSymEncap

// ///////////////////////////////////////////
// Inline members

// MatrixDenseSymEncap

inline
MatrixDenseSymEncap::MatrixDenseSymEncap( const MatrixSymOpGetGMSSym&  mat_get )
  :mat_get_(mat_get)
  ,sym_gms_view_(mat_get_.get_sym_gms_view())
{}

inline
MatrixDenseSymEncap::MatrixDenseSymEncap( const MatrixSymOp& mat )
  :mat_get_(Teuchos::dyn_cast<const MatrixSymOpGetGMSSym>(mat))
  ,sym_gms_view_(mat_get_.get_sym_gms_view())
{}

inline
MatrixDenseSymEncap::~MatrixDenseSymEncap()
{
  mat_get_.free_sym_gms_view(&sym_gms_view_);
}

inline
const DenseLinAlgPack::DMatrixSliceSym MatrixDenseSymEncap::operator()() const
{
  return sym_gms_view_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SYM_WITH_OP_GET_GMS_SYM_H
