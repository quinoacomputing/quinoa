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

#ifndef MATRIX_IDENT_CONCAT_STD_H
#define MATRIX_IDENT_CONCAT_STD_H

#include "ConstrainedOptPack_MatrixIdentConcat.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Teuchos_RCP.hpp"

namespace ConstrainedOptPack {

/** \brief Concrete implementation class for a matrix vertically concatonated with an identity matrix.
 *
 * Represents an interface for a matrix that represents:
 \verbatim
 
 M = [ alpha*op(D) ]    (TOP)
     [     I       ]

 or

 M = [     I       ]
     [ alpha*op(D) ]   (BOTTOM)
 \endverbatim
 *
 * This subclass allows a client to set the representation matrix \c D.
 */
class MatrixIdentConcatStd : public MatrixIdentConcat {
public:

  /** @name Public types */
  //@{
  /** \brief . */
  enum ETopBottom { TOP, BOTTOM };
  /** \brief . */
  typedef Teuchos::RCP<const MatrixOp> D_ptr_t;
  //@}

  /** @name Constructors/initializers. */
  //@{
  /** \brief Constructs to uninitialized.
   */
  MatrixIdentConcatStd();
  /** \brief Setup with a matrix object.
   *
   * @param  top_or_bottom
   *                 [in] If <tt>TOP</tt> then <tt>M = [ alpha*op(D); I ]</tt> and if <tt>BOTTOM</tt> then
   *                  <tt>M = [ I; alpha*op(D) ]</tt>
   * @param  alpha   [in] Scalar that multiplies \c op(D)
   * @param  D_ptr   [in] Smart pointer to a matrix object that represents \c D.  The matrix object pointed to must
   *                 not be altered until \c this object is no longer in use or <tt>this->set_uninitialized()</tt>
   *                 has been called.
   * @param  D_trans [in] Determines if <tt>op(D) = D</tt> (\c no_trans#) or <tt>op(D) = D'</tt> (\c trans).
   *
   * Preconditions:<ul>
   * <li> <tt>D.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>space_cols->dim() == rows(op(D)) + cols(op(D))</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>space_rows->dim() == cols(op(D))</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>space_cols->sub_space(D_rng)->is_compatible(op(D).space_cols())</tt>
   *      (throw <tt>std::invalid_argument</tt>) See \c D_rng defined below
   * <li> <tt>space_rows->is_compatible(op(D).space_rows())</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->D_ptr().get() == D_ptr.get()</tt>
   * <li> <tt>&this->D()          == this->D_ptr().get()</tt>
   * <li> <tt>this->D_trans()     == D_trans</tt>
   * <li> <tt>this->alpha()       == alpha</tt>
   * <li> <tt>this->rows()        == rows(op(D)) + cols(op(D))</tt>
   * <li> <tt>this->cols()        == cols(op(D))</tt>
   * <li> <tt>&this->space_cols() == space_cols.get()</tt>
   * <li> <tt>&this->space_rows() == space_rows.get()</tt>
   * <li> [<tt>top_or_bottom == TOP</tt>]    <tt>this->D_rng() = [1,rows(op(D))]</tt>
   * <li> [<tt>top_or_bottom == TOP</tt>]    <tt>this->I_rng() = [rows(op(D))+1,rows(op(D))+cols(op(D))]</tt>
   * <li> [<tt>top_or_bottom == BOTTOM</tt>] <tt>this->D_rng() = [cols(op(D))+1,rows(op(D))+cols(op(D))]</tt>
   * <li> [<tt>top_or_bottom == BOTTOM</tt>] <tt>this->I_rng() = [1,cols(op(D))]</tt>
   * </ul>
   */
  virtual void initialize(
    const VectorSpace::space_ptr_t&    space_cols
    ,const VectorSpace::space_ptr_t&   space_rows
    ,ETopBottom                        top_or_bottom
    ,value_type                        alpha
    ,const D_ptr_t                     &D_ptr
    ,BLAS_Cpp::Transp                  D_trans
    );
  /** \brief Set the matrix to uninitialized.
   *
   * Postconditions:<ul>
   * <li> <tt>this->space_cols()</tt> throws an exception
   * <li> <tt>this->space_rows()</tt> throws an exception
   * <li> <tt>this->D_ptr().get() == NULL</tt>
   * <li> <tt>&this->D()</tt> throws an exception
   * <li> <tt>this->D_trans() == no_trans</tt>
   * <li> <tt>this->alpha() == 0.0</tt>
   * <li> <tt>this->rows() == 0</tt>
   * <li> <tt>this->cols() == 0</tt>
   * <li> [<tt>top_or_bottom == TOP</tt>]    <tt>this->D_rng() = [1,rows(op(D))]</tt>
   * <li> [<tt>top_or_bottom == TOP</tt>]    <tt>this->I_rng() = [rows(op(D))+1,rows(op(D))+cols(op(D))]</tt>
   * <li> [<tt>top_or_bottom == BOTTOM</tt>] <tt>this->D_rng() = [cols(op(D))+1,rows(op(D))+cols(op(D))]</tt>
   * <li> [<tt>top_or_bottom == BOTTOM</tt>] <tt>this->I_rng() = [1,cols(op(D))]</tt>
   * </ul>
   */
  virtual void set_uninitialized();
  /** \brief Return the smart reference counted point to the \c D matrix.
   *
   * If the matrix object \c D is owned exclusively by \c this matrix object
   * then <tt>this->D_ptr().count() == 1</tt> on return.
   */
  virtual const D_ptr_t& D_ptr() const;
  //@}

  /** @name Overridden form MatrixIdentConcat */
  //@{
  /** \brief . */
  Range1D D_rng() const;
  /** \brief . */
  Range1D I_rng() const;
  /** \brief . */
  value_type alpha() const;
  /** \brief . */
  const MatrixOp& D() const;
  /** \brief . */
  BLAS_Cpp::Transp D_trans() const;
  //@}

  /** @name Overridden from MatrixOp */
  //@{
  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  const VectorSpace& space_rows() const;
  /** \brief The default just performs a shallow copy and just copies
   * the underlying smart reference counted pointer.  If other
   * behavior is desired then this method must be overridden.
   */
  MatrixOp& operator=(const MatrixOp& m);
  //@}

private:

#ifdef DOXYGEN_COMPILE
  AbstractLinAlgPack::VectorSpace    *space_cols;
  AbstractLinAlgPack::VectorSpace    *space_rows;
  AbstractLinAlgPack::MatrixOp   *D;
  RangePack::Range1D                 D_rng;
  RangePack::Range1D                 I_rng;
#else
  VectorSpace::space_ptr_t  space_cols_;
  VectorSpace::space_ptr_t  space_rows_;
  D_ptr_t           D_ptr_;
  Range1D           D_rng_;
  Range1D           I_rng_;
#endif
  value_type        alpha_;
  BLAS_Cpp::Transp  D_trans_;
 
  //
  void assert_initialized() const;

}; // end class MatrixIndentConcatStd

} // end namespace ConstrainedOptPack

#endif // MATRIX_IDENT_CONCAT_STD_H
