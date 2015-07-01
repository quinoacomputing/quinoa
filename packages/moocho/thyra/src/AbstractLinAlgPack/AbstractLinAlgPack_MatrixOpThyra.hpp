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

#ifndef ALAP_MATRIX_OP_Thyra_HPP
#define ALAP_MATRIX_OP_Thyra_HPP

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace AbstractLinAlgPack {

/** \brief <tt>MatrixOp</tt> adapter subclass for <tt>Thyra::LinearOpBase</tt>.
 */
class MatrixOpThyra : virtual public MatrixOp {
public:

  /** @name Constructors / Initializers */
  //@{

  /** \brief Construct to uninitialized.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_linear_op().get() == NULL</tt>
   * <li><tt>this->rows() == 0</tt>
   * <li><tt>this->cols() == 0</tt>
   * </ul>
   */
  MatrixOpThyra();
  /** \brief Calls <tt>this->initialize()</tt>.
   */
  MatrixOpThyra(
    const Teuchos::RCP<const Thyra::LinearOpBase<value_type> >   &thyra_linear_op
    ,BLAS_Cpp::Transp                                                    thyra_linear_op_trans = BLAS_Cpp::no_trans
    );
  /** \brief Initalize given a smart pointer to a <tt>Thyra::LinearOpBase</tt> object.
   *
   * @param  thyra_linear_op  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
   *
   * Preconditioins:<ul>
   * <li><tt>thyra_linear_op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li><tt>thyra_linear_op->opSupported(Thyra::NOTRANS) && thyra_linear_op->opSupported(Thyra::TRANS)
   *     (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_linear_op().get() == thyra_linear_op.get()</tt>
   * <li><tt>this->thyra_linear_op_trans() == thyra_linear_op_trans</tt>
   * <li><tt>this->rows() == thyra_linear_op->range()->dim()</tt>
   * <li><tt>this->cols() == thyra_linear_op->domain()->dim()</tt>
   * </ul>
   */
  virtual void initialize(
    const Teuchos::RCP<const Thyra::LinearOpBase<value_type> >   &thyra_linear_op
    ,BLAS_Cpp::Transp                                                    thyra_linear_op_trans = BLAS_Cpp::no_trans
    );
  /** \brief Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorBase</tt> object.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_linear_op().get() == NULL</tt>
   * </ul>
   */
  Teuchos::RCP<const Thyra::LinearOpBase<value_type> > set_uninitialized();
  /** \brief Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorBase</tt> object.
   *
   * If <tt>this->thyra_linear_op().count() == 1</tt>, then <tt>this</tt>
   * has sole ownership of the <tt>*this->thyra_linear_op()</tt> object.
   */
  const Teuchos::RCP<const Thyra::LinearOpBase<value_type> >& thyra_linear_op() const;
  /** \brief . */
  BLAS_Cpp::Transp thyra_linear_op_trans() const;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  const VectorSpace& space_rows() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  mat_mut_ptr_t clone();
  /** \brief . */
  MatrixOp& operator=(const MatrixOp& mwo_rhs);
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2, value_type beta
    ) const;
  /// Works for <tt>MultiVectorMutableThyra</tt> arguments
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ,value_type beta
    ) const;

  //@}

private:
  
#ifdef DOXYGEN_COMPILE
  Thyra::LinearOpBase<value_type>                                *thyra_linear_op;
#else
  Teuchos::RCP<const Thyra::LinearOpBase<value_type> >   thyra_linear_op_;
#endif
  BLAS_Cpp::Transp                                               thyra_linear_op_trans_;
  VectorSpaceThyra                                               space_cols_;
  VectorSpaceThyra                                               space_rows_;

}; // end class MatrixOpThyra

// //////////////////////////////////////////////////
// Inlined functions

inline
const Teuchos::RCP<const Thyra::LinearOpBase<value_type> >&
MatrixOpThyra::thyra_linear_op() const
{
  return thyra_linear_op_;
}

inline
BLAS_Cpp::Transp MatrixOpThyra::thyra_linear_op_trans() const
{
  return thyra_linear_op_trans_;
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MATRIX_OP_Thyra_HPP
