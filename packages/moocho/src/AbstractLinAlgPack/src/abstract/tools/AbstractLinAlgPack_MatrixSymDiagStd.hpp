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

#ifndef MATRIX_SYM_DIAGONAL_STD_H
#define MATRIX_SYM_DIAGONAL_STD_H

#include "AbstractLinAlgPack_MatrixSymInitDiag.hpp"
#include "AbstractLinAlgPack_MatrixSymDiag.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Simple diagonal matrix class.
 *
 * ToDo: Implement clone_mswons() and deal with this->diag_ptr().count() > 1
 * by cloning vector if told to.  This allows lazy evaluation of the clone_mswons()
 * method.
 */
class MatrixSymDiagStd
  : public virtual MatrixSymInitDiag
  , public virtual MatrixSymDiag
{
public:

  /** \brief PostMod class to use with <tt>MemMngPack::AbstractFactorStd</tt>.
   */
  class PostMod {
  public:
    PostMod(VectorSpace::space_ptr_t vectorSpace)
      : vectorSpace_(vectorSpace) {}

    void initialize(MatrixSymDiagStd* matrix) const
        { matrix->initialize(vectorSpace_->create_member()); }
         
  private:
    VectorSpace::space_ptr_t vectorSpace_;

  }; // end PostMod


  /** @name Constructors/initalizers */
  //@{

  /// Calls <tt>this->initialize()</tt>.
  MatrixSymDiagStd(
    const VectorSpace::vec_mut_ptr_t& diag   = Teuchos::null
    ,bool                             unique = true
    );

  /** \brief Initialize given the diagonal vector (or no vector at all).
   *
   * @param  diag   [in] Vector to be used for the diagonal.  If <tt>diag.get() == NULL</tt>
   *                then \c this will be uninitialized.
   * @param  unique [in] Determines if the underlying \c diag vector is guaranteed to be
   *                unique and not shared.
   */
  void initialize(
    const VectorSpace::vec_mut_ptr_t& diag
    ,bool                             unique = true
    );

  //@}

  /** @name Access */
  //@{

  /** \brief Give non-const access to the diagonal vector.
   *
   * Preconditions:<ul>
   * <li> <tt>this->diag_ptr().get() != NULL</tt> (throw <tt>???</tt>)
   * </ul>
   *
   * ToDo: Finish documentation!
   */
  VectorMutable& diag();
  /** \brief . */
  const VectorSpace::vec_mut_ptr_t& diag_ptr() const;
  /** \brief . */
  bool unique() const;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /// Returns 0 if not initalized (this->diag() == NULL).
  size_type rows() const;
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
  MatrixOp& operator=(const MatrixOp& mwo_rhs);
  /** \brief Add to a mutable matrix lhs.
   *
   * Preconditions:<ul>
   * <li> #dynamic_cast<MultiVectorMutable*>(m_lhs) != NULL#.
   * </ul>
   */
  bool Mp_StM(MatrixOp* g_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const;
  /** \brief . */
  void Vp_StMtV(VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const Vector& v_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief Implements the symmetric rank-k update for all diagonal matrix lhs
   *
   * @return Returns <tt>true</tt> if <tt>dynamic_cast<MatrixSymDiagStd>(sym_lhs) != NULL</tt>.
   * Otherwise, returns false.
   */
  bool syrk(
    BLAS_Cpp::Transp   M_trans
    ,value_type        alpha
    ,value_type        beta
    ,MatrixSymOp   *sym_lhs
    ) const;

  //@}

  /** Overridden from MatrixOpNonsing */
  //@{

  /** \brief . */
  void V_InvMtV(VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    , const Vector& v_rhs2) const;
  /** \brief . */
  void V_InvMtV(VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2) const;

  //@}

  /** @name Overridden from MatrixSymInitDiag */
  //@{

  /** \brief . */
  void init_identity( const VectorSpace& space_diag, value_type alpha );
  /** \brief . */
  void init_diagonal( const Vector& diag );

  //@}

  /** @name Overridden from MatrixSymDiag */
  //@{

  /** \brief . */
  const Vector& diag() const;

  //@}

private:

  VectorSpace::vec_mut_ptr_t     diag_;
  bool                           unique_;

  void copy_unique();

}; // end class MatrixSymDiagStd

// ////////////////////////////////////////
// Inline members

inline
bool MatrixSymDiagStd::unique() const
{
  return unique_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SYM_DIAGONAL_STD_H
