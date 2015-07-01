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

#ifndef MATRIX_DECOMP_RANGE_ORTHOG_H
#define MATRIX_DECOMP_RANGE_ORTHOG_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"

namespace ConstrainedOptPack {

/** \brief Matrix subclass for variable reduction orthogonal matrix </tt>R = Gc(:,con_decomp)'*Y</tt>.
 *
 * This matrix class is used to represent the matrix:
 \verbatim

 R = C*(I + D*D')

 inv(R) =  inv(I + D*D') * inv(C)

        = (I - D * inv(I + D'*D) * D') * inv(C)
                       \______/
                          S
 \endverbatim
 * Above, the expresion for <tt>inv(R)</tt> is derived using the Sherman-Morrison-Woodbury formula.
 * The nonsingular matrix <tt>S = I + D'*D</tt> is setup by the client, along with the basis matrix
 * \c C and the direct sensitivity matrix \c D.
 */
class MatrixDecompRangeOrthog
  : public AbstractLinAlgPack::MatrixOpNonsing // full path needed for doxygen
{
public:

  /** @name Public types */
  //@{
#ifndef DOXYGEN_COMPILE
  /** \brief . */
  typedef Teuchos::RCP<const MatrixOpNonsing>     C_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const MatrixOp>                D_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const MatrixSymOpNonsing>  S_ptr_t;
#endif
  //@}

  /** @name Constructors/initializers */
  //@{

  /** \brief Constructs uninitialized.
   *
   * Postconditions:<ul>
   * <li> Same as for <tt>this->set_uninitialized()</tt>.
   * </ul>
   */
  MatrixDecompRangeOrthog();

  /// Calls <tt>this->initialize()</tt>.
  MatrixDecompRangeOrthog(
    const C_ptr_t   &C_ptr
    ,const D_ptr_t  &D_ptr
    ,const S_ptr_t  &S_ptr
    );

  /** \brief Initialize the matrix object.
   *
   * @param  C_ptr  [in]
   * @param  D_ptr  [in]
   * @param  S_ptr  [in]
   *
   * Preconditions:<ul>
   * <li> <tt>C_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>D_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>S_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>C_ptr->space_rows().is_compatible(D_ptr->space_cols()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>S_ptr->space_cols().is_compatible(D_ptr->space_rows()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->C_ptr().get() == C_ptr.get()</tt>
   * <li> <tt>this->D_ptr().get() == D_ptr.get()</tt>
   * <li> <tt>this->S_ptr().get() == S_ptr.get()</tt>
   * <li> <tt>this->space_cols().is_compatible(C_ptr->space_cols())</tt>
   * <li> <tt>this->space_rows().is_compatible(C_ptr->space_rows())</tt>
   * </ul>
   */
  void initialize(
    const C_ptr_t   &C_ptr
    ,const D_ptr_t  &D_ptr
    ,const S_ptr_t  &S_ptr
    );

  /** \brief Make uninitialized.
   *
   * Postconditions:<ul>
   * <li> <tt>this->C_ptr().get() == NULL</tt>
   * <li> <tt>this->D_ptr().get() == NULL</tt>
   * <li> <tt>this->S_ptr().get() == NULL</tt>
   * <li> <tt>this->rows() == 0</tt>
   * <li> <tt>this->cols() == 0</tt>
   * </ul>
   */
  void set_uninitialized();

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  const C_ptr_t& C_ptr() const;
  /** \brief . */
  const D_ptr_t& D_ptr() const;
  /** \brief . */
  const S_ptr_t& S_ptr() const;

  //@}
  
  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  const VectorSpace& space_rows() const;
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const Vector& v_rhs2, value_type beta ) const;

  //@}

  /** @name Overridden from MatrixOpNonsing */
  //@{

  /** \brief . */
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    , const Vector& v_rhs2 ) const;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  AbstractLinAlgPack::MatrixOpNonsing      *C;
  AbstractLinAlgPack::MatrixOp                 *D;
  AbstractLinAlgPack::MatrixSymOpNonsing   *S;
#else
  C_ptr_t    C_ptr_;
  D_ptr_t    D_ptr_;
  S_ptr_t    S_ptr_;
#endif

  //
  void assert_initialized(const char func_name[]) const;

  // not defined and not to be called
  MatrixDecompRangeOrthog(const MatrixDecompRangeOrthog&);
  MatrixDecompRangeOrthog& operator=(const MatrixDecompRangeOrthog&);

}; // end class MatrixDecompRangeOrthog

// /////////////////////////////////
// Inline members

inline
const MatrixDecompRangeOrthog::C_ptr_t&
MatrixDecompRangeOrthog::C_ptr() const
{
  return C_ptr_;
}

inline
const MatrixDecompRangeOrthog::D_ptr_t&
MatrixDecompRangeOrthog::D_ptr() const
{
  return D_ptr_;
}

inline
const MatrixDecompRangeOrthog::S_ptr_t&
MatrixDecompRangeOrthog::S_ptr() const
{
  return S_ptr_;
}

} // end namespace ConstrainedOptPack

#endif // MATRIX_DECOMP_RANGE_ORTHOG_H
