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

#ifndef ALAP_MATRIX_SYM_IDENTITY_H
#define ALAP_MATRIX_SYM_IDENTITY_H

#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Matrix subclass for a scaled identity matrix.
 *
 * More operations will be overridden as they are needed by various applications.
 */
class MatrixSymIdent : virtual public MatrixSymOpNonsing {
public:
  
  /** @name Constructors/initializers */
  //@{

  /// Calls <tt>this->initialize()</tt>.
  MatrixSymIdent(
    const VectorSpace::space_ptr_t&          vec_space = Teuchos::null
    ,const value_type                        scale     = 1.0
    );

  /** \brief . */
  void initialize(
    const VectorSpace::space_ptr_t&          vec_space
    ,const value_type                        scale       = 1.0
  );

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  value_type scale() const;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /// Returns 0 if not initalized.
  size_type rows() const;
  /// Returns <tt>this->rows()</tt>
  size_type nz() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2, value_type beta ) const;

  //@}

  /** @name Overridden from MatrixNonsing */
  //@{

  /** \brief . */
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2 ) const;

  //@}

private:

  VectorSpace::space_ptr_t  vec_space_;
  value_type                scale_;
  
}; // end class MatrixSymIdent

// ///////////////////////////////////////////
// Inline members

inline
value_type MatrixSymIdent::scale() const
{
  return scale_;
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MATRIX_SYM_IDENTITY_H
