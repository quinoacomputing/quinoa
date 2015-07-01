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

#ifndef COP_MATRIX_SYM_IDENTITY_SERIAL_H
#define COP_MATRIX_SYM_IDENTITY_SERIAL_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsingSerial.hpp"

namespace ConstrainedOptPack {

/** \brief Matrix class for a serial scaled identity matrix.
 *
 * More operations will be overridden as they are needed by various applications.
 */
class MatrixSymIdentitySerial
  : virtual public AbstractLinAlgPack::MatrixSymOpNonsingSerial // doxygen needs full name
  , virtual public AbstractLinAlgPack::MatrixExtractInvCholFactor
{
public:
  
    /** @name Constructors/initalizers */
  //@{

  /// Calls <tt>this->initalize()</tt>
  MatrixSymIdentitySerial( size_type size = 1, value_type scale = 1.0 );

  /** \brief . */
  void initialize( size_type size, value_type scale );

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  value_type scale() const;

  //@}

  /** Overridden from MatrixBase */
  //@{

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type nz() const;

  //@}

  /** Overridden from MatrixOp */
  //@{

  /** \brief . */
  std::ostream& output(std::ostream& out) const;

  //@}

  /** Overridden from MatrixOpSerial */
  //@{

  /** \brief . */
  void Vp_StMtV(
    DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const DVectorSlice& vs_rhs2, value_type beta) const;

  //@}

  /** @name Overridden from MatrixNonsingSerial */
  //@{

  /** \brief . */
  void V_InvMtV(
    DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1,const DVectorSlice& vs_rhs2 ) const;

  //@}

  /** @name Overridden from MatrixSymNonsingSerial */
  //@{

  /** \brief . */
  void M_StMtInvMtM(
    DMatrixSliceSym* sym_gms_lhs, value_type alpha
    ,const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
    ,EMatrixDummyArg
    ) const;

  //@}

  /** @name Overridden from MatrixExtractInvCholFactor */
  //@{

  /** \brief . */
  void extract_inv_chol( DMatrixSliceTriEle* InvChol ) const;

  //@}

private:

  size_type   size_;
  value_type  scale_;

  
}; // end class MatrixSymIdentitySerial

// //////////////////////////////////////
// Inline members

inline
value_type MatrixSymIdentitySerial::scale() const
{
  return scale_;
}

} // end namespace ConstrainedOptPack

#endif // COP_MATRIX_SYM_IDENTITY_SERIAL_H
