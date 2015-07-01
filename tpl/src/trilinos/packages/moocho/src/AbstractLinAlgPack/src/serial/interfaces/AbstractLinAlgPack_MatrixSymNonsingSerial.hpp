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

#ifndef SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H
#define SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H

#include "AbstractLinAlgPack_MatrixNonsingSerial.hpp"
#include "AbstractLinAlgPack_MatrixSymNonsing.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all serial polymorphic symmetrix nonsingular matrices that
 * can be used to solve for linear systems relatively efficiently.
 *
 * The methods of this interface should not be called directly but instead through
 * the \ref MatrixSymNonsingularSerial_funcs "provided nonmember functions".
 */
class MatrixSymNonsingSerial
  : virtual public MatrixNonsingSerial
  , virtual public AbstractLinAlgPack::MatrixSymNonsing // doxygen needs full name
{
public:

  /** \brief . */
  using MatrixSymNonsing::M_StMtInvMtM;

  /** @name Level-3 */
  //@{

  /** \brief sym_gms_lhs = alpha * op(mwo) * inv(M) * op(mwo)'.
    *
    * The default implementation is based on the operation M_StInvMtM(...)
    * assuming that this \c M is a symmetric matrix.  For an efficient implementation
    * (for this = L*L' for instance) the subclass may want to override this function.
    */
  virtual void M_StMtInvMtM(
    DMatrixSliceSym* sym_gms_lhs, value_type alpha
    ,const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
    ,EMatrixDummyArg
    ) const;

  //@}

  /** @name Overridden from MatrixSymNonsing */
  //@{

  void M_StMtInvMtM(
    MatrixSymOp* sym_lhs, value_type alpha
    ,const MatrixOp& mwo, BLAS_Cpp::Transp mwo_trans
    ,EMatrixDummyArg
    ) const;

  //@}

};	// end class MatrixSymNonsingSerial

/** \defgroup MatrixSymNonsingularSerial_funcs MatrixSymNonsingSerial nonmember inline functions.
 *
 * These nonmember functions allow operations to be called on \c MatrixSymNonsingSerial objects
 * in similar manner to those in \c DenseLinAlgPack.
 */
//@{

inline
/// sym_gms_lhs = alpha * op(mwo) * inv(mswof) * op(mwo)'
void M_StMtInvMtM(
  DMatrixSliceSym* sym_gms_lhs, value_type alpha
  ,const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
  ,const MatrixSymNonsingSerial& mswons
  ,MatrixSymNonsingSerial::EMatrixDummyArg mwo_rhs
   )
{
  mswons.M_StMtInvMtM(sym_gms_lhs,alpha,mwo,mwo_trans,mwo_rhs);
}

//@}

} // end namespace AbstractLinAlgPack

#endif	// SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H
