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

#ifndef MATRIX_IDENT_CONCAT_H
#define MATRIX_IDENT_CONCAT_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"

namespace ConstrainedOptPack {

/** \brief Matrix class for a matrix vertically concatonated with an identity matrix {abstract}.
 *
 * Represents an interface for a matrix that represents:
 \verbatim

 M = [ alpha*op(D) ]
     [      I      ]
 where:
     D_rng = [1,rows(op(D))]
     I_rng = [rows(op(D))+1,rows(op(D))+cols(op(D))]

 or

 M = [      I      ]
     [ alpha*op(D) ]
 where:
     D_rng = [cols(op(D))+1,rows(op(D))+cols(op(D))]
     I_rng = [1,cols(op(D))]
 \endverbatim
 * and \c I is a <tt>op(D).cols() x op(D).cols()</tt> indentity matrix and
 * the full matrix \c M is of order <tt>(op(D).rows() + op(D).cols()) x op(D).cols()</tt>.
 */
class MatrixIdentConcat
  : virtual public AbstractLinAlgPack::MatrixOp
{
public:

  /** @name Access to representation.
   */
  //@{
  /** \brief . */
  virtual Range1D D_rng() const = 0;
  /** \brief . */
  virtual Range1D I_rng() const = 0;
  /** \brief . */
  virtual value_type alpha() const = 0;
  /** \brief . */
  virtual const MatrixOp& D() const = 0;
  /** \brief . */
  virtual BLAS_Cpp::Transp D_trans() const = 0;
  //@}

  /** @name Overridden from MatrixBase */
  //@{
  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;
  //@}

  /** @name Overridden from MatrixOp */
  //@{
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& vs_rhs2, value_type beta
    ) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2, value_type beta
    ) const;
  //@}

}; // end class MatrixIdentConcat

} // end namespace ConstrainedOptPack

#endif // MATRIX_IDENT_CONCAT_H
