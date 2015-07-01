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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_BASE_H
#define ABSTRACT_LINALG_PACK_MATRIX_BASE_H

#include <stdexcept>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Base class for all polymorphic matrices.
  */
class MatrixBase {
public:

  /// Thrown if matrices are incompatible
  class IncompatibleMatrices : public std::logic_error
  {public: IncompatibleMatrices(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Virtual destructor
  virtual ~MatrixBase() {}

  /** @name Vector spaces for the columns and rows of the matrix */
  //@{

  /// Vector space for vectors that are compatible with the columns of the matrix.
  virtual const VectorSpace& space_cols() const = 0;

  /// Vector space for vectors that are compatible with the rows of the matrix.
  virtual const VectorSpace& space_rows() const = 0;

  //@}

  /** @name Dimensionality */
  //@{

  /** \brief Return the number of rows in the matrix.
   *
   * The default implementation returns <tt>space_cols().dim()</tt>.
   */
  virtual size_type rows() const;

  /** \brief Return the number of columns in the matrix.
   *
   * The default implementation returns <tt>space_rows().dim()</tt>.
   */
  virtual size_type cols() const;

  /** \brief Return the number of nonzero elements in the matrix.
   *
   * The default is to just assume it is dense and to return
   * <tt>rows() * cols()</tt>.
   */
  virtual size_type nz() const;

  //@}

};	// end class MatrixBase

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_BASE_H
