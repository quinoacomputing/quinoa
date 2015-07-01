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

#ifndef MATRIX_SCALING_STRATEGY_H
#define MATRIX_SCALING_STRATEGY_H

namespace AbstractLinAlgPack {

/** \brief Abstract interface for sparse matrix scaling strategies
 *
 * ToDo: Finish documentation!
 */
class MatrixScaling_Strategy {
public:

  /** \brief . */
  virtual ~MatrixScaling_Strategy() {}

  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /// Scale the matrix and save the scalings for later use for rhs and lhs.
  virtual void scale_matrix(
    index_type m, index_type n, index_type nz
    ,const index_type row_i[], const index_type col_j[]
    ,bool new_matrix, value_type A[]
    ) = 0;
  
  /// Scale the rhs vector
  virtual void scale_rhs( BLAS_Cpp::Transp trans, value_type b[] ) const = 0;
  
  /// Scale the lhs vector
  virtual void scale_lhs( BLAS_Cpp::Transp trans, value_type x[] ) const = 0;
  
  //@}

}; // end class MatrixScaling_Strategy

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SCALING_STRATEGY_H
