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

#ifndef MATRIX_HESSIAN_SUPER_BASIC_INIT_DIAGONAL_H
#define MATRIX_HESSIAN_SUPER_BASIC_INIT_DIAGONAL_H

#include <vector>

#include "ConstrainedOptPack_MatrixHessianSuperBasic.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymInitDiag.hpp"

namespace ConstrainedOptPack {

/** \brief Matrix class that adds the ability to initialize to a diagonal
 * to a MatrixHessainSuperBasic object.
 *
 * Essentially, the matrix #B_RR# must support the #MatrixSymInitDiag#
 * interface.
 */

class MatrixHessianSuperBasicInitDiagonal
  : public virtual MatrixHessianSuperBasic
  , public virtual MatrixSymInitDiag
  {
public:

  /** \brief Constructs to uninitialized.
   */
  MatrixHessianSuperBasicInitDiagonal();

  /** @name Overridden from MatrixHessianSuperBasic.
   **/
  //@{

  /** \brief Initialize the matrix and require B_RR to support #MatrixSymInitDiag#.
   *
   * Preconditions:\begin{itemize}
   * \item #dynamic_cast<MatrixSymInitDiag*>(const_cast<MatrixSymWithOpFactorized*>(B_RR_ptr.get())) != NULL#
   * \end{itemize}
   *
   * This overridden function is a little bit of a hack but it is better than the
   * alternatives that I could think of.  What is important is that the #MatrixSymInitDiag#
   * interface for this class will be supported as long as this function executes
   * succesfully on initialization.  This is far better than waiting until later when
   * dynamic casting would be done and failed.  We need to catch any mistakes right away.
   *
   * In this manner we will have given up some compile-time checking but the run-time check will
   * be very good.
   */
  void initialize(
    size_type            n
    ,size_type           n_R
    ,const size_type     i_x_free[]
    ,const size_type     i_x_fixed[]
    ,const EBounds       bnd_fixed[]
    ,const B_RR_ptr_t&   B_RR_ptr
    ,const B_RX_ptr_t&   B_RX_ptr
    ,BLAS_Cpp::Transp    B_RX_trans
    ,const B_XX_ptr_t&   B_XX_ptr
    );
  
  //@}

  /** @name Overridden from MatrixSymInitDiag.
   *
   * These function call the corresponding functions on
   * #B_RR_ptr()# and then call the equivalent of:
   \begin{verbatim}
   MatrixHessianSuperBasic::initialize(
       this->B_RR_ptr()->rows(),this->B_RR_ptr()->cols()
     ,NULL,NULL,NULL
     ,this->B_RR_ptr()
     ,this->B_RX_ptr(),this->B_RX_trans()
     ,this->B_XX_ptr()
    );
  \end{verbatim}
   **/
  //@{

  /** \brief . */
  void init_identity( size_type n, value_type alpha );
  /** \brief . */
  void init_diagonal( const DVectorSlice& diag );

  //@}

private:

  // ///////////////////////////////////
  // Private data members

  MatrixSymInitDiag   *B_RR_init_;  // will be non-null for any valid initalization

  // //////////////////////////
  // Private member functions
  
}; // end class MatrixHessianSuperBasicInitDiagonal

} // end namespace ConstrainedOptPack

#endif // MATRIX_HESSIAN_SUPER_BASIC_INIT_DIAGONAL_H
