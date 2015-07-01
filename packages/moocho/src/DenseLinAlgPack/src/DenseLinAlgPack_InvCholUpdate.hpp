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

#ifndef INV_CHOL_UPDATE_H
#define INV_CHOL_UPDATE_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

/** \brief . */
/* * Perform an update of the Cholesky factor of a symetric matrix (SM = L * L^T)
  * and while maintaining the triangular structure.
  *
  * This function permforms an update of L^T or L^-1 depending on whether
  * we are updating the cholesky factor or its inverse.  This function
  * implements algorithm A3.4.1 in Dennis and Schabel.  The update is:
  * (J_new^T = L_old^T + u * v^T) where J_new is rotated back to triangular form.
  *
  * Preconditions: <ul>
  * <li> UpTriM.rows() == UpTriM.cols() == u.size() == v.size() (throw std::length_error)
  * </ul>
  */
void update_chol_factor(DMatrixSlice* UpTriM, DVectorSlice* u
  , const DVectorSlice& v);

/** \brief . */
/* * Perform a jacobi rotation or a matrix about row i.
  *
  * Preconditions: <ul>
  * <li> UpTriM.rows() == UpTriM.cols() (throw std::length_error)
  * </ul>
  */
void jacobi_rotate(DMatrixSlice* UpTriM, size_type row_i, value_type alpha
  , value_type beta); 

}	// end namespace DenseLinAlgPack

#endif	// INV_CHOL_UPDATE_H
