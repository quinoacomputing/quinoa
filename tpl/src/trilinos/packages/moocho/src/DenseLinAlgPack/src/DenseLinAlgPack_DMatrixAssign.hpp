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
//
// These are assignment functions that used to be in GenMatrixOp but because of
// some circular dependency problems their declarations where moved into a seprate
// file.

#ifndef GEN_MATRIX_ASSIGN_H
#define GEN_MATRIX_ASSIGN_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

// ///////////////////////////////////////////////////////////////////////////////////
/* * @name {\bf DMatrix/DMatrixSlice Assignment Fucntions}.
  *
  * These are functions for assigning one matrix to another.  The rhs matrices
  * can be transposed (op(rhs) = rhs) or non-transposed (op(rhs) = trans(rhs).
  * The assignment operators for
  * DMatrix and DMatrixSlice use these functions to implement their
  * assignment operator functions.  The assignmet functions for triangular
  * matrices are ment to be called from the assignment operator functions
  * for wrapper classes for triangular (upper and lower) and sysmetric 
  * (upper and lower storage).
  *
  * These assignment functions check for aliasing and assignment to self
  * and create any needed temporaries if needed.
  */
// @{

/// gm_lhs = alpha (elementwise)
void assign(DMatrix* gm_lhs, value_type alpha);

/// gm_lhs = op(gms_rhs)
void assign(DMatrix* gm_lhs, const DMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs);

/// gms_lhs = alpha (elementwise)
void assign(DMatrixSlice* gms_lhs, value_type alpha);

/// gms_lhs = op(gms_rhs)
void assign(DMatrixSlice* gms_lhs, const DMatrixSlice& gms_rhs, BLAS_Cpp::Transp trans_rhs);

/// tri_ele_gms_lhs = alpha (elementwise)
void assign(DMatrixSliceTriEle* tri_gms_lhs, value_type alpha);

/// tri_ele_gms_lhs = tri_ele_gms_rhs
void assign(DMatrixSliceTriEle* tri_ele_gms_lhs, const DMatrixSliceTriEle& tri_ele_gms_rhs);

//		end Assignment Fucntions
// @}

}	// end namespace DenseLinAlgPack

#endif	// GEN_MATRIX_ASSIGN_H
