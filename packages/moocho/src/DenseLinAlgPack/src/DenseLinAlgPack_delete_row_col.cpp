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

#include <assert.h>

#include "DenseLinAlgPack_delete_row_col.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"

void DenseLinAlgPack::delete_row_col( size_type kd, DMatrixSliceTriEle* tri_M )
{
  // Validate input
  TEUCHOS_TEST_FOR_EXCEPT( !(  tri_M  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  tri_M->rows()  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= kd && kd <= tri_M->rows()  ) );

  DMatrixSlice   M = tri_M->gms();
  const size_type  n = M.rows();

  if( tri_M->uplo() == BLAS_Cpp::lower ) {
    // Move M31 up one row at a time
    if( 1 < kd && kd < n ) {
      Range1D rng(1,kd-1);
      for( size_type i = kd; i < n; ++i )
        M.row(i)(rng) = M.row(i+1)(rng);
    }
    // Move M33 up and to the left one column at a time
    if( kd < n ) {
      for( size_type i = kd; i < n; ++i )
        M.col(i)(i,n-1) = M.col(i+1)(i+1,n);
    }
  }
  else if(  tri_M->uplo() == BLAS_Cpp::upper ) {
    // Move M13 left one column at a time.
    if( 1 < kd && kd < n ) {
      Range1D rng(1,kd-1);
      for( size_type j = kd; j < n; ++j )
        M.col(j)(rng) = M.col(j+1)(rng);
    }
    // Move the updated U33 up and left one column at a time.
    if( kd < n ) {
      for( size_type j = kd; j < n; ++j )
        M.col(j)(kd,j) = M.col(j+1)(kd+1,j+1);
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true); // Invalid input
  }
}
