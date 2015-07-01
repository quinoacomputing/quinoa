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

#include "AbstractLinAlgPack_MatrixExtractSparseElements.hpp"

namespace AbstractLinAlgPack {

// Overridden from MatrixConvertToSparseFortranCompatible

index_type
MatrixExtractSparseElements::num_nonzeros(
  EExtractRegion        extract_region
  ,EElementUniqueness   element_uniqueness
  ) const
{
  index_type dl,du;
  get_dl_du(extract_region,&dl,&du);
  return count_nonzeros(
    element_uniqueness,NULL,NULL,Range1D(1,rows()),Range1D(1,cols()),dl,du);
}

void MatrixExtractSparseElements::coor_extract_nonzeros(
  EExtractRegion                extract_region
  ,EElementUniqueness           element_uniqueness
  ,const index_type             len_Aval
  ,value_type                   Aval[]
  ,const index_type             len_Aij
  ,index_type                   Arow[]
  ,index_type                   Acol[]
  ,const index_type             row_offset
  ,const index_type             col_offset
   ) const
{
  index_type dl,du;
  get_dl_du(extract_region, &dl, &du);
  coor_extract_nonzeros(
    element_uniqueness
    ,NULL,NULL,Range1D(1,rows()),Range1D(1,cols()),dl,du
    ,1.0
    ,len_Aval,Aval,len_Aij,Arow,Acol,row_offset,col_offset);
}

// private

void MatrixExtractSparseElements::get_dl_du(
  EExtractRegion extract_region, index_type* dl, index_type* du
  ) const
{
  const size_type
    rows = this->rows(),
    cols = this->cols();
  switch(extract_region) {
    case EXTRACT_FULL_MATRIX:
      *dl = -(rows-1);
      *du = +(cols-1);
      break;
    case EXTRACT_UPPER_TRIANGULAR:
      *dl = 0;
      *du = +(cols-1);
      break;
    case EXTRACT_LOWER_TRIANGULAR:
      *dl = -(rows-1);
      *du = 0;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
      break;
  }
}


}	// end namespace AbstractLinAlgPack 
