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

#include <iomanip>

#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace DenseLinAlgPack {

// ////////////////////////////////////////////////////////////////////////////////
// DMatrixSlice

DVectorSlice DMatrixSlice::p_diag(difference_type k) const {
  if(k > 0) {
    validate_col_subscript(k+1);
    // upper diagonal (k > 0)
    return DVectorSlice( const_cast<value_type*>(col_ptr(1)) + k*max_rows()
      , cols()-k > rows() ? rows() : cols()-k, max_rows()+1 );
  }
  // lower diagonal (k < 0) or center diagonal (k = 0)
  validate_row_subscript(-k+1);
  return DVectorSlice( const_cast<value_type*>(col_ptr(1)) - k
    , rows()+k > cols() ? cols() : rows()+k, max_rows()+1 );
}

EOverLap DMatrixSlice::overlap(const DMatrixSlice& gms) const
{
  typedef DMatrixSlice::size_type size_type;
  
  const DVectorSlice::value_type
    *raw_ptr1 = col_ptr(1),
    *raw_ptr2 = gms.col_ptr(1);

  if( !raw_ptr1 || !raw_ptr2 )
    return NO_OVERLAP;		// If one of the views is unbound there can't be any overlap

  DVectorSlice::size_type
    max_rows1 = max_rows(),
    max_rows2 = gms.max_rows(),
    rows1 = rows(),
    rows2 = gms.rows(),
    cols1 = cols(),
    cols2 = gms.cols();

  // Establish a frame of reference where raw_ptr1 < raw_ptr2
  if(raw_ptr1 > raw_ptr2) {
    std::swap(raw_ptr1,raw_ptr2);
    std::swap(max_rows1,max_rows2);
    std::swap(rows1,rows2);
    std::swap(cols1,cols2);
  }

  if( raw_ptr2 > (raw_ptr1 + (cols1 - 1) * max_rows1 + (rows1 - 1)) ) {
    return NO_OVERLAP; // can't be any overlap
  }

  DVectorSlice::size_type
    start1 = 0,
    start2 = raw_ptr2 - raw_ptr1;

  if(start1 == start2 && max_rows1 == max_rows2 && rows1 == rows2 && cols1 == cols2)
    return SAME_MEM;
  if(start1 + (rows1 - 1) + (cols1 - 1) * max_rows1 < start2)
    return NO_OVERLAP;	// start2 is past the last element in m1 so no overlap.
  // There may be some overlap so determine if start2 lays in the region of m1.
  // Determine the number of elements in v that start2 is past the start of a
  // column of m1.  If start2 was the first element in one of m1's cols
  // row_i would be 1, and if it was just before for first element of the next
  // column of m1 then row_i would equal to max_rows1.
  size_type row_i = (start2 - start1 + 1) % max_rows1;
  if(row_i <= rows1)
    return SOME_OVERLAP; // start2 is in a column of m1 so there is some overlap
  // Determine how many rows in the original matrix are below the last row in m1
  size_type lower_rows = max_rows1 - (start1 % max_rows1 + rows1);
  if(row_i < rows1 + lower_rows)
    return NO_OVERLAP; // m2 lays below m1 in the original matrix
  // If you get here start2 lays above m1 in original matix so if m2 has enough
  // rows then the lower rows of m2 will overlap the upper rows of m1.
  if(row_i + rows2 - 1 <= max_rows1)
    return NO_OVERLAP; // m2 completely lays above m1
  return SOME_OVERLAP; // Some lower rows of m2 extend into m1 
}

#ifdef LINALGPACK_CHECK_RANGE
void DMatrixSlice::validate_row_subscript(size_type i) const
{
  if( i > rows() || !i )
    throw std::out_of_range( "DMatrixSlice::validate_row_subscript(i) :"
                  "row index i is out of bounds"				);
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void DMatrixSlice::validate_col_subscript(size_type j) const
{
  if( j > cols() || !j )
    throw std::out_of_range( "DMatrixSlice::validate_col_subscript(j) :"
                  "column index j is out of bounds"			);
}
#endif

#ifdef LINALGPACK_CHECK_SLICE_SETUP
void DMatrixSlice::validate_setup(size_type size) const
{
  if( !ptr_ && !rows() && !cols() && !max_rows() )
      return; // an unsized matrix slice is ok.
  if( (rows() - 1) + (cols() - 1) * max_rows() + 1 > size )
    throw std::out_of_range( "DMatrixSlice::validate_setup() : "
                  " DMatrixSlice constructed that goes past end of array" );
}
#endif

// /////////////////////////////////////////////////////////////////////////////////
// DMatrix

DVectorSlice DMatrix::p_diag(difference_type k) const {	
  if(k > 0) {
    validate_col_subscript(k+1);
    // upper diagonal (k > 0)
    return DVectorSlice( const_cast<value_type*>(col_ptr(1)) + k * rows()
      , cols()-k > rows() ? rows() : cols()-k, rows()+1 );
  }
  // lower diagonal (k < 0) or center diagonal (k = 0)
  validate_row_subscript(-k+1);
  return DVectorSlice( const_cast<value_type*>(col_ptr(1)) - k
    , rows()+k > cols() ? cols() : rows()+k, rows()+1 );
}

EOverLap DMatrix::overlap(const DMatrixSlice& gms) const {
  return (*this)().overlap(gms);
}

#ifdef LINALGPACK_CHECK_RANGE
void DMatrix::validate_row_subscript(size_type i) const {
  if( i > rows() || !i )
    throw std::out_of_range("DMatrix::validate_row_subscript(i) : row index out of bounds");
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void DMatrix::validate_col_subscript(size_type j) const {
  if( j > cols() || !j )
    throw std::out_of_range("DMatrix::validate_col_subscript(j) : column index out of bounds");
}
#endif

}	// end namespace DenseLinAlgPack

// ///////////////////////////////////////////////////////////////////////////////
// Non-member funcitons

void DenseLinAlgPack::assert_gms_sizes(const DMatrixSlice& gms1, BLAS_Cpp::Transp trans1
  , const DMatrixSlice& gms2, BLAS_Cpp::Transp trans2)
{
  if	(
      (trans1 == trans2) ? 
        gms1.rows() == gms2.rows() && gms1.cols() == gms2.cols() 
        : gms1.rows() == gms2.cols() && gms1.cols() == gms2.rows()
    )
    return; // compatible sizes so exit
  // not compatible sizes.
  throw std::length_error("Matrix sizes are not the compatible");
}
