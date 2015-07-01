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

#include <vector>

#include "ConstrainedOptPack_initialize_Q_R_Q_X.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"

void ConstrainedOptPack::initialize_Q_R_Q_X(
  size_type            n_R
  ,size_type           n_X
  ,const size_type     i_x_free[]
  ,const size_type     i_x_fixed[]
  ,bool                test_setup
  ,size_type           Q_R_row_i[]
  ,size_type           Q_R_col_j[]
  ,GenPermMatrixSlice  *Q_R
  ,size_type           Q_X_row_i[]
  ,size_type           Q_X_col_j[]
  ,GenPermMatrixSlice  *Q_X
  )
{
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  const size_type
    n = n_R + n_X;
  //
  //                   /  i_x_free[lR-1]    ,if lR = l_x_XR[i] > 0
  // Setup l_x_XR[i] = |  not set yet       ,if l_x_XR[i] == 0
  //                   \  i_x_fixed[lX-1]   ,if lX = l_x_XR[i] < 0
  //
  typedef std::vector<long int> l_x_XR_t;   // ToDo: use temporary workspace
  l_x_XR_t l_x_XR(n);
  std::fill( l_x_XR.begin(), l_x_XR.end(), 0 );
  // Setup the fixed portion of l_x_XR[] first.
  bool i_x_fixed_is_sorted = true;
  if(test_setup) {
    size_type last_set = 0;
    const size_type *i_x_X = i_x_fixed;
    for( long int lX = 1; lX <= n_X; ++lX, ++i_x_X ) {
      if( *i_x_X < 1 || *i_x_X > n ) {
        std::ostringstream omsg;
        omsg
          << "initialize_Q_R_Q_X(...) : Error, "
          << "i_x_fixed[" << (lX-1) << "] = "
          << (*i_x_X) << " is out of bounds";
        throw std::invalid_argument( omsg.str() );
      }
      const long int l = l_x_XR[*i_x_X-1];
      if( l != 0 ) {
        std::ostringstream omsg;
        omsg
          << "initialize_Q_R_Q_X(...) : Error, "
          << "Duplicate entries for i_x_fixed[" << (lX-1) << "] = "
          << (*i_x_X) << " == i_x_fixed[" << (-l-1) << "]";
        throw std::invalid_argument( omsg.str() );
      }
      l_x_XR[*i_x_X-1] = -lX;
      if( *i_x_X < last_set )
        i_x_fixed_is_sorted = false;
      last_set = *i_x_X;
    }
  }
  else {
    size_type last_set = 0;
    const size_type *i_x_X = i_x_fixed;
    for( size_type lX = 1; lX <= n_X; ++lX, ++i_x_X ) {
      l_x_XR[*i_x_X-1] = -lX;
      if( *i_x_X < last_set )
        i_x_fixed_is_sorted = false;
      last_set = *i_x_X;
    }
  }
  // Now setup the free portion of l_x_XR[]
  bool i_x_free_is_sorted = true;
  if( i_x_free == NULL ) {
    long int
      *l_x_XR_itr = &l_x_XR[0];
    size_type lR = 0;
    for( size_type i = 1; i <= n; ++i, ++l_x_XR_itr ) {
      if(*l_x_XR_itr == 0) {
        ++lR;
        *l_x_XR_itr = lR;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPT( !( lR == n_R ) );
  }
  else {
    if(test_setup) {
      size_type last_set = 0;
      const size_type *i_x_R = i_x_free;
      for( size_type lR = 1; lR <= n_R; ++lR, ++i_x_R ) {
        if( *i_x_R < 1 || *i_x_R > n ) {
          std::ostringstream omsg;
          omsg
            << "initialize_Q_R_Q_X(...) : Error, "
            << "i_x_free[" << (lR-1) << "] = "
            << (*i_x_R) << " is out of bounds";
          throw std::invalid_argument( omsg.str() );
        }
        const long int l = l_x_XR[*i_x_R-1];
        if( l != 0 ) {
          std::ostringstream omsg;
          omsg
            << "initialize_Q_R_Q_X(...) : Error, "
            << "Duplicate entries for i_x_free[" << (lR-1) << "] = "
            << (*i_x_R) << " == "
            << (l < 0 ? "i_x_fixed[" : "i_x_free[")
            << (l < 0 ? -l : l) << "]";
          throw std::invalid_argument( omsg.str() );
        }
        l_x_XR[*i_x_R-1] = lR;
        if( *i_x_R < last_set )
          i_x_free_is_sorted = false;
        last_set = *i_x_R;
      }
    }
    else {
      size_type last_set = 0;
      const size_type *i_x_R = i_x_free;
      for( size_type lR = 1; lR <= n_R; ++lR, ++i_x_R ) {
        l_x_XR[*i_x_R-1] = lR;
        if( *i_x_R < last_set )
          i_x_free_is_sorted = false;
        last_set = *i_x_R;
      }
    }
  }
  //
  // Setup Q_R and Q_X
  //
  if( i_x_free_is_sorted && Q_R_row_i == NULL && l_x_XR[n_R-1] == n_R ) {
    //
    // Here Q_R = [ I; 0 ] so lets set it
    //
    Q_R->initialize(n,n_R,GenPermMatrixSlice::IDENTITY_MATRIX);
    // Now we just have to set Q_X which may or may not be sorted
    const long int
      *l_x_X = &l_x_XR[n_R-1] + 1; // Just in case n_X == 0
    size_type
      *row_i_itr = Q_X_row_i,
      *col_j_itr = Q_X_col_j;
    for( size_type iX = n_R+1; iX <= n; ++iX, ++l_x_X, ++row_i_itr, ++col_j_itr ) {
      *row_i_itr = iX;         // This is sorted of course
      *col_j_itr = -(*l_x_X);  // This might be sorted
      TEUCHOS_TEST_FOR_EXCEPT( !(  *l_x_X < 0  ) );
    }					
    Q_X->initialize(
      n,n_X,n_X,0,0,i_x_fixed_is_sorted?GPMSIP::BY_ROW_AND_COL:GPMSIP::BY_ROW
      ,Q_X_row_i,Q_X_col_j,test_setup);
  }
  else {
    //
    // Here i_x_free[] and i_x_fixed[] may be sorted by Q_R != [ I; 0 ]
    //
    if( n_X > 0 && Q_R_row_i == NULL )
      throw std::invalid_argument(
        "initialize_Q_R_Q_X(...) : Error, "
        "Q_R_row_i != NULL since Q_R != [ I; 0 ]" );
    const long int
      *l_x = &l_x_XR[0];
    size_type
      *Q_R_row_i_itr = Q_R_row_i,   // Okay if NULL and n_R == 0
      *Q_R_col_j_itr = Q_R_col_j,
      *Q_X_row_i_itr = Q_X_row_i,   // Okay if NULL and n_X == 0
      *Q_X_col_j_itr = Q_X_col_j;
    for( size_type i = 1; i <= n; ++i, ++l_x ) {
      const long int l = *l_x;
      TEUCHOS_TEST_FOR_EXCEPT( !(  l != 0  ) );
      if( l > 0 ) {
        *Q_R_row_i_itr++ = i;
        *Q_R_col_j_itr++ = l;
      }
      else {
        *Q_X_row_i_itr++ = i;
        *Q_X_col_j_itr++ = -l;
      }
    }					
    TEUCHOS_TEST_FOR_EXCEPT( !(  Q_R_row_i_itr - Q_R_row_i == n_R  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  Q_X_row_i_itr - Q_X_row_i == n_X  ) );
    Q_R->initialize(
      n,n_R,n_R,0,0,i_x_free_is_sorted?GPMSIP::BY_ROW_AND_COL:GPMSIP::BY_ROW
      ,Q_R_row_i,Q_R_col_j,test_setup);
    Q_X->initialize(
      n,n_X,n_X,0,0,i_x_fixed_is_sorted?GPMSIP::BY_ROW_AND_COL:GPMSIP::BY_ROW
      ,Q_X_row_i,Q_X_col_j,test_setup);
  }
}
