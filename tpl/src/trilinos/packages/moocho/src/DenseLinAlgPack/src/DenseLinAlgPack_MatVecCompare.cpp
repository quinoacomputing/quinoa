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

#include <limits>
#include <algorithm>

#include "DenseLinAlgPack_MatVecCompare.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"

namespace {
//
using DenseLinAlgPack::value_type;
//
using DenseLinAlgPack::sqrt_eps;
//
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
//
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }
//
inline
bool _comp(value_type val1, value_type val2)
{
  const value_type denom = my_max( ::fabs(val1), 1.0 ); // compare relative errors.
  return ::fabs(val1 - val2) / denom < sqrt_eps;
} 

}	// end namespace

bool DenseLinAlgPack::comp(const DVectorSlice& vs1, const DVectorSlice& vs2) {
  DVectorSlice::const_iterator
    vs1_itr = vs1.begin(),
    vs2_itr = vs2.begin();
  for(; vs1_itr != vs1.end() && vs2_itr != vs2.end(); ++vs1_itr, ++vs2_itr)
    if( !_comp(*vs1_itr,*vs2_itr) ) return false;
  return true;
}

bool DenseLinAlgPack::comp(const DVectorSlice& vs, value_type alpha) {
  DVectorSlice::const_iterator vs_itr = vs.begin();
  for(; vs_itr != vs.end(); ++vs_itr)
    if( !_comp(*vs_itr,alpha) ) return false;
  return true;
}

bool DenseLinAlgPack::comp(const DMatrixSlice& gms1, BLAS_Cpp::Transp trans1
  , const DMatrixSlice& gms2, BLAS_Cpp::Transp trans2)
{
  for(size_type i = 1; i < my_min(gms1.cols(),gms2.cols()); ++i)
    if( !comp( col(gms1,trans1,i) , col( gms2, trans2, i ) ) ) return false;
  return true;
}

bool DenseLinAlgPack::comp(const DMatrixSlice& gms, value_type alpha)
{
  for(size_type i = 1; i < gms.cols(); ++i)
    if( !comp( gms.col(i) , alpha ) ) return false;
  return true;
}

bool DenseLinAlgPack::comp(const DMatrixSliceTriEle& tri_gms1, const DMatrixSliceTriEle& tri_gms2)
{
  using BLAS_Cpp::bool_to_trans;
  BLAS_Cpp::Transp
    trans1 = bool_to_trans(tri_gms1.uplo() != BLAS_Cpp::lower),
    trans2 = bool_to_trans(tri_gms2.uplo() != BLAS_Cpp::lower);

  const size_type n = tri_gms1.rows();	// same as cols()
  for(size_type i = 1; i < n; ++i) {
    if( !comp( col(tri_gms1.gms(),trans1,i)(i,n), col(tri_gms2.gms(),trans2,i)(i,n) ) )
      return false;
  }
  return true;
}

bool DenseLinAlgPack::comp(const DMatrixSliceTriEle& tri_gms1, value_type alpha)
{
  using BLAS_Cpp::bool_to_trans;
  BLAS_Cpp::Transp
    trans1 = bool_to_trans(tri_gms1.uplo() != BLAS_Cpp::lower);

  const size_type n = tri_gms1.rows();	// same as cols()
  for(size_type i = 1; i < n; ++i) {
    if( !comp( col(tri_gms1.gms(),trans1,i)(i,n), alpha ) )
      return false;
  }
  return true;
}

bool DenseLinAlgPack::comp_less(const DVectorSlice& vs, value_type alpha)
{
  DVectorSlice::const_iterator vs_itr = vs.begin();
  const value_type denom = my_max( ::fabs(alpha), 1.0 );
  for(; vs_itr != vs.end(); ++vs_itr)
    if( *vs_itr > alpha ) return false;
  return true;
}
