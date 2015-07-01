#if 0

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

#include "AbstractLinAlgPack_sparse_bounds_diff.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

void AbstractLinAlgPack::imp_sparse_bnd_diff(
    int						sign
  , const SpVectorSlice		&sv
  , BLAS_Cpp::Uplo			uplo
  , const DVectorSlice			&v
  , DVectorSlice				*r
  )
{
  DenseLinAlgPack::Vp_V_assert_sizes(r->size(),sv.size());
  DenseLinAlgPack::VopV_assert_sizes(sv.size(),v.size());

  typedef DenseLinAlgPack::value_type value_type;
  const value_type
    inf = std::numeric_limits<value_type>::max();
  *r = ( uplo == BLAS_Cpp::upper ? inf : -inf );
  const AbstractLinAlgPack::SpVectorSlice::difference_type o = sv.offset();
  for( AbstractLinAlgPack::SpVectorSlice::const_iterator itr = sv.begin();
      itr != sv.end(); ++itr )
  {
    (*r)(itr->indice() + o) = itr->value();
  }
  DenseLinAlgPack::Vp_StV( r, -1.0, v );
  if( sign < 0 )
    DenseLinAlgPack::Vt_S( r, -1.0 );
}

#endif // 0
