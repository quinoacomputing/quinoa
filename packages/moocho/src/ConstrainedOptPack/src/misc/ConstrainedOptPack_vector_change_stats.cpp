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

#include "ConstrainedOptPack_vector_change_stats.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "Teuchos_ScalarTraits.hpp"

void ConstrainedOptPack::vector_change_stats(
    const DVectorSlice& x, const DVectorSlice& d
  , value_type* max_term, size_type* max_k
  , value_type* min_term, size_type* min_k
  , value_type* av_term )
{
  typedef Teuchos::ScalarTraits<value_type> ST;
  DenseLinAlgPack::VopV_assert_sizes( x.dim(), d.dim() );
  const value_type
    min_num		= std::numeric_limits<value_type>::min(),
    inf			= std::numeric_limits<value_type>::max();
  // Initialize statistics
  *max_term	= 0.0;
  *max_k		= 0;
  *min_term	= inf;
  *min_k		= 0;
  *av_term	= 0.0;
  // Compute statistics
  DVectorSlice::const_iterator
    x_itr	= x.begin(),
    d_itr	= d.begin();
  for( size_type i = 1; x_itr != x.end(); ++i, ++d_itr, ++x_itr ) {
    // Compute this ratio and make sure we are not dividing by zero.
    // We only care about ratios less than 1.0 and also deal
    // with the case that both x(i) and d(i) are zero (in which
    // case the ratio should be zero since x(i) == x(i) + d(i)).
    const value_type
      term = ST::magnitude(*d_itr) / ( 1 + ST::magnitude(*x_itr) );
    if( term > *max_term ) {
      *max_term	= term;
      *max_k		= i;
    }
    if( term < *min_term ) {
      *min_term	= term;
      *min_k		= i;
    }
    *av_term += term;
  }
  *av_term = *av_term / x.dim();
}
