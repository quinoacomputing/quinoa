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

#include "ConstrainedOptPack_ComputeMinMult.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

namespace {
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }
} // end namespace

ConstrainedOptPack::value_type
ConstrainedOptPack ::min_abs( const DVectorSlice& mu )
{
  if( !mu.dim() )
    return 0.0;
  value_type min = std::fabs(mu(1));
  for( DVectorSlice::const_iterator itr = mu.begin() + 1; itr != mu.end(); )
    min = my_min( min, std::fabs(*itr++) );
  return min;
}

ConstrainedOptPack::value_type
ConstrainedOptPack ::min_abs( const SpVectorSlice& mu )
{
  if( !mu.dim() )
    return 0.0;
  if( !mu.nz() )
    return 0.0;
  value_type min = std::fabs(mu.begin()->value());
  for( SpVectorSlice::const_iterator itr = mu.begin() + 1; itr != mu.end(); ++itr )
    min = my_min( min, std::fabs(itr->value()) );
  return min;
}
