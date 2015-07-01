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

#include "AbstractLinAlgPack_sparse_bounds.hpp"

/** \brief Count the number of sparse bounds where at least one bound is
  * finite.
  */
AbstractLinAlgPack::size_type
AbstractLinAlgPack::num_bounds( const SpVectorSlice& bl, const SpVectorSlice& bu )
{
  SpVectorSlice::const_iterator
    bl_itr			= bl.begin(),
    bl_itr_end		= bl.end(),
    bu_itr			= bu.begin(),
    bu_itr_end		= bu.end();
  size_type num_bounds = 0;
  while( bl_itr != bl_itr_end || bu_itr != bu_itr_end ) {
    if( ( bl_itr != bl_itr_end )
      && ( bu_itr == bu_itr_end || bl_itr->indice() + bl.offset() < bu_itr->indice() + bu.offset() ) )
    {
      // Only the lower bound is finite
      ++bl_itr;
    }
    else if( ( bu_itr != bu_itr_end )
      && ( bl_itr == bl_itr_end || bu_itr->indice() + bu.offset() < bl_itr->indice() + bl.offset()) )
    {
      // Only the upper bound is finite
      ++bu_itr;
    }
    else if(bl_itr->indice() == bu_itr->indice()) {
      // Both bounds exist.
      ++bl_itr;
      ++bu_itr; 
    }
    ++num_bounds;
  }
  return num_bounds;
}


#endif // 0
