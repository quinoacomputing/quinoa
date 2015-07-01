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

#ifndef VECTOR_CHANGE_STATS_H
#define VECTOR_CHANGE_STATS_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Compute statistics for change in a vector.
  *
  * Given two vectors x and d where we wish to generate statistics
  * for the update x+d this function computes the following
  * quantitines:
  *
  * max( |d(i)|/(1+|x(i)|), i=1...n )
  *   => #max_k# = k, #max_term# = |d(k)|/(1+|x(k)|) <= 1\\
  * #min( |d(i)|/(1+|x(i)|), i=1...n ) => #min_k# = k, #min_term# = |d(k)|/(1+|x(k)|)#\\
  * #average( |d(i)|/(1+|x(i)|), i=1...10 )# => #av_term#\\
  * 
  * The purpose of generating these satistics is to determine
  * by how much x+d differs from x.
  *
  * If |d(i)|/|x(i)| < mach_eps with x(i) > 0 then we know that d(i) will
  * be lost when added to x(i) so x(i) + d(i) == x(i).
  *
  */
void vector_change_stats( const DVectorSlice& x, const DVectorSlice& d
  , value_type* max_term, size_type* max_k
  , value_type* min_term, size_type* min_k
  , value_type* av_term );

}	// end namespace ConstrainedOptPack

#endif	// VECTOR_CHANGE_STATS_H
