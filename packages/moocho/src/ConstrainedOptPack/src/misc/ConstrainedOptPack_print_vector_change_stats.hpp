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

#ifndef PRINT_VECTOR_CHANGE_STATS_H
#define PRINT_VECTOR_CHANGE_STATS_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Compute statistics for change in a vector and output to a stream.
  *
  * Calls the function vector_change_stats(x,d,max_term,max_k,min_term,min_k
  ,av_term) then produces the following output to the given stream.
  *
  * max(|d(i)|/(1+|x(i)|)) => |d(max_k)|/(1+|x(max_k)|) = max_term \\
  * min(|d(i)|/(1+|x(i)|)) => |d(min_k)|/(1+|x(min_k)|) = min_term \\ 
  * average(|d(i)|/(1+|x(i)|)) = av_term \\
  *
  * Note that above the names x and d are replaced with their input
  * names #x_name# and #d_name# and max_term, max_k etc.
  * are replaced with their computed values.
  */
void print_vector_change_stats( const DVectorSlice& x, const char x_name[]
  , const DVectorSlice& d, const char d_name[], std::ostream& out );

}	// end namespace ConstrainedOptPack

#endif	// PRINT_VECTOR_CHANGE_STATS_H
