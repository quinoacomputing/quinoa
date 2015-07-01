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

#ifndef ACTIVE_SET_CHANGE_H
#define ACTIVE_SET_CHANGE_H

#include <iosfwd>

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Calculate the change in the active set and output change
  * if asked to.
  *
  * ToDo: Add more description of the output you get.
  *
  * @param	nu_k	[in] Multipliers for variable bounds for iteration k.
  * @param	num_km1	[in] Multipliers for variable bounds for iteration k-1
  * @param	olevel	[in] Specifies the output level.  We have:\\
  *				PRINT_NOTHING : No output is sent to out\\
  *				PRINT_ALGORITHM_STEPS :
  *					Just the number of additions
  *					and deletions to the active set and the total number
  *					of active constraints is output.\\ 
  *				PRINT_ACTIVE_SET : Enumerates
  *					which variable were added and dropped from the active set.\\
  * @param	num_adds
  *                 [out] Gives the total number of variables fixed at a bound
  *					added to the active set.
  * @param	num_drops
  *                 [out] Gives the total number of variables freed from a 
  *					bound and dropped from the active set.
  * @param	num_adds_indep
  *                 [out] Gives the number of independent variables fixed at a bound
  *					added to the active set.
  * @param	num_drops_indep
  *                 [out] Gives the number of independent variables freed from a 
  *					bound and dropped from the active set.
  * @param	out	[O] Target for output.
  */
void active_set_change(
  const SpVectorSlice& nu_k, const SpVectorSlice& nu_km1, Range1D indep
  ,EJournalOutputLevel olevel, std::ostream* out
  ,size_type* num_adds, size_type* num_drops
  ,size_type* num_active_indep, size_type* num_adds_indep, size_type* num_drops_indep
  ); 

}	// end namespace MoochoPack

#endif	// ACTIVE_SET_CHANGE_H
