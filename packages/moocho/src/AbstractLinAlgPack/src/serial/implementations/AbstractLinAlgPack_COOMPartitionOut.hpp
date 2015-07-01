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
// Output stream operators for Partition<> and TransposedPartition<>

#ifndef COOM_PARTITION_OUT_H
#define COOM_PARTITION_OUT_H

#include "AbstractLinAlgPack_COOMatrixTmplOutFunc.hpp"

namespace AbstractLinAlgPack {

/** \brief Partition<> output stream operator.
  *
  * This operator function calls the function output_COOM(os,part,0).
  */
template <class T_Indice, class T_Value>
inline std::ostream& operator<<(std::ostream& os
  , const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& part) {
  return output_COOM(os,part,0);
}

/** \brief TransposedPartition<> output stream operator.
  *
  * This operator function calls the function output_COOM(os,trans_part,0).
  */
template <class T_Indice, class T_Value>
inline std::ostream& operator<<(std::ostream& os
  , const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& trans_part)
{
  return output_COOM(os,trans_part,0);	
}

}	// end namespace AbstractLinAlgPack

#endif // VECTOROUT_H
