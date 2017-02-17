
// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_util_parallel_GenerateParallelUniqueIDs_hpp
#define stk_util_parallel_GenerateParallelUniqueIDs_hpp

#include "stk_util/parallel/ParallelVectorConcat.hpp" 
#include "stk_util/parallel/Parallel.hpp" 
#include "stk_util/parallel/ParallelIndexGapFinder.hpp"
#include <vector> 
#include <algorithm> 
#include <stdexcept>     
#include <string>                       // for string
#include <sstream>   
#include "mpi.h"  
#include <assert.h>

namespace stk {




  //--------------------------------------------------------------------------------------------
  //
  //  Given an input set of existing ids, generate the requested number of new ids.  Subject to the following
  //  requirements:
  //
  //  1) The new ids will not overlap any existing ids, this will hold globally in parallel
  //  2) All new ids will be less than or equal to the provided maxAllowedId
  //
  //  Note this if this routine is called on different parallel decompositions it may return different id sets both in 
  //  order and in content.
  //
  //  Arguments:
  //  existingIds    : input  : List of ids already in use.  Note, the id list may be unique to a processor 
  //                            or different processor lists may overlap
  //  numNewIdsLocal : input  : Number of new ids requested on this processor
  //  maxAllowedId   : input  : Max id value allowed to be returned.  If a valid set of ids cannot be found
  //                            below maxAllowedId routine will return a failure code.
  //
  //  Returns:  Vector of ids created for the local processor.  Will generally throw if encountering any internal 
  //            errors.  Note, the error throws are NOT guaranteed to be parallel consistent and are usually unrecoverable.
  //
  //  Assumptions and Usage Guidelines
  //    This algorithm assumes current existingIds are relatively dense packed.  The specific requirement for
  //    success is 'maxAllowedId - global_max(existingIds) > orderArray.size()'
  //

  std::vector<uint64_t> generate_parallel_unique_ids(const uint64_t maxAllowedId,
                                                     const std::vector<uint64_t>& existingIds,
                                                     uint64_t numNewIdsLocal,
                                                     ParallelMachine comm);
}

#endif

