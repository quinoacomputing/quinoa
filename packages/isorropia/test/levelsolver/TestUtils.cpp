/*
//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER
*/

#include "TestUtils.h"

int parseNumThreads(const char *arg, std::vector<int> &numThreads) 
{
  int maxT, minT, modT;
  if (sscanf(arg,"%d:%d+%d",&minT,&maxT,&modT) == 3) {
    for (int nt=minT; nt<=maxT; nt += modT) numThreads.push_back(nt);
  }
  else if (sscanf(arg,"%d:%d*%d",&minT,&maxT,&modT) == 3) {
    for (int nt=minT; nt<=maxT; nt *= modT) numThreads.push_back(nt);
  }
  else if (sscanf(arg,"%d",&minT) == 1) {
    if (minT > 0) numThreads.push_back(minT);
    arg = strchr(arg,',');
    while (arg != NULL) {
      if (sscanf(arg+1,"%d", &minT) != 1) {
        break;
      }
      if (minT > 0) numThreads.push_back(minT);
      arg = strchr(arg+1,',');
    }
  }
  else {
    return -1;
  }
  return 0;
}
