/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */


#include "SundanceBasicVertexView.hpp"
#include "Teuchos_Utils.hpp"

using namespace Teuchos;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;

string VertexView::toString() const
{
  int* ptr = *base_ +  offset_*length_;
	string rtn="{";
	for (int i=0; i<length_; i++) 
		{
			rtn += Teuchos::toString(ptr[i]);
			if (i < length_-1) rtn += ", ";
		}
	rtn += "}";
	return rtn;
}

/*
 * Return a hash code for the vertex set. 
 */
int VertexView::hashCode() const
{
  int rtn = 0;
      int* p = *base_ + offset_*length_;

      for (int i=0; i<length_; i++)
        {
          rtn += p[i];
        }

      return rtn;

#ifdef BLARF
      // fails with sign int
  int* p = *base_ + offset_*length_;

  unsigned char* key = reinterpret_cast<unsigned char*>(p);
  unsigned int key_len = length_ * sizeof(int);

  /* Jenkins hash */

  unsigned int hash = 0;
  int i;
  
  for (i = 0; i < key_len; i++) {
    hash += key[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  return hash;
#endif
}
