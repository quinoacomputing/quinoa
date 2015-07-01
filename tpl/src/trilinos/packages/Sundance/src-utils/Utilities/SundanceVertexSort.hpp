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

#ifndef SUNDANCE_VERTEX_SORT_H
#define SUNDANCE_VERTEX_SORT_H


#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"


namespace Sundance
{
using Teuchos::Array;

/** */
template <class T> inline
void insertionSort(Teuchos::Array<T>& A)
{
  int N = A.size();
  for (int i=1; i<N; i++)
  {
    T val = A[i];
    int j = i-1;
    bool done = false;
    while ( !done )
    {
      if (A[j] > val)
      {
        A[j+1]=A[j];
        j=j-1;
        if ( j<0 ) done = true;
      }
      else done = true;
      A[j+1]=val;
    }
  }
}
 
/* Sort, returning the key associated with the permutation */
void vertexSort(Array<int>& verts, int* key);

/** Return the permutation that produces the specified key */
void getKeyedPerm(int key, Array<int>& digits);

/** Compute base^N */
int iPow(int base, int n);

/** */
int exFacetIndexToUFCFacetIndex(int meshDim, int permKey,
  int exFacetID);

/** */
int ufcFacetIndexToExFacetIndex(int meshDim, int ufcFacetID);

/** */
int exVertPosToUFCVertPos(int meshDim, int permKey, int exVertPos);

/** */
int mapExSideToMissingVertex(int dim, int exFaceID);

/** */
Array<int> exSideVertPos(int dim, int f);
/** */
Array<int> ufcSideVertPos(int dim, int f);

/** */
Array<int> ufcSide(int f, const Array<int>& verts);

/** */
Array<int> exSide(int f, const Array<int>& verts);


}


#endif
