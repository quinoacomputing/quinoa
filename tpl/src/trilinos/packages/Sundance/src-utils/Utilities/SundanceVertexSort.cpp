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

#include "SundanceVertexSort.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"

namespace Sundance
{
using Teuchos::Array;
using Teuchos::tuple;
using std::endl;


int mapExSideToMissingVertex(int dim, int exFaceID)
{
  static Array<int> map3D = tuple(2, 0, 1, 3);
  static Array<int> map2D = tuple(2, 0, 1);
  switch(dim)
  {
    case 3:
      return map3D[exFaceID];
    case 2:
      return map2D[exFaceID];
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(dim != 3, std::runtime_error,
        "dimension " << dim << " not supported "
        "in mapExSideToMissingVertex()");
  }
  return -1; // -Wall;
}

Array<int> exSideVertPos(int dim, int sideIndex)
{
  static Array<Array<int> > faceVerts = tuple<Array<int> >
    (
      tuple(0,1,3),
      tuple(1,2,3),
      tuple(0,2,3),
      tuple(0,1,2)
      );
  static Array<Array<int> > edgeVerts = tuple<Array<int> >
    (
      tuple(0,1),
      tuple(1,2),
      tuple(0,2)
      );
  
  if (dim==2) return edgeVerts[sideIndex];
  else if (dim==3) return faceVerts[sideIndex];
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return faceVerts[0]; // -Wall
}

Array<int> ufcSideVertPos(int dim, int f)
{
  static Array<Array<int> > faceVerts = tuple<Array<int> >
    (
      tuple(1,2,3),
      tuple(0,2,3),
      tuple(0,1,3),
      tuple(0,1,2)
      );
  static Array<Array<int> > edgeVerts = tuple<Array<int> >
    (
      tuple(1,2),
      tuple(0,2),
      tuple(0,1)
      );
  
  if (dim==2) return edgeVerts[f];
  else if (dim==3) return faceVerts[f];
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return faceVerts[0]; // -Wall

}

Array<int> ufcSide(int f, const Array<int>& verts)
{
  Array<int> rtn;
  for (int i=0; i<verts.size(); i++)
  {
    if (f==i) continue;
    rtn.append(verts[i]);
  }
  return rtn;
}

Array<int> exSide(int f, const Array<int>& verts)
{
  Array<int> rtn;
  Array<int> exPos = exSideVertPos(verts.size()-1, f);
  for (int i=0; i<exPos.size(); i++)
  {
    rtn.append(verts[exPos[i]]);
  }
  return rtn;
}

void vertexSort(Array<int>& verts, int* key)
{
  Array<std::pair<int,int> > A(verts.size());
  for (int i=0; i<verts.size(); i++) 
  {
    A[i].first = verts[i];
    A[i].second = i;
  }

  insertionSort(A);

  *key = 0;
  int base = verts.size();
  int p = 1;
  for (int i=verts.size()-1; i>=0; i--) 
  {
    verts[i] = A[i].first;
    *key += p * A[i].second;
    p *= base;
  }
}

void getKeyedPerm(int key, Array<int>& digits)
{
  int r = key;
  int B = digits.size(); 
  for (int b=B-1; b>=0; b--)
  {
    int B_toThe_b = iPow(B,b);
    digits[B-b-1] = r/B_toThe_b;
    r = r - digits[B-b-1]*B_toThe_b;
  }
} 

/** Compute base^N */
int iPow(int base, int n)
{
  int rtn = 1;
  for (int i=0; i<n; i++) rtn = rtn*base;
  return rtn;
}


int exVertPosToUFCVertPos(int meshDim, int permKey, int exVertPos)
{
  Array<int> digits(meshDim+1);
  getKeyedPerm(permKey, digits);
  for (int i=0; i<digits.size(); i++) 
  {
    if (digits[i] == exVertPos) return i;
  }
  /* loop should never be completed w/o return */
  Out::os() << "vert=" << exVertPos << ", key=" << permKey 
            << ", digits = " << digits << endl;
  return -1; // indicates an error. Should be checked by client
}

int exFacetIndexToUFCFacetIndex(int meshDim, int permKey,
  int exFacetID)
{
  int missingExVert = mapExSideToMissingVertex(meshDim, exFacetID);
  int missingUFCVert = exVertPosToUFCVertPos(meshDim, permKey, missingExVert);
  return missingUFCVert; // returns -1 if error
}


int ufcFacetIndexToExFacetIndex(int meshDim, int ufcFacetIndex)
{
  static Array<int> map2D = tuple(1,2,0);
  static Array<int> map3D = tuple(1,2,0,3);

  if (meshDim==2)
  {
    return map2D[ufcFacetIndex];
  }
  else if (meshDim==3)
  {
    return map3D[ufcFacetIndex];
  }
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPT( meshDim!=3 && meshDim!=2 );
  }
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return -1; // -Wall
}

}



