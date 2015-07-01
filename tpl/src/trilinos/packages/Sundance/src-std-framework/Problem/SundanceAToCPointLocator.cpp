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

#include "SundanceAToCPointLocator.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceGeomUtils.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <queue>

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





static Time& pointLocatorCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("point locator ctor"); 
  return *rtn;
}

AToCPointLocator::AToCPointLocator(const Mesh& mesh, 
                                   const CellFilter& subdomain,
                                   const std::vector<int>& nx)
  : dim_(mesh.spatialDim()),
    mesh_(mesh),
    nFacets_(mesh.numFacets(dim_, 0, 0)),
    nx_(nx),
    low_(nx.size(), 1.0/ScalarTraits<double>::sfmin()),
    high_(nx.size(), -1.0/ScalarTraits<double>::sfmin()),
    dx_(nx.size()),
    table_(),
    subdomain_(subdomain),
    neighborSet_()
{
  TimeMonitor timer(pointLocatorCtorTimer());
  
  /* allocate the neighbor set table */
  neighborSet_.resize(mesh.numCells(dim_));

  /* first pass to find bounding box */
  CellSet cells = subdomain.getCells(mesh);
  
  for (CellIterator i = cells.begin(); i!= cells.end(); i++)
    {
      int cellLID = *i;
      Array<int> facetLIDs;
      Array<int> facetOri;
      mesh.getFacetArray(dim_, cellLID, 0, facetLIDs, facetOri);
      for (int f=0; f<facetLIDs.size(); f++)
        {
          Point x = mesh.nodePosition(facetLIDs[f]);
          for (int d=0; d<dim_; d++)
            {
              if (x[d] < low_[d]) low_[d] = x[d];
              if (x[d] > high_[d]) high_[d] = x[d];
            }
        }
    }

  /* fudge the bounding box */
  for (int d=0; d<dim_; d++)
    {
      low_[d] -= 0.01 * (high_[d] - low_[d]);
      high_[d] += 0.01 * (high_[d] - low_[d]);
    }

  /* second pass to put cells in lookup table */

  int s = 1;
  for (int d=0; d<dim_; d++) 
    {
      dx_[d] = (high_[d] - low_[d])/nx_[d];
      s *= nx_[d];
    }


  table_ = rcp(new Array<int>(s, -1));


  Array<int> lowIndex;
  Array<int> highIndex;
  for (CellIterator i = cells.begin(); i!= cells.end(); i++)
    {
      int cellLID = *i;
      getGridRange(mesh, dim_, cellLID, lowIndex, highIndex);
      if (dim_==2)
        {
          for (int ix=lowIndex[0]; ix<=highIndex[0]; ix++)
            {
              for (int iy=lowIndex[1]; iy<=highIndex[1]; iy++)
                {
                  int index = nx_[1]*ix + iy;
                  (*table_)[index] = cellLID;
                }
            }
        }
      else
        {
          TEUCHOS_TEST_FOR_EXCEPT(true);
        }
    }
}

int AToCPointLocator::getGridIndex(const double* x) const 
{
  int index = 0;
  for (int d=0; d<dim_; d++) 
    {
      double r = (x[d] - low_[d])/dx_[d];
      int ix = (int) floor(r);
      index = index*nx_[d] + ix;
    }

  return index;
}


void AToCPointLocator::getGridRange(const Mesh& mesh, int cellDim, int cellLID,
                                    Array<int>& lowIndex, Array<int>& highIndex) const
{
  Array<int> facetLIDs;
  Array<int> facetOri;
  mesh.getFacetArray(cellDim, cellLID, 0, facetLIDs, facetOri);

  lowIndex.resize(cellDim);
  highIndex.resize(cellDim);
  for (int d=0; d<cellDim; d++) 
    {
      highIndex[d] = -1;
      lowIndex[d] = 1000000; /* This magic number should be much bigger than any 
                              * reasonable index size in any one dimension */
    }

  for (int f=0; f<facetLIDs.size(); f++)
    {
      Point x = mesh.nodePosition(facetLIDs[f]);
      for (int d=0; d<cellDim; d++)
        {
          double r = (x[d] - low_[d])/dx_[d];
          int ix = (int) floor(r);
          if (ix < lowIndex[d]) lowIndex[d] = ix;
          if (ix > highIndex[d]) highIndex[d] = ix;
        }
    }
}
                 

void AToCPointLocator::fillMaximalNeighbors(int cellLID, 
                                              const int* facetLID) const
{
  if (neighborSet_[cellLID].get()==0)
    {
      neighborSet_[cellLID] = rcp(new Set<int>());

      for (int f=0; f<nFacets_; f++)
        {
          Array<int> cofacets;
          mesh_.getCofacets(0, facetLID[f], dim_, cofacets);
          for (int c=0; c<cofacets.size(); c++)
            {
              if (cofacets[c] != cellLID) neighborSet_[cellLID]->put(cofacets[c]);
            }
        }
    }
}


bool AToCPointLocator::cellContainsPoint(int cellLID, 
                                           const double* x, 
                                           const int* facetLID) const
{
  if (dim_==2)
    {
      const double* A = mesh_.nodePositionView(facetLID[0]);
      const double* B = mesh_.nodePositionView(facetLID[1]);
      const double* C = mesh_.nodePositionView(facetLID[2]);
      /* first determine whether the three points of the triangle
       * are in ccw or cw order. */
      double sign = orient2D(A, B, C);
      if (sign > 0.0)
        {
          double s1 = orient2D(A, B, x);
          double s2 = orient2D(B, C, x);
          double s3 = orient2D(C, A, x);
          if (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) return true;
          return false;
        }
      else
        {
          double s1 = orient2D(A, C, x);
          double s2 = orient2D(C, B, x);
          double s3 = orient2D(B, A, x);
          if (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) return true;
          return false;
        }
    }
  else if (dim_==3)
    {
      TEUCHOS_TEST_FOR_EXCEPT(true);
      return false; // -Wall
    }
  else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(dim_<=0 || dim_>3, std::runtime_error,
                         "invalid point dimension " << dim_);
      return false; // -Wall
    }
}

bool AToCPointLocator::cellContainsPoint(int cellLID, 
                                         const double* x, 
                                         const int* facetLID,
                                         double* localCoords) const
{
  if (dim_==2)
    {
      const double* A = mesh_.nodePositionView(facetLID[0]);
      const double* B = mesh_.nodePositionView(facetLID[1]);
      const double* C = mesh_.nodePositionView(facetLID[2]);
      /* first determine whether the three points of the triangle
       * are in ccw or cw order. */
      double sign = orient2D(A, B, C);
      if (sign > 0.0)
        {
          double s1 = orient2D(A, B, x);
          double s2 = orient2D(B, C, x);
          double s3 = orient2D(C, A, x);
          if (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) 
            {
              double bax = B[0] - A[0];
              double bay = B[1] - A[1];
              double cax = C[0] - A[0];
              double cay = C[1] - A[1];
              double delta = bax*cay - bay*cax;
              
              double xax = x[0] - A[0];
              double xay = x[1] - A[1];

              localCoords[0] = (cay*xax - cax*xay)/delta;
              localCoords[1] = (bax*xay - bay*xax)/delta;
              return true;
            }
          return false;
        }
      else
        {
          double s1 = orient2D(A, C, x);
          double s2 = orient2D(C, B, x);
          double s3 = orient2D(B, A, x);
          if (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) 
            {
              /* swap B and C if cell is CW oriented */
              std::cout << "swapping!" << std::endl;
              double bax = C[0] - A[0];
              double bay = C[1] - A[1];
              double cax = B[0] - A[0];
              double cay = B[1] - A[1];
              double delta = bax*cay - bay*cax;
              
              double xax = x[0] - A[0];
              double xay = x[1] - A[1];

              localCoords[0] = (cay*xax - cax*xay)/delta;
              localCoords[1] = (bax*xay - bay*xax)/delta;
              return true;
            }
          return false;
        }
    }
  else if (dim_==3)
    {
      TEUCHOS_TEST_FOR_EXCEPT(true);
      return false; // -Wall
    }
  else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(dim_<=0 || dim_>3, std::runtime_error,
                         "invalid point dimension " << dim_);
      return false; // -Wall
    }
}

int AToCPointLocator::findEnclosingCell(int initialGuessLID,
                                        const double* x) const
{
  std::queue<int> Q;
  Set<int> repeats;
  int d = mesh_.spatialDim();
  Array<int> facets(mesh_.numFacets(d, 0, 0));
  Array<int> facetOr(facets.size());

  Q.push(initialGuessLID);

  while (!Q.empty())
    {
      int next = Q.front();
      Q.pop();
      if (repeats.contains(next)) continue;
      
//      const int* facets = mesh_.elemZeroFacetView(next);
      mesh_.getFacetArray(d, next, 0, facets, facetOr);
      if (cellContainsPoint(next, x, &(facets[0]))) return next;
      repeats.put(next);
      
      fillMaximalNeighbors(next, &(facets[0]));
      
      for (Set<int>::const_iterator 
             i=neighborSet_[next]->begin(); i!=neighborSet_[next]->end(); i++)
        {
          Q.push(*i);
        }
    }
  return -1; // no containing cell found
}




int AToCPointLocator::findEnclosingCell(int initialGuessLID,
                                        const double* x,
                                        double* xLocal) const
{
  std::queue<int> Q;
  Set<int> repeats;
  int d = mesh_.spatialDim();
  Array<int> facets(mesh_.numFacets(d, 0, 0));
  Array<int> facetOr(facets.size());

  Q.push(initialGuessLID);

  while (!Q.empty())
    {
      int next = Q.front();
      Q.pop();
      if (repeats.contains(next)) continue;
      
//      const int* facets = mesh_.elemZeroFacetView(next);
      mesh_.getFacetArray(d, next, 0, facets, facetOr);
      if (cellContainsPoint(next, x, &(facets[0]), xLocal)) return next;
      repeats.put(next);
      
      fillMaximalNeighbors(next, &(facets[0]));
      
      for (Set<int>::const_iterator 
             i=neighborSet_[next]->begin(); i!=neighborSet_[next]->end(); i++)
        {
          Q.push(*i);
        }
    }
  return -1; // no containing cell found
}


Point AToCPointLocator::makePoint(int dim, const double* x) 
{
  if (dim==1) return Point(x[0]);
  else if (dim==2) return Point(x[0], x[1]);
  else return Point(x[0], x[1], x[2]);
}
