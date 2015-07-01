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

#include "SundanceMeshBase.hpp"
#include "SundanceMesh.hpp"
#include "PlayaExceptions.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;



MeshBase::MeshBase(int dim, const MPIComm& comm, 
  const MeshEntityOrder& order) 
  : dim_(dim), 
    comm_(comm),
    order_(order),
    reorderer_(Mesh::defaultReorderer().createInstance(this)),
    validWeights_(true),
    specialWeights_(),
    curvePoints_Are_Valid_(true),
    nrCurvesForIntegral_(0),
    curvePoints_() ,
    curveDerivative_(),
    curveNormal_(),
    curveID_to_ArrayIndex_()
{ specialWeights_.resize(dim_);
  curvePoints_.resize(0);
  curveDerivative_.resize(0);
  curveNormal_.resize(0);
}

const int* MeshBase::elemZeroFacetView(int maxCellLID) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "MeshBase::elemZeroFacetView() is deprecated");
  return (const int*) 0;
}


Point MeshBase::centroid(int cellDim, int cellLID) const
{
  if (cellDim==0) return nodePosition(cellLID);
  int dummy;
  Point x = nodePosition(facetLID(cellDim, cellLID, 0, 0, dummy));
  int nf = numFacets(cellDim, cellLID, 0);
  for (int f=1; f<nf; f++) 
    x += nodePosition(facetLID(cellDim, cellLID, 0, f, dummy));
  return x / ((double) nf);
}

void MeshBase::outwardNormals(
  const Array<int>& cellLIDs,
  Array<Point>& outwardNormals
  ) const 
{
  int D = spatialDim();
  outwardNormals.resize(cellLIDs.size());
  for (int c=0; c<cellLIDs.size(); c++)
  {
    int f=-1;
    TEUCHOS_TEST_FOR_EXCEPTION(numMaxCofacets(D-1, cellLIDs[c]) > 1, 
      std::runtime_error,
      "cell #" << cellLIDs[c] << " is not a boundary cell");
    int maxLID = maxCofacetLID(D-1, cellLIDs[c], 0, f);
    Point cInterior = centroid(D, maxLID);
    Point cBdry = centroid(D-1, cellLIDs[c]);
    Point q = cBdry - cInterior;
    Point s;
    if (D==1) 
    {
      s = Point(1.0);
    }
    else if (D==2)
    {
      Point A = nodePosition(facetLID(D-1, cellLIDs[c], 0, 0, f));
      Point B = nodePosition(facetLID(D-1, cellLIDs[c], 0, 1, f));
      Point t = B - A;
      s = Point(-t[1], t[0]);
    }
    else 
    {
      Point A = nodePosition(facetLID(D-1, cellLIDs[c], 0, 0, f));
      Point B = nodePosition(facetLID(D-1, cellLIDs[c], 0, 1, f));
      Point C = nodePosition(facetLID(D-1, cellLIDs[c], 0, 2, f));
      s = cross(B-A, C-A);
    }
    if (q * s > 0.0)
    {
      outwardNormals[c] = s/::sqrt(s*s);
    }
    else
    {
      outwardNormals[c] = -s/::sqrt(s*s);
    }
  }
}


void  MeshBase::tangentsToEdges(
  const Array<int>& cellLIDs,
  Array<Point>& tangentVectors
  ) const 
{
  TEUCHOS_TEST_FOR_EXCEPT(spatialDim() <= 1);

  tangentVectors.resize(cellLIDs.size());

  for (int c=0; c<cellLIDs.size(); c++)
  {
    int fOrient=1;
    Point A = nodePosition(facetLID(1, cellLIDs[c], 0, 0, fOrient));
    Point B = nodePosition(facetLID(1, cellLIDs[c], 0, 1, fOrient));
    Point t = B - A;
    tangentVectors[c] = t/(sqrt(t*t));
  }
}




void MeshBase::getFacetArray(int cellDim, int cellLID, int facetDim, 
  Array<int>& facetLIDs,
  Array<int>& facetOrientations) const
{
  int nf = numFacets(cellDim, cellLID, facetDim);
  facetLIDs.resize(nf);
  facetOrientations.resize(nf);
  for (int f=0; f<nf; f++) 
  {
    facetLIDs[f] = facetLID(cellDim, cellLID, facetDim, f, 
      facetOrientations[f]);
  }
}


void MeshBase::getLabels(int cellDim, const Array<int>& cellLID, 
  Array<int>& labels) const
{
  labels.resize(cellLID.size());
  for (int i=0; i<cellLID.size(); i++) labels[i] = label(cellDim, cellLID[i]);
}


void MeshBase::getMaxCofacetLIDs(
  const Array<int>& cellLIDs,
  MaximalCofacetBatch& cofacets) const
{
  TEUCHOS_TEST_FOR_EXCEPT(cellLIDs.size()==0);
  int d = spatialDim() - 1;
  int nc = numMaxCofacets(d, cellLIDs[0]);  

  cofacets.reset(cellLIDs.size(), nc);

  for (int i=0; i<cellLIDs.size(); i++) 
  {
    int f1;
    int cofacet1 = maxCofacetLID(d, cellLIDs[i], 0, f1);
    if (nc==1) cofacets.addSingleCofacet(i, cofacet1, f1);
    else
    {
      int f2;
      int cofacet2 = maxCofacetLID(d, cellLIDs[i], 1, f2);
      cofacets.addTwoCofacets(i, cofacet1, f1, cofacet2, f2);
    }
  }
}

void MeshBase::getLIDsForLabel(int cellDim, int labelToCheck, 
  Array<int>& cellLIDs) const
{
  cellLIDs.resize(0);
  int nc = numCells(cellDim); 

  for (int c=0; c<nc; c++)
  {
    int lc = label(cellDim, c);
    if (lc==labelToCheck) cellLIDs.append(c);
  }
}

Set<int> MeshBase::getAllLabelsForDimension(int cellDim) const
{
  int nc = numCells(cellDim); 
  Set<int> rtn;

  for (int c=0; c<nc; c++)
  {
    rtn.put(label(cellDim, c));
  }
  return rtn;
}

// ===================== storing special weights for specfic cells  ======================

bool MeshBase::hasSpecialWeight(int dim, int cellLID) const {
	if (specialWeights_[dim-1].containsKey(cellLID))
        return ((specialWeights_[dim-1].get(cellLID).size() > 0));
	else
		return false;
}

void MeshBase::setSpecialWeight(int dim, int cellLID, Array<double>& w) const {
	specialWeights_[dim-1].put(cellLID,w);
}

void MeshBase::getSpecialWeight(int dim, int cellLID, Array<double>& w) const {
	w = specialWeights_[dim-1].get(cellLID);
}

/** deletes all special weights so those have to be recreated*/
void MeshBase::flushSpecialWeights() const {
	// delete all the special weights for all cells
	for (int d = 0 ; d < dim_ ; d++){
		// reset the map to a new and empty one
		specialWeights_[d] = Sundance::Map< int , Array<double> >();
	}
	validWeights_ = true;
}

// ===================== storing curve intersection/quadrature points ======================

void MeshBase::flushCurvePoints() const {
	// all the curve points should be deleted
	// first reset the Map between curve ID and index
	curveID_to_ArrayIndex_ = Sundance::Map< int , int >();
	nrCurvesForIntegral_ = 0;
	curvePoints_Are_Valid_ = true;
	for (int c = 0 ; c < curvePoints_.size() ; c++ ){
		// here delete each Map which existed before
		curvePoints_[c] = Sundance::Map< int , Array<Point> >();
		curveDerivative_[c] = Sundance::Map< int , Array<Point> >();
		curveNormal_[c] = Sundance::Map< int , Array<Point> >();
	}
}

bool MeshBase::hasCurvePoints(int maxCellLID , int curveID) const {
   	if (curvePoints_[mapCurveID_to_Index(curveID)].containsKey(maxCellLID))
		return ( curvePoints_[mapCurveID_to_Index(curveID)].get(maxCellLID).size() > 0 );
   	else
   		return false;
}

void MeshBase::setCurvePoints(int maxCellLID, int curveID ,
		Array<Point>& points , Array<Point>& derivs , Array<Point>& normals) const {

	  Tabs tabs;
	  int verbo = 0;
	  SUNDANCE_MSG3(verbo, tabs << "MeshBase::setCurvePoints , nr:" << nrCurvesForIntegral_);
   	  curvePoints_[mapCurveID_to_Index(curveID)].put( maxCellLID , points );
   	  curveDerivative_[mapCurveID_to_Index(curveID)].put( maxCellLID , derivs );
   	  curveNormal_[mapCurveID_to_Index(curveID)].put( maxCellLID , normals );

}

void MeshBase::getCurvePoints(int maxCellLID, int curveID ,
		Array<Point>& points , Array<Point>& derivs , Array<Point>& normals) const {

	  Tabs tabs;
	  int verbo = 0;
	  SUNDANCE_MSG3(verbo, tabs << "MeshBase::getCurvePoints , nr:" << nrCurvesForIntegral_);
   	  points = curvePoints_[mapCurveID_to_Index(curveID)].get( maxCellLID );
   	  derivs = curveDerivative_[mapCurveID_to_Index(curveID)].get( maxCellLID );
   	  normals = curveNormal_[mapCurveID_to_Index(curveID)].get( maxCellLID );

}

int MeshBase::mapCurveID_to_Index(int curveID) const {

	 Tabs tabs;
	 int verbo = 0;

	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index curveID:" << curveID);
     if (curveID_to_ArrayIndex_.containsKey(curveID)){
    	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index value found for curveID:" << curveID << " ret:" << curveID_to_ArrayIndex_.get(curveID));
       	 return curveID_to_ArrayIndex_.get(curveID);
     } else {
    	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index create new :" << nrCurvesForIntegral_);
       	 curveID_to_ArrayIndex_.put( curveID , nrCurvesForIntegral_ );
       	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index , increment ");
       	 nrCurvesForIntegral_++;
    	 curvePoints_.resize(nrCurvesForIntegral_);
    	 curveDerivative_.resize(nrCurvesForIntegral_);
    	 curveNormal_.resize(nrCurvesForIntegral_);
    	 SUNDANCE_MSG3(verbo, tabs << "MeshBase::mapCurveID_to_Index create new :" << nrCurvesForIntegral_);
       	 return nrCurvesForIntegral_-1;
     }

}
