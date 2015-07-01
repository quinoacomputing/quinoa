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

#include "SundancePeriodicSingleCellMesh1D.hpp"

#include "SundanceCellJacobianBatch.hpp"
#include "SundanceMaximalCofacetBatch.hpp"
#include "SundanceDebug.hpp"
#include "SundanceOut.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceObjectWithVerbosity.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace std;
using Playa::MPIComm;
using Playa::MPIContainerComm;


PeriodicSingleCellMesh1D::PeriodicSingleCellMesh1D(double xMin, double xMax)
  : MeshBase(1, MPIComm::self(), ExodusMeshOrder),
    xMin_(xMin),
    xMax_(xMax),
    vert_(0),
    labels_(2)
{
  labels_[0] = 0;
  labels_[1] = 0;
}

int PeriodicSingleCellMesh1D::numCells(int cellDim) const
{
  switch(cellDim)
  {
    case 0 :
      return 1;
    case 1:
      return 1;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return -1; // -Wall
}


Point PeriodicSingleCellMesh1D::nodePosition(int i) const 
{
  TEUCHOS_TEST_FOR_EXCEPT(i != 0);
  return xMin_;
}


const double* PeriodicSingleCellMesh1D::nodePositionView(int i) const 
{
  TEUCHOS_TEST_FOR_EXCEPT(i != 0);
  return &(xMin_);
}

void PeriodicSingleCellMesh1D::getJacobians(int cellDim, const Array<int>& cellLID,
    CellJacobianBatch& jBatch) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
    "cellDim=" << cellDim 
    << " is not in expected range [0, " << spatialDim()
    << "]");

  jBatch.resize(cellLID.size(), spatialDim(), cellDim);

  int nCells = cellLID.size();

  TEUCHOS_TEST_FOR_EXCEPT(nCells != 1);

  if (cellDim==0)
  {
    for (int i=0; i<nCells; i++)
    {
      double* detJ = jBatch.detJ(i);
      *detJ = 1.0;
    }
  }
  else
  {
    for (int i=0; i<nCells; i++)
    {
      double* J = jBatch.jVals(i);
      J[0] = fabs(xMin_-xMax_);
    }
  }
}


void PeriodicSingleCellMesh1D::getCellDiameters(int cellDim, const Array<int>& cellLID,
  Array<double>& cellDiameters) const
{
  cellDiameters.resize(1);
  
  TEUCHOS_TEST_FOR_EXCEPT(cellDim != 1);

  cellDiameters[0] = ::fabs(xMax_-xMin_);
}

void PeriodicSingleCellMesh1D::pushForward(int cellDim, const Array<int>& cellLID,
    const Array<Point>& refQuadPts,
    Array<Point>& physQuadPts) const
{
  TEUCHOS_TEST_FOR_EXCEPT(cellDim < 0 || cellDim > 1);

  TEUCHOS_TEST_FOR_EXCEPT(cellLID.size() > 1);

  if (cellDim==1)
  {
    if (physQuadPts.size() > 0) physQuadPts.resize(0);
    physQuadPts.reserve(refQuadPts.size() * cellLID.size());
    
    for (int i=0; i<cellLID.size(); i++)
    {
      double h = xMax_ - xMin_;
      for (int q=0; q<refQuadPts.size(); q++)
      {
        physQuadPts.append(xMin_ + refQuadPts[q][0] * h);
      }
    }
  }
  else
  {
    for (int i=0; i<cellLID.size(); i++)
    {
      physQuadPts.append(xMin_);
    }
  }
}

int PeriodicSingleCellMesh1D::numFacets(int cellDim, int cellLID,
    int facetDim) const
{
  TEUCHOS_TEST_FOR_EXCEPT(cellLID != 0);
  if (cellDim == 1 && facetDim==0) return 1;
  return 0;
}


    
int PeriodicSingleCellMesh1D::facetLID(int cellDim, int cellLID,
  int facetDim, int facetIndex,
  int& facetOrientation) const
{
  TEUCHOS_TEST_FOR_EXCEPT(cellLID < 0 || cellLID >= 1);

  TEUCHOS_TEST_FOR_EXCEPT(cellDim != 1);
  TEUCHOS_TEST_FOR_EXCEPT(facetDim != 0);
  TEUCHOS_TEST_FOR_EXCEPT(facetIndex < 0);
  TEUCHOS_TEST_FOR_EXCEPT(facetIndex > 1);

  return vert_;
}

void PeriodicSingleCellMesh1D::getFacetLIDs(int cellDim,
    const Array<int>& cellLID,
    int facetDim,
    Array<int>& facetLID,
    Array<int>& facetSign) const
{
  facetLID.resize(2*cellLID.size());
  facetSign.resize(2*cellLID.size());

  for (int i=0; i<cellLID.size(); i++) 
  {
    facetLID[2*i] = this->facetLID(cellDim, cellLID[i], facetDim, 0, facetSign[2*i]);
    facetLID[2*i+1] = this->facetLID(cellDim, cellLID[i], facetDim, 1, facetSign[2*i]);
  }
}


const int* PeriodicSingleCellMesh1D::elemZeroFacetView(int cellLID) const
{
  TEUCHOS_TEST_FOR_EXCEPT(cellLID != 0);
  return &vert_;
}


int PeriodicSingleCellMesh1D::numMaxCofacets(int cellDim, int cellLID) const
{
  TEUCHOS_TEST_FOR_EXCEPT(cellDim != 0);
  return 1;
}

int PeriodicSingleCellMesh1D::maxCofacetLID(int cellDim, int cellLID,
    int cofacetIndex,
    int& facetIndex) const
{
  TEUCHOS_TEST_FOR_EXCEPT(cellDim != 0 || cellLID != 0);
  facetIndex = 0;
  return 0;
}

void PeriodicSingleCellMesh1D::getMaxCofacetLIDs(const Array<int>& cellLIDs,
    MaximalCofacetBatch& cofacets) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

void PeriodicSingleCellMesh1D::getCofacets(int cellDim, int cellLID,
    int cofacetDim, Array<int>& cofacetLIDs) const
{
  TEUCHOS_TEST_FOR_EXCEPT(cellDim != 0);
  TEUCHOS_TEST_FOR_EXCEPT(cofacetDim != 1);

  cofacetLIDs.resize(1);
  cofacetLIDs[0] = 0;
}

int PeriodicSingleCellMesh1D::mapGIDToLID(int cellDim, int globalIndex) const
{
  return globalIndex;
}

bool PeriodicSingleCellMesh1D::hasGID(int cellDim, int globalIndex) const
{
  return globalIndex==0;
}

int PeriodicSingleCellMesh1D::mapLIDToGID(int cellDim, int localIndex) const
{
  return localIndex;
}

CellType PeriodicSingleCellMesh1D::cellType(int cellDim) const
{
  if (cellDim==0) return PointCell;
  else if (cellDim==1) return LineCell;
  else return NullCell;
}

int PeriodicSingleCellMesh1D::label(int cellDim, int cellLID) const
{
  return labels_[cellDim];
}

void PeriodicSingleCellMesh1D::getLabels(int cellDim, const Array<int>& cellLID,
  Array<int>& labels) const
{
  labels.resize(cellLID.size());
  for (int i=0; i<cellLID.size(); i++)
  {
    labels[i] = labels_[cellDim];
  }
}

Set<int> PeriodicSingleCellMesh1D::getAllLabelsForDimension(int cellDim) const
{
  Set<int> rtn;

  rtn.put(labels_[cellDim]);
    
  return rtn;
}

void PeriodicSingleCellMesh1D::setLabel(int cellDim, int cellLID, int label)
{
  labels_[cellDim] = label;
}


void PeriodicSingleCellMesh1D::getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const
{
  cellLIDs.resize(0);
  if (label != labels_[cellDim]) cellLIDs.append(0);
}
