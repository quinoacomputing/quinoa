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
/*
 * SundanceHNMesh2D.cpp
 *
 *  Created on: April 30, 2009
 *      Author: benk
 */

#include "SundanceHNMesh2D.hpp"

#include "SundanceMeshType.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceMaximalCofacetBatch.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceDebug.hpp"
#include "SundanceOut.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceCollectiveExceptionCheck.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace std;
using Playa::MPIComm;
using Playa::MPIContainerComm;

int HNMesh2D::offs_Points_x_[4] = {0, 1, 0, 1};

int HNMesh2D::offs_Points_y_[4] = {0, 0, 1, 1};

int HNMesh2D::edge_Points_localIndex[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };


HNMesh2D::HNMesh2D(int dim, const MPIComm& comm ,
	    const MeshEntityOrder& order)
: MeshBase(dim, comm , order),_dimension(dim), _comm(comm)
{
	setVerb(0);
	// get the number of processors
	nrProc_ = MPIComm::world().getNProc();
	myRank_ = MPIComm::world().getRank();
	//------ Point storage ----
	points_.resize(0);
    nrElem_.resize(3,0);
	nrElemOwned_.resize(3,0);
	//----- Facets -----
    cellsPoints_.resize(0);
    cellsEdges_.resize(0);
    isCellOut_.resize(0);
	edgePoints_.resize(0);
	edgeVertex_.resize(0);
	// ----- MaxCofacets ----
	edgeMaxCoF_.resize(0);
    pointMaxCoF_.resize(0);
    //------ Element (processor) ownership -----
	elementOwner_.resize(3); elementOwner_[0].resize(0); elementOwner_[1].resize(0); elementOwner_[2].resize(0);
	// ------------- different boundary ---------
    edgeIsProcBonudary_.resize(0);
    edgeIsMeshDomainBonudary_.resize(0);
    //---- hierarchical storage -----
    indexInParent_.resize(0);
    parentCellLID_.resize(0);
	cellLevel_.resize(0);
	isCellLeaf_.resize(0);
	// ---- "hanging" info storage ---
	isPointHanging_.resize(0);
	isEdgeHanging_.resize(0);
	// ---- hanging element and refinement (temporary) storage ---
	hangElmStore_.resize(2);
	hangElmStore_[0] = Hashtable< int, Array<int> >();
	hangElmStore_[1] = Hashtable< int, Array<int> >();
	refineCell_.resize(0);
    // set the leaf counter to zero
	nrCellLeafGID_ = 0; nrEdgeLeafGID_ = 0; nrVertexLeafGID_ = 0;
	nrVertexLeafLID_ = 0; nrCellLeafLID_ = 0; nrEdgeLeafLID_ = 0;
}

int HNMesh2D::numCells(int dim) const  {
	SUNDANCE_MSG3(verb(),"HNMesh2D::numCells(int dim):   dim:" << dim );
	switch (dim){
	case 0: return nrVertexLeafLID_;
	case 1: return nrEdgeLeafLID_;
	case 2: return nrCellLeafLID_;
	}
	return 0;
}

Point HNMesh2D::nodePosition(int i) const {
	SUNDANCE_MSG3(verb(),"HNMesh2D::nodePosition(int i)   i:"<< i);
    // point do not need leaf indexing
	return points_[vertexLeafToLIDMapping_[i]];
}

const double* HNMesh2D::nodePositionView(int i) const {
	SUNDANCE_MSG3(verb(),"HNMesh2D::nodePositionView(int i)   i:" << i);
	//SUNDANCE_VERB_HIGH("nodePosition(int i)");
	return &(points_[vertexLeafToLIDMapping_[i]][0]);;
}

void HNMesh2D::getJacobians(int cellDim, const Array<int>& cellLID,
                          CellJacobianBatch& jBatch) const
{
	  SUNDANCE_MSG3(verb(),"HNMesh2D::getJacobians  cellDim:"<<cellDim<<" _x:"<<_ofs_x<<" _y:"<<_ofs_y);
	  SUNDANCE_VERB_HIGH("getJacobians()");
	  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	  int nCells = cellLID.size();
	  int LID;
	  Point pnt(0.0,0.0);
	  jBatch.resize(cellLID.size(), spatialDim(), cellDim);
	  if (cellDim < spatialDim()) // they need the Jacobian of a lower dinemsional element
	  {
		   for (int i=0; i<nCells; i++)
		    {
		      double* detJ = jBatch.detJ(i);
		      switch(cellDim)
		      {
		        case 0: *detJ = 1.0;
		          break;
		        case 1:
				  LID = edgeLeafToLIDMapping_[cellLID[i]];
			      pnt = (points_[edgePoints_[LID][1]] - points_[edgePoints_[LID][0]]);
		          *detJ = sqrt(pnt * pnt); // the length of the edge
		        break;
		        default:
		          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "  "cellDim=" << cellDim << " in HNMesh2D::getJacobians()");
		      }
		    }
	  }else{ // they request the complete Jacoby matrix for this bunch of elements
		    //Array<double> J(cellDim*cellDim);
		    SUNDANCE_VERB_HIGH("cellDim == spatialDim()");
		    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
		    {
		      double* J = jBatch.jVals(i);
		      switch(cellDim)
		      {
		        case 2:
				  LID = cellLeafToLIDMapping_[cellLID[i]];
				  // Jacobi for unstructured quad this will not work, but because of linear Jacoby we only have structured quads
		          J[0] = points_[cellsPoints_[LID][1]][0] - points_[cellsPoints_[LID][0]][0];
		          J[1] = 0.0;
		          J[2] = 0.0;
		          J[3] = points_[cellsPoints_[LID][2]][1] - points_[cellsPoints_[LID][0]][1];
			      SUNDANCE_MSG3(verb() , "HNMesh2D::getJacobians X:" << J[0] << " Y:" << J[3] );
		        break;
		        default:
		          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value " "cellDim=" << cellDim << " in HNMesh2D::getJacobians()");
		      }
		    }
	  }
}

void HNMesh2D::getCellDiameters(int cellDim, const Array<int>& cellLID,
                              Array<double>& cellDiameters) const {
	 TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	 SUNDANCE_VERB_HIGH("getCellDiameters()");
	  cellDiameters.resize(cellLID.size());
	  Point pnt(0.0,0.0);
	  int LID;
	  if (cellDim < spatialDim())
	  {
		SUNDANCE_MSG3(verb(),"HNMesh2D::getCellDiameters(), cellDim < spatialDim() ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
	        case 0:
	             cellDiameters[i] = 1.0;
	          break;
	        case 1:  //length of the edge
				  LID = edgeLeafToLIDMapping_[cellLID[i]];
			      pnt = (points_[edgePoints_[LID][1]] - points_[edgePoints_[LID][0]]);
			      cellDiameters[i] = sqrt(pnt * pnt); // the length of the edge
	        break;
	        default:
	        	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "  "cellDim=" << cellDim << " in HNMesh2D::getCellDiameters()");
	      }
	    }
	  }
	  else
	  {
		SUNDANCE_MSG3(verb(),"HNMesh2D::getCellDiameters(), cellDim == spatialDim() ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
	        case 2:
	          LID = cellLeafToLIDMapping_[cellLID[i]];
	          pnt = points_[cellsPoints_[LID][3]] - points_[cellsPoints_[LID][0]];
	          cellDiameters[i] = sqrt(pnt * pnt);
	        break;
	        default:
	          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
	           "cellDim=" << cellDim  << " in HNMesh2D::getCellDiameters()");
	      }
	    }
	  }
}

void HNMesh2D::pushForward(int cellDim, const Array<int>& cellLID,
                         const Array<Point>& refQuadPts,
                         Array<Point>& physQuadPts) const {

	  SUNDANCE_MSG3(verb(),"HNMesh2D::pushForward cellDim:"<<cellDim);
	  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");

	  int nQuad = refQuadPts.size();
	  Point pnt( 0.0 , 0.0 );
	  Point pnt1( 0.0 , 0.0 );

	  if (physQuadPts.size() > 0) physQuadPts.resize(0);
	  physQuadPts.reserve(cellLID.size() * refQuadPts.size());
	  for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	  {
	    switch(cellDim)
	    {
	      case 0: // integrate one point
	         physQuadPts.append(pnt);
	        break;
	      case 1:{ // integrate on one line
			     int LID = edgeLeafToLIDMapping_[cellLID[i]];
		         pnt = points_[edgePoints_[LID][0]];
		         pnt1 = points_[edgePoints_[LID][1]] - points_[edgePoints_[LID][0]];
	             for (int q=0; q<nQuad; q++) {
	               physQuadPts.append(pnt + (pnt1)*refQuadPts[q][0]);
	        	}
	      break;}
	      case 2:{
	             int LID = cellLeafToLIDMapping_[cellLID[i]];
	             pnt = points_[cellsPoints_[LID][0]];
	             // this works only for structured, but we only work on structured quads
	             pnt1 = points_[cellsPoints_[LID][3]] - points_[cellsPoints_[LID][0]];
		         for (int q=0; q<nQuad; q++) {
		          	  physQuadPts.append( pnt + Point( refQuadPts[q][0] * pnt1[0] , refQuadPts[q][1] * pnt1[1] ) );
		         }
	      break;}
	      default:
	    	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value " "in HNMesh2D::getJacobians()");
	    }
	  }
}

int HNMesh2D::ownerProcID(int cellDim, int cellLID) const  {
	 int ID = -1;
	 if (cellDim == 0) ID = vertexLeafToLIDMapping_[cellLID];
     if (cellDim == 1) ID = edgeLeafToLIDMapping_[cellLID];
     if (cellDim == 2) ID = cellLeafToLIDMapping_[cellLID];
     SUNDANCE_MSG3(verb() , " HNMesh2D::ownerProcID ,cellDim:" << cellDim << ", cellLID:"
    		 << cellLID <<" ,ID:" << ID << ", ret:"<< elementOwner_[cellDim][ID] );
	 return elementOwner_[cellDim][ID];
}


int HNMesh2D::numFacets(int cellDim, int cellLID,
                      int facetDim) const  {
	//SUNDANCE_VERB_HIGH("numFacets()");
	if (cellDim==1) { // 1 dimension
         return 2; //one line has 2 points
    }
    else if (cellDim==2) { // 2 dimensions
         return 4; //one quad has 4 edges and 4 points
    }
	return -1;
}

int HNMesh2D::facetLID(int cellDim, int cellLID,
                     int facetDim, int facetIndex,
                     int& facetOrientation) const  {
	// todo: check weather facet orientation is right
	facetOrientation = 1;
	SUNDANCE_MSG3(verb(),"HNMesh2D::facetLID  cellDim:"<<cellDim<<", cellLID:"<<cellLID<<", facetDim:"<<facetDim<< ", facetIndex:" << facetIndex);
	int rnt = -1 , LID=-1 , tmp=-1;
	if (facetDim == 0){ // return the Number/ID of a Vertex
		if (cellDim == 2 ){
		    LID = cellLeafToLIDMapping_[cellLID];
		    rnt = cellsPoints_[LID][facetIndex]; tmp = rnt;
		    rnt = vertexLIDToLeafMapping_[rnt];
	    }
	    else if ((cellDim==1)){
	        LID = edgeLeafToLIDMapping_[cellLID];
	        rnt = edgePoints_[LID][facetIndex]; tmp = rnt;
	        rnt = vertexLIDToLeafMapping_[rnt];
	    }
	} else if (facetDim == 1){
	        LID = cellLeafToLIDMapping_[cellLID];
	        rnt = cellsEdges_[LID][facetIndex]; tmp = rnt;
			rnt = edgeLIDToLeafMapping_[rnt];
	       }
	SUNDANCE_MSG3(verb()," RET = " << rnt << ", LID:" << LID << ", tmp:" <<  tmp);
	return rnt;
}


void HNMesh2D::getFacetLIDs(int cellDim,
                          const Array<int>& cellLID,
                          int facetDim,
                          Array<int>& facetLID,
                          Array<int>& facetSign) const {

	SUNDANCE_MSG3(verb(),"HNMesh2D::getFacetLIDs()  cellDim:"<<cellDim<<"  cellLID.size():"<<cellLID.size()<<"  facetDim:" <<facetDim);
    int LID = 0 , cLID , facetOrientation ;
    int ptr = 0;

    int nf = numFacets(cellDim, cellLID[0], facetDim);
    facetLID.resize(cellLID.size() * nf);
    facetSign.resize(cellLID.size() * nf);
    // At this moment we just use the previous function
	for (unsigned int i = 0 ; i < (unsigned int)cellLID.size() ; i++){
		  cLID = cellLID[i];
	      for (int f=0; f<nf; f++, ptr++) {
	    	  // we use this->facetLID caz facetLID is already used as variable
			  LID = this->facetLID( cellDim, cLID, facetDim, f , facetOrientation);
	          facetLID[ptr] = LID;
	          facetSign[ptr] = facetOrientation;
	      }
	}
	SUNDANCE_MSG3(verb() ,"HNMesh2D::getFacetLIDs()  DONE. ");
}


const int* HNMesh2D::elemZeroFacetView(int cellLID) const {
    int LID = cellLeafToLIDMapping_[cellLID];
    SUNDANCE_MSG3(verb() , "HNMesh2D::elemZeroFacetView ");
	return (const int*)(&cellsPoints_[LID]);
}


int HNMesh2D::numMaxCofacets(int cellDim, int cellLID) const  {
	//SUNDANCE_VERB_HIGH("numMaxCofacets()");
	SUNDANCE_MSG3(verb() , "HNMesh2D::numMaxCofacets():  cellDim:"<<cellDim<<" cellLID:"<<cellLID );
	int rnt = -1;
	if (cellDim==0) { // 1 dimension
		int LID = vertexLeafToLIDMapping_[cellLID];
        int sum = 0;
        for (int i = 0 ; i < 4 ; i++)
        	if ( (pointMaxCoF_[LID][i] >= 0) && ( cellLIDToLeafMapping_[pointMaxCoF_[LID][i]] >= 0) )
        		sum++;
        // return the value, how many cells has this point, on the leaf level
        rnt = sum;
    }
    else if (cellDim==1) { // 2 dimensions
        int LID = edgeLeafToLIDMapping_[cellLID];
		if ((edgeMaxCoF_[LID][0] >= 0) && ( cellLIDToLeafMapping_[edgeMaxCoF_[LID][0]] >= 0) &&
			(edgeMaxCoF_[LID][1] >= 0) && ( cellLIDToLeafMapping_[edgeMaxCoF_[LID][1]] >= 0))
			rnt = 2;
		else
			rnt = 1;
    }
	SUNDANCE_MSG3(verb() ," RET = " << rnt );
	return rnt;
}


int HNMesh2D::maxCofacetLID(int cellDim, int cellLID,
                       int cofacetIndex,
                       int& facetIndex) const  {

	SUNDANCE_MSG3(verb() ,"HNMesh2D::maxCofacetLID() cellDim:"<<cellDim<<" cellLID:"<<cellLID<<" cofacetIndex:"<<cofacetIndex<< " facetIndex:"
			<< facetIndex);
	int rnt =-1;
	if (cellDim==0) { // 1 dimension
		//facetIndex = cofacetIndex;
		int actCoFacetIndex = 0;
    	int LID = vertexLeafToLIDMapping_[cellLID];
		for (int ii = 0 ; ii < 4 ; ii++){
			// take this maxCoFacet only if that exist and is inside
			if ( (pointMaxCoF_[LID][ii] >= 0) && (cellLIDToLeafMapping_[pointMaxCoF_[LID][ii]] >= 0) ){
				if ( actCoFacetIndex < cofacetIndex ){
					actCoFacetIndex++;
				}else{
					facetIndex = ii;
					rnt = pointMaxCoF_[LID][ii];
					break;
				}
			}
		}
    }
    else if (cellDim==1) { // 2 dimensions
    	int orientation = 0;
    	int addFakt = 0;
    	int maxCoFacet = 0;
        int LID = edgeLeafToLIDMapping_[cellLID];
        // this works only in structured mesh case, but we will only use structured quads because of the Jacoby
		if ( edgeVertex_[LID] == 0 ){
			orientation = 0;   addFakt = 2;
		}else{
			orientation = 1;   addFakt = 1;
		}
		SUNDANCE_MSG3(verb() ," HNMesh2D::maxCofacetLID() , orientation:" <<orientation<< " addFakt:" << addFakt );
        // return the index in the vector, which later will be corrected later
		int actCoFacetIndex = 0;
		for (int ii = 0 ; ii < 2 ; ii++){
			// take this maxCoFacet only if that exist and is inside
			if ( (edgeMaxCoF_[LID][ii] >= 0) && (cellLIDToLeafMapping_[edgeMaxCoF_[LID][ii]] >= 0) ){
				if ( actCoFacetIndex < cofacetIndex ){
					actCoFacetIndex++;
				}else{
					facetIndex = 1-ii;
					maxCoFacet = edgeMaxCoF_[LID][ii];
					break;
				}
			}
		}
		SUNDANCE_MSG3(verb() ,"HNMesh2D::maxCofacetLID() , facetIndex:" << facetIndex <<", _edgeMaxCoFacet[0]:" << edgeMaxCoF_[LID][0] <<
				"_edgeMaxCoFacet[1]:" <<  edgeMaxCoF_[LID][1] );
		// calculate the correct facetIndex of the edge in the cell of cofacetIndex
		if ( orientation == 0 ){
			facetIndex = facetIndex + facetIndex*addFakt; // this should be eighter 0 or 3
		} else {
			facetIndex = facetIndex + addFakt; // this should be eighter 1 or 2
		}
		SUNDANCE_MSG3(verb() ," maxCoFacet = " << maxCoFacet );
		rnt = ( maxCoFacet );
    }
	// transform back to leaf indexing
	rnt = cellLIDToLeafMapping_[ rnt ];
	SUNDANCE_MSG3(verb() ," RET = " << rnt << ",  facetIndex:" << facetIndex);
	return rnt;
}

void HNMesh2D::getCofacets(int cellDim, int cellLID,
                 int cofacetDim, Array<int>& cofacetLIDs) const {
	// Nothing to do
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh2D::getCofacets() not implemented");
}


void HNMesh2D::getMaxCofacetLIDs(const Array<int>& cellLIDs,
  MaximalCofacetBatch& cofacets) const {
	// nothing to do here
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh2D::getMaxCofacetLIDs() not implemented");
}


int HNMesh2D::mapGIDToLID(int cellDim, int globalIndex) const  {
	//SUNDANCE_VERB_HIGH("mapGIDToLID()");
	switch (cellDim){
	case 0:{
	         int ID = vertexLeafToGIDMapping_[globalIndex];
		     SUNDANCE_MSG3(verb() , " HNMesh2D::mapGIDToLID 0 , globalIndex:" << globalIndex << " ,ID:" << ID << ", ret:"<< vertexLIDToLeafMapping_[ID]);
	         return vertexLIDToLeafMapping_[ID];
		    break;}
	case 1:{
		     int ID = edgeLeafToGIDMapping_[globalIndex];
		     SUNDANCE_MSG3(verb() , " HNMesh2D::mapGIDToLID 1 , globalIndex:" << globalIndex << " ,ID:" << ID << ", ret:"<< edgeLIDToLeafMapping_[ID]);
		     return edgeLIDToLeafMapping_[ID];
		    break;}
	case 2:{
             int ID = cellLeafToGIDMapping_[globalIndex];
             SUNDANCE_MSG3(verb() , " HNMesh2D::mapGIDToLID 2 , globalIndex:" << globalIndex << " ,ID:" << ID << ", ret:"<< cellLIDToLeafMapping_[ID]);
             return cellLIDToLeafMapping_[ID];
		    break;}
	}
	return -1; //Wall
}


bool HNMesh2D::hasGID(int cellDim, int globalIndex) const {
	//SUNDANCE_VERB_HIGH("hasGID()");
	// we should always have all GIDs
	return true;
}


int HNMesh2D::mapLIDToGID(int cellDim, int localIndex) const  {
	//SUNDANCE_VERB_HIGH("mapLIDToGID()");
	switch (cellDim){
	case 0:{
	         int ID = vertexLeafToLIDMapping_[localIndex];
		     SUNDANCE_MSG3(verb() , " HNMesh2D::mapLIDToGID 0 , localIndex:" << localIndex << " ,ID:" << ID << ", ret:"<< vertexGIDToLeafMapping_[ID]);
	         return vertexGIDToLeafMapping_[ID];
		    break;}
	case 1:{
		     int ID = edgeLeafToLIDMapping_[localIndex];
		     SUNDANCE_MSG3(verb() , " HNMesh2D::mapLIDToGID 1 , localIndex:" << localIndex << " ,ID:" << ID << ", ret:"<< edgeGIDToLeafMapping_[ID]);
		     return edgeGIDToLeafMapping_[ID];
		    break;}
	case 2:{
             int ID = cellLeafToLIDMapping_[localIndex];
             SUNDANCE_MSG3(verb() , " HNMesh2D::mapLIDToGID 2 , localIndex:" << localIndex << " ,ID:" << ID << ", ret:"<< cellGIDToLeafMapping_[ID]);
             return cellGIDToLeafMapping_[ID];
		    break;}
	}
	return -1; //Wall
}


CellType HNMesh2D::cellType(int cellDim) const  {
	 switch(cellDim)
	  {
	    case 0:  return PointCell;
	    case 1:  return LineCell;
	    case 2:  return QuadCell;
	    default:
	      return NullCell; // -Wall
	  }
}


int HNMesh2D::label(int cellDim, int cellLID) const {
   // not used
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh2D::label() not implemented yet");
   return 0;
}


void HNMesh2D::getLabels(int cellDim, const Array<int>& cellLID,
		Array<int>& labels) const {
   // not used
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh2D::getLabels() not implemented yet");
}

Set<int> HNMesh2D::getAllLabelsForDimension(int cellDim) const {
   Set<int>                 rtn;
   // not used
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh2D::getAllLabelsForDimension() not implemented yet");
   return rtn;
}

void HNMesh2D::getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const {
    // not used
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh2D::getLIDsForLabel() not implemented yet");
}

void HNMesh2D::setLabel(int cellDim, int cellLID, int label) {
   // not used
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh2D::setLabel() not implemented yet");
}


void HNMesh2D::assignIntermediateCellGIDs(int cellDim) {
	// The GIDs are assigned
	//SUNDANCE_VERB_HIGH("assignIntermediateCellGIDs()");
}


bool HNMesh2D::hasIntermediateGIDs(int dim) const {
	//SUNDANCE_VERB_HIGH("hasIntermediateGIDs()");
	// the mesh always has intermediate cells
	return true; // true means they have been synchronized ... not used now
}


// =============================== HANGING NODE FUNCTIONS ==========================
bool HNMesh2D::isElementHangingNode(int cellDim , int cellLID) const {
	SUNDANCE_MSG3(verb() ,"HNMesh2D::isElementHangingNode  cellDim:"<<cellDim<<" LID:"<< cellLID);
	if (cellDim==0) { // 1 dimension
    	int LID = vertexLeafToLIDMapping_[cellLID];
		return (isPointHanging_[LID]);
    }
    else if (cellDim==1) { // 2 dimensions
    	int LID = edgeLeafToLIDMapping_[cellLID];
        return (isEdgeHanging_[LID]);
    } else {
	return false; //Wall
    }
	return false; //Wall
}

int HNMesh2D::indexInParent(int maxCellLID) const {
	// this is just a mapping from the HNMesh2D child numbering to the DoFMapHN child numbering in 3D
	int mapMyChilderIndexToStandardDoFMAP[9] = {0,5,6,1,4,7,2,3,8};
	int LID = cellLeafToLIDMapping_[maxCellLID];
	int indexInPar = indexInParent_[LID];
	//map to the trisection standard which is used to the DoFMaps STANDARD NUMBERING(PEANO CURVE) !
	return mapMyChilderIndexToStandardDoFMAP[indexInPar];
}

void HNMesh2D::returnParentFacets(  int childCellLID , int dimFacets ,
		                         Array<int> &facetsLIDs , int &parentCellLIDs ) const {
	int LID = cellLeafToLIDMapping_[childCellLID];
	parentCellLIDs = parentCellLID_[LID];
	facetsLIDs.resize(4);
	SUNDANCE_MSG3(verb() , "HNMesh2D::returnParentFacets  childCellLID:"<<childCellLID<<" dimFacets:"<<dimFacets<<
			"  parentCellLIDs:"<< parentCellLIDs);
	// this is the same for edges and for points
	facetsLIDs[0] = facetLID_tree( 2 , parentCellLIDs ,  dimFacets , 0 );
	facetsLIDs[1] = facetLID_tree( 2 , parentCellLIDs ,  dimFacets , 1 );
	facetsLIDs[2] = facetLID_tree( 2 , parentCellLIDs ,  dimFacets , 2 );
	facetsLIDs[3] = facetLID_tree( 2 , parentCellLIDs ,  dimFacets , 3 );
	// map parent cell ID back to leaf indexing
	parentCellLIDs = cellLIDToLeafMapping_[parentCellLIDs];
}

// only used in determining the parents
int HNMesh2D::facetLID_tree(int cellDim, int cellLID,
                     int facetDim, int facetIndex) const{
    int rnt = -1;
	if (facetDim == 0){ // return the Number/ID of a Vertex
	     rnt = cellsPoints_[cellLID][facetIndex];
	     rnt = vertexLIDToLeafMapping_[rnt];
	     // rnt must be greater than 0
	} else if (facetDim == 1){
    	 rnt = cellsEdges_[cellLID][facetIndex];
	     rnt = edgeLIDToLeafMapping_[rnt];
	     // rnt must be greater than 0
	}
	SUNDANCE_MSG3(verb() , "HNMesh2D::facetLID_tree cellDim:"<<cellDim<<", cellLID:"<<cellLID<<", facetDim:"<<facetDim<<
			", facetIndex:"<<facetIndex<<" RET = " << rnt );
	return rnt;
}

// =========================== MESH CREATION ========================================

/** adds one vertex to the mesh */
void HNMesh2D::addVertex(int vertexLID , int ownerProc , bool isHanging ,
		 double coordx , double coordy , const Array<int> &maxCoF){
  // add only when the LID is new
  // WE ASSUME THAT THE "vertexLID" will come in increasing manner
  if (points_.size() <= vertexLID){
	  TEUCHOS_TEST_FOR_EXCEPTION(vertexLID != nrElem_[0] , std::logic_error ,"HNMesh2D::addVertex " <<
			 " vertexLID:" << vertexLID << " nrElem_[0]:" << nrElem_[0] );
     Point pt(coordx,coordy);
     points_.append( pt );
     pointMaxCoF_.append( maxCoF );
     isPointHanging_.append( isHanging );
     elementOwner_[0].append( ownerProc );
     SUNDANCE_MSG3(verb() , "HNMesh2D::addVertex: " << nrElem_[0] << " P:" << pt << " ,  maxCoF:" << maxCoF);
     SUNDANCE_MSG3(verb() , " ownerProc:" << ownerProc << " , isHanging:" << isHanging);
     nrElem_[0]++;
  }
}

/** adds one edge to the mesh */
void HNMesh2D::addEdge(int edgeLID , int ownerProc , bool isHanging , int edgeVertex ,
		               bool isProcBoundary , bool isMeshBoundary ,
		               const Array<int> &vertexLIDs , const Array<int> &maxCoF){
	  // add only when the edgeLID is new
	  if (edgePoints_.size() <= edgeLID ){
		  TEUCHOS_TEST_FOR_EXCEPTION(edgeLID != nrElem_[1], std::logic_error, "HNMesh2D::addEdge edgeLID != nrElem_[1]");
		 edgePoints_.append( vertexLIDs );
		 edgeVertex_.append( edgeVertex );
		 edgeMaxCoF_.append( maxCoF );
		 edgeIsProcBonudary_.append( isProcBoundary );
		 edgeIsMeshDomainBonudary_.append( isMeshBoundary );
		 isEdgeHanging_.append(isHanging);
	     elementOwner_[1].append( ownerProc );
	     SUNDANCE_MSG3(verb() , "HNMesh2D::addEdge: " << nrElem_[1] << " vertexLIDs:" << vertexLIDs << " ,  maxCoF:" << maxCoF );
	     SUNDANCE_MSG3(verb() , "          ownerProc:" << ownerProc << ", isHanging:" << isHanging << ", edgeVertex:" << edgeVertex <<
	    		 ", isProcBoundary:" << isProcBoundary << ", isMeshBoundary:" << isMeshBoundary);
	     nrElem_[1]++;
	  }
}

/** adds one cell(2D) to the mesh */
void HNMesh2D::addCell(int cellLID , int ownerProc ,
		               int indexInParent , int parentCellLID , int level ,
		               const Array<int> &edgeLIDs , const Array<int> &vertexLIDs){
	  // add only when the edgeLID is new
	  if (cellsPoints_.size() <= cellLID ) {
		  TEUCHOS_TEST_FOR_EXCEPTION(cellLID != nrElem_[2], std::logic_error, "HNMesh2D::cellLID cellLID != nrElem_[2]");
		 cellsPoints_.append( vertexLIDs );
		 cellsEdges_.append( edgeLIDs );
	     indexInParent_.append( indexInParent );
	     parentCellLID_.append( parentCellLID );
	     cellLevel_.append( level );
	     isCellLeaf_.append( true );
	     refineCell_.append( 0 );
	     cellsChildren_.append( tuple(1) );
	     elementOwner_[2].append( ownerProc );
	     SUNDANCE_MSG3(verb() , "HNMesh2D::addCell: " << nrElem_[2] << " vertexLIDs:" << vertexLIDs << " edgeLIDs:" << edgeLIDs);
	     SUNDANCE_MSG3(verb() , "HNMesh2D::addCell p0:" << points_[vertexLIDs[0]] );
	     SUNDANCE_MSG3(verb() , "HNMesh2D::addCell p1:" << points_[vertexLIDs[1]] );
	     SUNDANCE_MSG3(verb() , "HNMesh2D::addCell p2:" << points_[vertexLIDs[2]] );
	     SUNDANCE_MSG3(verb() , "HNMesh2D::addCell p3:" << points_[vertexLIDs[3]] );
		 // calculate if the cell is complete outside the mesh domain
		 isCellOut_.append( !( meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[0]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[1]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[2]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[3]]) ) );
		 SUNDANCE_MSG3(verb() , "HNMesh2D::addCell IN DOMAIN:" <<  isCellOut_[nrElem_[2]] );
	     nrElem_[2]++;
	  }
}

/** creates one regular mesh without refinement. With a different function the
 * refinement can start later , independently from this function. <br>
 * The structure of this mesh also supports unstructured storage of the cells,
 * so we might create unstructured mesh and later refine in the same way */
void HNMesh2D::createMesh(
                      double position_x,
			          double position_y,
			          double offset_x,
			          double offset_y,
			          int resolution_x,
			          int resolution_y,
			          const RefinementClass& refineClass ,
			          const MeshDomainDef& meshDomain
){

	setVerb(0);

	// initialize object fields
	_pos_x = position_x; _pos_y = position_y;
	_ofs_x = offset_x; _ofs_y = offset_y;
	_res_x = resolution_x; _res_y = resolution_y;
	refineClass_ = refineClass;
	meshDomain_ = meshDomain;

	// create coarsest mesh
	createCoarseMesh();

	// loop as long there is no refinement
    bool doRefinement = true;
    while (doRefinement){
    	doRefinement = oneRefinementIteration();
    }

	// calculate global IDs and create leaf Numbering
    //createLeafNumbering();
    createLeafNumbering_sophisticated();

}

void HNMesh2D::createCoarseMesh(){
	// estimate load for parallel case,
	// assign cells to each processors, based on the load and the domain
	// we assign cells to processors, (y,x) (optimal for vertical channel flow)
	// the estimation is done in each processor, globally but then only the local mesh will be build
	int nrCoarseCell = _res_x*_res_y;
	int nrCoarsePoints = (_res_x+1)*(_res_y+1);
	int nrCoarseEdge = (_res_x+1)*_res_y + _res_x*(_res_y+1);
	Array<int> coarseProcessCell( nrCoarseCell , -1 );
	Array<int> coarseCellOwner( nrCoarseCell , -1 );
	Array<int> coarsePointOwner( nrCoarsePoints , -1 );
	Array<int> coarseEdgeOwner( nrCoarseEdge , -1 );
	Array<int> coarseCellLID( nrCoarseCell , -1 );
	Array<int> coarsePointLID( nrCoarsePoints , -1 );
	Array<int> coarseEdgeLID( nrCoarseEdge , -1 );
	Array<int> coarseCellsLoad( nrCoarseCell , 0 );
	int totalLoad = 0;

	SUNDANCE_MSG3(verb() , "HNMesh2D::createMesh nrCoarseCell:" << nrCoarseCell << " nrCoarsePoints:" << nrCoarsePoints
			          << " nrCoarseEdge:" << nrCoarseEdge << " nrProc_:" << nrProc_ << " myRank_:" << myRank_);
	TEUCHOS_TEST_FOR_EXCEPTION( nrCoarseCell < nrProc_ , std::logic_error," HNMesh2D::createMesh nrCoarseCell < nrProc_ ");
	// now always divide as a flow channel , no resolution driven division

    // calculate total load and load per coarse cell
	double h[2];
	h[0] = _ofs_x/(double)_res_x; h[1] = _ofs_y/(double)_res_y;
	Point pos(h[0],h[1]);
	Point res(h[0],h[1]);
	// estimate total estimated load of the mesh
	for (int i=0; i < nrCoarseCell; i++){
		// midpoint of the cell
		pos[0] = _pos_x + (double)(i / _res_y)*h[0] + 0.5*h[0];
		pos[1] = _pos_y + (double)(i % _res_y)*h[1] + 0.5*h[1];
		// todo: take the domain in consideration (take the 4 points) (when cells are turned off)
		coarseCellsLoad[i] = refineClass_.estimateRefinementLevel( pos , res );
		totalLoad += coarseCellsLoad[i];
	}
	// calculate average load per cell
	double loadPerProc = (double)totalLoad / (double)nrProc_;
	int actualProc=0;
	totalLoad = 0;
	Array<int>  ind(2); // vertex is mesh boundary
	// offsets for vertex and edge index
	int vertexOffs[4] = {0 , _res_y + 1 , 1 ,  _res_y + 2};
	int edgeOffs[4] = {_res_y , 0 , 2*_res_y + 1 , _res_y + 1};

	// assign owners to the cells, edges, vertexes , greedy method
	SUNDANCE_MSG3(verb() , "Processor asign, loadPerProc:" << loadPerProc );
	double diff_load = 0.0;
	for (int i=0; i < nrCoarseCell; i++){
		ind[0] = (i / _res_y);
		ind[1] = (i % _res_y);
		int vertexInd = (_res_y+1)*ind[0] + ind[1];
		int edgeInd = (2*_res_y+1)*ind[0] + ind[1];
		//SUNDANCE_MSG3(verb() , "Cell ID:" << i << " vertexInd:" <<vertexInd << " edgeInd:" << edgeInd );
		//SUNDANCE_MSG3(verb() , "Cell, actual index" << ind  );
		// assign ownership for vertexes
		for (int jj = 0 ; jj < 4 ; jj++){
			if (coarsePointOwner[vertexInd+vertexOffs[jj]] < 0){
				coarsePointOwner[vertexInd+vertexOffs[jj]] = actualProc;
				SUNDANCE_MSG3(verb() , "Vertex CPU assign " << vertexInd+vertexOffs[jj] << " ->" << actualProc );
			}
		}
		// assign ownership for edges
		for (int jj = 0 ; jj < 4 ; jj++){
			if (coarseEdgeOwner[edgeInd+edgeOffs[jj]] < 0){
				coarseEdgeOwner[edgeInd+edgeOffs[jj]] = actualProc;
				//SUNDANCE_MSG3(verb() , "Edge CPU assign " << edgeInd+edgeOffs[jj] << " ->" << actualProc );
			}
		}
		// assign ownership for the cell
		coarseCellOwner[i] = actualProc;
		totalLoad += coarseCellsLoad[i];
		SUNDANCE_MSG3(verb() , "Cell CPU assign " << i << " ->" << actualProc <<
				", totalLoad:" << totalLoad << " loadPerProc:" << loadPerProc);
		// the rounding of the load estimator is in favor to late earlier
		if (((double)totalLoad >= (loadPerProc - 1e-8 - diff_load)) && ( actualProc < nrProc_-1 )){
			SUNDANCE_MSG3(verb() , "Increase CPU , totalLoad:" << totalLoad << " loadPerProc:" << loadPerProc );
			// compensate the load difference for the next CPU
			diff_load = totalLoad - loadPerProc;
			actualProc = actualProc + 1;
			totalLoad = 0;
		}
	}

	// next step is to see which cell we will process (we also have to add the remote cells)
	for (int i=0; i < nrCoarseCell; i++){
		ind[0] = (i / _res_y);
		ind[1] = (i % _res_y);
		coarseProcessCell[i] = 1;
	}

	// now go trough all cells which have to be added to this processor
	SUNDANCE_MSG3(verb() , " Process Cells:" << nrCoarseCell );
	for (int i=0; i < nrCoarseCell; i++){
		ind[0] = (i / _res_y);
		ind[1] = (i % _res_y);
		pos[0] = _pos_x + (double)(ind[0])*h[0];
		pos[1] = _pos_y + (double)(ind[1])*h[1];
		SUNDANCE_MSG3(verb() , "PCell ID:" << i << " coarseProcessCell[i]:" <<coarseProcessCell[i] << " pos:" << pos );
		SUNDANCE_MSG3(verb() , "PCell, actual index" << ind  );
		// this condition is so that remote cells will be added
		if (coarseProcessCell[i] > 0) {
			int vertexInd = (_res_y+1)*ind[0] + ind[1];
			int edgeInd = (2*_res_y+1)*ind[0] + ind[1];
			int cellLID = coarseCellLID[i];
			Array<int> vertexLID(4,-1) , vertexMaxCoF(4,-1) , vLID(4,-1);;
			Array<int> edgeLID(4,-1) , edgeVert(2,-1) , edgeMaxCoef(2,-1);
			//SUNDANCE_MSG3(verb() , "Cell ID:" << i << " vertexInd:" <<vertexInd << " edgeInd:" << edgeInd << " cellLID:" << cellLID);
			//SUNDANCE_MSG3(verb() , "Cell, actual index" << ind << " pos:" << pos);
			// assign new cellID if necessary
			if (coarseCellLID[i] < 0 ){
				coarseCellLID[i] = nrElem_[2];
				cellLID = coarseCellLID[i];
			}
			// add all Vertexes , ignoring neighbor cells (maxcofacets)
            for (int jj = 0 ; jj < 4 ; jj++){
            	vLID[jj] = coarsePointLID[vertexInd+vertexOffs[jj]];
            	//SUNDANCE_MSG3(verb() , "Vertex  vertexInd:" << vertexInd);
            	//SUNDANCE_MSG3(verb() , "Vertex  vertexOffs[jj]:" << vertexOffs[jj]);
            	//SUNDANCE_MSG3(verb() , "Vertex  index:" << vertexInd+vertexOffs[jj]);
            	if (coarsePointLID[vertexInd+vertexOffs[jj]] < 0){
            		coarsePointLID[vertexInd+vertexOffs[jj]] = nrElem_[0];
            	}
            	vLID[jj] = coarsePointLID[vertexInd+vertexOffs[jj]];
                // add vertex with -1 maxCOfacets
            	//SUNDANCE_MSG3(verb() , "Vertex  X:" << ((double)offs_Points_x_[jj])*h[0] << "  Y:" << pos[1] + ((double)offs_Points_y_[jj])*h[1]);
            	//SUNDANCE_MSG3(verb() , "Vertex  vLID[jj]:" << vLID[jj] << "  , coarsePointOwner[vertexInd+vertexOffs[jj]]" << coarsePointOwner[vertexInd+vertexOffs[jj]] );
            	addVertex( vLID[jj] , coarsePointOwner[vertexInd+vertexOffs[jj]] , false ,
            			   pos[0] + ((double)offs_Points_x_[jj])*h[0] , pos[1] + ((double)offs_Points_y_[jj])*h[1] ,
            	           vertexMaxCoF );
            }
			// add all Edges , ignoring neighbor cells (maxcofacets)
            for (int jj = 0 ; jj < 4 ; jj++){
            	edgeLID[jj] = coarseEdgeLID[edgeInd+edgeOffs[jj]];
            	if (coarseEdgeLID[edgeInd+edgeOffs[jj]] < 0){
            		coarseEdgeLID[edgeInd+edgeOffs[jj]] = nrElem_[1];
            	}
            	edgeLID[jj] = coarseEdgeLID[edgeInd+edgeOffs[jj]];
            	edgeVert[0] = vLID[edge_Points_localIndex[jj][0]];
            	edgeVert[1] = vLID[edge_Points_localIndex[jj][1]];
                // add edge with -1 maxCOfacets
            	addEdge( edgeLID[jj] , coarseEdgeOwner[edgeInd+edgeOffs[jj]] , false , (jj==0 || jj==3) ? 0 : 1 ,
            		            false , false , // these two flags later need to be updated
            			        edgeVert , edgeMaxCoef );
            }
			// add the Cell
            addCell( cellLID , coarseCellOwner[i] , 0 , cellLID , 0 , edgeLID , vLID);
		}
	} // --- end from for loop

	// next is maxCoFacet and boundary info update for vertexes and edges
	SUNDANCE_MSG3(verb() , " Process maxCofacets:" << nrCoarseCell );
	for (int i=0; i < nrCoarseCell; i++){
		// this condition is so that remote cells will be added
		if (coarseProcessCell[i] > 0){
			Array<int> vLID(4,-1) , eLID(4,-1) , maxVertexCoF(4,-1);
			ind[0] = (i / _res_y);
			ind[1] = (i % _res_y);
			int vertexInd = (_res_y+1)*ind[0] + ind[1];
			int edgeInd = (2*_res_y+1)*ind[0] + ind[1];
			int cellLID = coarseCellLID[i];
			SUNDANCE_MSG3(verb() , "MCell cellLID:" << cellLID << " coarseProcessCell[i]:" <<coarseProcessCell[i] );
			//SUNDANCE_MSG3(verb() , "MCell ID:" << i << " vertexInd:" <<vertexInd << " edgeInd:" << edgeInd );
			//SUNDANCE_MSG3(verb() , "MCell, actual index" << ind  );
			// vertex maxCoFac
            for (int jj = 0 ; jj < 4 ; jj++){
            	vLID[jj] = coarsePointLID[vertexInd+vertexOffs[jj]];
            	// if all elements are added then this is OK
            	pointMaxCoF_[vLID[jj]][jj] = cellLID;
            	SUNDANCE_MSG3(verb() , "Vertex MaxCoFacet vLID[jj]:" << vLID[jj] << " jj:" << jj << " cellLID:" << cellLID );
            }
            // edge maxCoFac
            for (int jj = 0 ; jj < 4 ; jj++){
            	eLID[jj] = coarseEdgeLID[edgeInd+edgeOffs[jj]];
            	// if all elements are added then this is OK
            	edgeMaxCoF_[eLID[jj]][(3-jj)/2] = cellLID;
            	SUNDANCE_MSG3(verb() , "Edge MaxCoFacet eLID[jj]:" << eLID[jj] << " (3-jj)/2:" << (3-jj)/2 << " cellLID:" << cellLID );
            	// update boundary info of the edge
            	int orientation = ( (jj==0) || (jj==3)) ? 1 : 0;

            	if (orientation == 0){ // vertical edge
            	   // mesh boundary
                   if ( (ind[0] == 0) && (jj==1))   edgeIsMeshDomainBonudary_[eLID[jj]] = true;
                   if ( (ind[0] == _res_x - 1 ) && (jj==2))   edgeIsMeshDomainBonudary_[eLID[jj]] = true;
                   // process boundary
                   if ( edgeIsMeshDomainBonudary_[eLID[jj]] != true){
                	   // now here is dummy implementation
                       edgeIsProcBonudary_[eLID[jj]] = false; //(bool)(edgeOwner != edgeOwnerN);
                   }
            	}else{ // horizontal edge
            	   // mesh boundary
                   if ( (ind[1] == 0) && (jj==0))  { edgeIsMeshDomainBonudary_[eLID[jj]] = true; }
                   if ( (ind[1] == _res_y - 1 ) && (jj==3))  { edgeIsMeshDomainBonudary_[eLID[jj]] = true; }
                   // process boundary
                   if ( edgeIsMeshDomainBonudary_[eLID[jj]] != true){
                	   // now there is only the dummy implementation
                       edgeIsProcBonudary_[eLID[jj]] = false; //(bool)(edgeOwner != edgeOwnerN);
                   }
                   SUNDANCE_MSG3(verb() , "Neighb. H  Cell DONE ");
            	}
            }
		}
	}
	// basic rectangular mesh is build
}

// -----------
bool HNMesh2D::oneRefinementIteration(){

	Array<int> ind(2);
    int nrActualCell = nrElem_[2];
    bool rtn = false;
    SUNDANCE_MSG3(verb() , " HNMesh2D::oneRefinementIteration, start one refinement iteration cycle ");
    // we iterate only over the existing cells (not the ones which will be later created)
	for (int i=0 ; i < nrActualCell ; i++){

		ind[0] = (i / _res_y);
		ind[1] = (i % _res_y);

		// cell is owned by the current processor, and is leaf and is inside the mesh domain
	    SUNDANCE_MSG3(verb() , " Test cell " << i << ", elementOwner_[2][i]:" << elementOwner_[2][i] <<
	    		               ", isCellLeaf_[i]:" << isCellLeaf_[i] << ", out:" << (!isCellOut_[i]));

		if ( (isCellLeaf_[i] == true) )
			//(elementOwner_[2][i] == myRank_) && (isCellLeaf_[i] == true) )
			//( (elementOwner_[2][i] == myRank_) && (isCellLeaf_[i] == true) && (!isCellOut_[i]))
			//  mark neighbors for refinements if because of the hanging nodes a refinement is not possible, regardless of the IN or OUT of Mesh domain
		{
			// check if refinement is needed and possible, none of the vertexes can be hanging
			Array<int>& cellsEdges = cellsEdges_[i];
            bool doRefined = true;
            bool refFunction = false;
            for (int jj = 0 ; jj < cellsEdges.size() ; jj++){
            	SUNDANCE_MSG3(verb() , " eLID: " << cellsEdges[jj] << ", isHanging:" << isEdgeHanging_[cellsEdges[jj]]);
            	doRefined = ( doRefined && ((isEdgeHanging_[cellsEdges[jj]]) == false) );
            }
            // --- if refinement is possible
            SUNDANCE_MSG3(verb() , " is possible to refine cell nr: " << i << ", doRefined:" << doRefined);
            // call refinement function
			Array<int>& cellsVertex = cellsPoints_[i];
			Point h = points_[cellsVertex[3]] - points_[cellsVertex[0]];
			Point p2 = points_[cellsVertex[0]] + 0.5*h;
			refFunction = ((refineCell_[i] == 1) || refineClass_.refine( cellLevel_[i] , p2 , h ));

            // decide if we refine this cell
            SUNDANCE_MSG3( verb() , " execute refinement on cell nr: " << i << ", refFunction:" << refFunction);
            if (doRefined && refFunction)
            {
            	// -- call the function to refine the actual cell ---
            	refineCell(i);
            	rtn = true;
            	refineCell_[i] = 0;
            }
            else
            {
            	// Cell can not be refined
                // we mark neighbor cells, based on "refFunction"
            	if (refFunction){
            		//SUNDANCE_MSG3( verb() , " HNMesh2D::oneRefinementIteration mark neighbor cells ");
                    for (int jj = 0 ; jj < cellsEdges.size() ; jj++){
                    	if (isEdgeHanging_[cellsEdges[jj]]){
                    		//SUNDANCE_MSG3( verb() , " HNMesh2D::oneRefinementIteration isEdgeHanging_[cellsEdges[jj]] == true , " << jj);
                    		// get the parent cell
                    		int cLevel = cellLevel_[i];
                    		int parentID = parentCellLID_[i];
                    		int parentCEdge = cellsEdges_[parentID][jj];
                    		int refCell = -1;
                    		// look for maxCoFacets
                    		if ( (jj == 0) || ( jj == 1) )
                    			refCell = edgeMaxCoF_[parentCEdge][0];
                    		else
                    			refCell = edgeMaxCoF_[parentCEdge][1];
                    		// refCell should be refined and mark for refinement
                    		//SUNDANCE_MSG3( verb() , " cLevel:" << cLevel << " parentID:" << parentID << " parentCEdge:" << parentCEdge
                    		//		<< " refCell:" << refCell );
                    		//SUNDANCE_MSG3( verb() ," edgeMaxCoF_[parentCEdge]:" << edgeMaxCoF_[parentCEdge]);
                            if ( (refCell >= 0) && (cellLevel_[refCell] < cLevel) && (isCellLeaf_[refCell])){

                            	refineCell_[refCell] = 1;
                            	rtn = true;
                            	//SUNDANCE_MSG3( verb() , " HNMesh2D::oneRefinementIteration refineCell_[refCell] = 1 , " << refCell
                            	//		<< ", cellLevel_[refCell]:" << cellLevel_[refCell]);
                            }
                    	}
                    }
            	}
            }
		}
	}

	//  -> communicate with neighbor processor
	//     not needed if the mesh exist globally

	SUNDANCE_MSG3(verb() , " HNMesh2D::oneRefinementIteration DONE with one refinement iteration");
	// return if there was refinement or attempt to refine
	return rtn;
}

// -------- refine cell cellID (assuming all conditions are satisfied)--
void HNMesh2D::refineCell(int cellLID){
    // data initialization
	Array<int>& cellsEdges = cellsEdges_[cellLID];
	Array<int>& cellsVertex = cellsPoints_[cellLID];
	int cellOwner = elementOwner_[2][cellLID];
	Point p = points_[cellsVertex[0]];
	Point h = points_[cellsVertex[3]] - points_[cellsVertex[0]];
    // the X and the Y coordinates of the newly
	double vertex_X[16] = { 0.0 , 0.0 , 0.0 , 0.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 ,
			                2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 1.0 , 1.0 , 1.0 , 1.0 };
	double vertex_Y[16] = { 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 ,
			                0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 };
    // if we have pre stored  hanging edge , we might just load them and these are the indexes (in 3X3 refinement)
	int hanginEdgeI[4][3] = { {3,10,17} , {0,1,2} , {21,22,23} , {6,13,20} };
	// existing vertexes and their index in the 3X3 cell refinement
	int hanginVertexI[4][2] = { {4,8} , {1,2} , {13,14} , {7,11} };
	// orientation of the newly (or existing) edges
	int edgeVertex[24] = { 1,1,1 , 0,0,0,0 , 1,1,1  ,  0,0,0,0  ,  1,1,1  ,  0,0,0,0  ,  1,1,1};
	// the index in the vertex for one edge where the existing elements are stored
	int VertexIndex[4] = { 0 , 1 , 1 , 0 };
	// for each edge (where the existing vertexes are stored)
	int VertexI[4] = { 0 , 0 , 1 , 2 };

    // vertex info , owner of the new vertexes
	int vertexOwner[16] = { elementOwner_[0][cellsVertex[0]] , elementOwner_[1][cellsEdges[1]] ,
			                elementOwner_[1][cellsEdges[1]]  , elementOwner_[0][cellsVertex[2]] ,
			                elementOwner_[1][cellsEdges[0]]  , cellOwner , cellOwner , elementOwner_[1][cellsEdges[3]] ,
			                elementOwner_[1][cellsEdges[0]]  , cellOwner , cellOwner , elementOwner_[1][cellsEdges[3]] ,
			                elementOwner_[0][cellsVertex[1]] , elementOwner_[1][cellsEdges[2]] ,
    			            elementOwner_[1][cellsEdges[2]]  , elementOwner_[0][cellsVertex[3]]  };
	// if vertex is hanging
	bool vertexIsHanging[16] = { false , !edgeIsMeshDomainBonudary_[cellsEdges[1]] , !edgeIsMeshDomainBonudary_[cellsEdges[1]], false ,
			                     !edgeIsMeshDomainBonudary_[cellsEdges[0]] , false , false , !edgeIsMeshDomainBonudary_[cellsEdges[3]] ,
			                     !edgeIsMeshDomainBonudary_[cellsEdges[0]] , false , false , !edgeIsMeshDomainBonudary_[cellsEdges[3]] ,
	                             false , !edgeIsMeshDomainBonudary_[cellsEdges[2]] , !edgeIsMeshDomainBonudary_[cellsEdges[2]], false };
	// edge info
	int edgeOwner[24] = { elementOwner_[1][cellsEdges[1]] , elementOwner_[1][cellsEdges[1]] , elementOwner_[1][cellsEdges[1]] ,
			              elementOwner_[1][cellsEdges[0]] , cellOwner , cellOwner , elementOwner_[1][cellsEdges[3]] ,
			              cellOwner , cellOwner , cellOwner ,
			              elementOwner_[1][cellsEdges[0]] , cellOwner , cellOwner , elementOwner_[1][cellsEdges[3]] ,
			              cellOwner , cellOwner , cellOwner ,
			              elementOwner_[1][cellsEdges[0]] , cellOwner , cellOwner , elementOwner_[1][cellsEdges[3]] ,
			              elementOwner_[1][cellsEdges[2]] , elementOwner_[1][cellsEdges[2]] , elementOwner_[1][cellsEdges[2]]  };
	// edge hanging
	bool edgeIsHanging[24] = { !edgeIsMeshDomainBonudary_[cellsEdges[1]] , !edgeIsMeshDomainBonudary_[cellsEdges[1]] , !edgeIsMeshDomainBonudary_[cellsEdges[1]] ,
			                   !edgeIsMeshDomainBonudary_[cellsEdges[0]] , false , false , !edgeIsMeshDomainBonudary_[cellsEdges[3]] ,
	                           false , false , false ,
	                           !edgeIsMeshDomainBonudary_[cellsEdges[0]] , false , false , !edgeIsMeshDomainBonudary_[cellsEdges[3]] ,
	                           false , false , false ,
	                           !edgeIsMeshDomainBonudary_[cellsEdges[0]] , false , false , !edgeIsMeshDomainBonudary_[cellsEdges[3]] ,
	                           !edgeIsMeshDomainBonudary_[cellsEdges[2]] , !edgeIsMeshDomainBonudary_[cellsEdges[2]] , !edgeIsMeshDomainBonudary_[cellsEdges[2]] };
	// vertexes which form the new (or existing) edges
	int edgeVertexes[24][2] = { {0,1} , {1,2} , {2,3} , {0,4} , {1,5} , {2,6} , {3,7} , {4,5} , {5,6} , {6,7} ,
			                  {4,8} , {5,9} , {6,10} , {7,11} , {8,9} , {9,10} , {10,11} ,
			                  {8,12} , {9,13} , {10,14} , {11,15} , {12,13} , {13,14} , {14,15}};
	// if edge is processor boundary
	bool edgePBnd[24] = { edgeIsProcBonudary_[cellsEdges[1]] , edgeIsProcBonudary_[cellsEdges[1]] , edgeIsProcBonudary_[cellsEdges[1]] ,
			              edgeIsProcBonudary_[cellsEdges[0]] , false , false , edgeIsProcBonudary_[cellsEdges[3]] ,
	                      false , false , false ,
	                      edgeIsProcBonudary_[cellsEdges[0]] , false , false , edgeIsProcBonudary_[cellsEdges[3]] ,
	                      false , false , false ,
	                      edgeIsProcBonudary_[cellsEdges[0]] , false , false , edgeIsProcBonudary_[cellsEdges[3]] ,
	                      edgeIsProcBonudary_[cellsEdges[2]] , edgeIsProcBonudary_[cellsEdges[2]] , edgeIsProcBonudary_[cellsEdges[2]] };
	// if edge is mesh boundary
	bool edgeMBnd[24] = { edgeIsMeshDomainBonudary_[cellsEdges[1]] , edgeIsMeshDomainBonudary_[cellsEdges[1]] , edgeIsMeshDomainBonudary_[cellsEdges[1]] ,
			              edgeIsMeshDomainBonudary_[cellsEdges[0]] , false , false , edgeIsMeshDomainBonudary_[cellsEdges[3]] ,
						  false , false , false ,
						  edgeIsMeshDomainBonudary_[cellsEdges[0]] , false , false , edgeIsMeshDomainBonudary_[cellsEdges[3]] ,
                          false , false , false ,
                          edgeIsMeshDomainBonudary_[cellsEdges[0]] , false , false , edgeIsMeshDomainBonudary_[cellsEdges[3]] ,
                          edgeIsMeshDomainBonudary_[cellsEdges[2]] , edgeIsMeshDomainBonudary_[cellsEdges[2]] , edgeIsMeshDomainBonudary_[cellsEdges[2]] };
	// vertex index in the 3X3 context which is needed for the cell creation
	int refinedCellsVertexes[9][4] = { {0,4,1,5} , {1,5,2,6} , {2,6,3,7} ,
			                  {4,8,5,9} , {5,9,6,10} , {6,10,7,11} ,
			                  {8,12,9,13} , {9,13,10,14} , {10,14,11,15} };
	// edge context in the 3X3 context, for cell creation
	int refinedCellsEdges[9][4] =   { {3,0,7,4} , {4,1,8,5} , {5,2,9,6} ,
			                  {10,7,14,11} , {11,8,15,12} , {12,9,16,13} ,
			                  {17,14,21,18} , {18,15,22,19} , {19,16,23,20} };
	// the index and offset to store the hanging nodes
	int hVertexIndex[16]       = { -1 , 0 , 1 , -1 , 0 , -1 , -1 , 0 , 1 , -1 , -1 , 1 , -1 , 0 , 1 , -1 };
	int hEdgeIndex[24]         = { 0 , 1 , 2 , 0 , -1 , -1 , 0 , -1 , -1 , -1 ,  1 , -1 , -1 , 1 , -1 , -1 , -1 , 2 , -1 , -1 , 2 , 0 , 1 , 2 };
	int hVertexVertexIndex[16] = { -1 , 1 , 1 , -1 , 0 , -1 , -1 , 3 , 0 , -1 , -1 , 3 , -1 , 2 , 2 , -1 };
	int hEdgeVertexIndex[24]   = { 1 , 1 , 1 , 0 , -1 , -1 , 3 , -1 , -1 , -1 ,  0 , -1 , -1 , 3 , -1 , -1 , -1 , 0 , -1 , -1 , 3 , 2 , 2 , 2 };

	// cell index in the 3X3 for vertex maxCoFacets
	int vertexMaxCoFacet[16][4] = { {0,-1,-1,-1} , {1,-1,0,-1} , {2 ,-1,1,-1} , {-1,-1,2,-1} ,
			                        {3,0,-1,-1} , {4,1,3,0}    , {5 ,2,4,1} , {-1,-1,5,2}    ,
			                        {6,3,-1,-1} , {7,4,6,3}    , {8 ,5,7,4} , {-1,-1,8,5}    ,
			                        {-1,6,-1,-1} , {-1,7,-1,6} , {-1,8,-1,7} , {-1,-1,-1,8}  };
	// there will be 9 newly created cells, these should be maxCoFactes in the vertexes and the edges
	int edgeMaxCoFacet[24][2] = { {-1,0} , {-1,1} , {-1,2} ,
			                      {-1,0} , {0,1} , {1,2} , {2,-1} ,
			                       {0,3} , {1,4} , {2,5} ,
			                      {-1,3} , {3,4} , {4,5} , {5,-1} ,
			                       {3,6} , {4,7} , {5,8} ,
			                      {-1,6} , {6,7} , {7,8} , {8,-1} ,
			                      {6,-1} , {7,-1} , {8,-1} };
	// we might need to copy the MaxCoFacet from the parent edge
	int indexEdgeParentMaxCoFacet[24] = { 1 , 1, 1 , 0 , -1 , -1 , 3 , -1 , -1 , - 1, 0 , -1 , -1 , 3 , -1 , -1 , -1 , 0 ,-1,-1,3, 2,2,2 };
	// in which direction to look
	int offsEdgeParentMaxCoFacet[24] = { 0 , 0, 0 , 0 , -1 , -1 , 1 , -1 , -1 , - 1, 0 , -1 , -1 , 1 , -1 , -1 , -1 , 0 ,-1,-1,1, 1,1,1 };

    // hanging elements information
	Array< Array<int> > hangingInfo(4);
	hangingInfo[0].resize(0); hangingInfo[1].resize(0); hangingInfo[2].resize(0); hangingInfo[3].resize(0);

    // the created vertex, edge and cell LIDs
    Array<int> vertexLIDs(16,-1);
    Array<int> edgeLIDs(24,-1);
    Array<int> cellLIDs(9,-1);

    SUNDANCE_MSG3(verb() , " ========== HNMesh2D::refineCell, cellLID:" << cellLID << " ==========");
    //SUNDANCE_MSG3( verb() , " cellsVertex:" << cellsVertex );
    // the 4 "corner" points are already created (according to refinement numbering!)
    vertexLIDs[0] = cellsVertex[0]; vertexLIDs[3] = cellsVertex[2];
    vertexLIDs[12] = cellsVertex[1]; vertexLIDs[15] = cellsVertex[3];

    // look for existing elements
    for (int v = 0 ; v < 4 ; v++){
    	// look for already stored elements
    	//SUNDANCE_MSG3( verb() , " Test vertex VertexIndex[v]:" << VertexIndex[v] << ", cellsVertex[VertexI[v]]:" << cellsVertex[VertexI[v]] );
    	//SUNDANCE_MSG3( verb() , " hangElmStore_[VertexIndex[v]].size()" << hangElmStore_[VertexIndex[v]].size() );
    	if (hangElmStore_[VertexIndex[v]].containsKey(cellsVertex[VertexI[v]])){
           const Array<int>& hnInfo = hangElmStore_[VertexIndex[v]].get(cellsVertex[VertexI[v]]);
           // this is the convention how we store the hanging info
           // these elements are not hanging after this
           //SUNDANCE_MSG3( verb() , " Found key: " << cellsVertex[VertexI[v]] << ", size:" << hnInfo.size() << " , array:" << hnInfo );
           vertexLIDs[hanginVertexI[v][0]] = hnInfo[0];  isPointHanging_[hnInfo[0]] = false;
           vertexLIDs[hanginVertexI[v][1]] = hnInfo[1];  isPointHanging_[hnInfo[1]] = false;
           edgeLIDs[hanginEdgeI[v][0]] = hnInfo[2];
           isEdgeHanging_[hnInfo[2]] = false;
           edgeLIDs[hanginEdgeI[v][1]] = hnInfo[3];
           isEdgeHanging_[hnInfo[3]] = false;
           edgeLIDs[hanginEdgeI[v][2]] = hnInfo[4];  isEdgeHanging_[hnInfo[4]] = false;
           // remove this from the tmp storage
           hangElmStore_[VertexIndex[v]].remove(cellsVertex[VertexI[v]]);
    	}
    }
    SUNDANCE_MSG3(verb() , " refineCell bef. vertexLIDs:" << vertexLIDs );
    SUNDANCE_MSG3(verb() , " refineCell bef. edgeLIDs:" << edgeLIDs );

    // create all vertexes, those which are not created
    for (int v = 0 ; v < 16 ; v++){
    	// create vertex if necessary
    	if (vertexLIDs[v] < 0){
    		// this means vertex is not created we create now
    		Array<int> maxCoFacets(4,-1);
    		// add the hanging node
    		int vLID = nrElem_[0];
    	    addVertex( nrElem_[0] , vertexOwner[v] , vertexIsHanging[v] ,
    	  		         p[0] + vertex_X[v]*h[0] , p[1] + vertex_Y[v]*h[1] , maxCoFacets );
    	    vertexLIDs[v] = vLID;
	    	// if vertex is hanging add to hangElmStore_
    	    if (vertexIsHanging[v]){
    	    	Array<int>& hinf = hangingInfo[hVertexVertexIndex[v]];
    	    	if (hinf.size() < 1) hinf.resize(5,-1);
    	    	hinf[hVertexIndex[v]] = vLID;
    	    	// for hanging Vertexes we might add the parent neighbor cell twice as maxCoFacet
    	    	if ( (v == 1) || ( v==2 )) {  pointMaxCoF_[vertexLIDs[v]][1] = edgeMaxCoF_[cellsEdges[1]][0];
    	    	                              pointMaxCoF_[vertexLIDs[v]][3] = edgeMaxCoF_[cellsEdges[1]][0]; }
    	    	if ( (v == 4) || ( v==8 )) {  pointMaxCoF_[vertexLIDs[v]][2] = edgeMaxCoF_[cellsEdges[0]][0];
    	    	                              pointMaxCoF_[vertexLIDs[v]][3] = edgeMaxCoF_[cellsEdges[0]][0]; }
    	    	if ( (v == 7) || ( v==11 )) {  pointMaxCoF_[vertexLIDs[v]][0] = edgeMaxCoF_[cellsEdges[3]][1];
    	    	                               pointMaxCoF_[vertexLIDs[v]][1] = edgeMaxCoF_[cellsEdges[3]][1]; }
    	    	if ( (v == 13) || ( v==14 )) {  pointMaxCoF_[vertexLIDs[v]][0] = edgeMaxCoF_[cellsEdges[2]][1];
    	    	                                pointMaxCoF_[vertexLIDs[v]][2] = edgeMaxCoF_[cellsEdges[2]][1]; }
    	    	//SUNDANCE_MSG3(verb() , " ParentEdgeMaxCofacet 0 v:" << v << "," << edgeMaxCoF_[cellsEdges[0]]);
    	    	//SUNDANCE_MSG3(verb() , " ParentEdgeMaxCofacet 1 v:" << v << "," << edgeMaxCoF_[cellsEdges[1]]);
    	    	//SUNDANCE_MSG3(verb() , " ParentEdgeMaxCofacet 2 v:" << v << "," << edgeMaxCoF_[cellsEdges[2]]);
    	    	//SUNDANCE_MSG3(verb() , " ParentEdgeMaxCofacet 3 v:" << v << "," << edgeMaxCoF_[cellsEdges[3]]);
    	    	//SUNDANCE_MSG3(verb() , " vertexLIDs[v]:" << vertexLIDs[v] << ", Maxcof:" << pointMaxCoF_[vertexLIDs[v]] );
    	    	//pointMaxCoF_[vertexLIDs[v]][0] = -1;
    	    }
    	}
    	// update maxCoFacet info (no additional information is needed)
    	for (int hh = 0 ; hh < 4 ; hh++){
    		SUNDANCE_MSG3(verb() , " vertexMaxCoFacet[v][hh]:" << vertexMaxCoFacet[v][hh] << " , vertexLIDs[v]:"<< vertexLIDs[v] );
    		if (vertexMaxCoFacet[v][hh] >= 0)
    			pointMaxCoF_[vertexLIDs[v]][hh] = nrElem_[2] + vertexMaxCoFacet[v][hh];
    	}
    	SUNDANCE_MSG3(verb() , " MaxCoFac point:" << vertexLIDs[v] << " , MaxCoFac:"<< pointMaxCoF_[vertexLIDs[v]] );
    }
    SUNDANCE_MSG3(verb() , " refineCell after Vertex creation vertexLIDs:" << vertexLIDs );

    // create all edges , those which are not created
    for (int e = 0 ; e < 24 ; e++){
    	// create edge if necessary
    	if (edgeLIDs[e] < 0){
    		Array<int> maxCoFacets(2,-1);
    		Array<int> edgeVertexLIDs(2,-1);
    		edgeVertexLIDs[0] = vertexLIDs[edgeVertexes[e][0]];
    		edgeVertexLIDs[1] = vertexLIDs[edgeVertexes[e][1]];
    		int eLID = nrElem_[1];
    		addEdge( nrElem_[1] , edgeOwner[e] , edgeIsHanging[e] , edgeVertex[e] ,
    				edgePBnd[e] , edgeMBnd[e] ,  edgeVertexLIDs , maxCoFacets );
    		edgeLIDs[e] = eLID;
			//if edge is hanging add to hangElmStore_
    		if (edgeIsHanging[e]){
    	    	Array<int>& hinf = hangingInfo[hEdgeVertexIndex[e]];
    	    	if (hinf.size() < 1) hinf.resize(5,-1);
    	    	SUNDANCE_MSG3(verb() , " HNI e:" << e << ", hEdgeIndex[e]:" << hEdgeIndex[e] );
    	    	hinf[2+hEdgeIndex[e]] = eLID;
    	    	// assign MaxCoefs only when is not Mesh boundary
    	    	if ( edgeMBnd[e] == false){
        	    	// add maxCoFs from parent edge, so that hanging edges will have also 2 maxCoFs (for boundary)
    	    		edgeMaxCoF_[edgeLIDs[e]][offsEdgeParentMaxCoFacet[e]] =
    	    				edgeMaxCoF_[cellsEdges[indexEdgeParentMaxCoFacet[e]]][offsEdgeParentMaxCoFacet[e]];
    	    	}
    		}
    	}
    	// update maxCoFacet info for this edge
    	for (int hh = 0 ; hh < 2 ; hh++){
    		//SUNDANCE_MSG3(verb() , " edgeMaxCoFacet[e][hh]:" << edgeMaxCoFacet[e][hh] << " , edgeLIDs[e]:"<< edgeLIDs[e] );
    		if (edgeMaxCoFacet[e][hh] >= 0)
    			edgeMaxCoF_[edgeLIDs[e]][hh] = nrElem_[2] + edgeMaxCoFacet[e][hh];
    	}
    	//SUNDANCE_MSG3(verb() , " MaxCoFac edge:" << edgeLIDs[e] << " , MaxCoFac:"<< edgeMaxCoF_[edgeLIDs[e]] );
    }
    SUNDANCE_MSG3(verb() , " refineCell after Edges creation edgeLIDs:" << edgeLIDs );

    // add all 9 cells
    for (int c = 0 ; c < 9 ; c++){
    	Array<int> eLIDs(4,-1);
    	Array<int> vLIDs(4,-1);
    	vLIDs[0] = vertexLIDs[refinedCellsVertexes[c][0]];  vLIDs[1] = vertexLIDs[refinedCellsVertexes[c][1]];
    	vLIDs[2] = vertexLIDs[refinedCellsVertexes[c][2]];  vLIDs[3] = vertexLIDs[refinedCellsVertexes[c][3]];
    	eLIDs[0] = edgeLIDs[refinedCellsEdges[c][0]]; eLIDs[1] = edgeLIDs[refinedCellsEdges[c][1]];
    	eLIDs[2] = edgeLIDs[refinedCellsEdges[c][2]]; eLIDs[3] = edgeLIDs[refinedCellsEdges[c][3]];
    	// add cell
    	cellLIDs[c] = nrElem_[2];
        addCell( nrElem_[2] , cellOwner , c , cellLID , cellLevel_[cellLID] + 1 , eLIDs , vLIDs );
    }

    SUNDANCE_MSG3(verb() , " ---- Adding hanging node information ----- " );

    if ( (hangElmStore_[0].containsKey(cellsPoints_[cellLID][0]) == false) && ( hangingInfo[0].size() > 0 ))
    {
    	hangElmStore_[0].put(cellsPoints_[cellLID][0],hangingInfo[0]);
    	//SUNDANCE_MSG3( verb() , " Store array 1 : " << hangingInfo[0] );
    }
    if ( (hangElmStore_[1].containsKey(cellsPoints_[cellLID][0]) == false) && ( hangingInfo[1].size() > 0 ))
    {
    	hangElmStore_[1].put(cellsPoints_[cellLID][0],hangingInfo[1]);
    	//SUNDANCE_MSG3( verb() , " Store array 2 : " << hangingInfo[1] );
    }
    if ( (hangElmStore_[1].containsKey(cellsPoints_[cellLID][1]) == false) && ( hangingInfo[2].size() > 0 ))
    {
    	hangElmStore_[1].put(cellsPoints_[cellLID][1],hangingInfo[2]);
    	//SUNDANCE_MSG3( verb() , " Store array 3 : " << hangingInfo[2] );
    }
    if ( (hangElmStore_[0].containsKey(cellsPoints_[cellLID][2]) == false) && ( hangingInfo[3].size() > 0 ))
    {
    	hangElmStore_[0].put(cellsPoints_[cellLID][2],hangingInfo[3]);
    	//SUNDANCE_MSG3( verb() , " Store array 4 : " << hangingInfo[3] );
    }

    SUNDANCE_MSG3(verb() , " ---- Refinement cell DONE , mark cell as not leaf and store children LID ----");
    isCellLeaf_[cellLID] = false;
    // store the children of the parent cell
    cellsChildren_[cellLID] = cellLIDs;
}

// -----------
void HNMesh2D::createLeafNumbering(){

	// set all leaf numbers to -1

	// - iterate trough the mesh and in the leaf cells , distribute leaf numbering
	// , detect if one cell is leaf ()
	// , have a tree similar tree traversal ... todo: later

	SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering nrPoint:" << nrElem_[0] << " , nrEdge:" << nrElem_[1] << ", nrCell:" << nrElem_[2]);
	// we resize the leafID - > global
	vertexGIDToLeafMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexGIDToLeafMapping_[dd] = -1;
	vertexLeafToGIDMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLeafToGIDMapping_[dd] = -1;

	edgeGIDToLeafMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeGIDToLeafMapping_[dd] = -1;
	edgeLeafToGIDMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLeafToGIDMapping_[dd] = -1;

	cellGIDToLeafMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) cellGIDToLeafMapping_[dd] = -1;
	cellLeafToGIDMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) cellLeafToGIDMapping_[dd] = -1;

	nrVertexLeafGID_ = 0; nrCellLeafGID_ = 0; nrEdgeLeafGID_ = 0;

	nrVertexLeafLID_ = 0; nrCellLeafLID_ = 0; nrEdgeLeafLID_ = 0;
	vertexLIDToLeafMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLIDToLeafMapping_[dd] = -1;
	vertexLeafToLIDMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLeafToLIDMapping_[dd] = -1;

	edgeLIDToLeafMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLIDToLeafMapping_[dd] = -1;
	edgeLeafToLIDMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLeafToLIDMapping_[dd] = -1;

	cellLIDToLeafMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) cellLIDToLeafMapping_[dd] = -1;
	cellLeafToLIDMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) cellLeafToLIDMapping_[dd] = -1;

	SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering , start assigning leaf numbers");

	for (int ind = 0 ; ind < nrElem_[0] ; ind++){
		//SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering point TID :" << ind << " maxCoFac: " << pointMaxCoF_[ind]);
	}

	// look for those leaf cells which points have a cell which maxCoFacet owner = myRank_
	// only those will have an LID
	Array<bool> hasCellLID(nrElem_[2],false);

	for (int ind = 0 ; ind < nrElem_[2] ; ind++){
		Array<int>& vertexIDs = cellsPoints_[ind];
		hasCellLID[ind] = false;
		for (int v = 0 ; v < 4 ; v++){
			Array<int>& maxCoFacet = pointMaxCoF_[vertexIDs[v]];
			hasCellLID[ind] =  ( hasCellLID[ind]
					|| ( (maxCoFacet[0] >= 0) && (elementOwner_[2][maxCoFacet[0]] == myRank_) )
                    || ( (maxCoFacet[1] >= 0) && (elementOwner_[2][maxCoFacet[1]] == myRank_) )
                    || ( (maxCoFacet[2] >= 0) && (elementOwner_[2][maxCoFacet[2]] == myRank_) )
                    || ( (maxCoFacet[3] >= 0) && (elementOwner_[2][maxCoFacet[3]] == myRank_) ) ) ;
			//SUNDANCE_MSG3(verb() , " Point ID"<< vertexIDs[v] << ", MacCoFacet:" << maxCoFacet);

			// add cells with hanging nodes which have contribution to element which are owned by this processor
			// if vertex is hanging look into the parent cell at the same index and if the owner is myRank_ then add
			// to the cells which should be processed
			if ( (hasCellLID[ind] == false) && (isPointHanging_[vertexIDs[v]] == true)){
				int parentID = parentCellLID_[ind];
				Array<int>& parentVertexIDs = cellsPoints_[parentID];
				hasCellLID[ind] = hasCellLID[ind] || (elementOwner_[0][parentVertexIDs[v]] == myRank_);
			}
		}
		SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering Cell ID :" << ind << " should be LID: " << hasCellLID[ind] <<
				" ,isCellLeaf_[ind]:" << isCellLeaf_[ind]);
	}

	//  treat special case, so that each hanging element has its parents
	// if we add one cell check hanging face, then add the maxCoF from the parent face if is leaf
	// if this is not successful then do the same thing for edges
	// - from each hanging edge there should be at least one cell on this processor which contains that parent edge !
	bool check_Ghost_cells = true;
	while (check_Ghost_cells){
		check_Ghost_cells = false;
	    for (int ind = 0 ; ind < nrElem_[2] ; ind++){
		   if ( (hasCellLID[ind] == true) && (elementOwner_[2][ind] != myRank_ ) ){
			  // check edges
			  Array<int>& edgeIDs = cellsEdges_[ind];
			  // we have this if only
		      for (int ii = 0 ; ii < 4 ; ii++ ){
				  // if the face is hanging and does not belong to me
				  if (isEdgeHanging_[edgeIDs[ii]] && ( elementOwner_[1][edgeIDs[ii]] != myRank_)){
                    // get parent cells same face
					int parentCell = parentCellLID_[ind];
					Array<int>& parentEdgesIDs = cellsEdges_[parentCell];
					for (int f = 0 ; f < 2 ; f++)
					if ( ( edgeMaxCoF_[parentEdgesIDs[ii]][f] >= 0 ) &&
						 ( elementOwner_[2][ edgeMaxCoF_[parentEdgesIDs[ii]][f] ] != myRank_ ) &&
						 ( hasCellLID[edgeMaxCoF_[parentEdgesIDs[ii]][f]] == false)   &&
						 ( isCellLeaf_[edgeMaxCoF_[parentEdgesIDs[ii]][f]] )
					   ){
						hasCellLID[edgeMaxCoF_[parentEdgesIDs[ii]][f]] = true;
						check_Ghost_cells = true;
					}
				  }
			   } // from loop
		   }
	    }
	}

	// we also have to list the cells which are not owned by the processor
	for (int ind = 0 ; ind < nrElem_[2] ; ind++)
	{
		 // GID numbering
		 // if cell is leaf and if is inside the computational domain
         if ( (isCellLeaf_[ind] == true) && (!isCellOut_[ind]) )
         {
        	 Array<int>& vertexIDs = cellsPoints_[ind];
           	 for (int v = 0; v < 4 ; v++)
             {
           		SUNDANCE_MSG3(verb() , " createLeafGIDNumbering  vertexIDs[v]:" << vertexIDs[v] );
            	if (vertexGIDToLeafMapping_[vertexIDs[v]] < 0)
            	{
            	   SUNDANCE_MSG3(verb() , " createLeafGIDNumbering -> VertexID:" << vertexIDs[v] << " , nrVertexLeafGID_:" << nrVertexLeafGID_ );
            	   vertexLeafToGIDMapping_[nrVertexLeafGID_] = vertexIDs[v];
            	   vertexGIDToLeafMapping_[vertexIDs[v]] = nrVertexLeafGID_;
            	   nrVertexLeafGID_++;
            	}
             }
        	 Array<int>& edgeIDs = cellsEdges_[ind];
        	 // for each edge check weather it already has a leaf index, if not create one
        	 for (int e = 0; e < 4 ; e++)
        	 {
        		 //SUNDANCE_MSG3(verb() , " createLeafNumbering  edgeLIDs[e]:" << edgeLIDs[e] );
        		 if (edgeGIDToLeafMapping_[edgeIDs[e]] < 0)
        		 {
        			 SUNDANCE_MSG3(verb() , " createLeafGIDNumbering -> edgeID:" << edgeIDs[e] << " , nrEdgeLeafGID_:" << nrEdgeLeafGID_ );
        			 //SUNDANCE_MSG3(verb() , " MaxCoFacet:" << edgeMaxCoF_[edgeLIDs[e]] << " edgeVertex:" << edgeVertex_[edgeLIDs[e]]);
        			 edgeLeafToGIDMapping_[nrEdgeLeafGID_] = edgeIDs[e];
        			 edgeGIDToLeafMapping_[edgeIDs[e]] = nrEdgeLeafGID_;
        			 nrEdgeLeafGID_++;
        		 }
        	 }
        	 // create leaf index for the leaf cell
			 SUNDANCE_MSG3(verb() , " createLeafGIDNumbering CELL cellID:" << ind << " , nrCellLeafGID_:" << nrCellLeafGID_ );
        	 cellLeafToGIDMapping_[nrCellLeafGID_] = ind;
        	 cellGIDToLeafMapping_[ind] = nrCellLeafGID_;
        	 nrCellLeafGID_++;

        	 // LID numbering
        	 // create leaf LID numbering , if this cell needs to be processed
        	 if (hasCellLID[ind]){
        		 // vertex
              	 for (int v = 0; v < 4 ; v++)
                 {
                	if (vertexLIDToLeafMapping_[vertexIDs[v]] < 0)
                	{
                	   SUNDANCE_MSG3(verb() , " createLeafLIDNumbering -> VertexID:" << vertexIDs[v] << " , nrVertexLeafLID_:" << nrVertexLeafLID_ );
                	   vertexLeafToLIDMapping_[nrVertexLeafLID_] = vertexIDs[v];
                	   vertexLIDToLeafMapping_[vertexIDs[v]] = nrVertexLeafLID_;
                	   nrVertexLeafLID_++;
                	}
                 }
            	 // for each edge check weather it already has a leaf index, if not create one
            	 for (int e = 0; e < 4 ; e++)
            	 {
            		 if (edgeLIDToLeafMapping_[edgeIDs[e]] < 0)
            		 {
            			 SUNDANCE_MSG3(verb() , " createLeafLIDNumbering -> edgeID:" << edgeIDs[e] << " , nrEdgeLeafLID_:" << nrEdgeLeafLID_ );
            			 edgeLeafToLIDMapping_[nrEdgeLeafLID_] = edgeIDs[e];
            			 edgeLIDToLeafMapping_[edgeIDs[e]] = nrEdgeLeafLID_;
            			 nrEdgeLeafLID_++;
            		 }
            	 }
            	 // create leaf index for the leaf cell
            	 SUNDANCE_MSG3(verb() , " createLeafLIDNumbering CELL cellID:" << ind << " , nrCellLeafLID_:" << nrCellLeafLID_ );
            	 cellLeafToLIDMapping_[nrCellLeafLID_] = ind;
            	 cellLIDToLeafMapping_[ind] = nrCellLeafLID_;
            	 nrCellLeafLID_++;
        	 }
         }
	}
	SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering , DONE");
}


// ====================================== OTHER LEAF NUMBERING ALGORITHM ==================

int HNMesh2D::estimateCellLoad(int ID){
	int rtn = 0;
	if (isCellLeaf_[ID]){
		if (!isCellOut_[ID]) rtn = 1;
	} else {
		// for each child call recursivly the function
		for (int r = 0 ; r < (int)cellsChildren_[ID].size() ; r++){
			rtn = rtn + estimateCellLoad(cellsChildren_[ID][r]);
		}
	}
    return rtn;
}

/** mark the cells and its facets for one processor*/
void HNMesh2D::markCellsAndFacets(int cellID , int procID){
	// mark the cell and the facets
	if (elementOwner_[2][cellID] < 0)  { elementOwner_[2][cellID] = procID; }
	//SUNDANCE_MSG3(verb() , "mark cell: " << cellID );
	if (elementOwner_[1][cellsEdges_[cellID][0]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][0]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][1]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][1]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][2]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][2]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][3]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][3]] = procID;}
	//SUNDANCE_MSG3(verb() , " ,mark edge: " << cellsEdges_[cellID][0] << " ,mark edge: " << cellsEdges_[cellID][1]
	//                      << " ,mark edge: " << cellsEdges_[cellID][2] << " ,mark edge: " << cellsEdges_[cellID][3]);
	if (elementOwner_[0][cellsPoints_[cellID][0]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][0]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][1]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][1]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][2]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][2]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][3]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][3]] = procID;}
	//SUNDANCE_MSG3(verb() , " ,mark point: " << cellsPoints_[cellID][0] << " ,mark point: " << cellsPoints_[cellID][1]
	//                      << " ,mark point: " << cellsPoints_[cellID][2] << " ,mark point: " << cellsPoints_[cellID][3] );
	if (!isCellLeaf_[cellID]){
		// for each child cell do it recursively
		for (int r = 0 ; r < (int)cellsChildren_[cellID].size() ; r++){
			markCellsAndFacets(cellsChildren_[cellID][r] , procID);
		}
	}
}

void HNMesh2D::createLeafNumbering_sophisticated(){

	// this array shows which cell will belong to this processor
	Array<bool> hasCellLID(nrElem_[2],false);
	double total_load = 0.0;
	int nrCoarseCell = _res_x * _res_y;
	Array<int> coarseCellLoad( _res_x * _res_y , 1 );

	// the principle for load is that each cell is one unit load
	// count the total number of cells which are inside the computational domain and are leaf cells
	// make a space filling curve traversal and assign each cell to one processor
	// on the coarser level make a Z-curve traversal, and there for each cell make a recursive traversal
	// distribute only the coarsest cells, since the tree traversal is not continuous
	// "elementOwner_" has to be changed!!!

	for (int ind = 0 ; ind < nrElem_[2] ; ind++){
        if (ind < nrCoarseCell) {
        	// estimate cells load
        	coarseCellLoad[ind] = estimateCellLoad(ind);
        }
		if ((isCellLeaf_[ind] == true) && (!isCellOut_[ind]) )
		{ total_load = total_load + 1 ; }
	}

	SUNDANCE_MSG3(verb() , "total_load = " << total_load << " , nrCell = " << nrElem_[2]);

	// generate the space filling curve traversal for a given level and unit square
	// and assign the coarsest cells to processors
	int levelM = ::ceil( ::fmax( ::log2(_res_x) , ::log2(_res_y ) ) );
	//int unitN = (int)::pow(2, levelM );
	Array<int> vectX1(4), vectY1(4), vectX2(4), vectY2(4);
	vectX1[0] = 0; vectX1[1] = (int)::pow(2,levelM-1); vectX1[2] = 0; vectX1[3] = (int)::pow(2,levelM-1);
	vectY1[0] = 0; vectY1[1] = 0; vectY1[2] = (int)::pow(2,levelM-1); vectY1[3] = (int)::pow(2,levelM-1);
	vectX2[0] = 0; vectX2[1] = (int)::pow(2,levelM-1); vectX2[2] = 0; vectX2[3] = (int)::pow(2,levelM-1);
	vectY2[0] = 0; vectY2[1] = 0; vectY2[2] = (int)::pow(2,levelM-1); vectY2[3] = (int)::pow(2,levelM-1);
	int addX[4] = { 0 , 1 , 0 , 1};
	int addY[4] = { 0 , 0 , 1 , 1};
	Array<int> *inX = &vectX1 , *inY = &vectY1 , *outX = &vectX2 , *outY = &vectY2 , *tmpVectP;
	int levelActual = levelM - 2;
	// this method generates the index for a unit square Z-curve traversal
	while (levelActual >= 0){
		outX->resize( 4 * inX->size() );
		outY->resize( 4 * inY->size() );
		int cI = 0 , addO = (int)::pow(2,levelActual);
		SUNDANCE_MSG3(verb() , " outX->size():" << outX->size() << ", levelActual:" << levelActual << " , addO:" << addO);
		// here create the 4 recursive cells
		for (int ce = 0 ; ce < inX->size() ; ce++){
			(*outX)[cI+0] = (*inX)[ce] + addO*addX[0];
			(*outX)[cI+1] = (*inX)[ce] + addO*addX[1];
			(*outX)[cI+2] = (*inX)[ce] + addO*addX[2];
			(*outX)[cI+3] = (*inX)[ce] + addO*addX[3];
			(*outY)[cI+0] = (*inY)[ce] + addO*addY[0];
			(*outY)[cI+1] = (*inY)[ce] + addO*addY[1];
			(*outY)[cI+2] = (*inY)[ce] + addO*addY[2];
			(*outY)[cI+3] = (*inY)[ce] + addO*addY[3];
			cI = cI + 4;
		}
		SUNDANCE_MSG3(verb() , " EX: " << (*outX)[0] << " , " << (*outX)[1] << " , " << (*outX)[2]);
		SUNDANCE_MSG3(verb() , " EY: " << (*outY)[0] << " , " << (*outY)[1] << " , " << (*outY)[2]);
		// decrease the level
		levelActual = levelActual - 1;
		tmpVectP = inX; inX = outX; outX = tmpVectP;
		tmpVectP = inY; inY = outY; outY = tmpVectP;
	}
	// switch the vectors back once we are finished
	tmpVectP = inX; inX = outX; outX = tmpVectP;
	tmpVectP = inY; inY = outY; outY = tmpVectP;

	// unmark the cells owners
	for (int tmp = 0 ; tmp < nrElem_[0] ; tmp++ ){ elementOwner_[0][tmp] = -1; }
	for (int tmp = 0 ; tmp < nrElem_[1] ; tmp++ ){ elementOwner_[1][tmp] = -1; }
	for (int tmp = 0 ; tmp < nrElem_[2] ; tmp++ ){ elementOwner_[2][tmp] = -1; }

	//mark the cells, vertex and edge to which cell they belong, recursively for each cell
	int coarseCellID , actProcID = 0 , actualLoad = 0;
	double loadPerProc = (double)total_load / (double)nrProc_ , diff_load = 0.0;
	for (int ind = 0 ; ind < outX->size() ; ind++){
		// first test the combinaiton if this is in the range
        if ( ((*outX)[ind] < _res_x) && ((*outY)[ind] < _res_y) ){
        	// !!!! --- here is very important that we compute the right index
        	coarseCellID = ((*outX)[ind])*_res_y + ((*outY)[ind]);
        	SUNDANCE_MSG3(verb(),"Z-curve trav. ind:" << ind << " , coarseCellID:" << coarseCellID << " , indX:" << (*outX)[ind] << " , indY:" << (*outY)[ind]);
        	//the level of this cell with the ID should be zero
        	TEUCHOS_TEST_FOR_EXCEPTION( cellLevel_[coarseCellID] > 0 , std::logic_error, " coarseCellID:" << coarseCellID << " has level:" << cellLevel_[coarseCellID] );
        	markCellsAndFacets( coarseCellID , actProcID);
        	actualLoad = actualLoad + coarseCellLoad[coarseCellID];
        	// increment the processor if necessary
    		if (((double)actualLoad >= (loadPerProc - 1e-8 - diff_load)) && ( actProcID < nrProc_-1 )){
    			SUNDANCE_MSG3(verb() , "Increase CPU , actualLoad:" << actualLoad << " loadPerProc:" << loadPerProc );
    			// compensate the load difference for the next CPU
    			diff_load = actualLoad - loadPerProc;
    			actProcID = actProcID + 1;
    			actualLoad = 0;
    		}
        }
	}

	// unmark the cells owners
	SUNDANCE_MSG3(verb()," nrElem_[0]:" << nrElem_[0] << " , nrElem_[1]:" << nrElem_[1] << " , nrElem_[2]" << nrElem_[2]);
	for (int tmp = 0 ; tmp < nrElem_[0] ; tmp++ ){ TEUCHOS_TEST_FOR_EXCEPTION( elementOwner_[0][tmp] < 0 , std::logic_error, " 0 tmp:" << tmp); }
	for (int tmp = 0 ; tmp < nrElem_[1] ; tmp++ ){ TEUCHOS_TEST_FOR_EXCEPTION( elementOwner_[1][tmp] < 0 , std::logic_error, " 1 tmp:" << tmp); }
	for (int tmp = 0 ; tmp < nrElem_[2] ; tmp++ ){ TEUCHOS_TEST_FOR_EXCEPTION( elementOwner_[2][tmp] < 0 , std::logic_error, " 2 tmp:" << tmp); }

// ==== what comes here is a code duplication from the method above ===========
	// set all leaf numbers to -1
	// - iterate trough the mesh and in the leaf cells , distribute leaf numbering
	// , detect if one cell is leaf ()
	// , have a tree similar tree traversal ... todo: later

	SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering nrPoint:" << nrElem_[0] << " , nrEdge:" << nrElem_[1] << ", nrCell:" << nrElem_[2]);
	// we resize the leafID - > global
	vertexGIDToLeafMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexGIDToLeafMapping_[dd] = -1;
	vertexLeafToGIDMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLeafToGIDMapping_[dd] = -1;

	edgeGIDToLeafMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeGIDToLeafMapping_[dd] = -1;
	edgeLeafToGIDMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLeafToGIDMapping_[dd] = -1;

	cellGIDToLeafMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) cellGIDToLeafMapping_[dd] = -1;
	cellLeafToGIDMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) cellLeafToGIDMapping_[dd] = -1;

	nrVertexLeafGID_ = 0; nrCellLeafGID_ = 0; nrEdgeLeafGID_ = 0;

	nrVertexLeafLID_ = 0; nrCellLeafLID_ = 0; nrEdgeLeafLID_ = 0;
	vertexLIDToLeafMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLIDToLeafMapping_[dd] = -1;
	vertexLeafToLIDMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLeafToLIDMapping_[dd] = -1;

	edgeLIDToLeafMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLIDToLeafMapping_[dd] = -1;
	edgeLeafToLIDMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLeafToLIDMapping_[dd] = -1;

	cellLIDToLeafMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) cellLIDToLeafMapping_[dd] = -1;
	cellLeafToLIDMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) cellLeafToLIDMapping_[dd] = -1;

	SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering , start assigning leaf numbers");

	for (int ind = 0 ; ind < nrElem_[2] ; ind++){
		Array<int>& vertexIDs = cellsPoints_[ind];
		hasCellLID[ind] = false;
		for (int v = 0 ; v < 4 ; v++){
			Array<int>& maxCoFacet = pointMaxCoF_[vertexIDs[v]];
			hasCellLID[ind] =  ( hasCellLID[ind]
					|| ( (maxCoFacet[0] >= 0) && (elementOwner_[2][maxCoFacet[0]] == myRank_) )
                    || ( (maxCoFacet[1] >= 0) && (elementOwner_[2][maxCoFacet[1]] == myRank_) )
                    || ( (maxCoFacet[2] >= 0) && (elementOwner_[2][maxCoFacet[2]] == myRank_) )
                    || ( (maxCoFacet[3] >= 0) && (elementOwner_[2][maxCoFacet[3]] == myRank_) ) ) ;
			//SUNDANCE_MSG3(verb() , " Point ID"<< vertexIDs[v] << ", MacCoFacet:" << maxCoFacet);

			// add cells with hanging nodes which have contribution to element which are owned by this processor
			// if vertex is hanging look into the parent cell at the same index and if the owner is myRank_ then add
			// to the cells which should be processed
			if ( (hasCellLID[ind] == false) && (isPointHanging_[vertexIDs[v]] == true)){
				int parentID = parentCellLID_[ind];
				Array<int>& parentVertexIDs = cellsPoints_[parentID];
				hasCellLID[ind] = hasCellLID[ind] || (elementOwner_[0][parentVertexIDs[v]] == myRank_);
			}
		}
		SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering Cell ID :" << ind << " should be LID: " << hasCellLID[ind] <<
				" ,isCellLeaf_[ind]:" << isCellLeaf_[ind]);
	}

	//  treat special case, so that each hanging element has its parents
	// if we add one cell check hanging face, then add the maxCoF from the parent face if is leaf
	// if this is not successful then do the same thing for edges
	// - from each hanging edge there should be at least one cell on this processor which contains that parent edge !
	bool check_Ghost_cells = true;
	while (check_Ghost_cells){
		check_Ghost_cells = false;
	    for (int ind = 0 ; ind < nrElem_[2] ; ind++){
		   if ( (hasCellLID[ind] == true) && (elementOwner_[2][ind] != myRank_ ) ){
			  // check edges
			  Array<int>& edgeIDs = cellsEdges_[ind];
			  // we have this if only
		      for (int ii = 0 ; ii < 4 ; ii++ ){
				  // if the face is hanging and does not belong to me
				  if (isEdgeHanging_[edgeIDs[ii]] && ( elementOwner_[1][edgeIDs[ii]] != myRank_)){
                    // get parent cells same face
					int parentCell = parentCellLID_[ind];
					Array<int>& parentEdgesIDs = cellsEdges_[parentCell];
					for (int f = 0 ; f < 2 ; f++)
					if ( ( edgeMaxCoF_[parentEdgesIDs[ii]][f] >= 0 ) &&
						 ( elementOwner_[2][ edgeMaxCoF_[parentEdgesIDs[ii]][f] ] != myRank_ ) &&
						 ( hasCellLID[edgeMaxCoF_[parentEdgesIDs[ii]][f]] == false)   &&
						 ( isCellLeaf_[edgeMaxCoF_[parentEdgesIDs[ii]][f]] )
					   ){
						hasCellLID[edgeMaxCoF_[parentEdgesIDs[ii]][f]] = true;
						check_Ghost_cells = true;
					}
				  }
			   } // from loop
		   }
	    }
	}

	// we also have to list the cells which are not owned by the processor
	for (int ind = 0 ; ind < nrElem_[2] ; ind++)
	{
		 // GID numbering
		 // if cell is leaf and if is inside the computational domain
         if ( (isCellLeaf_[ind] == true) && (!isCellOut_[ind]) )
         {
        	 Array<int>& vertexIDs = cellsPoints_[ind];
           	 for (int v = 0; v < 4 ; v++)
             {
           		SUNDANCE_MSG3(verb() , " createLeafGIDNumbering  vertexIDs[v]:" << vertexIDs[v] );
            	if (vertexGIDToLeafMapping_[vertexIDs[v]] < 0)
            	{
            	   SUNDANCE_MSG3(verb() , " createLeafGIDNumbering -> VertexID:" << vertexIDs[v] << " , nrVertexLeafGID_:" << nrVertexLeafGID_ );
            	   vertexLeafToGIDMapping_[nrVertexLeafGID_] = vertexIDs[v];
            	   vertexGIDToLeafMapping_[vertexIDs[v]] = nrVertexLeafGID_;
            	   nrVertexLeafGID_++;
            	}
             }
        	 Array<int>& edgeIDs = cellsEdges_[ind];
        	 // for each edge check weather it already has a leaf index, if not create one
        	 for (int e = 0; e < 4 ; e++)
        	 {
        		 //SUNDANCE_MSG3(verb() , " createLeafNumbering  edgeLIDs[e]:" << edgeLIDs[e] );
        		 if (edgeGIDToLeafMapping_[edgeIDs[e]] < 0)
        		 {
        			 SUNDANCE_MSG3(verb() , " createLeafGIDNumbering -> edgeID:" << edgeIDs[e] << " , nrEdgeLeafGID_:" << nrEdgeLeafGID_ );
        			 //SUNDANCE_MSG3(verb() , " MaxCoFacet:" << edgeMaxCoF_[edgeLIDs[e]] << " edgeVertex:" << edgeVertex_[edgeLIDs[e]]);
        			 edgeLeafToGIDMapping_[nrEdgeLeafGID_] = edgeIDs[e];
        			 edgeGIDToLeafMapping_[edgeIDs[e]] = nrEdgeLeafGID_;
        			 nrEdgeLeafGID_++;
        		 }
        	 }
        	 // create leaf index for the leaf cell
			 SUNDANCE_MSG3(verb() , " createLeafGIDNumbering CELL cellID:" << ind << " , nrCellLeafGID_:" << nrCellLeafGID_ );
        	 cellLeafToGIDMapping_[nrCellLeafGID_] = ind;
        	 cellGIDToLeafMapping_[ind] = nrCellLeafGID_;
        	 nrCellLeafGID_++;

        	 // LID numbering
        	 // create leaf LID numbering , if this cell needs to be processed
        	 if (hasCellLID[ind]){
        		 // vertex
              	 for (int v = 0; v < 4 ; v++)
                 {
                	if (vertexLIDToLeafMapping_[vertexIDs[v]] < 0)
                	{
                	   SUNDANCE_MSG3(verb() , " createLeafLIDNumbering -> VertexID:" << vertexIDs[v] << " , nrVertexLeafLID_:" << nrVertexLeafLID_ );
                	   vertexLeafToLIDMapping_[nrVertexLeafLID_] = vertexIDs[v];
                	   vertexLIDToLeafMapping_[vertexIDs[v]] = nrVertexLeafLID_;
                	   nrVertexLeafLID_++;
                	}
                 }
            	 // for each edge check weather it already has a leaf index, if not create one
            	 for (int e = 0; e < 4 ; e++)
            	 {
            		 if (edgeLIDToLeafMapping_[edgeIDs[e]] < 0)
            		 {
            			 SUNDANCE_MSG3(verb() , " createLeafLIDNumbering -> edgeID:" << edgeIDs[e] << " , nrEdgeLeafLID_:" << nrEdgeLeafLID_ );
            			 edgeLeafToLIDMapping_[nrEdgeLeafLID_] = edgeIDs[e];
            			 edgeLIDToLeafMapping_[edgeIDs[e]] = nrEdgeLeafLID_;
            			 nrEdgeLeafLID_++;
            		 }
            	 }
            	 // create leaf index for the leaf cell
            	 SUNDANCE_MSG3(verb() , " createLeafLIDNumbering CELL cellID:" << ind << " , nrCellLeafLID_:" << nrCellLeafLID_ );
            	 cellLeafToLIDMapping_[nrCellLeafLID_] = ind;
            	 cellLIDToLeafMapping_[ind] = nrCellLeafLID_;
            	 nrCellLeafLID_++;
        	 }
         }
	}
	SUNDANCE_MSG3(verb() , "HNMesh2D::createLeafNumbering , DONE");
}
