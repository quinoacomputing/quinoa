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
 * SundanceHNMesh3D.cpp
 *
 *  Created on: May 30, 2009
 *      Author: benk
 */

#include "SundanceHNMesh3D.hpp"

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

int HNMesh3D::offs_Points_x_[8] = {0, 1, 0, 1 , 0 , 1 , 0 , 1};

int HNMesh3D::offs_Points_y_[8] = {0, 0, 1, 1 , 0 , 0 , 1 , 1};

int HNMesh3D::offs_Points_z_[8] = {0, 0, 0, 0 , 1 , 1 , 1 , 1 };

int HNMesh3D::edge_Points_localIndex[12][2] = { {0,1} , {0,2} , {0,4} , {1,3} , {1,5} , {2,3} , {2,6} , {3,7} ,
		                                        {4,5} , {4,6} , {5,7} , {6,7} };

int HNMesh3D::edge_Orientation[12] = { 0, 1, 2, 1, 2, 0, 2, 2, 0, 1, 1, 0 };
int HNMesh3D::edge_MaxCofacetIndex[3][4] = { {0,5,8,11},{1,3,9,10},{2,4,6,7} };
int HNMesh3D::edge_MaxCof[12] = { 0,0,0, 1, 1, 1, 2, 3, 2, 2, 3, 3 };

int HNMesh3D::face_Points_localIndex[6][4] = { {0,1,2,3} , {0,1,4,5} , {0,2,4,6} , {1,3,5,7} , {2,3,6,7} , {4,5,6,7} };

int HNMesh3D::face_Edges_localIndex[6][4]= { {0,1,3,5} , {0,2,4,8} , {1,2,6,9} , {3,4,7,10}, {5,6,7,11}, {8,9,10,11}};

int HNMesh3D::face_Orientation[6] = { 0,1,2,2,1,0 };
int HNMesh3D::face_MaxCofacetIndex[3][2] = { {0,5},{1,4},{2,3}};
int HNMesh3D::face_MaxCof[6] = { 0,0,0,1,1,1 };

// -----------------------------------
int HNMesh3D::vInd[8];
int HNMesh3D::eInd[12];
int HNMesh3D::fInd[6];

// the X and the Y coordinates of the newly
double HNMesh3D::vertex_X[64] =
                      { 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 ,
                        0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 ,
                        0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 ,
                        0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 ,
                        0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 ,
                        0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 ,
                        0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 ,
                        0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 , 0.0 , 1.0/3.0 , 2.0/3.0 , 1.0 };

double HNMesh3D::vertex_Y[64] =
                      { 0.0 , 0.0 , 0.0 , 0.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 ,
                        2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
                        0.0 , 0.0 , 0.0 , 0.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 ,
                        2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
                        0.0 , 0.0 , 0.0 , 0.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 ,
                        2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
                        0.0 , 0.0 , 0.0 , 0.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 ,
                        2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 1.0 , 1.0 , 1.0 , 1.0 };
double HNMesh3D::vertex_Z[64] =
                      { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0  , 0.0 , 0.0 ,
		                0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0  , 0.0 , 0.0 ,
                        1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 ,
                        1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 , 1.0/3.0 ,
                        2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 ,
                        2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 , 2.0/3.0 ,
                        1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0  , 1.0 , 1.0 ,
                        1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0  , 1.0 , 1.0 };
// face index is above 20
int HNMesh3D::vertexToParentEdge[64]  =
                      { -1,  0,  0, -1,  1, 20, 20,  3,  1, 20, 20,  3, -1,  5,  5, -1,  2, 21, 21,  4,
                        22, -1, -1, 23, 22, -1, -1, 23,  6, 24, 24,  7,  2, 21, 21,  4, 22, -1, -1, 23,
                        22, -1, -1, 23,  6, 24, 24,  7, -1,  8,  8, -1,  9, 25, 25, 10,  9, 25, 25, 10,
                        -1, 11, 11, -1  };
//
int HNMesh3D::vertexInParentIndex[64]  =
                      { -1,  0,  1, -1,  0, 20, 21,  0,  1, 22, 23,  1, -1,  0,  1, -1,  0, 20, 21,  0,
                    	20, -1, -1, 20, 21, -1, -1, 21,  0, 20, 21,  0,  1, 22, 23,  1, 22, -1, -1, 22,
                    	23, -1, -1, 23,  1, 22, 23,  1, -1,  0,  1, -1,  0, 20, 21,  0,  1, 22, 23,  1,
                        -1,  0,  1, -1  };
//
int HNMesh3D::edgeToParentEdge[144]  =
                      {  0,  0,  0,  1, 20, 20,  3, 20, 20, 20,  1, 20, 20,  3, 20, 20, 20,  1, 20, 20,
                         3,  5,  5,  5,  2, 21, 21,  4, 22, -1, -1, 23, 22, -1, -1, 23,  6, 24, 24,  7,
                        21, 21, 21, 22, -1, -1, 23, -1, -1, -1, 22, -1, -1, 23, -1, -1, -1, 22, -1, -1,
                        23, 24, 24, 24,  2, 21, 21,  4, 22, -1, -1, 23, 22, -1, -1, 23,  6, 24, 24,  7,
                        21, 21, 21, 22, -1, -1, 23, -1, -1, -1, 22, -1, -1, 23, -1, -1, -1, 22, -1, -1,
                        23, 24, 24, 24,  2, 21, 21,  4, 22, -1, -1, 23, 22, -1, -1, 23,  6, 24, 24,  7,
                         8,  8,  8,  9, 25, 25, 10, 25, 25, 25,  9, 25, 25, 10, 25, 25, 25,  9, 25, 25,
                        10, 11, 11, 11  };
//
int HNMesh3D::edgeInParentIndex[144]  =
                      {  0,  1,  2,  0, 20, 21,  0, 22, 23, 24,  1, 25, 26,  1, 27, 28, 29,  2, 30, 31,
                    	 2,  0,  1,  2,  0, 20, 21,  0, 20, -1, -1, 20, 21, -1, -1, 21,  0, 20, 21,  0,
                    	22, 23, 24, 22, -1, -1, 22, -1, -1, -1, 23, -1, -1, 23, -1, -1, -1, 24, -1, -1,
                    	24, 22, 23, 24,  1, 25, 26,  1, 25, -1, -1, 25, 26, -1, -1, 26,  1, 25, 26,  1,
                        27, 28, 29, 27, -1, -1, 27, -1, -1, -1, 28, -1, -1, 28, -1, -1, -1, 29, -1, -1,
                    	29, 27, 28, 29,  2, 30, 31,  2, 30, -1, -1, 30, 31, -1, -1, 31,  2, 30, 31,  2,
                    	 0,  1,  2,  0, 20, 21,  0, 22, 23, 24,  1, 25, 26,  1, 27, 28, 29,  2, 30, 31,
                         2,  0,  1,  2  };
//
int HNMesh3D::faceToParentFace[108]  =
                      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  2, -1, -1,  3, -1, -1, -1,  2,
                        -1, -1,  3, -1, -1, -1,  2, -1, -1,  3,  4,  4,  4, -1, -1, -1, -1, -1, -1, -1,
                        -1, -1,  1,  1,  1,  2, -1, -1,  3, -1, -1, -1,  2, -1, -1,  3, -1, -1, -1,  2,
                        -1, -1,  3,  4,  4,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  2, -1,
                        -1,  3, -1, -1, -1,  2, -1, -1,  3, -1, -1, -1,  2, -1, -1,  3,  4,  4,  4,  5,
                         5,  5,  5,  5,  5,  5,  5,  5  };
//
int HNMesh3D::faceInParentIndex[108]  =
                      {  0,  1,  2,  3,  4,  5,  6,  7,  8,  0,  1,  2,  0, -1, -1,  0, -1, -1, -1,  1,
                    	-1, -1,  1, -1, -1, -1,  2, -1, -1,  2,  0,  1,  2, -1, -1, -1, -1, -1, -1, -1,
                    	-1, -1,  3,  4,  5,  3, -1, -1,  3, -1, -1, -1,  4, -1, -1,  4, -1, -1, -1,  5,
                    	-1, -1,  5,  3,  4,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1,  6,  7,  8,  6, -1,
                        -1,  6, -1, -1, -1,  7, -1, -1,  7, -1, -1, -1,  8, -1, -1,  8,  6,  7,  8,  0,
                    	 1,  2,  3,  4,  5,  6,  7,  8 };


HNMesh3D::HNMesh3D(int dim, const MPIComm& comm ,
	    const MeshEntityOrder& order)
: MeshBase(dim, comm , order), _comm(comm)
{
	setVerb(0);

	// get the number of processors
	nrProc_ = MPIComm::world().getNProc();
	myRank_ = MPIComm::world().getRank();
	//------ Point storage ----
	points_.resize(0);
    nrElem_.resize(4,0);
	nrElemOwned_.resize(4,0);
	//----- Facets -----
    cellsPoints_.resize(0);
    cellsEdges_.resize(0);
    cellsFaces_.resize(0);
    isCellOut_.resize(0);
    faceEdges_.resize(0);
    facePoints_.resize(0);
	edgePoints_.resize(0);
	edgeOrientation_.resize(0);
	faceOrientation_.resize(0);
	// ----- MaxCofacets ----
	faceMaxCoF_.resize(0);
	edgeMaxCoF_.resize(0);
    pointMaxCoF_.resize(0);
    //------ Element (processor) ownership -----
	elementOwner_.resize(4); elementOwner_[0].resize(0); elementOwner_[1].resize(0); elementOwner_[2].resize(0); elementOwner_[3].resize(0);
    //---- hierarchical storage -----
    indexInParent_.resize(0);
    parentCellLID_.resize(0);
	cellLevel_.resize(0);
	isCellLeaf_.resize(0);
	// ---- "hanging" info storage ---
	isPointHanging_.resize(0);
	isEdgeHanging_.resize(0);
	// ---- hanging element and refinement (temporary) storage ---
	edgeHangingElmStore_ = Hashtable< int, Array<int> >();
	hangingAccessCount_.resize(0);
	faceHangingElmStore_ = Hashtable< int, Array<int> >();
	refineCell_.resize(0);
    // set the leaf counter to zero
	nrCellLeafGID_ = 0; nrEdgeLeafGID_ = 0; nrFaceLeafGID_ = 0; nrVertexLeafGID_ = 0;
	nrVertexLeafLID_ = 0; nrCellLeafLID_ = 0; nrFaceLeafLID_ = 0; nrEdgeLeafLID_ = 0;
}

int HNMesh3D::numCells(int dim) const  {
	SUNDANCE_MSG3(verb(),"HNMesh3D::numCells(int dim):   dim:" << dim );
	switch (dim){
	case 0: return nrVertexLeafLID_;
	case 1: return nrEdgeLeafLID_;
	case 2: return nrFaceLeafLID_;
	case 3: return nrCellLeafLID_;
	}
	return 0;
}

Point HNMesh3D::nodePosition(int i) const {
	SUNDANCE_MSG3(verb(),"HNMesh3D::nodePosition(int i)   i:"<< i);
    // point do not need leaf indexing
	return points_[vertexLeafToLIDMapping_[i]];
}

const double* HNMesh3D::nodePositionView(int i) const {
	SUNDANCE_MSG3(verb(),"HNMesh3D::nodePositionView(int i)   i:" << i);
	//SUNDANCE_VERB_HIGH("nodePosition(int i)");
	return &(points_[vertexLeafToLIDMapping_[i]][0]);;
}

void HNMesh3D::getJacobians(int cellDim, const Array<int>& cellLID,
                          CellJacobianBatch& jBatch) const
{

	  SUNDANCE_MSG3(verb(),"HNMesh3D::getJacobians  cellDim:"<<cellDim<<" _x:"<<_ofs_x<<" _y:"<<_ofs_y<<" _z:"<<_ofs_z);
	  SUNDANCE_VERB_HIGH("getJacobians()");
	  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	  int nCells = cellLID.size();
	  int LID;
	  Point pnt(0.0,0.0,0.0);
	  jBatch.resize(cellLID.size(), spatialDim(), cellDim);
	  if (cellDim < spatialDim()) // they need the Jacobian of a lower dinemsional element
	  {
		   for (int i=0; i<nCells; i++)
		    {
		      double* detJ = jBatch.detJ(i);
		      switch(cellDim)
		      {
		        case 0:{ *detJ = 1.0;
		          break;}
		        case 1:{
				  LID = edgeLeafToLIDMapping_[cellLID[i]];
			      pnt = (points_[edgePoints_[LID][1]] - points_[edgePoints_[LID][0]]);
		          *detJ = sqrt(pnt * pnt); // the length of the edge
		        break;}
		        case 2:{
		          LID = faceLeafToLIDMapping_[cellLID[i]];
		          int a = facePoints_[LID][0];
		          int b = facePoints_[LID][1];
		          int c = facePoints_[LID][2];
		          const Point& pa = points_[a];
		          const Point& pb = points_[b];
		          const Point& pc = points_[c];
			      Point directedArea = cross( pb - pa , pc - pa );
		          *detJ = sqrt(directedArea * directedArea); // the area of the face
		        break;}
		        default:
		          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "  "cellDim=" << cellDim << " in HNMesh3D::getJacobians()");
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
		        case 3:{
				  LID = cellLeafToLIDMapping_[cellLID[i]];
				  // Jacobi for unstructured brick this will not work, but because of linear Jacoby we only have structured brick
		          J[0] = points_[cellsPoints_[LID][1]][0] - points_[cellsPoints_[LID][0]][0];
		          J[1] = 0.0; J[2] = 0.0; J[3] = 0.0;
		          J[4] = points_[cellsPoints_[LID][2]][1] - points_[cellsPoints_[LID][0]][1];
		          J[5] = 0.0; J[6] = 0.0; J[7] = 0.0;
		          J[8] = points_[cellsPoints_[LID][4]][2] - points_[cellsPoints_[LID][0]][2]; // the Jacobi of the brick cell
			      SUNDANCE_MSG3(verb() , "HNMesh3D::getJacobians LID:" << LID << " X:" << J[0] << " Y:" << J[4] << " Z:" << J[8]);
			      //SUNDANCE_MSG3(verb() , "HNMesh3D::getJacobians P0:" << points_[cellsPoints_[LID][0]] );
			      //SUNDANCE_MSG3(verb() , "HNMesh3D::getJacobians P1:" << points_[cellsPoints_[LID][1]] );
			      //SUNDANCE_MSG3(verb() , "HNMesh3D::getJacobians P2:" << points_[cellsPoints_[LID][2]] );
			      //SUNDANCE_MSG3(verb() , "HNMesh3D::getJacobians P4:" << points_[cellsPoints_[LID][4]] );
		        break;}
		        default:
		          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value " "cellDim=" << cellDim << " in HNMesh3D::getJacobians()");
		      }
		    }
	  }
}

void HNMesh3D::getCellDiameters(int cellDim, const Array<int>& cellLID,
                              Array<double>& cellDiameters) const {

	 TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	 SUNDANCE_VERB_HIGH("getCellDiameters()");
	  cellDiameters.resize(cellLID.size());
	  Point pnt(0.0 , 0.0 , 0.0 );
	  int LID;
	  if (cellDim < spatialDim())
	  {
		SUNDANCE_MSG3(verb(),"HNMesh3D::getCellDiameters(), cellDim < spatialDim() ");
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
	        case 2:  //length of the edge
				  LID = faceLeafToLIDMapping_[cellLID[i]];
			      pnt = (points_[facePoints_[LID][3]] - points_[facePoints_[LID][0]]);
			      cellDiameters[i] = sqrt(pnt * pnt); // the diameter of the face
	        break;
	        default:
	        	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "  "cellDim=" << cellDim << " in HNMesh3D::getCellDiameters()");
	      }
	    }
	  }
	  else
	  {
		SUNDANCE_MSG3(verb(),"HNMesh3D::getCellDiameters(), cellDim == spatialDim() ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
	        case 3:
	          LID = cellLeafToLIDMapping_[cellLID[i]];
	          pnt = points_[cellsPoints_[LID][7]] - points_[cellsPoints_[LID][0]];
	          cellDiameters[i] = sqrt(pnt * pnt);
	        break;
	        default:
	          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
	           "cellDim=" << cellDim  << " in HNMesh3D::getCellDiameters()");
	      }
	    }
	  }
}

void HNMesh3D::pushForward(int cellDim, const Array<int>& cellLID,
                         const Array<Point>& refQuadPts,
                         Array<Point>& physQuadPts) const {

	  SUNDANCE_MSG3(verb(),"HNMesh3D::pushForward cellDim:"<<cellDim);
	  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");

	  int nQuad = refQuadPts.size();
	  Point pnt( 0.0 , 0.0 , 0.0);
	  Point pnt1( 0.0 , 0.0 , 0.0);

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
	             int LID = faceLeafToLIDMapping_[cellLID[i]];
	             pnt = points_[facePoints_[LID][0]];
	             // this works only for structured, but we only work on structured quads
	             pnt1 = points_[facePoints_[LID][3]] - points_[facePoints_[LID][0]];
		         for (int q=0; q<nQuad; q++) {
		          	  physQuadPts.append( pnt + Point( refQuadPts[q][0] * pnt1[0] , refQuadPts[q][1] * pnt1[1] , refQuadPts[q][2] * pnt1[2]) );
		         }
		         break;}
	      case 3:{
	             int LID = cellLeafToLIDMapping_[cellLID[i]];
	             pnt = points_[cellsPoints_[LID][0]];
	             // this works only for structured, but we only work on structured quads
	             pnt1 = points_[cellsPoints_[LID][7]] - points_[cellsPoints_[LID][0]];
		         for (int q=0; q<nQuad; q++) {
		          	  physQuadPts.append( pnt + Point( refQuadPts[q][0] * pnt1[0] , refQuadPts[q][1] * pnt1[1] , refQuadPts[q][2] * pnt1[2] ) );
		         }
	      break;}
	      default:
	    	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value " "in HNMesh3D::getJacobians()");
	    }
	  }
}

int HNMesh3D::ownerProcID(int cellDim, int cellLID) const  {
	 int ID = -1;
	 if (cellDim == 0) ID = vertexLeafToLIDMapping_[cellLID];
     if (cellDim == 1) ID = edgeLeafToLIDMapping_[cellLID];
     if (cellDim == 2) ID = faceLeafToLIDMapping_[cellLID];
     if (cellDim == 3) ID = cellLeafToLIDMapping_[cellLID];
     SUNDANCE_MSG3(verb() , " HNMesh3D::ownerProcID ,cellDim:" << cellDim << ", cellLID:"
    		 << cellLID <<" ,ID:" << ID << ", ret:"<< elementOwner_[cellDim][ID] );
	 return elementOwner_[cellDim][ID];
}


int HNMesh3D::numFacets(int cellDim, int cellLID,
                      int facetDim) const  {
	//SUNDANCE_VERB_HIGH("numFacets()");
	if (cellDim==1) { // 1 dimension
         return 2; //one line has 2 points
    }
    else if (cellDim==2) { // 2 dimensions
         return 4; //one quad has 4 edges and 4 points
    }
    else if (cellDim==3) { // brick cell
    	if (facetDim == 0) return 8;
    	if (facetDim == 1) return 12;
    	if (facetDim == 2) return 6;
    }
	return -1;
}

int HNMesh3D::facetLID(int cellDim, int cellLID,
                     int facetDim, int facetIndex,
                     int& facetOrientation) const  {

	// todo: check weather facet orientation is right
	facetOrientation = 1;
	SUNDANCE_MSG3(verb(),"HNMesh3D::facetLID  cellDim:"<<cellDim<<", cellLID:"<<cellLID<<", facetDim:"<<facetDim<< ", facetIndex:" << facetIndex);
	int rnt = -1 , LID=-1 , tmp=-1;
	if (facetDim == 0){ // return the Number/ID of a Vertex
		if (cellDim == 3 ){
		    LID = cellLeafToLIDMapping_[cellLID];
		    rnt = cellsPoints_[LID][facetIndex]; tmp = rnt;
		    rnt = vertexLIDToLeafMapping_[rnt];
	    }
	    else if ((cellDim==2)){
		    LID = faceLeafToLIDMapping_[cellLID];
		    rnt = facePoints_[LID][facetIndex]; tmp = rnt;
		    rnt = vertexLIDToLeafMapping_[rnt];
	    }
	    else if ((cellDim==1)){
	        LID = edgeLeafToLIDMapping_[cellLID];
	        rnt = edgePoints_[LID][facetIndex]; tmp = rnt;
	        rnt = vertexLIDToLeafMapping_[rnt];
	    }
	} else if (facetDim == 1){
		if (cellDim == 3 ){
	        LID = cellLeafToLIDMapping_[cellLID];
	        rnt = cellsEdges_[LID][facetIndex]; tmp = rnt;
			rnt = edgeLIDToLeafMapping_[rnt];
	    } else if ((cellDim==2)){
	        LID = faceLeafToLIDMapping_[cellLID];
	        rnt = faceEdges_[LID][facetIndex]; tmp = rnt;
			rnt = edgeLIDToLeafMapping_[rnt];
	    }
	} else if (facetDim == 2){
		if (cellDim == 3 ){
	        LID = cellLeafToLIDMapping_[cellLID];
	        rnt = cellsFaces_[LID][facetIndex]; tmp = rnt;
			rnt = faceLIDToLeafMapping_[rnt];
	    }
	}
	SUNDANCE_MSG3(verb()," RET = " << rnt << ", LID:" << LID << ", tmp:" <<  tmp);
	return rnt;
}


void HNMesh3D::getFacetLIDs(int cellDim,
                          const Array<int>& cellLID,
                          int facetDim,
                          Array<int>& facetLID,
                          Array<int>& facetSign) const {

	SUNDANCE_MSG3(verb(),"HNMesh3D::getFacetLIDs()  cellDim:"<<cellDim<<"  cellLID.size():"<<cellLID.size()<<"  facetDim:" <<facetDim);
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
	SUNDANCE_MSG3(verb() ,"HNMesh3D::getFacetLIDs()  DONE. ");
}


const int* HNMesh3D::elemZeroFacetView(int cellLID) const {
    int LID = cellLeafToLIDMapping_[cellLID];
    SUNDANCE_MSG3(verb() , "HNMesh3D::elemZeroFacetView ");
	return (const int*)(&cellsPoints_[LID]);
}


int HNMesh3D::numMaxCofacets(int cellDim, int cellLID) const  {
	//SUNDANCE_VERB_HIGH("numMaxCofacets()");
	SUNDANCE_MSG3(verb() , "HNMesh3D::numMaxCofacets():  cellDim:"<<cellDim<<" cellLID:"<<cellLID );
	int rnt = -1;

	if (cellDim==0) { // point MaxCoFacet
		int LID = vertexLeafToLIDMapping_[cellLID];
        int sum = 0;
        SUNDANCE_MSG3(verb() ," pointMaxCoF_[LID] = " << pointMaxCoF_[LID] );
        for (int i = 0 ; i < 8 ; i++)
        	if ( (pointMaxCoF_[LID][i] >= 0) && ( cellLIDToLeafMapping_[pointMaxCoF_[LID][i]] >= 0) )
        		sum++;
        // return the value, how many cells has this point, on the leaf level
        rnt = sum;
    }
    else if (cellDim==1) { // edge MaxCoFacet
        int LID = edgeLeafToLIDMapping_[cellLID];
        int sum = 0;
        SUNDANCE_MSG3(verb() ," edgeMaxCoF_[LID] = " << edgeMaxCoF_[LID] );
        for (int i = 0 ; i < 4 ; i++)
        	if ( (edgeMaxCoF_[LID][i] >= 0) && ( cellLIDToLeafMapping_[edgeMaxCoF_[LID][i]] >= 0) )
        		sum++;
        // return the value, how many cells has this point, on the leaf level
        rnt = sum;
    }
    else if (cellDim==2) { // face MaxCoFacet
        int LID = faceLeafToLIDMapping_[cellLID];
        int sum = 0;
        SUNDANCE_MSG3(verb() ," faceMaxCoF_[LID] = " << faceMaxCoF_[LID] );
        for (int i = 0 ; i < 2 ; i++)
        	if ( (faceMaxCoF_[LID][i] >= 0) && ( cellLIDToLeafMapping_[faceMaxCoF_[LID][i]] >= 0) )
        		sum++;
        // return the value, how many cells has this point, on the leaf level
        rnt = sum;
    }
	SUNDANCE_MSG3(verb() ," RET = " << rnt );
	return rnt;
}


int HNMesh3D::maxCofacetLID(int cellDim, int cellLID,
                       int cofacetIndex,
                       int& facetIndex) const  {

	SUNDANCE_MSG3(verb() ,"HNMesh3D::maxCofacetLID() cellDim:"<<cellDim<<" cellLID:"<<cellLID<<" cofacetIndex:"<<cofacetIndex<< " facetIndex:"
			<< facetIndex);
	int rnt =-1;

	if (cellDim==0) { // 0 dimension
		//facetIndex = cofacetIndex;
		int actCoFacetIndex = 0;
    	int LID = vertexLeafToLIDMapping_[cellLID];
		for (int ii = 0 ; ii < 8 ; ii++){
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
    else if (cellDim==1) { // 1 dimensions
    	int maxCoFacet = 0;
        int LID = edgeLeafToLIDMapping_[cellLID];
    	int orientation = edgeOrientation_[LID];
		SUNDANCE_MSG3(verb() ," HNMesh3D::maxCofacetLID() 1 , orientation:" << orientation );
        // return the index in the vector, which later will be corrected later
		int actCoFacetIndex = 0;
		for (int ii = 0 ; ii < 4 ; ii++){
			// take this maxCoFacet only if that exist and is inside
			if ( (edgeMaxCoF_[LID][ii] >= 0) && (cellLIDToLeafMapping_[edgeMaxCoF_[LID][ii]] >= 0) ){
				if ( actCoFacetIndex < cofacetIndex ){
					actCoFacetIndex++;
				}else{
					facetIndex = ii;
					maxCoFacet = edgeMaxCoF_[LID][ii];
					break;
				}
			}
		}
		// calculate the correct facetIndex of the edge in the cell of cofacetIndex
		facetIndex = edge_MaxCofacetIndex[orientation][facetIndex];
		SUNDANCE_MSG3(verb() ,"HNMesh3D::maxCofacetLID() 1 , facetIndex:" << facetIndex << " maxCoFacet = " << maxCoFacet );
		rnt = ( maxCoFacet );
    }
    else if (cellDim==2) { // 2 dimensions
    	int maxCoFacet = 0;
        int LID = faceLeafToLIDMapping_[cellLID];
    	int orientation = faceOrientation_[LID];
		SUNDANCE_MSG3(verb() ," HNMesh3D::maxCofacetLID() 2 , orientation:" << orientation );
        // return the index in the vector, which later will be corrected later
		int actCoFacetIndex = 0;
		for (int ii = 0 ; ii < 2 ; ii++){
			// take this maxCoFacet only if that exist and is inside
			if ( (faceMaxCoF_[LID][ii] >= 0) && (cellLIDToLeafMapping_[faceMaxCoF_[LID][ii]] >= 0) ){
				if ( actCoFacetIndex < cofacetIndex ){
					actCoFacetIndex++;
				}else{
					facetIndex = ii;
					maxCoFacet = faceMaxCoF_[LID][ii];
					break;
				}
			}
		}
		// calculate the correct facetIndex of the edge in the cell of cofacetIndex
		facetIndex = face_MaxCofacetIndex[orientation][facetIndex];
		SUNDANCE_MSG3(verb() ,"HNMesh3D::maxCofacetLID() 2 , facetIndex:" << facetIndex << " maxCoFacet = " << maxCoFacet );
		rnt = ( maxCoFacet );
    }
	// transform back to leaf indexing
	rnt = cellLIDToLeafMapping_[ rnt ];

	SUNDANCE_MSG3(verb() ," RET = " << rnt << ",  facetIndex:" << facetIndex);
	return rnt;
}

void HNMesh3D::getCofacets(int cellDim, int cellLID,
                 int cofacetDim, Array<int>& cofacetLIDs) const {
	// Nothing to do
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh3D::getCofacets() not implemented");
}


void HNMesh3D::getMaxCofacetLIDs(const Array<int>& cellLIDs,
  MaximalCofacetBatch& cofacets) const {
	// nothing to do here
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh3D::getMaxCofacetLIDs() not implemented");
}


int HNMesh3D::mapGIDToLID(int cellDim, int globalIndex) const  {
	//SUNDANCE_VERB_HIGH("mapGIDToLID()");
	switch (cellDim){
	case 0:{
	         int ID = vertexLeafToGIDMapping_[globalIndex];
		     SUNDANCE_MSG3(verb() , " HNMesh3D::mapGIDToLID 0 , globalIndex:" << globalIndex << " ,ID:" << ID << ", ret:"<< vertexLIDToLeafMapping_[ID]);
	         return vertexLIDToLeafMapping_[ID];
		    break;}
	case 1:{
		     int ID = edgeLeafToGIDMapping_[globalIndex];
		     SUNDANCE_MSG3(verb() , " HNMesh3D::mapGIDToLID 1 , globalIndex:" << globalIndex << " ,ID:" << ID << ", ret:"<< edgeLIDToLeafMapping_[ID]);
		     return edgeLIDToLeafMapping_[ID];
		    break;}
	case 2:{
             int ID = faceLeafToGIDMapping_[globalIndex];
             SUNDANCE_MSG3(verb() , " HNMesh3D::mapGIDToLID 2 , globalIndex:" << globalIndex << " ,ID:" << ID << ", ret:"<< faceLIDToLeafMapping_[ID]);
             return faceLIDToLeafMapping_[ID];
		    break;}
	case 3:{
             int ID = cellLeafToGIDMapping_[globalIndex];
             SUNDANCE_MSG3(verb() , " HNMesh3D::mapGIDToLID 3 , globalIndex:" << globalIndex << " ,ID:" << ID << ", ret:"<< cellLIDToLeafMapping_[ID]);
             return cellLIDToLeafMapping_[ID];
		    break;}
	}
	return -1; //Wall
}


bool HNMesh3D::hasGID(int cellDim, int globalIndex) const {
	//SUNDANCE_VERB_HIGH("hasGID()");
	// we should always have all GIDs
	return true;
}


int HNMesh3D::mapLIDToGID(int cellDim, int localIndex) const  {
	//SUNDANCE_VERB_HIGH("mapLIDToGID()");
	switch (cellDim){
	case 0:{
	         int ID = vertexLeafToLIDMapping_[localIndex];
		     SUNDANCE_MSG3(verb() , " HNMesh3D::mapLIDToGID 0 , localIndex:" << localIndex << " ,ID:" << ID << ", ret:"<< vertexGIDToLeafMapping_[ID]);
	         return vertexGIDToLeafMapping_[ID];
		    break;}
	case 1:{
		     int ID = edgeLeafToLIDMapping_[localIndex];
		     SUNDANCE_MSG3(verb() , " HNMesh3D::mapLIDToGID 1 , localIndex:" << localIndex << " ,ID:" << ID << ", ret:"<< edgeGIDToLeafMapping_[ID]);
		     return edgeGIDToLeafMapping_[ID];
		    break;}
	case 2:{
             int ID = faceLeafToLIDMapping_[localIndex];
             SUNDANCE_MSG3(verb() , " HNMesh3D::mapLIDToGID 2 , localIndex:" << localIndex << " ,ID:" << ID << ", ret:"<< faceGIDToLeafMapping_[ID]);
             return faceGIDToLeafMapping_[ID];
		    break;}
	case 3:{
             int ID = cellLeafToLIDMapping_[localIndex];
             SUNDANCE_MSG3(verb() , " HNMesh3D::mapLIDToGID 3 , localIndex:" << localIndex << " ,ID:" << ID << ", ret:"<< cellGIDToLeafMapping_[ID]);
             return cellGIDToLeafMapping_[ID];
		    break;}
	}
	return -1; //Wall
}


CellType HNMesh3D::cellType(int cellDim) const  {
	 switch(cellDim)
	  {
	    case 0:  return PointCell;
	    case 1:  return LineCell;
	    case 2:  return QuadCell;
	    case 3:  return BrickCell;
	    default:
	      return NullCell; // -Wall
	  }
}


int HNMesh3D::label(int cellDim, int cellLID) const {
   // not used
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh3D::label() not implemented yet");
   return 0;
}


void HNMesh3D::getLabels(int cellDim, const Array<int>& cellLID,
		Array<int>& labels) const {
   // not used
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh3D::getLabels() not implemented yet");
}

Set<int> HNMesh3D::getAllLabelsForDimension(int cellDim) const {
   Set<int>                 rtn;
   // not used
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh3D::getAllLabelsForDimension() not implemented yet");
   return rtn;
}

void HNMesh3D::getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const {
    // not used
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh3D::getLIDsForLabel() not implemented yet");
}

void HNMesh3D::setLabel(int cellDim, int cellLID, int label) {
   // not used
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," HNMesh3D::setLabel() not implemented yet");
}


void HNMesh3D::assignIntermediateCellGIDs(int cellDim) {
	// The GIDs are assigned
}


bool HNMesh3D::hasIntermediateGIDs(int dim) const {
	// the mesh always has intermediate cells
	return true; // true means they have been synchronized ... not used now
}


// =============================== HANGING NODE FUNCTIONS ==========================
bool HNMesh3D::isElementHangingNode(int cellDim , int cellLID) const {
	SUNDANCE_MSG3(verb() ,"HNMesh3D::isElementHangingNode  cellDim:"<<cellDim<<" LID:"<< cellLID);
	if (cellDim==0) { // 1 dimension
    	int LID = vertexLeafToLIDMapping_[cellLID];
		return (isPointHanging_[LID]);
    }
    else if (cellDim==1) { // 2 dimensions
    	int LID = edgeLeafToLIDMapping_[cellLID];
        return (isEdgeHanging_[LID]);
    }
    else if (cellDim==2)
    {
    	int LID = faceLeafToLIDMapping_[cellLID];
        return (isFaceHanging_[LID] );
        // todo:
        //|| isEdgeHanging_[faceEdges_[LID][0]] || isEdgeHanging_[faceEdges_[LID][1]]
        //|| isEdgeHanging_[faceEdges_[LID][2]] || isEdgeHanging_[faceEdges_[LID][3]] );
    }
	return false; //Wall
}

int HNMesh3D::indexInParent(int maxCellLID) const {
	int ID = cellLeafToLIDMapping_[maxCellLID];
	int indexInPar = indexInParent_[ID];
	return indexInPar;

}

void HNMesh3D::returnParentFacets(  int childCellLID , int dimFacets ,
		                         Array<int> &facetsLIDs , int &parentCellLIDs ) const {
	int LID = cellLeafToLIDMapping_[childCellLID];
	parentCellLIDs = parentCellLID_[LID];

	SUNDANCE_MSG3( verb() , "HNMesh3D::returnParentFacets  childCellLID:"<<childCellLID<<" dimFacets:"<<dimFacets<<
			"  parentCellLIDs:"<< parentCellLIDs);
	//SUNDANCE_MSG3( verb() , "LID:"<<LID<<" parentCellLID_:"<<parentCellLID_);
	// this is the same for edges and for points
	if (dimFacets == 0){
		facetsLIDs.resize(8);
		for (int kuku = 0; kuku < 8 ; kuku++)
		facetsLIDs[kuku] = facetLID_tree( 3 , parentCellLIDs ,  dimFacets , kuku );
	}
	else if (dimFacets == 1){
		facetsLIDs.resize(12);
		for (int kuku = 0; kuku < 12 ; kuku++)
		facetsLIDs[kuku] = facetLID_tree( 3 , parentCellLIDs ,  dimFacets , kuku );
	}
	else if (dimFacets == 2){
		facetsLIDs.resize(6);
		for (int kuku = 0; kuku < 6 ; kuku++)
		facetsLIDs[kuku] = facetLID_tree( 3 , parentCellLIDs ,  dimFacets , kuku );
	}
	// map parent cell ID back to leaf indexing
	//parentCellLIDs = cellLIDToLeafMapping_[parentCellLIDs];

}

// only used in determining the parents
int HNMesh3D::facetLID_tree(int cellDim, int cellLID,
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
	} else if (facetDim == 2){
   	     rnt = cellsFaces_[cellLID][facetIndex];
	     rnt = faceLIDToLeafMapping_[rnt];
	     // rnt must be greater than 0
	}
	SUNDANCE_MSG3(verb() , "HNMesh3D::facetLID_tree cellDim:"<<cellDim<<", cellLID:"<<cellLID<<", facetDim:"<<facetDim<<
			", facetIndex:"<<facetIndex<<" RET = " << rnt );
	return rnt;
}

// =========================== MESH CREATION ========================================

/** adds one vertex to the mesh */
void HNMesh3D::addVertex(int vertexLID , int ownerProc , bool isHanging ,
		 double coordx , double coordy , double coordz , const Array<int> &maxCoF){
  // add only when the LID is new
  if (points_.size() <= vertexLID){
	  TEUCHOS_TEST_FOR_EXCEPTION(vertexLID != nrElem_[0] , std::logic_error ,"HNMesh3D::addVertex " <<
			 " vertexLID:" << vertexLID << " nrElem_[0]:" << nrElem_[0] );
     Point pt(coordx, coordy, coordz );
     points_.append( pt );
     pointMaxCoF_.append( maxCoF );
     isPointHanging_.append( isHanging );
     elementOwner_[0].append( (short int)ownerProc );
     SUNDANCE_MSG3(verb() , "HNMesh3D::addVertex: " << nrElem_[0] << " P:" << pt << " ,  maxCoF:" << maxCoF);
     SUNDANCE_MSG3(verb() , " ownerProc:" << ownerProc << " , isHanging:" << isHanging);
     nrElem_[0]++;
  }
}

/** adds one edge to the mesh */
void HNMesh3D::addEdge(int edgeLID , int ownerProc , bool isHanging , int edgeOrientation ,
        const Array<int> &vertexLIDs , const Array<int> &maxCoF){
	  // add only when the edgeLID is new
	  SUNDANCE_MSG3(verb() , "HNMesh3D -- addEdge: " << edgeLID << " nrElem_[1]: " << nrElem_[1] << " edgePoints_.size():" << edgePoints_.size() );
	  if (edgePoints_.size() <= edgeLID ){
		  TEUCHOS_TEST_FOR_EXCEPTION(edgeLID != nrElem_[1], std::logic_error, "HNMesh3D::addEdge edgeLID != nrElem_[1]");
		 edgePoints_.append( vertexLIDs );
		 edgeOrientation_.append( (short int)edgeOrientation );
		 edgeMaxCoF_.append( maxCoF );
		 isEdgeHanging_.append(isHanging);
		 hangingAccessCount_.append( (short int)0);
	     elementOwner_[1].append( (short int)ownerProc );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addEdge: " << nrElem_[1] << " vertexLIDs:" << vertexLIDs << " ,  maxCoF:" << maxCoF );
	     SUNDANCE_MSG3(verb() , "          ownerProc:" << ownerProc << ", isHanging:" << isHanging << ", edgeOrientation:" << edgeOrientation );
	     nrElem_[1]++;
	  }
}

void HNMesh3D::addFace(int faceLID , int ownerProc , bool isHanging , int faceOrientation ,
		        const Array<int> &vertexLIDs , const Array<int> &edgeLIDs ,
		        const Array<int> &maxCoF){

	  // add only when the edgeLID is new
	  if (facePoints_.size() <= faceLID ){
		 TEUCHOS_TEST_FOR_EXCEPTION(faceLID != nrElem_[2], std::logic_error, "HNMesh3D::addFace faceLID != nrElem_[2]");
		 facePoints_.append( vertexLIDs );
		 faceEdges_.append( edgeLIDs );
		 faceOrientation_.append( (short int)faceOrientation );
		 faceMaxCoF_.append( maxCoF );
		 isFaceHanging_.append(isHanging);
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addFace: " << nrElem_[2] << " vertexLIDs:" << vertexLIDs << " ,  maxCoF:" << maxCoF );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addFace ,  edgeLIDs:" << edgeLIDs );
	     SUNDANCE_MSG3(verb() , "          ownerProc:" << ownerProc << ", isHanging:" << isHanging << ", faceOrientation:" << faceOrientation);
	     elementOwner_[2].append( (short int)ownerProc );
	     nrElem_[2]++;
	  }
}

/** adds one cell(3D) to the mesh */
void HNMesh3D::addCell(int cellLID , int ownerProc ,
        int indexInParent , int parentCellLID , int level ,
        const Array<int> &faceLIDs , const Array<int> &edgeLIDs ,
        const Array<int> &vertexLIDs)
{
	  // add only when the edgeLID is new
	  if (cellsPoints_.size() <= cellLID ) {
		 TEUCHOS_TEST_FOR_EXCEPTION(cellLID != nrElem_[3], std::logic_error, "HNMesh3D::cellLID cellLID != nrElem_[3]");
		 cellsFaces_.append( faceLIDs );
		 cellsEdges_.append( edgeLIDs );
		 cellsPoints_.append( vertexLIDs );
	     indexInParent_.append( indexInParent );
	     parentCellLID_.append( parentCellLID );
	     cellLevel_.append( level );
	     isCellLeaf_.append( true );
	     refineCell_.append( 0 );
	     cellsChildren_.append( tuple(1) );
	     elementOwner_[3].append( (short int)ownerProc );
		 // calculate if the cell is complete outside the mesh domain
		 isCellOut_.append( !( meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[0]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[1]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[2]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[3]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[4]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[5]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[6]]) ||
				               meshDomain_.isInsideComputationalDomain(points_[vertexLIDs[7]]) ) );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell: " << nrElem_[3] <<
	    		 " vertexLIDs:" << vertexLIDs << " edgeLIDs:" << edgeLIDs << " faceLIDs:" << faceLIDs);
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell p0:" << points_[vertexLIDs[0]] );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell p1:" << points_[vertexLIDs[1]] );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell p2:" << points_[vertexLIDs[2]] );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell p3:" << points_[vertexLIDs[3]] );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell p4:" << points_[vertexLIDs[4]] );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell p5:" << points_[vertexLIDs[5]] );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell p6:" << points_[vertexLIDs[6]] );
	     SUNDANCE_MSG3(verb() , "HNMesh3D::addCell p7:" << points_[vertexLIDs[7]] );
		 SUNDANCE_MSG3(verb() , "HNMesh3D::addCell IN DOMAIN:" <<  isCellOut_[nrElem_[3]] );
	     nrElem_[3]++;
	  }
}

/** creates one regular mesh without refinement. With a different function the
 * refinement can start later , independently from this function. <br>
 * The structure of this mesh also supports unstructured storage of the cells,
 * so we might create unstructured mesh and later refine in the same way */
void HNMesh3D::createMesh(
                      double position_x,
			          double position_y,
			          double position_z,
			          double offset_x,
			          double offset_y,
			          double offset_z,
			          int resolution_x,
			          int resolution_y,
			          int resolution_z,
			          const RefinementClass& refineClass ,
			          const MeshDomainDef& meshDomain
){

	setVerb(0);

	// initialize object fields
	_pos_x = position_x; _pos_y = position_y; _pos_z = position_z;
	_ofs_x = offset_x; _ofs_y = offset_y;  _ofs_z = offset_z;
	_res_x = resolution_x; _res_y = resolution_y; _res_z = resolution_z;
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

void HNMesh3D::updateLocalCoarseNumbering(int ix, int iy , int iz , int Nx , int Ny){
/* ==== vertex indexing =====*/
	vInd[0] = (Nx+1)*(Ny+1)*iz + (Nx+1)*iy + ix;
	vInd[1] = vInd[0] + 1;
	vInd[2] = vInd[0] + Nx+1;
	vInd[3] = vInd[0] + Nx+2;
	vInd[4] = vInd[0] + (Nx+1)*(Ny+1);
	vInd[5] = vInd[4] + 1;
	vInd[6] = vInd[4] + Nx+1;
	vInd[7] = vInd[4] + Nx+2;
/*  ===== edge indexing ====== */
	eInd[0] = (Nx*(Ny+1) + Ny*(Nx+1) + (Nx+1)*(Ny+1))*iz + (Nx+1+Nx)*iy + ix;
	eInd[1] = eInd[0] + Nx;
	eInd[3] = eInd[0] + Nx+1;
	eInd[5] = eInd[0] + 2*Nx+1;
	int ed5 = (2*Nx+1)*(iy+1)+ix;
	eInd[2] = eInd[5]+(Nx*(Ny+1)+Ny*(Nx+1)  - ed5) + (Nx+1)*iy+ix;
	eInd[4] = eInd[2]+1;
	eInd[6] = eInd[2]+Nx+1;
	eInd[7] = eInd[2]+Nx+2;
	int ed7 = (Nx+1)*(iy+1)+ix+1;
	eInd[8] = eInd[7] + ( (Nx+1)*(Ny+1) - ed7) + (2*Nx+1)*iy+ix;
	eInd[9] = eInd[8] + Nx;
	eInd[10] = eInd[8] + Nx+1;
	eInd[11] = eInd[8] + 2*Nx+1;
/* ======== face numbering ======== */
	fInd[0] = (Nx*Ny+Nx*(Ny+1)+ Ny*(Nx+1))*iz+Nx*iy+ix;
	fInd[1] = fInd[0] + (Nx*Ny - (Nx*iy+ix)) + (iy*(2*Nx+1) + ix);
	fInd[2] = fInd[1] + Nx;
	fInd[3] = fInd[2] + 1;
	fInd[4] = fInd[3] + Nx;
	fInd[5] = fInd[4] + ((Nx+1)*Ny + (Ny+1)*Nx - (2*Nx+1)*(iy+1) - (ix)  + Nx*iy + ix);
}

void HNMesh3D::createCoarseMesh(){

	// estimate load for parallel case,
	// assign cells to each processors, based on the load and the domain
	// we assign cells to processors, (y,x) (optimal for vertical channel flow)
	// the estimation is done in each processor, globally but then only the local mesh will be build
	int nrCoarseCell = _res_x * _res_y * _res_z;
	int nrCoarsePoints = (_res_x+1)*(_res_y+1)*(_res_z+1);
	int nrCoarseEdge = (_res_x+1)*(_res_y+1)*_res_z + _res_x*(_res_y+1)*(_res_z+1) + (_res_x+1)*_res_y*(_res_z+1);
	int nrCoarseFace = (_res_x+1)*_res_y*_res_z + _res_x*(_res_y+1)*_res_z+ + _res_x*_res_y*(_res_z+1);
	Array<int> coarseCellOwner( nrCoarseCell , -1 );
	Array<int> coarsePointOwner( nrCoarsePoints , -1 );
	Array<int> coarseEdgeOwner( nrCoarseEdge , -1 );
	Array<int> coarseFaceOwner( nrCoarseFace , -1 );
	Array<int> coarseCellLID( nrCoarseCell , -1 );
	Array<int> coarsePointLID( nrCoarsePoints , -1 );
	Array<int> coarseEdgeLID( nrCoarseEdge , -1 );
	Array<int> coarseFaceLID( nrCoarseFace , -1 );
	Array<int> coarseCellsLoad( nrCoarseCell , 0 );
	int totalLoad = 0;

	SUNDANCE_MSG3(verb() , "HNMesh3D::createMesh nrCoarseCell:" << nrCoarseCell << " nrCoarsePoints:" << nrCoarsePoints
			          << " nrCoarseEdge:" << nrCoarseEdge <<  " nrCoarseFace:" << nrCoarseFace << " nrProc_:" << nrProc_ << " myRank_:" << myRank_);
	TEUCHOS_TEST_FOR_EXCEPTION( nrCoarseCell < nrProc_ , std::logic_error," HNMesh3D::createMesh nrCoarseCell < nrProc_ ");
	// now always divide as a flow channel , no resolution driven division

    // calculate total load and load per coarse cell
	double h[3];
	Array<int>  ind(3);
	h[0] = _ofs_x/(double)_res_x; h[1] = _ofs_y/(double)_res_y; h[2] = _ofs_z/(double)_res_z;
	Point pos(h[0],h[1],h[2]);
	Point res(h[0],h[1],h[2]);
	// estimate total estimated load of the mesh
	for (int i=0; i < nrCoarseCell; i++){
		// midpoint of the cell
		ind[0] = ((i % (_res_x*_res_y)) % _res_x);
		ind[1] = ((i % (_res_x*_res_y)) / _res_x);
		ind[2] = (i / (_res_x*_res_y));
		pos[0] = _pos_x + (double)(ind[0])*h[0] + 0.5*h[0];
		pos[1] = _pos_y + (double)(ind[1])*h[1] + 0.5*h[1];
		pos[2] = _pos_z + (double)(ind[2])*h[2] + 0.5*h[2];
		// todo: take the domain in consideration (take the 8 points) (when cells are turned off)
		coarseCellsLoad[i] = refineClass_.estimateRefinementLevel( pos , res );
		totalLoad += coarseCellsLoad[i];
	}

	// calculate average load per cell
	double loadPerProc = (double)totalLoad / (double)nrProc_;
	int actualProc=0;
	totalLoad = 0;
	// assign owners to the cells, edges, vertexes , greedy method
	SUNDANCE_MSG3(verb() , "Processor asign, loadPerProc:" << loadPerProc );
	double diff_load = 0.0;
	for (int i=0; i < nrCoarseCell; i++){
		ind[0] = ((i % (_res_x*_res_y)) % _res_x);
		ind[1] = ((i % (_res_x*_res_y)) / _res_x);
		ind[2] = (i / (_res_x*_res_y));
		// call the function to update
		updateLocalCoarseNumbering( ind[0] , ind[1] , ind[2] , _res_x , _res_y);
		//SUNDANCE_MSG3(verb() , "Cell ID:" << i << " vertexInd:" <<vertexInd << " edgeInd:" << edgeInd );
		//SUNDANCE_MSG3(verb() , "Cell, actual index" << ind  );
		// assign ownership for vertexes
		for (int jj = 0 ; jj < 8 ; jj++){
			if (coarsePointOwner[vInd[jj]] < 0){
				coarsePointOwner[vInd[jj]] = actualProc;
				SUNDANCE_MSG3(verb() , "Vertex CPU assign " << vInd[jj] << " ->" << actualProc );
			}
		}
		// assign ownership for edges
		for (int jj = 0 ; jj < 12 ; jj++){
			if (coarseEdgeOwner[eInd[jj]] < 0){
				coarseEdgeOwner[eInd[jj]] = actualProc;
				SUNDANCE_MSG3(verb() , "Edge CPU assign " << eInd[jj] << " ->" << actualProc );
			}
		}
		// assign ownership for faces
		for (int jj = 0 ; jj < 6 ; jj++){
			if (coarseFaceOwner[fInd[jj]] < 0){
				coarseFaceOwner[fInd[jj]] = actualProc;
				SUNDANCE_MSG3(verb() , "Face CPU assign " << fInd[jj] << " ->" << actualProc );
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

	// now go trough all cells which have to be added to this processor
	SUNDANCE_MSG3(verb() , " Process Cells:" << nrCoarseCell );
	for (int i=0; i < nrCoarseCell; i++)
	{
		ind[0] = ((i % (_res_x*_res_y)) % _res_x);
		ind[1] = ((i % (_res_x*_res_y)) / _res_x);
		ind[2] = (i / (_res_x*_res_y));
		// calculate the local index
		updateLocalCoarseNumbering( ind[0] , ind[1] , ind[2] , _res_x , _res_y );
		pos[0] = _pos_x + (double)(ind[0])*h[0];
		pos[1] = _pos_y + (double)(ind[1])*h[1];
		pos[2] = _pos_z + (double)(ind[2])*h[2];
		SUNDANCE_MSG3(verb() , "PCell ID:" << i <<" pos:"<<pos<<" _res_x:"<<_res_x<<" _res_y:"<<_res_y<<" _res_z:"<<_res_z);
		SUNDANCE_MSG3(verb() , "PCell, actual index" << ind  );
		// this condition is so that remote cells will be added
		int cellLID = coarseCellLID[i];
		Array<int> vLID(8,-1) , vertexMaxCoF(8,-1) ;
		Array<int> edgeLID(12,-1) , edgeVert(2,-1) , edgeMaxCoef(4,-1);
		Array<int> faceLID(6,-1) , faceVert(4,-1) , faceEdge(4,-1) , faceMaxCoef(2,-1);
		//SUNDANCE_MSG3(verb() , "Cell ID:" << i << " vertexInd:" <<vertexInd << " edgeInd:" << edgeInd << " cellLID:" << cellLID);
		//SUNDANCE_MSG3(verb() , "Cell, actual index" << ind << " pos:" << pos);
		// assign new cellID if necessary
		if (coarseCellLID[i] < 0 ){
			coarseCellLID[i] = nrElem_[3];
			cellLID = coarseCellLID[i];
		}
		// add all Vertexes , ignoring neighbor cells (maxcofacets)
        for (int jj = 0 ; jj < 8 ; jj++){
           	if (coarsePointLID[vInd[jj]] < 0){
           		coarsePointLID[vInd[jj]] = nrElem_[0];
           	}
           	vLID[jj] = coarsePointLID[vInd[jj]];
            // add vertex with -1 maxCOfacets
          	//SUNDANCE_MSG3(verb() , "Vertex  X:" << ((double)offs_Points_x_[jj])*h[0] << "  Y:" << pos[1] + ((double)offs_Points_y_[jj])*h[1]);
           	//SUNDANCE_MSG3(verb() , "Vertex  vLID[jj]:" << vLID[jj] << "  , coarsePointOwner[vertexInd+vertexOffs[jj]]" << coarsePointOwner[vertexInd+vertexOffs[jj]] );
           	addVertex( vLID[jj] , coarsePointOwner[vInd[jj]] , false ,
           			   pos[0] + ((double)offs_Points_x_[jj])*h[0] , pos[1] + ((double)offs_Points_y_[jj])*h[1] ,
           			   pos[2] + ((double)offs_Points_z_[jj])*h[2] ,
           	           vertexMaxCoF );
        }
		// add all Edges , ignoring neighbor cells (maxcofacets)
        for (int jj = 0 ; jj < 12 ; jj++){
           	if (coarseEdgeLID[eInd[jj]] < 0){
           		coarseEdgeLID[eInd[jj]] = nrElem_[1];
           	}
           	edgeLID[jj] = coarseEdgeLID[eInd[jj]];
           	edgeVert[0] = vLID[edge_Points_localIndex[jj][0]];
           	edgeVert[1] = vLID[edge_Points_localIndex[jj][1]];
           	SUNDANCE_MSG3(verb() , "Edge local Index:" << eInd[jj] );
            // add edge with -1 maxCOfacets
           	addEdge( edgeLID[jj] , coarseEdgeOwner[eInd[jj]] , false , edge_Orientation[jj] , edgeVert , edgeMaxCoef );
        }
		// add all Faces , ignoring neighbor cells (maxcofacets)
        for (int jj = 0 ; jj < 6 ; jj++){
           	if (coarseFaceLID[fInd[jj]] < 0){
           		coarseFaceLID[fInd[jj]] = nrElem_[2];
           	}
           	faceLID[jj] = coarseFaceLID[fInd[jj]];
           	//
           	faceVert[0] = vLID[face_Points_localIndex[jj][0]];
           	faceVert[1] = vLID[face_Points_localIndex[jj][1]];
           	faceVert[2] = vLID[face_Points_localIndex[jj][2]];
           	faceVert[3] = vLID[face_Points_localIndex[jj][3]];
           	faceEdge[0] = edgeLID[face_Edges_localIndex[jj][0]];
           	faceEdge[1] = edgeLID[face_Edges_localIndex[jj][1]];
           	faceEdge[2] = edgeLID[face_Edges_localIndex[jj][2]];
           	faceEdge[3] = edgeLID[face_Edges_localIndex[jj][3]];
            // add face with -1 maxCOfacets
           	addFace( faceLID[jj] , coarseFaceOwner[fInd[jj]] , false , face_Orientation[jj] ,
           		         faceVert , faceEdge , faceMaxCoef );
        }
		// add the Cell
        addCell( cellLID , coarseCellOwner[i] , 0 , cellLID , 0 , faceLID , edgeLID , vLID);
	} // --- end from for loop

	// next is maxCoFacet and boundary info update for vertexes and edges
	SUNDANCE_MSG3(verb() , " Process maxCofacets:" << nrCoarseCell );
	for (int i=0; i < nrCoarseCell; i++){
		Array<int> vLID(8,-1) , eLID(12,-1) ,fLID(6,-1) ;
		ind[0] = ((i % (_res_x*_res_y)) % _res_x);
		ind[1] = ((i % (_res_x*_res_y)) / _res_x);
		ind[2] = (i / (_res_x*_res_y));
		// calculate the local index
		updateLocalCoarseNumbering( ind[0] , ind[1] , ind[2] , _res_x , _res_y );
		int cellLID = coarseCellLID[i];
		SUNDANCE_MSG3(verb() , " MaxCoFs in Cell cellLID:" << cellLID  );
		//SUNDANCE_MSG3(verb() , "MaxCoFs in Cell ID:" << i << " vertexInd:" <<vertexInd << " edgeInd:" << edgeInd );
		//SUNDANCE_MSG3(verb() , "MaxCoFs, actual index" << ind  );
		// vertex maxCoFac
        for (int jj = 0 ; jj < 8 ; jj++){
           	vLID[jj] = coarsePointLID[vInd[jj]];
           	// if all elements are added then this is OK
           	pointMaxCoF_[vLID[jj]][jj] = cellLID;
           	SUNDANCE_MSG3(verb() , "Vertex MaxCoFacet vLID[jj]:" << vLID[jj] << " jj:" << jj << " cellLID:" << cellLID );
        }
        // edge maxCoFac
        for (int jj = 0 ; jj < 12 ; jj++){
           	eLID[jj] = coarseEdgeLID[eInd[jj]];
           	// if all elements are added then this is OK
           	edgeMaxCoF_[eLID[jj]][edge_MaxCof[jj]] = cellLID;
           	SUNDANCE_MSG3(verb() , "Edge MaxCoFacet eLID[jj]:" << eLID[jj] << " edge_MaxCof[jj]:" << edge_MaxCof[jj] << " cellLID:" << cellLID );
        }
        // face maxCoFac
        for (int jj = 0 ; jj < 6 ; jj++){
        	fLID[jj] = coarseFaceLID[fInd[jj]];
           	// if all elements are added then this is OK
           	faceMaxCoF_[fLID[jj]][face_MaxCof[jj]] = cellLID;
           	SUNDANCE_MSG3(verb() , "Face MaxCoFacet fLID[jj]:" << fLID[jj] << " face_MaxCof[jj]:" << face_MaxCof[jj] << " cellLID:" << cellLID );
        }
	}
	// basic rectangular mesh is build
}

// -----------
bool HNMesh3D::oneRefinementIteration(){

    int nrActualCell = nrElem_[3];
    bool rtn = false;
    SUNDANCE_MSG3(verb() , " HNMesh3D::oneRefinementIteration, start one refinement iteration cycle ");
    // we iterate only over the existing cells (not the ones which will be later created)
	for (int i=0 ; i < nrActualCell ; i++){
		// cell is owned by the current processor, and is leaf and is inside the mesh domain
	    SUNDANCE_MSG3(verb() , " Test cell " << i << ", elementOwner_[3][i]:" << elementOwner_[3][i] <<
	    		               ", isCellLeaf_[i]:" << isCellLeaf_[i] << ", out:" << (!isCellOut_[i]));

		if ( (isCellLeaf_[i] == true) )
			//  mark neighbors for refinements if because of the hanging nodes a refinement is not possible, regardless of the IN or OUT of Mesh domain
		{
			// check if refinement is needed and possible, none of the edges can be hanging (we could check also vertexes and faces as well)
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
			Point h = points_[cellsVertex[7]] - points_[cellsVertex[0]];
			Point p2 = points_[cellsVertex[0]] + 0.5*h;
			SUNDANCE_MSG3(verb() , " points_[cellsVertex[0]]: " << points_[cellsVertex[0]]);
			SUNDANCE_MSG3(verb() , " points_[cellsVertex[7]]: " << points_[cellsVertex[7]]);
			refFunction = ( (refineCell_[i] == 1) || refineClass_.refine( cellLevel_[i] , p2 , h ) );

            // decide if we refine this cell
			//SUNDANCE_OUT(cellLevel_[i] < 1 , " execute refinement on cell nr: " << i << ", refFunction:" << refFunction << " , p0:"
			//		<< points_[cellsVertex[0]] << " , h=" << h);
            SUNDANCE_MSG3(verb() , " execute refinement on cell nr: " << i << ", refFunction:" << refFunction);
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
                // we might use only edges
            	if (refFunction) {
                    // now take all hanging edge neighbors
                    for (int jj = 0 ; jj < cellsEdges.size() ; jj++)
                    	if (isEdgeHanging_[cellsEdges[jj]]){
                    		// get the parent cell
                    		int cLevel = cellLevel_[i];
                    		int parentID = parentCellLID_[i];
                    		int parentCEdge = cellsEdges_[parentID][jj];
                    		int refCell = -1;
                    		// take all maxCoFacets
                    		for (int coF = 0 ; coF < 4 ; coF++){
                    			refCell = edgeMaxCoF_[parentCEdge][coF];
								// refCell should be refined and mark for refinement
								if ( (refCell >= 0) && (cellLevel_[refCell] < cLevel) && (isCellLeaf_[refCell])){
									// todo: in some particular case this leads to error, (debug that error)
									// when in one points 2 different level meet (so we do a conservative refinement in 3D)
									refineCell_[refCell] = 1;
									rtn = true;
									//SUNDANCE_MSG3( verb() , " HNMesh3D::oneRefinementIteration refineCell_[refCell] = 1" << refCell);
								}
                    		}
                    	}
            	}
            }
		}
	}

	SUNDANCE_MSG3(verb() , " HNMesh3D::oneRefinementIteration DONE with one refinement iteration");
	// return if there was refinement or attempt to refine
	return rtn;
}

// -------- refine cell cellID (assuming all conditions are satisfied)--
void HNMesh3D::refineCell(int cellLID){

    // data initialization
	int cellOwner = elementOwner_[3][cellLID];
	Point p = points_[cellsPoints_[cellLID][0]];
	Point h = points_[cellsPoints_[cellLID][7]] - points_[cellsPoints_[cellLID][0]];

    // the created vertex, edge and cell LIDs
    Array<int> vertexLIDs(64,-1);
    // some of the vertexes already exist
    vertexLIDs[0] = cellsPoints_[cellLID][0];  vertexLIDs[3] = cellsPoints_[cellLID][1];
    vertexLIDs[12] = cellsPoints_[cellLID][2]; vertexLIDs[15] = cellsPoints_[cellLID][3];
    vertexLIDs[48] = cellsPoints_[cellLID][4]; vertexLIDs[51] = cellsPoints_[cellLID][5];
    vertexLIDs[60] = cellsPoints_[cellLID][6]; vertexLIDs[63] = cellsPoints_[cellLID][7];
    Array<int> edgeLIDs(144,-1);
    Array<int> faceLIDs(108,-1);
    // the cell IDs can be determined in advance
    Array<int> cellLIDs(27,-1);
    for (int c = 0; c < 27 ; cellLIDs[c] = nrElem_[3] + c , c++);

    // the parent cell is no leaf any more
	isCellLeaf_[cellLID] = false;

    SUNDANCE_MSG3(verb() , " ========== HNMesh3D::refineCell, cellLID:" << cellLID << " ==========");

    Array<int> ind(3,-1);
    for (int childCell = 0 ; childCell < 27 ; childCell++){
        Array<int> currCellV(8,-1);
        Array<int> currCellE(12,-1);
        Array<int> currCellF(6,-1);

    	// calculate index
    	ind[0] = ( (childCell % 9) % 3);
    	ind[1] = ( (childCell % 9) / 3);
    	ind[2] = ( childCell / 9);
    	// calculate all index
    	updateLocalCoarseNumbering( ind[0] , ind[1] , ind[2] , 3 , 3);

    	// ------------ VERTEX ----------------
    	for (int v = 0; v < 8 ; v++){
    		int  vertexOwner = cellOwner;
    		bool vertexIsHanging = false;
    		bool useEdge = true;
    		int pIndex = 0 , pOffset = 0 , pID = 0;
    		// look for existing vertex
    		if ( vertexToParentEdge[vInd[v]] >= 0){
    		  if ( vertexToParentEdge[vInd[v]] > 19){
    			useEdge = false;
    			pIndex = vertexToParentEdge[vInd[v]] - 20;
    			pID = cellsFaces_[cellLID][ pIndex ];
    			pOffset = vertexInParentIndex[vInd[v]] - 20;
    			SUNDANCE_MSG3(verb() , "Vertex -> Face vInd[v]:" << vInd[v] << " pIndex: " << pIndex << " pID:" << pID << " pOffset:" << pOffset);
    			SUNDANCE_MSG3(verb() , " cellsFaces:" << cellsFaces_[cellLID] << " elementOwner_[2][pID]:" << elementOwner_[2][pID]);
    			vertexOwner = elementOwner_[2][pID];
    			vertexIsHanging = ( numMaxCofacets_ID( 2 , pID ) > 1);
    		  } else {
    			pIndex = vertexToParentEdge[vInd[v]];
    			pID = cellsEdges_[cellLID][ pIndex ];
    			pOffset = vertexInParentIndex[vInd[v]];
    			SUNDANCE_MSG3(verb() , "Vertex -> Edge vInd[v]:" << vInd[v] << " pIndex: " << pIndex << " pID:" << pID << " pOffset:" << pOffset);
    			SUNDANCE_MSG3(verb() , " cellsEdges:" << cellsEdges_[cellLID] << " elementOwner_[1][pID]:" << elementOwner_[1][pID]);
    			vertexOwner = elementOwner_[1][pID];
    			vertexIsHanging = ( numMaxCofacets_ID( 1 , pID ) > 1);
    		  }
    		}
    		// see if the vertex already exists
    		if ( (vertexLIDs[vInd[v]] < 0) && (vertexToParentEdge[vInd[v]] >= 0) ){
    		       vertexLIDs[vInd[v]] = getHangingElement(0, useEdge, pID , pOffset );
    		}
			SUNDANCE_MSG3(verb() , "Vertex vInd[v]:" << vInd[v] << " pIndex: " << pIndex << " pID:" << pID << " pOffset:" << pOffset);
            // create vertex if is not already created
    		if ( vertexLIDs[vInd[v]] < 0){
    			Array<int> maxCoF(8,-1);
    			addVertex( nrElem_[0] , vertexOwner , vertexIsHanging ,
    					   p[0] + h[0]*vertex_X[vInd[v]] , p[1] + h[1]*vertex_Y[vInd[v]] , p[2] + h[2]*vertex_Z[vInd[v]] , maxCoF);
    			vertexLIDs[vInd[v]] = nrElem_[0] - 1;
    			// add to parent elem if this is hanging
    			if (vertexIsHanging)
    			   addHangingElement(0 , vertexLIDs[vInd[v]] , useEdge , pID , pOffset);
    		}
    		currCellV[v] = vertexLIDs[vInd[v]];
			// set MaxCoFacet
    		SUNDANCE_MSG3(verb() , " vertexLIDs[vInd[v]]: " << vertexLIDs[vInd[v]] );
    		pointMaxCoF_[ vertexLIDs[vInd[v]] ][v] = cellLIDs[childCell];
    		// set maxCofacet also in other directions  ... not necessary
    	}

    	// ------------ EDGES ----------------
    	for (int e = 0; e < 12 ; e++){
    		int  edgeOwner = cellOwner;
    		bool edgeIsHanging = false;
    		bool useEdge = true;
    		int pIndex = 0 , pOffset = 0 , pID = 0;
    		// look for existing edge
    		if ( edgeToParentEdge[eInd[e]] >= 0){
    		  if ( edgeToParentEdge[eInd[e]] > 19){
    			useEdge = false;
    			pIndex = edgeToParentEdge[eInd[e]] - 20;
    			pID = cellsFaces_[cellLID][ pIndex ];
    			pOffset = edgeInParentIndex[eInd[e]] - 20;
    			SUNDANCE_MSG3(verb() , "Edge -> Face eInd[e]:" << eInd[e] <<" pIndex: " << pIndex << " pID:" << pID << " pOffset:" << pOffset);
    			SUNDANCE_MSG3(verb() , " cellsFaces:" << cellsFaces_[cellLID] );
    			edgeOwner = elementOwner_[2][pID];
    			edgeIsHanging = ( numMaxCofacets_ID( 2 , pID ) > 1);
    		  } else {
    			pIndex = edgeToParentEdge[eInd[e]];
    			pID = cellsEdges_[cellLID][ pIndex ];
    			pOffset = edgeInParentIndex[eInd[e]];
    			SUNDANCE_MSG3(verb() , "Edge -> Edge eInd[e]:" << eInd[e] <<" pIndex: " << pIndex << " pID:" << pID << " pOffset:" << pOffset);
    			//SUNDANCE_MSG3(verb() , " cellsEdges:" << cellsEdges_[cellLID] );
    			edgeOwner = elementOwner_[1][pID];
    			edgeIsHanging = ( numMaxCofacets_ID( 1 , pID ) > 1);
    		  }
    		}
    		// see if the edge already exists
    		if ( (edgeLIDs[eInd[e]] < 0) && ( edgeToParentEdge[eInd[e]] >= 0) ){
    			edgeLIDs[eInd[e]] = getHangingElement( 1 , useEdge, pID , pOffset );
    		}
    		SUNDANCE_MSG3(verb() , "Edge eInd[e]:" << eInd[e] <<" pIndex: " << pIndex << " pID:" << pID << " pOffset:" << pOffset);
    		if ( edgeLIDs[eInd[e]] < 0){
    			Array<int> maxCoF(4,-1);
    			Array<int> edgeVertexs(2,-1);
    			edgeVertexs[0] = currCellV[edge_Points_localIndex[e][0]]; edgeVertexs[1] = currCellV[edge_Points_localIndex[e][1]];
        		// create edges for the current cell
    			addEdge( nrElem_[1] , edgeOwner , edgeIsHanging , edge_Orientation[e] , edgeVertexs , maxCoF );
    			edgeLIDs[eInd[e]] = nrElem_[1] - 1;
    			// add to parent elem if this is hanging
    			if (edgeIsHanging)
    				addHangingElement(1 , edgeLIDs[eInd[e]] , useEdge , pID , pOffset);
    		}
    		currCellE[e] = edgeLIDs[eInd[e]];
    		// set MaXCoFacet
    		edgeMaxCoF_[edgeLIDs[eInd[e]]][edge_MaxCof[e]] = cellLIDs[childCell];
    		SUNDANCE_MSG3(verb() , " edgeLIDs[eInd[e]]: " << edgeLIDs[eInd[e]] );
    		// set maxCofacet also in other directions  ... not necessary
    	}

    	// ------------ FACES ----------------
    	for (int f = 0; f < 6 ; f++){
    		int  faceOwner = cellOwner;
    		bool faceIsHanging = false;
    		bool useEdge = false;
    		int pIndex = 0 , pOffset = 0 , pID = 0;
    		// look for existing face
    		if ( faceToParentFace[fInd[f]] >= 0 ){
    		  pIndex = faceToParentFace[fInd[f]];
    		  pID = cellsFaces_[cellLID][ pIndex ];
    		  pOffset = faceInParentIndex[fInd[f]];
    		  SUNDANCE_MSG3(verb() , " cellsFaces:" << cellsFaces_[cellLID] );
    		  faceOwner = elementOwner_[2][pID];
    		  faceIsHanging = ( numMaxCofacets_ID( 2 , pID ) > 1);
    		}
    		// see if the face already exists
    		if ( (faceLIDs[fInd[f]] < 0) && ( faceToParentFace[fInd[f]] >= 0 )){
    			faceLIDs[fInd[f]] = getHangingElement( 2 , useEdge, pID , pOffset );
    		}
  		    SUNDANCE_MSG3(verb() , "Face -> Face pIndex: " << pIndex << " pID:" << pID << " pOffset:" << pOffset);
    		// add face if necessary
    		if ( faceLIDs[fInd[f]] < 0){
    			Array<int> maxCoF(2,-1);
    			Array<int> faceVertexs(4,-1);
    			faceVertexs[0] = currCellV[face_Points_localIndex[f][0]]; faceVertexs[1] = currCellV[face_Points_localIndex[f][1]];
    			faceVertexs[2] = currCellV[face_Points_localIndex[f][2]]; faceVertexs[3] = currCellV[face_Points_localIndex[f][3]];
    			Array<int> faceEdges(4,-1);
    			faceEdges[0] = currCellE[face_Edges_localIndex[f][0]]; faceEdges[1] = currCellE[face_Edges_localIndex[f][1]];
    			faceEdges[2] = currCellE[face_Edges_localIndex[f][2]]; faceEdges[3] = currCellE[face_Edges_localIndex[f][3]];
        		//create all faces for the current cell
    			addFace( nrElem_[2] , faceOwner , faceIsHanging , face_Orientation[f] , faceVertexs , faceEdges , maxCoF );
    			faceLIDs[fInd[f]] = nrElem_[2] - 1;
    			// add to parent elem if this is hanging
    			if (faceIsHanging)
    				addHangingElement(2 , faceLIDs[fInd[f]] , useEdge , pID , pOffset);
    		}
    		// set MaxCoFacet
    		currCellF[f] = faceLIDs[fInd[f]];
    		faceMaxCoF_[ currCellF[f] ][face_MaxCof[f]] = cellLIDs[childCell];
    		SUNDANCE_MSG3(verb() , " faceLIDs[fInd[f]]: " << faceLIDs[fInd[f]] );
    		// set maxCofacet also in other directions  ... not necessary
    	}

    	// ------------ CELL ----------------
    	addCell( nrElem_[3] , cellOwner , childCell , cellLID ,
    			cellLevel_[cellLID] + 1 , currCellF , currCellE , currCellV );
    	// append this child to the father list
    	cellsChildren_[cellLID].append(cellLIDs[childCell]);
    }
}

void HNMesh3D::addHangingElement(int cellDim, int cellID ,bool useEdge,
		int parentID , int parentOffset){
	// store in edge
	if (useEdge){
		if (!edgeHangingElmStore_.containsKey(parentID)){
			Array<int> hnInfo(6,-1);
			hnInfo[5] = numMaxCofacets_ID( 1 , parentID );
			hangingAccessCount_[parentID] = hnInfo[5];
			edgeHangingElmStore_.put( parentID , hnInfo );
		}
		const Array<int>& hnInfo = edgeHangingElmStore_.get(parentID);
		Array<int> tmp_vect(hnInfo.size());
		for (int ii = 0; ii < hnInfo.size() ; ii++) tmp_vect[ii] = hnInfo[ii];
        int realOffset = 0;
		// point in edge
		if (cellDim == 0){ realOffset = 0;}
		// edge in edge
		else{ realOffset = 2; }
		// store hanging element
  	    SUNDANCE_MSG3(verb() ,"HNMesh3D::addHangingElement EDGE  realOffset:" << realOffset << " parentOffset:" << parentOffset);
  	    SUNDANCE_MSG3(verb() ,"HNMesh3D::addHangingElement EDGE i:" << realOffset+parentOffset <<" hnInfo:" << hnInfo);
		tmp_vect[realOffset+parentOffset] = cellID;
		// store the new vector (with the updated information)
		edgeHangingElmStore_.remove( parentID );
		edgeHangingElmStore_.put( parentID , tmp_vect );
	}
	// store in face
	else{
		if (!faceHangingElmStore_.containsKey(parentID)){
			Array<int> hnInfo(25,-1);
			faceHangingElmStore_.put( parentID , hnInfo );
		}
		const Array<int>& hnInfo = faceHangingElmStore_.get(parentID);
		Array<int> tmp_vect(hnInfo.size());
        for (int ii = 0; ii < hnInfo.size() ; ii++) tmp_vect[ii] = hnInfo[ii];
        int realOffset = 0;
        // vertex in face
		if (cellDim == 0) { realOffset = 0; }
		// edge in face
		else if (cellDim == 1) { realOffset = 4; }
		// face in face
		else { realOffset = 16; }
		// store hanging element
  	    SUNDANCE_MSG3(verb() ,"HNMesh3D::addHangingElement FACE  realOffset:" << realOffset << " parentOffset:" << parentOffset);
  	    SUNDANCE_MSG3(verb() ,"HNMesh3D::addHangingElement FACE i:" << realOffset+parentOffset <<" hnInfo:" << hnInfo);
		tmp_vect[realOffset+parentOffset] = cellID;
		// store the new vector with the updated information
		faceHangingElmStore_.remove( parentID );
		faceHangingElmStore_.put( parentID , tmp_vect );
	}
}

int HNMesh3D::getHangingElement(int cellDim, bool useEdge, int parentID , int parentOffset){
	int rtn = -1;
	int realOffset = 0;
	if (useEdge){
		// get point from edge
		if (cellDim == 0){ realOffset = 0;}
		// get edge from edge
		else{ realOffset = 2; }
        if (edgeHangingElmStore_.containsKey(parentID)){
     	   const Array<int>& hnInfo = edgeHangingElmStore_.get(parentID);
     	   rtn = hnInfo[realOffset+parentOffset];
     	   SUNDANCE_MSG3(verb() ,"HNMesh3D::getHangingElement EDGE rtn = " << rtn << " realOffset:" << realOffset << " parentOffset:" << parentOffset);
     	   SUNDANCE_MSG3(verb() ,"HNMesh3D::getHangingElement EDGE i:" << realOffset+parentOffset <<" hnInfo:" << hnInfo);
     	   // mark stored elements as not hanging
     	   if ( rtn >= 0){
     		   if ((cellDim == 0) && (parentOffset == 1)) {
     			  if (hangingAccessCount_[parentID] <= 2){
     			  isPointHanging_[hnInfo[0]] = false; isPointHanging_[hnInfo[1]] = false;
     			  } else {
     				 //the edge will we accessed at the last time so that has to decrement
     			  }
     		   }
			   if ((cellDim == 1) && (parentOffset == 2)) {
	     		  if (hangingAccessCount_[parentID] <= 2){
	     			isEdgeHanging_[hnInfo[2]] = false; isEdgeHanging_[hnInfo[3]] = false;
	     			isEdgeHanging_[hnInfo[4]] = false;
	     		  } else {
	     			 hangingAccessCount_[parentID]--;
	     		  }
			   }
     	   }
        }
	}
	else
	{
		// get point from face
		if (cellDim == 0) { realOffset = 0; }
		// get edge from face
		else if (cellDim == 1) { realOffset = 4; }
		// get face from face
		else { realOffset = 16; }
        if (faceHangingElmStore_.containsKey(parentID)){
     	   const Array<int>& hnInfo = faceHangingElmStore_.get(parentID);
     	   rtn = hnInfo[realOffset+parentOffset];
     	   SUNDANCE_MSG3(verb() ,"HNMesh3D::getHangingElement FACE rtn = " << rtn << " realOffset:" << realOffset << " parentOffset:" << parentOffset);
     	   SUNDANCE_MSG3(verb() ,"HNMesh3D::getHangingElement FACE i:" << realOffset+parentOffset <<" hnInfo:" << hnInfo);
     	   // mark element as not hanging
     	   if ( rtn >= 0){
     		   if (cellDim == 0) {  isPointHanging_[rtn] = false; }
			   else if (cellDim == 1) { isEdgeHanging_[rtn] = false; }
			   else { isFaceHanging_[rtn] = false; }
     	   }
        }
	}
	return rtn; // Wall
}

int HNMesh3D::numMaxCofacets_ID(int cellDim, int cellID)
{
	int rtn = -1;
    if (cellDim==1) { // edge MaxCoFacet
        int LID = cellID;
        int sum = 0;
        SUNDANCE_MSG3(verb() ,"HNMesh3D::numMaxCofacets_ID ID:" << LID << " edgeMaxCoF_[LID] = " << edgeMaxCoF_[LID] );
        for (int i = 0 ; i < 4 ; i++){
        	if (edgeMaxCoF_[LID][i] >= 0)
        		sum++;
        }
        // return the value, how many cells has this point, on the leaf level
        rtn = sum;
    }
    else if (cellDim==2) { // face MaxCoFacet
        int LID = cellID;
        int sum = 0;
        SUNDANCE_MSG3(verb() ,"HNMesh3D::numMaxCofacets_ID ID:" << LID << " faceMaxCoF_[LID] = " << faceMaxCoF_[LID] );
        for (int i = 0 ; i < 2 ; i++){
        	if (faceMaxCoF_[LID][i] >= 0)
        		sum++;
        }
        // return the value, how many cells has this point, on the leaf level
        rtn = sum;
    }
    return rtn;
}

// -----------
void HNMesh3D::createLeafNumbering(){

	// set all leaf numbers to -1

	// - iterate trough the mesh and in the leaf cells , distribute leaf numbering
	// , detect if one cell is leaf ()
	// , have a tree similar tree traversal ... todo: later

	SUNDANCE_MSG3(verb() , "HNMesh3D::createLeafNumbering nrPoint:" << nrElem_[0] << " , nrEdge:"
			<< nrElem_[1] << ", nrFace:" << nrElem_[2] << ", nrCell:" << nrElem_[3]);
	// we resize the leafID - > global
	vertexGIDToLeafMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexGIDToLeafMapping_[dd] = -1;
	vertexLeafToGIDMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLeafToGIDMapping_[dd] = -1;

	edgeGIDToLeafMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeGIDToLeafMapping_[dd] = -1;
	edgeLeafToGIDMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLeafToGIDMapping_[dd] = -1;

	faceGIDToLeafMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) faceGIDToLeafMapping_[dd] = -1;
	faceLeafToGIDMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) faceLeafToGIDMapping_[dd] = -1;

	cellGIDToLeafMapping_.resize(nrElem_[3],-1);
	for (int dd = 0 ; dd < nrElem_[3] ; dd++) cellGIDToLeafMapping_[dd] = -1;
	cellLeafToGIDMapping_.resize(nrElem_[3],-1);
	for (int dd = 0 ; dd < nrElem_[3] ; dd++) cellLeafToGIDMapping_[dd] = -1;

	nrVertexLeafGID_ = 0; nrCellLeafGID_ = 0; nrEdgeLeafGID_ = 0; nrFaceLeafGID_ = 0;

	nrVertexLeafLID_ = 0; nrCellLeafLID_ = 0; nrEdgeLeafLID_ = 0;
	vertexLIDToLeafMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLIDToLeafMapping_[dd] = -1;
	vertexLeafToLIDMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLeafToLIDMapping_[dd] = -1;

	edgeLIDToLeafMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLIDToLeafMapping_[dd] = -1;
	edgeLeafToLIDMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLeafToLIDMapping_[dd] = -1;

	faceLIDToLeafMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) faceLIDToLeafMapping_[dd] = -1;
	faceLeafToLIDMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) faceLeafToLIDMapping_[dd] = -1;

	cellLIDToLeafMapping_.resize(nrElem_[3],-1);
	for (int dd = 0 ; dd < nrElem_[3] ; dd++) cellLIDToLeafMapping_[dd] = -1;
	cellLeafToLIDMapping_.resize(nrElem_[3],-1);
	for (int dd = 0 ; dd < nrElem_[3] ; dd++) cellLeafToLIDMapping_[dd] = -1;

	SUNDANCE_MSG3(verb() , "HNMesh3D::createLeafNumbering , start assigning leaf numbers");

	// look for those leaf cells which points have a cell which maxCoFacet owner = myRank_
	// only those will have an LID
	Array<bool> hasCellLID(nrElem_[3],false);

	for (int ind = 0 ; ind < nrElem_[3] ; ind++){
		Array<int>& vertexIDs = cellsPoints_[ind];
		hasCellLID[ind] = false;
		for (int v = 0 ; v < 8 ; v++){
			Array<int>& maxCoFacet = pointMaxCoF_[vertexIDs[v]];
			hasCellLID[ind] =  ( hasCellLID[ind]
					|| ( (maxCoFacet[0] >= 0) && (elementOwner_[3][maxCoFacet[0]] == myRank_) )
                    || ( (maxCoFacet[1] >= 0) && (elementOwner_[3][maxCoFacet[1]] == myRank_) )
                    || ( (maxCoFacet[2] >= 0) && (elementOwner_[3][maxCoFacet[2]] == myRank_) )
                    || ( (maxCoFacet[3] >= 0) && (elementOwner_[3][maxCoFacet[3]] == myRank_) )
                    || ( (maxCoFacet[4] >= 0) && (elementOwner_[3][maxCoFacet[4]] == myRank_) )
                    || ( (maxCoFacet[5] >= 0) && (elementOwner_[3][maxCoFacet[5]] == myRank_) )
                    || ( (maxCoFacet[6] >= 0) && (elementOwner_[3][maxCoFacet[6]] == myRank_) )
                    || ( (maxCoFacet[7] >= 0) && (elementOwner_[3][maxCoFacet[7]] == myRank_) )) ;

			// add cells with hanging nodes which have contribution to element which are owned by this processor
			// if vertex is hanging look into the parent cell at the same index and if the owner is myRank_ then add
			// to the cells which should be processed
			if ( (hasCellLID[ind] == false) && (isPointHanging_[vertexIDs[v]] == true)){
				int parentID = parentCellLID_[ind];
				Array<int>& parentVertexIDs = cellsPoints_[parentID];
				hasCellLID[ind] = hasCellLID[ind] || (elementOwner_[0][parentVertexIDs[v]] == myRank_);
			}
		}
		SUNDANCE_MSG3(verb() , "HNMesh3D::createLeafNumbering Cell ID :" << ind << " should be LID: " << hasCellLID[ind] <<
				" ,isCellLeaf_[ind]:" << isCellLeaf_[ind]);
	}

	//  treat special case, so that each hanging element has its parents
	// if we add one cell check hanging face, then add the maxCoF from the parent face if is leaf
	// if this is not successful then do the same thing for edges
	// - from each hanging edge there should be at least one cell on this processor which contains that parent edge !
	bool check_Ghost_cells = true;
	while (check_Ghost_cells){
		check_Ghost_cells = false;
	    for (int ind = 0 ; ind < nrElem_[3] ; ind++){
		   if ( (hasCellLID[ind] == true) && (elementOwner_[3][ind] != myRank_ ) ){
			  bool lookforEdges = true;
			  // check faces
			  Array<int>& faceIDs = cellsFaces_[ind];
			  for (int ii = 0 ; ii < 6 ; ii++ ){
				  // if the face is hanging and does not belong to me
				  if (isFaceHanging_[faceIDs[ii]] && ( elementOwner_[2][faceIDs[ii]] != myRank_)){
                    // get parent cells same face
					int parentCell = parentCellLID_[ind];
					Array<int>& parentfaceIDs = cellsFaces_[parentCell];
					for (int f = 0 ; f < 2 ; f++)
					if ( ( faceMaxCoF_[parentfaceIDs[ii]][f] >= 0 ) &&
						 ( elementOwner_[3][ faceMaxCoF_[parentfaceIDs[ii]][f] ] != myRank_ ) &&
						 ( hasCellLID[faceMaxCoF_[parentfaceIDs[ii]][f]] == false)  &&
						 ( isCellLeaf_[faceMaxCoF_[parentfaceIDs[ii]][f]] ) ){
						hasCellLID[faceMaxCoF_[parentfaceIDs[ii]][f]] = true;
						check_Ghost_cells = true;
						lookforEdges = false;
					}
				  }
			  }
			  // check edges
			  Array<int>& edgeIDs = cellsEdges_[ind];
			  // we have this if only
			  if (lookforEdges){
			    for (int ii = 0 ; ii < 12 ; ii++ ){
				  // if the face is hanging and does not belong to me
				  if (isEdgeHanging_[edgeIDs[ii]] && ( elementOwner_[1][edgeIDs[ii]] != myRank_)){
                    // get parent cells same face
					int parentCell = parentCellLID_[ind];
					Array<int>& parentEdgesIDs = cellsEdges_[parentCell];
					for (int f = 0 ; f < 4 ; f++)
					if ( ( edgeMaxCoF_[parentEdgesIDs[ii]][f] >= 0 ) &&
						 ( elementOwner_[3][ edgeMaxCoF_[parentEdgesIDs[ii]][f] ] != myRank_ ) &&
						 ( hasCellLID[edgeMaxCoF_[parentEdgesIDs[ii]][f]] == false)   &&
						 ( isCellLeaf_[edgeMaxCoF_[parentEdgesIDs[ii]][f]] )
					   ){
						hasCellLID[edgeMaxCoF_[parentEdgesIDs[ii]][f]] = true;
						check_Ghost_cells = true;
					}
				  }
			    } // from loop
			  }// from lookforEdges IF
		   }
	    }
	}

	// we also have to list the cells which are not owned by the processor
	for (int ind = 0 ; ind < nrElem_[3] ; ind++)
	{
		 // --------- GID numbering -----------
		 // if cell is leaf and if is inside the computational domain
         if ( (isCellLeaf_[ind] == true) && (!isCellOut_[ind]) )
         {
        	 Array<int>& vertexIDs = cellsPoints_[ind];
           	 for (int v = 0; v < 8 ; v++)
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
        	 for (int e = 0; e < 12 ; e++)
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
        	 Array<int>& faceIDs = cellsFaces_[ind];
        	 // for each face check weather it already has a leaf index, if not create one
        	 for (int f = 0; f < 6 ; f++)
        	 {
        		 //SUNDANCE_MSG3(verb() , " createLeafNumbering  edgeLIDs[e]:" << edgeLIDs[e] );
        		 if (faceGIDToLeafMapping_[faceIDs[f]] < 0)
        		 {
        			 SUNDANCE_MSG3(verb() , " createLeafGIDNumbering -> faceID:" << faceIDs[f] << " , nrFaceLeafGID_:" << nrFaceLeafGID_ );
        			 //SUNDANCE_MSG3(verb() , " MaxCoFacet:" << faceMaxCoF_[faceLIDs[e]] << " edgeVertex:" << faceVertex_[edgeLIDs[e]]);
        			 faceLeafToGIDMapping_[nrFaceLeafGID_] = faceIDs[f];
        			 faceGIDToLeafMapping_[faceIDs[f]] = nrFaceLeafGID_;
        			 nrFaceLeafGID_++;
        		 }
        	 }
        	 // create leaf index for the leaf cell
			 SUNDANCE_MSG3(verb() , " createLeafGIDNumbering CELL cellID:" << ind << " , nrCellLeafGID_:" << nrCellLeafGID_ );
        	 cellLeafToGIDMapping_[nrCellLeafGID_] = ind;
        	 cellGIDToLeafMapping_[ind] = nrCellLeafGID_;
        	 nrCellLeafGID_++;

        	 // --------- LID numbering -----------
        	 // create leaf LID numbering , if this cell needs to be processed
        	 if (hasCellLID[ind]){
        		 // vertex
              	 for (int v = 0; v < 8 ; v++)
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
            	 for (int e = 0; e < 12 ; e++)
            	 {
            		 if (edgeLIDToLeafMapping_[edgeIDs[e]] < 0)
            		 {
            			 SUNDANCE_MSG3(verb() , " createLeafLIDNumbering -> edgeID:" << edgeIDs[e] << " , nrEdgeLeafLID_:" << nrEdgeLeafLID_ );
            			 edgeLeafToLIDMapping_[nrEdgeLeafLID_] = edgeIDs[e];
            			 edgeLIDToLeafMapping_[edgeIDs[e]] = nrEdgeLeafLID_;
            			 nrEdgeLeafLID_++;
            		 }
            	 }
            	 // face LID
            	 for (int f = 0; f < 6 ; f++)
            	 {
            		 if (faceLIDToLeafMapping_[faceIDs[f]] < 0)
            		 {
            			 SUNDANCE_MSG3(verb() , " createLeafLIDNumbering -> faceID:" << faceIDs[f] << " , nrFaceLeafLID_:" << nrFaceLeafLID_ );
            			 faceLeafToLIDMapping_[nrFaceLeafLID_] = faceIDs[f];
            			 faceLIDToLeafMapping_[faceIDs[f]] = nrFaceLeafLID_;
            			 nrFaceLeafLID_++;
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
	SUNDANCE_MSG3(verb() , " nrVertexLeafGID_:" << nrVertexLeafGID_ << " nrEdgeLeafGID_:" << nrEdgeLeafGID_
			<< " nrFaceLeafGID_:" << nrFaceLeafGID_ << " nrCellLeafGID_:" << nrCellLeafGID_ );
	SUNDANCE_MSG3(verb() , " nrVertexLeafLID_:" << nrVertexLeafLID_ << " nrEdgeLeafLID_:" << nrEdgeLeafLID_
			<< " nrFaceLeafLID_:" <<nrFaceLeafLID_ << " nrCellLeafLID_:" << nrCellLeafLID_);
	SUNDANCE_MSG3(verb() , " vertexLIDToLeafMapping_: " << vertexLIDToLeafMapping_);
	SUNDANCE_MSG3(verb() , "HNMesh3D::createLeafNumbering , DONE");
}


// ====================================== OTHER LEAF NUMBERING ALGORITHM ==================

int HNMesh3D::estimateCellLoad(int ID){
	int rtn = 0;
	//SUNDANCE_MSG3( 5 , "HNMesh3D::estimateCellLoad , ID:" << ID);
	if (isCellLeaf_[ID]){
		if (!isCellOut_[ID]) rtn = 1;
	} else {
		// for each child call recursivly the function
		for (int r = 0 ; r < (int)cellsChildren_[ID].size() ; r++){
			if (cellLevel_[ID] < cellLevel_[cellsChildren_[ID][r]]){
				rtn = rtn + estimateCellLoad(cellsChildren_[ID][r]);
			}
		}
	}
    return rtn;
}

/** mark the cells and its facets for one processor*/
void HNMesh3D::markCellsAndFacets(int cellID , int procID){
	// mark the cell and the facets
	if (elementOwner_[3][cellID] < 0)  { elementOwner_[3][cellID] = procID; }
	//SUNDANCE_MSG3(verb() , "mark cell: " << cellID );
	if (elementOwner_[2][cellsFaces_[cellID][0]] < 0 ) { elementOwner_[2][cellsFaces_[cellID][0]] = procID;}
	if (elementOwner_[2][cellsFaces_[cellID][1]] < 0 ) { elementOwner_[2][cellsFaces_[cellID][1]] = procID;}
	if (elementOwner_[2][cellsFaces_[cellID][2]] < 0 ) { elementOwner_[2][cellsFaces_[cellID][2]] = procID;}
	if (elementOwner_[2][cellsFaces_[cellID][3]] < 0 ) { elementOwner_[2][cellsFaces_[cellID][3]] = procID;}
	if (elementOwner_[2][cellsFaces_[cellID][4]] < 0 ) { elementOwner_[2][cellsFaces_[cellID][4]] = procID;}
	if (elementOwner_[2][cellsFaces_[cellID][5]] < 0 ) { elementOwner_[2][cellsFaces_[cellID][5]] = procID;}

	if (elementOwner_[1][cellsEdges_[cellID][0]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][0]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][1]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][1]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][2]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][2]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][3]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][3]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][4]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][4]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][5]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][5]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][6]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][6]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][7]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][7]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][8]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][8]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][9]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][9]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][10]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][10]] = procID;}
	if (elementOwner_[1][cellsEdges_[cellID][11]] < 0 ) { elementOwner_[1][cellsEdges_[cellID][11]] = procID;}

	if (elementOwner_[0][cellsPoints_[cellID][0]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][0]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][1]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][1]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][2]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][2]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][3]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][3]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][4]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][4]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][5]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][5]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][6]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][6]] = procID;}
	if (elementOwner_[0][cellsPoints_[cellID][7]] < 0 ) { elementOwner_[0][cellsPoints_[cellID][7]] = procID;}

	if (!isCellLeaf_[cellID]){
		// for each child cell do it recursively
		for (int r = 0 ; r < (int)cellsChildren_[cellID].size() ; r++){
			if (cellLevel_[cellID] < cellLevel_[cellsChildren_[cellID][r]]){
				markCellsAndFacets(cellsChildren_[cellID][r] , procID);
			}
		}
	}
}

void HNMesh3D::createLeafNumbering_sophisticated(){

	// this array shows which cell will belong to this processor
	Array<bool> hasCellLID(nrElem_[3],false);
	double total_load = 0.0;
	//int nrCoarseCell = _res_x * _res_y * _res_z;
	Array<int> coarseCellLoad( _res_x * _res_y * _res_z , 1 );

	// the principle for load is that each cell is one unit load
	// count the total number of cells which are inside the computational domain and are leaf cells
	// make a space filling curve traversal and assign each cell to one processor
	// on the coarser level make a Z-curve traversal, and there for each cell make a recursive traversal
	// distribute only the coarsest cells, since the tree traversal is not continuous
	// "elementOwner_" has to be changed!!!

	for (int ind = 0 ; ind < nrElem_[3] ; ind++){
        if (cellLevel_[ind] < 1) {
        	// estimate cells load
        	coarseCellLoad[ind] = estimateCellLoad(ind);
        }
		if ((isCellLeaf_[ind] == true) && (!isCellOut_[ind]) )
		{ total_load = total_load + 1 ; }
	}

	SUNDANCE_MSG3(verb() , "total_load = " << total_load << " , nrCell = " << nrElem_[3]);

	// generate the space filling curve traversal for a given level and unit square
	// and assign the coarsest cells to processors
	int levelM = ::ceil( ::fmax( ::fmax( ::log2(_res_x) , ::log2(_res_y ) ) , ::log2(_res_z ) ) );
	//int unitN = (int)::pow(2, levelM );
	Array<int> vectX1(8), vectY1(8), vectZ1(8), vectX2(8), vectY2(8), vectZ2(8);
	vectX1[0] = 0; vectX1[1] = (int)::pow(2,levelM-1); vectX1[2] = 0; vectX1[3] = (int)::pow(2,levelM-1);
	vectX1[4] = 0; vectX1[5] = (int)::pow(2,levelM-1); vectX1[6] = 0; vectX1[7] = (int)::pow(2,levelM-1);
	vectY1[0] = 0; vectY1[1] = 0; vectY1[2] = (int)::pow(2,levelM-1); vectY1[3] = (int)::pow(2,levelM-1);
	vectY1[4] = 0; vectY1[5] = 0; vectY1[6] = (int)::pow(2,levelM-1); vectY1[7] = (int)::pow(2,levelM-1);
	vectZ1[0] = 0; vectZ1[1] = 0; vectZ1[2] = 0; vectZ1[3] = 0;
	vectZ1[4] = (int)::pow(2,levelM-1); vectZ1[5] = (int)::pow(2,levelM-1); vectZ1[6] = (int)::pow(2,levelM-1); vectZ1[7] = (int)::pow(2,levelM-1);

	vectX2[0] = 0; vectX2[1] = (int)::pow(2,levelM-1); vectX2[2] = 0; vectX2[3] = (int)::pow(2,levelM-1);
	vectX2[4] = 0; vectX2[5] = (int)::pow(2,levelM-1); vectX2[6] = 0; vectX2[7] = (int)::pow(2,levelM-1);
	vectY2[0] = 0; vectY2[1] = 0; vectY2[2] = (int)::pow(2,levelM-1); vectY2[3] = (int)::pow(2,levelM-1);
	vectY2[4] = 0; vectY2[5] = 0; vectY2[6] = (int)::pow(2,levelM-1); vectY2[7] = (int)::pow(2,levelM-1);
	vectZ2[0] = 0; vectZ2[1] = 0; vectZ2[2] = 0; vectZ2[3] = 0;
	vectZ2[4] = (int)::pow(2,levelM-1); vectZ2[5] = (int)::pow(2,levelM-1); vectZ2[6] = (int)::pow(2,levelM-1); vectZ2[7] = (int)::pow(2,levelM-1);

	int addX[8] = { 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1};
	int addY[8] = { 0 , 0 , 1 , 1 , 0 , 0 , 1 , 1};
	int addZ[8] = { 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1};
	Array<int> *inX = &vectX1 , *inY = &vectY1 , *inZ = &vectZ1 ,*outX = &vectX2 , *outY = &vectY2 , *outZ = &vectZ2 , *tmpVectP;
	int levelActual = levelM - 2;
	// this method generates the index for a unit square Z-curve traversal
	while (levelActual >= 0){
		outX->resize( 8 * inX->size() );
		outY->resize( 8 * inY->size() );
		outZ->resize( 8 * inZ->size() );
		int cI = 0 , addO = (int)::pow(2,levelActual);
		SUNDANCE_MSG3(verb() , " outX->size():" << outX->size() << ", levelActual:" << levelActual << " , addO:" << addO);
		// here create the 8 recursive cells
		for (int ce = 0 ; ce < inX->size() ; ce++){
			(*outX)[cI+0] = (*inX)[ce] + addO*addX[0];  (*outY)[cI+0] = (*inY)[ce] + addO*addY[0];  (*outZ)[cI+0] = (*inZ)[ce] + addO*addZ[0];
			(*outX)[cI+1] = (*inX)[ce] + addO*addX[1];  (*outY)[cI+1] = (*inY)[ce] + addO*addY[1];  (*outZ)[cI+1] = (*inZ)[ce] + addO*addZ[1];
			(*outX)[cI+2] = (*inX)[ce] + addO*addX[2];  (*outY)[cI+2] = (*inY)[ce] + addO*addY[2];  (*outZ)[cI+2] = (*inZ)[ce] + addO*addZ[2];
			(*outX)[cI+3] = (*inX)[ce] + addO*addX[3];  (*outY)[cI+3] = (*inY)[ce] + addO*addY[3];  (*outZ)[cI+3] = (*inZ)[ce] + addO*addZ[3];
			(*outX)[cI+4] = (*inX)[ce] + addO*addX[4];  (*outY)[cI+4] = (*inY)[ce] + addO*addY[4];  (*outZ)[cI+4] = (*inZ)[ce] + addO*addZ[4];
			(*outX)[cI+5] = (*inX)[ce] + addO*addX[5];  (*outY)[cI+5] = (*inY)[ce] + addO*addY[5];  (*outZ)[cI+5] = (*inZ)[ce] + addO*addZ[5];
			(*outX)[cI+6] = (*inX)[ce] + addO*addX[6];  (*outY)[cI+6] = (*inY)[ce] + addO*addY[6];  (*outZ)[cI+6] = (*inZ)[ce] + addO*addZ[6];
			(*outX)[cI+7] = (*inX)[ce] + addO*addX[7];  (*outY)[cI+7] = (*inY)[ce] + addO*addY[7];  (*outZ)[cI+7] = (*inZ)[ce] + addO*addZ[7];
			cI = cI + 8;
		}
		SUNDANCE_MSG3(verb() , " EX: " << (*outX)[0] << " , " << (*outX)[1] << " , " << (*outX)[2]);
		SUNDANCE_MSG3(verb() , " EY: " << (*outY)[0] << " , " << (*outY)[1] << " , " << (*outY)[2]);
		SUNDANCE_MSG3(verb() , " EZ: " << (*outZ)[0] << " , " << (*outZ)[1] << " , " << (*outZ)[2]);
		// decrease the level
		levelActual = levelActual - 1;
		tmpVectP = inX; inX = outX; outX = tmpVectP;
		tmpVectP = inY; inY = outY; outY = tmpVectP;
		tmpVectP = inZ; inZ = outZ; outZ = tmpVectP;
	}
	// switch the vectors back once we are finished
	tmpVectP = inX; inX = outX; outX = tmpVectP;
	tmpVectP = inY; inY = outY; outY = tmpVectP;
	tmpVectP = inZ; inZ = outZ; outZ = tmpVectP;

	// unmark the cells owners
	for (int tmp = 0 ; tmp < nrElem_[0] ; tmp++ ){ elementOwner_[0][tmp] = -1; }
	for (int tmp = 0 ; tmp < nrElem_[1] ; tmp++ ){ elementOwner_[1][tmp] = -1; }
	for (int tmp = 0 ; tmp < nrElem_[2] ; tmp++ ){ elementOwner_[2][tmp] = -1; }
	for (int tmp = 0 ; tmp < nrElem_[3] ; tmp++ ){ elementOwner_[3][tmp] = -1; }

	//mark the cells, vertex and edge to which cell they belong, recursively for each cell
	int coarseCellID , actProcID = 0 , actualLoad = 0;
	double loadPerProc = (double)total_load / (double)nrProc_ , diff_load = 0.0;
	for (int ind = 0 ; ind < outX->size() ; ind++){
		// first test the combinaiton if this is in the range
        if ( ((*outX)[ind] < _res_x) && ((*outY)[ind] < _res_y) && ((*outZ)[ind] < _res_z) ){
        	// !!!! --- here is very important that we compute the right index
        	coarseCellID = ((*outZ)[ind])*_res_x*_res_y + ((*outY)[ind])*_res_x + ((*outX)[ind]) ;
        	SUNDANCE_MSG3(verb(),"Z-curve trav. ind:" << ind << " , coarseCellID:" << coarseCellID
        			<< " , indX:" << (*outX)[ind] << " , indY:" << (*outY)[ind] << " , indZ:" << (*outZ)[ind]);
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
	for (int tmp = 0 ; tmp < nrElem_[3] ; tmp++ ){ TEUCHOS_TEST_FOR_EXCEPTION( elementOwner_[3][tmp] < 0 , std::logic_error, " 3 tmp:" << tmp); }

// ==== what comes here is a code duplication from the method above ===========

	SUNDANCE_MSG3(verb() , "HNMesh3D::createLeafNumbering nrPoint:" << nrElem_[0] << " , nrEdge:"
			<< nrElem_[1] << ", nrFace:" << nrElem_[2] << ", nrCell:" << nrElem_[3]);
	// we resize the leafID - > global
	vertexGIDToLeafMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexGIDToLeafMapping_[dd] = -1;
	vertexLeafToGIDMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLeafToGIDMapping_[dd] = -1;

	edgeGIDToLeafMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeGIDToLeafMapping_[dd] = -1;
	edgeLeafToGIDMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLeafToGIDMapping_[dd] = -1;

	faceGIDToLeafMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) faceGIDToLeafMapping_[dd] = -1;
	faceLeafToGIDMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) faceLeafToGIDMapping_[dd] = -1;

	cellGIDToLeafMapping_.resize(nrElem_[3],-1);
	for (int dd = 0 ; dd < nrElem_[3] ; dd++) cellGIDToLeafMapping_[dd] = -1;
	cellLeafToGIDMapping_.resize(nrElem_[3],-1);
	for (int dd = 0 ; dd < nrElem_[3] ; dd++) cellLeafToGIDMapping_[dd] = -1;

	nrVertexLeafGID_ = 0; nrCellLeafGID_ = 0; nrEdgeLeafGID_ = 0; nrFaceLeafGID_ = 0;

	nrVertexLeafLID_ = 0; nrCellLeafLID_ = 0; nrEdgeLeafLID_ = 0;
	vertexLIDToLeafMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLIDToLeafMapping_[dd] = -1;
	vertexLeafToLIDMapping_.resize(nrElem_[0],-1);
	for (int dd = 0 ; dd < nrElem_[0] ; dd++) vertexLeafToLIDMapping_[dd] = -1;

	edgeLIDToLeafMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLIDToLeafMapping_[dd] = -1;
	edgeLeafToLIDMapping_.resize(nrElem_[1],-1);
	for (int dd = 0 ; dd < nrElem_[1] ; dd++) edgeLeafToLIDMapping_[dd] = -1;

	faceLIDToLeafMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) faceLIDToLeafMapping_[dd] = -1;
	faceLeafToLIDMapping_.resize(nrElem_[2],-1);
	for (int dd = 0 ; dd < nrElem_[2] ; dd++) faceLeafToLIDMapping_[dd] = -1;

	cellLIDToLeafMapping_.resize(nrElem_[3],-1);
	for (int dd = 0 ; dd < nrElem_[3] ; dd++) cellLIDToLeafMapping_[dd] = -1;
	cellLeafToLIDMapping_.resize(nrElem_[3],-1);
	for (int dd = 0 ; dd < nrElem_[3] ; dd++) cellLeafToLIDMapping_[dd] = -1;

	SUNDANCE_MSG3(verb() , "HNMesh3D::createLeafNumbering , start assigning leaf numbers");

	// look for those leaf cells which points have a cell which maxCoFacet owner = myRank_
	// only those will have an LID
	for (int ind = 0 ; ind < nrElem_[3] ; ind++){
		Array<int>& vertexIDs = cellsPoints_[ind];
		hasCellLID[ind] = false;
		for (int v = 0 ; v < 8 ; v++){
			Array<int>& maxCoFacet = pointMaxCoF_[vertexIDs[v]];
			hasCellLID[ind] =  ( hasCellLID[ind]
					|| ( (maxCoFacet[0] >= 0) && (elementOwner_[3][maxCoFacet[0]] == myRank_) )
                    || ( (maxCoFacet[1] >= 0) && (elementOwner_[3][maxCoFacet[1]] == myRank_) )
                    || ( (maxCoFacet[2] >= 0) && (elementOwner_[3][maxCoFacet[2]] == myRank_) )
                    || ( (maxCoFacet[3] >= 0) && (elementOwner_[3][maxCoFacet[3]] == myRank_) )
                    || ( (maxCoFacet[4] >= 0) && (elementOwner_[3][maxCoFacet[4]] == myRank_) )
                    || ( (maxCoFacet[5] >= 0) && (elementOwner_[3][maxCoFacet[5]] == myRank_) )
                    || ( (maxCoFacet[6] >= 0) && (elementOwner_[3][maxCoFacet[6]] == myRank_) )
                    || ( (maxCoFacet[7] >= 0) && (elementOwner_[3][maxCoFacet[7]] == myRank_) )) ;

			// add cells with hanging nodes which have contribution to element which are owned by this processor
			// if vertex is hanging look into the parent cell at the same index and if the owner is myRank_ then add
			// to the cells which should be processed
			if ( (hasCellLID[ind] == false) && (isPointHanging_[vertexIDs[v]] == true)){
				int parentID = parentCellLID_[ind];
				Array<int>& parentVertexIDs = cellsPoints_[parentID];
				hasCellLID[ind] = hasCellLID[ind] || (elementOwner_[0][parentVertexIDs[v]] == myRank_);
			}
		}
		SUNDANCE_MSG3(verb() , "HNMesh3D::createLeafNumbering Cell ID :" << ind << " should be LID: " << hasCellLID[ind] <<
				" ,isCellLeaf_[ind]:" << isCellLeaf_[ind]);
	}

	//  treat special case, so that each hanging element has its parents
	// if we add one cell check hanging face, then add the maxCoF from the parent face if is leaf
	// if this is not successful then do the same thing for edges
	// - from each hanging edge there should be at least one cell on this processor which contains that parent edge !
	bool check_Ghost_cells = true;
	while (check_Ghost_cells){
		check_Ghost_cells = false;
	    for (int ind = 0 ; ind < nrElem_[3] ; ind++){
		   if ( (hasCellLID[ind] == true) && (elementOwner_[3][ind] != myRank_ ) ){
			  bool lookforEdges = true;
			  // check faces
			  Array<int>& faceIDs = cellsFaces_[ind];
			  for (int ii = 0 ; ii < 6 ; ii++ ){
				  // if the face is hanging and does not belong to me
				  if (isFaceHanging_[faceIDs[ii]] && ( elementOwner_[2][faceIDs[ii]] != myRank_)){
                    // get parent cells same face
					int parentCell = parentCellLID_[ind];
					Array<int>& parentfaceIDs = cellsFaces_[parentCell];
					for (int f = 0 ; f < 2 ; f++)
					if ( ( faceMaxCoF_[parentfaceIDs[ii]][f] >= 0 ) &&
						 ( elementOwner_[3][ faceMaxCoF_[parentfaceIDs[ii]][f] ] != myRank_ ) &&
						 ( hasCellLID[faceMaxCoF_[parentfaceIDs[ii]][f]] == false)  &&
						 ( isCellLeaf_[faceMaxCoF_[parentfaceIDs[ii]][f]] ) ){
						hasCellLID[faceMaxCoF_[parentfaceIDs[ii]][f]] = true;
						check_Ghost_cells = true;
						lookforEdges = false;
					}
				  }
			  }
			  // check edges
			  Array<int>& edgeIDs = cellsEdges_[ind];
			  // we have this if only
			  if (lookforEdges){
			    for (int ii = 0 ; ii < 12 ; ii++ ){
				  // if the face is hanging and does not belong to me
				  if (isEdgeHanging_[edgeIDs[ii]] && ( elementOwner_[1][edgeIDs[ii]] != myRank_)){
                    // get parent cells same face
					int parentCell = parentCellLID_[ind];
					Array<int>& parentEdgesIDs = cellsEdges_[parentCell];
					for (int f = 0 ; f < 4 ; f++)
					if ( ( edgeMaxCoF_[parentEdgesIDs[ii]][f] >= 0 ) &&
						 ( elementOwner_[3][ edgeMaxCoF_[parentEdgesIDs[ii]][f] ] != myRank_ ) &&
						 ( hasCellLID[edgeMaxCoF_[parentEdgesIDs[ii]][f]] == false)   &&
						 ( isCellLeaf_[edgeMaxCoF_[parentEdgesIDs[ii]][f]] )
					   ){
						hasCellLID[edgeMaxCoF_[parentEdgesIDs[ii]][f]] = true;
						check_Ghost_cells = true;
					}
				  }
			    } // from loop
			  }// from lookforEdges IF
		   }
	    }
	}

	// we also have to list the cells which are not owned by the processor
	for (int ind = 0 ; ind < nrElem_[3] ; ind++)
	{
		 // --------- GID numbering -----------
		 // if cell is leaf and if is inside the computational domain
         if ( (isCellLeaf_[ind] == true) && (!isCellOut_[ind]) )
         {
        	 Array<int>& vertexIDs = cellsPoints_[ind];
           	 for (int v = 0; v < 8 ; v++)
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
        	 for (int e = 0; e < 12 ; e++)
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
        	 Array<int>& faceIDs = cellsFaces_[ind];
        	 // for each face check weather it already has a leaf index, if not create one
        	 for (int f = 0; f < 6 ; f++)
        	 {
        		 //SUNDANCE_MSG3(verb() , " createLeafNumbering  edgeLIDs[e]:" << edgeLIDs[e] );
        		 if (faceGIDToLeafMapping_[faceIDs[f]] < 0)
        		 {
        			 SUNDANCE_MSG3(verb() , " createLeafGIDNumbering -> faceID:" << faceIDs[f] << " , nrFaceLeafGID_:" << nrFaceLeafGID_ );
        			 //SUNDANCE_MSG3(verb() , " MaxCoFacet:" << faceMaxCoF_[faceLIDs[e]] << " edgeVertex:" << faceVertex_[edgeLIDs[e]]);
        			 faceLeafToGIDMapping_[nrFaceLeafGID_] = faceIDs[f];
        			 faceGIDToLeafMapping_[faceIDs[f]] = nrFaceLeafGID_;
        			 nrFaceLeafGID_++;
        		 }
        	 }
        	 // create leaf index for the leaf cell
			 SUNDANCE_MSG3(verb() , " createLeafGIDNumbering CELL cellID:" << ind << " , nrCellLeafGID_:" << nrCellLeafGID_ );
        	 cellLeafToGIDMapping_[nrCellLeafGID_] = ind;
        	 cellGIDToLeafMapping_[ind] = nrCellLeafGID_;
        	 nrCellLeafGID_++;

        	 // --------- LID numbering -----------
        	 // create leaf LID numbering , if this cell needs to be processed
        	 if (hasCellLID[ind]){
        		 // vertex
              	 for (int v = 0; v < 8 ; v++)
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
            	 for (int e = 0; e < 12 ; e++)
            	 {
            		 if (edgeLIDToLeafMapping_[edgeIDs[e]] < 0)
            		 {
            			 SUNDANCE_MSG3(verb() , " createLeafLIDNumbering -> edgeID:" << edgeIDs[e] << " , nrEdgeLeafLID_:" << nrEdgeLeafLID_ );
            			 edgeLeafToLIDMapping_[nrEdgeLeafLID_] = edgeIDs[e];
            			 edgeLIDToLeafMapping_[edgeIDs[e]] = nrEdgeLeafLID_;
            			 nrEdgeLeafLID_++;
            		 }
            	 }
            	 // face LID
            	 for (int f = 0; f < 6 ; f++)
            	 {
            		 if (faceLIDToLeafMapping_[faceIDs[f]] < 0)
            		 {
            			 SUNDANCE_MSG3(verb() , " createLeafLIDNumbering -> faceID:" << faceIDs[f] << " , nrFaceLeafLID_:" << nrFaceLeafLID_ );
            			 faceLeafToLIDMapping_[nrFaceLeafLID_] = faceIDs[f];
            			 faceLIDToLeafMapping_[faceIDs[f]] = nrFaceLeafLID_;
            			 nrFaceLeafLID_++;
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
	SUNDANCE_MSG3(verb() , " nrVertexLeafGID_:" << nrVertexLeafGID_ << " nrEdgeLeafGID_:" << nrEdgeLeafGID_
			<< " nrFaceLeafGID_:" << nrFaceLeafGID_ << " nrCellLeafGID_:" << nrCellLeafGID_ );
	SUNDANCE_MSG3(verb() , " nrVertexLeafLID_:" << nrVertexLeafLID_ << " nrEdgeLeafLID_:" << nrEdgeLeafLID_
			<< " nrFaceLeafLID_:" <<nrFaceLeafLID_ << " nrCellLeafLID_:" << nrCellLeafLID_);
	SUNDANCE_MSG3(verb() , " vertexLIDToLeafMapping_: " << vertexLIDToLeafMapping_);
	SUNDANCE_MSG3(verb() , "HNMesh3D::createLeafNumbering , DONE");
}
