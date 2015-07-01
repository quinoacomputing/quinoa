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
 * SundancePeanoMesh2D.cpp
 *
 *  Created on: Sep 8, 2009
 *      Author: benk
 */

#include "SundancePeanoMesh2D.hpp"

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

#ifdef HAVE_SUNDANCE_PEANO

using namespace Sundance;
using namespace Teuchos;
using Playa::MPIComm;
using Playa::MPIContainerComm;

//#define printf(msg)
//#define SUNDANCE_VERB_HIGH(msg) printf(msg);printf("\n");

Point PeanoMesh2D::returnPoint(0.0 , 0.0);

PeanoMesh2D::PeanoMesh2D(int dim, const MPIComm& comm ,
	    const MeshEntityOrder& order)
: MeshBase(dim, comm , order),_dimension(dim), _comm(comm)
 ,_peanoMesh(NULL)
{
	_uniqueResolution = 1.0;
}

void PeanoMesh2D::createMesh(
                      double position_x,
			          double position_y,
			          double offset_x,
			          double offset_y,
			          double resolution
){
      double position[2];
      double offset[2];
      double res[2];

      // setting the values of the ctor argument
      position[0] = position_x; position[1] = position_y;
      offset[0] = offset_x;     offset[1] = offset_y;
      res[0] = resolution;      res[1] = resolution;


      // this is a 2D case
      // call the ctor for the Peano mesh
      SUNDANCE_VERB_LOW(" create Peano Mesh ... ");
      _dimension = 2;
      // here we create the Peano grid
      _peanoMesh = new SundancePeanoInterface2D( position , offset , res );
      _uniqueResolution = _peanoMesh->returnResolution(0);

      SUNDANCE_VERB_LOW(" Peano Mesh created ... \n");
      //_peanoMesh->plotVTK("Peano2D");
      SUNDANCE_VERB_LOW(" After Plot ... \n");
}


PeanoMesh2D::~PeanoMesh2D() {
	//delete _peanoMesh;
}


int PeanoMesh2D::numCells(int dim) const  {
	//printf("PeanoMesh2D::numCells(int dim):%d   dim:%d \n",_peanoMesh->numCells(dim),dim);
	return _peanoMesh->numCells(dim);
}

Point PeanoMesh2D::nodePosition(int i) const {
	//SUNDANCE_VERB_HIGH("nodePosition(int i)");
	//printf("PeanoMesh2D::nodePosition(int i)   i:%d \n", i);
	double* coords;
	coords = _peanoMesh->nodePositionView(i);
	// set the values
	PeanoMesh2D::returnPoint[0] = coords[0];
	PeanoMesh2D::returnPoint[1] = coords[1];
	return PeanoMesh2D::returnPoint;
}

const double* PeanoMesh2D::nodePositionView(int i) const {
	//printf("PeanoMesh2D::nodePositionView(int i)   i:%d \n", i);
	//SUNDANCE_VERB_HIGH("nodePosition(int i)");
	nodePosition(i);
	return &(PeanoMesh2D::returnPoint[0]);
}

void PeanoMesh2D::getJacobians(int cellDim, const Array<int>& cellLID,
                          CellJacobianBatch& jBatch) const
{
	  //printf("cellDim:%d  _uniqueResolution:%f ",cellDim, _uniqueResolution);
	  SUNDANCE_VERB_HIGH("getJacobians()");
	  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	  int nCells = cellLID.size();
	  int tmp , tmp_index , tmp_index1;
	  Point pnt(0.0,0.0);
	  Point pnt1(0.0,0.0);

	  jBatch.resize(cellLID.size(), spatialDim(), cellDim);
	  if (cellDim < spatialDim()) // they need the Jacobian of a lower dinemsional element
	  {
		  //printf("PeanoMesh2D::getJacobians() cellDim < spatialDim() \n");
		   for (int i=0; i<nCells; i++)
		    {
		     //printf("PeanoMesh2D::getJacobian() cellDim < spatialDim() cellDim:%d , ret:%f \n",cellDim , _uniqueResolution);
		      double* detJ = jBatch.detJ(i);
		      switch(cellDim)
		      {
		        case 0: *detJ = 1.0;
		          break;
		        case 1:
			    	 tmp_index  = this->facetLID(cellDim,  cellLID[i] , 0 , 0 , tmp );
			    	 tmp_index1 = this->facetLID(cellDim,  cellLID[i] , 0 , 1 , tmp );
			    	 pnt = nodePosition(tmp_index);
			    	 pnt1 = nodePosition(tmp_index1);
			    	 pnt = pnt1 - pnt;
		             *detJ = sqrt(pnt * pnt); // the length of the edge
		        break;
		        default:
		          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
		            "cellDim=" << cellDim << " in PeanoMesh2D::getJacobians()");
		      }
		    }
	  }else{ // they request the complete Jacoby matrix for this bunch of elements
		    //Array<double> J(cellDim*cellDim);
		    SUNDANCE_VERB_HIGH("cellDim == spatialDim()");
		    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
		    {
			  //printf("PeanoMesh2D::getJacobian() cellDim == spatialDim() cellDim:%d , ret:%f \n",cellDim , _uniqueResolution);
		      double* J = jBatch.jVals(i);
		      switch(cellDim)
		      {
		        case 2:
		          J[0] = _peanoMesh->returnResolution(0);
		          J[1] = 0.0;  J[2] = 0.0;   //
		          J[3] = _peanoMesh->returnResolution(1); // the Jacobi of the quad,
		        break;
		        default:
		          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
		            "cellDim=" << cellDim
		            << " in SundancePeano2D::getJacobians()");
		      }
		    }
	  }
}

void PeanoMesh2D::getCellDiameters(int cellDim, const Array<int>& cellLID,
                              Array<double>& cellDiameters) const {
	 TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	 SUNDANCE_VERB_HIGH("getCellDiameters()");
	  cellDiameters.resize(cellLID.size());

	  int tmp , tmp_index , tmp_index1;
	  Point pnt(0.0,0.0);
	  Point pnt1(0.0,0.0);

	  if (cellDim < spatialDim())
	  {
		//printf("PeanoMesh2D::getCellDiameters(), cellDim < spatialDim() \n ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
	        case 0:
	             cellDiameters[i] = 1.0;
	          break;
	        case 1:  //length of the edge
		    	 tmp_index = this->facetLID(cellDim,  cellLID[i] , 0 , 0 , tmp );
		    	 tmp_index1= this->facetLID(cellDim,  cellLID[i] , 0 , 1 , tmp );
		    	 pnt = nodePosition(tmp_index);
		    	 pnt1 = nodePosition(tmp_index1);
		    	 pnt = pnt1 - pnt;
		    	 cellDiameters[i] = sqrt(pnt * pnt); // the length of the edge
	        break;
	        default:
	          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
	            "cellDim=" << cellDim << " in PeanoMesh2D::getCellDiameters()");
	      }
	    }
	  }
	  else
	  {
		//printf("PeanoMesh2D::getCellDiameters(), cellDim == spatialDim() \n ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
	        case 2:
	          cellDiameters[i] = (_peanoMesh->returnResolution(0)+_peanoMesh->returnResolution(1))/2.0;
	        break;
	        default:
	          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
	            "cellDim=" << cellDim
	            << " in PeanoMesh2D::getCellDiameters()");
	      }
	    }
	  }
}

void PeanoMesh2D::pushForward(int cellDim, const Array<int>& cellLID,
                         const Array<Point>& refQuadPts,
                         Array<Point>& physQuadPts) const {

	  //printf("PeanoMesh2D::pushForward cellDim:%d\n",cellDim);
	  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim
	    << " is not in expected range [0, " << spatialDim()
	    << "]");

	  int nQuad = refQuadPts.size();
	  double start_point[2] , end_point[2];

	  Array<double> J(cellDim*cellDim);

	  if (physQuadPts.size() > 0) physQuadPts.resize(0);
	  physQuadPts.reserve(cellLID.size() * refQuadPts.size());
	  for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	  {
	    int lid = cellLID[i];
	    _peanoMesh->pushForward( cellDim, lid ,start_point , end_point );
   	    Point pnt( start_point[0] , start_point[1] );
   	    Point pnt1( end_point[0] , end_point[1] );
	    switch(cellDim)
	    {
	      case 0: // integrate one point
	         physQuadPts.append(pnt);
	        break;
	      case 1:{ // integrate on one line
	         for (int q=0; q<nQuad; q++) {
	               physQuadPts.append(pnt + (pnt1-pnt)*refQuadPts[q][0]);
	        	}
	      break;}
	      case 2:{
		         for (int q=0; q<nQuad; q++) {
		          	  physQuadPts.append( pnt
		           	    + Point(refQuadPts[q][0]*_peanoMesh->returnResolution(0),
		           	    		refQuadPts[q][1]*_peanoMesh->returnResolution(1)));
		         }
	      break;}
	      default:
	        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
	          "in PeanoMesh2D::getJacobians()");
	    }
	  }
}

int PeanoMesh2D::ownerProcID(int cellDim, int cellLID) const  {
	 SUNDANCE_VERB_HIGH("ownerProcID()"); return 0; }


int PeanoMesh2D::numFacets(int cellDim, int cellLID,
                      int facetDim) const  {
	SUNDANCE_VERB_HIGH("numFacets()");
    return _peanoMesh->numFacets(cellDim, cellLID, facetDim);
}


int PeanoMesh2D::facetLID(int cellDim, int cellLID,
                     int facetDim, int facetIndex,
                     int& facetOrientation) const  {
    int LID;
    LID = _peanoMesh->facetLID( cellDim,cellLID, facetDim, facetIndex, facetOrientation);
  	//printf("PeanoMesh2D::facetLID  cellDim: %d , cellLID: %d , facetDim %d , facetIndex:%d  %d\n" , cellDim , cellLID , facetDim , facetIndex , LID );
	return LID;
}

void PeanoMesh2D::getFacetLIDs(int cellDim,
                          const Array<int>& cellLID,
                          int facetDim,
                          Array<int>& facetLID,
                          Array<int>& facetSign) const {
	  SUNDANCE_VERB_HIGH("getFacetLIDs()");
	  //printf("PeanoMesh2D::getFacetLIDs()  cellDim:%d  cellLID.size():%d  facetDim:%d\n" , cellDim, (int)cellLID.size() , facetDim);
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
			  //printf("LID:%d , cellDim:%d , cLID:%d , facetDim:%d , f:%d , facetOrientation:%d \n"
			  //	  ,LID , cellDim, cLID, facetDim, f , facetOrientation );
	          facetLID[ptr] = LID;
	          facetSign[ptr] = facetOrientation;
	      }
	  }
}

const int* PeanoMesh2D::elemZeroFacetView(int cellLID) const {
	return _peanoMesh->elemZeroFacetView(cellLID);
}

int PeanoMesh2D::numMaxCofacets(int cellDim, int cellLID) const  {
	  //SUNDANCE_VERB_HIGH("numMaxCofacets()");
      int coFacetCounter;
      coFacetCounter = _peanoMesh->numMaxCofacets( cellDim, cellLID);
	  //printf("numMaxCofacets:  cellDim:%d cellLID:%d ret:%d\n",cellDim, cellLID, coFacetCounter);
	  return coFacetCounter;
}

int PeanoMesh2D::maxCofacetLID(int cellDim, int cellLID,
                       int cofacetIndex,
                       int& facetIndex) const  {
	  int rtn;
	  rtn = _peanoMesh->maxCofacetLID(cellDim, cellLID, cofacetIndex, facetIndex);
	  //printf("maxCofacetLID() cellDim:%d,  cellLID:%d, cofacetIndex:%d , rtn:%d , facetIndex:%d\n",
		//	  cellDim,  cellLID, cofacetIndex , rtn , facetIndex);
	  return rtn;
}

void PeanoMesh2D::getMaxCofacetLIDs(const Array<int>& cellLIDs,
  MaximalCofacetBatch& cofacets) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," PeanoMesh2D::getMaxCofacetLIDs() not implemented yet");
	//TODO: Implement this, uses only in ExodusWriter::writeMesh
}


void PeanoMesh2D::getCofacets(int cellDim, int cellLID,
                 int cofacetDim, Array<int>& cofacetLIDs) const {
    int tmpVect[12] , nrCofacets;
    _peanoMesh->getCofacets( cellDim, cellLID, cofacetDim, &tmpVect[0], nrCofacets);
    cofacetLIDs.resize(nrCofacets);
    for (int ii = 0 ; ii < nrCofacets ; ii++ ) cofacetLIDs[ii] = tmpVect[ii];
}


int PeanoMesh2D::mapGIDToLID(int cellDim, int globalIndex) const  {
	SUNDANCE_VERB_HIGH("mapGIDToLID()");
	// in the serial implementation GID = LID
	// in the parallel version this should be done differently
	return globalIndex;
}

bool PeanoMesh2D::hasGID(int cellDim, int globalIndex) const {
	SUNDANCE_VERB_HIGH("hasGID()");
	// since currently we have a serial implementation , this is always true
	// in the parallel version this function has to be implemented differetly
	return true;
}

int PeanoMesh2D::mapLIDToGID(int cellDim, int localIndex) const  {
	SUNDANCE_VERB_HIGH("mapLIDToGID()");
	// at the current stage we have only serial implementation,
	// parallel implementation will(should) come soon
	return localIndex;
}

CellType PeanoMesh2D::cellType(int cellDim) const  {
	//printf("cellType() cellDim:%d\n",cellDim);
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

int PeanoMesh2D::label(int cellDim, int cellLID) const {
   return _peanoMesh->label( cellDim, cellLID);
}

void PeanoMesh2D::getLabels(int cellDim, const Array<int>& cellLID,
		Array<int>& labels) const {
    int tmpIndex;
    SUNDANCE_VERB_HIGH("getLabels()");
    // resize the array
	labels.resize(cellLID.size());

    for (tmpIndex = 0 ; tmpIndex < (int)cellLID.size() ; tmpIndex++){
    	labels[tmpIndex] = _peanoMesh->label( cellDim, cellLID[tmpIndex]);
    }
}

Set<int> PeanoMesh2D::getAllLabelsForDimension(int cellDim) const {
	  Set<int>                 rtn;
	  int                      tmpIndex;
	  SUNDANCE_VERB_HIGH("getAllLabelsForDimension()");

	  for (tmpIndex = 0 ; tmpIndex < _peanoMesh->numCells(cellDim) ; tmpIndex++){
		  rtn.put( _peanoMesh->label( cellDim, tmpIndex) );
	  }
	  return rtn;
}

void PeanoMesh2D::getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const {
    int                      tmpIndex , tmpLabel;
	SUNDANCE_VERB_HIGH("getLIDsForLabel()");
    for (tmpIndex = 0 ; tmpIndex < _peanoMesh->numCells(cellDim) ; tmpIndex++){
    	tmpLabel = this->label( cellDim , tmpIndex);
    	if (tmpLabel == label) cellLIDs.append( tmpIndex );
    }
}

void PeanoMesh2D::setLabel(int cellDim, int cellLID, int label) {
	_peanoMesh->setLabel(cellDim, cellLID, label);
}


void PeanoMesh2D::assignIntermediateCellGIDs(int cellDim) {
	SUNDANCE_VERB_HIGH("assignIntermediateCellGIDs()");
	// in this method we could do synchronization between processors, not usede now
}


bool PeanoMesh2D::hasIntermediateGIDs(int dim) const {
	SUNDANCE_VERB_HIGH("hasIntermediateGIDs()");
	return true; // true means they have been synchronized ... not used now
}


#endif
