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
 * SundancePeanoMesh3D.cpp
 *
 *  Created on: Sep 8, 2009
 *      Author: benk
 */

#include "SundancePeanoMesh3D.hpp"

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

Point PeanoMesh3D::returnPoint(0.0 , 0.0 , 0.0);

PeanoMesh3D::PeanoMesh3D(int dim, const MPIComm& comm,
	    const MeshEntityOrder& order )
: MeshBase(dim, comm , order),_dimension(dim), _comm(comm)
 ,_peanoMesh(NULL)
{
	_uniqueResolution = 1.0;
}

void PeanoMesh3D::createMesh(
                      double position_x,
			          double position_y,
			          double position_z,
			          double offset_x,
			          double offset_y,
			          double offset_z,
			          double resolution
){
      double position[3];
      double offset[3];
      double res[3];

      // setting the values of the ctor argument
      position[0] = position_x; position[1] = position_y;  position[2] = position_z;
      offset[0] = offset_x;     offset[1] = offset_y;      offset[2] = offset_z;
      res[0] = resolution;      res[1] = resolution;       res[2] = resolution;


      // this is a 2D case
      // call the ctor for the Peano mesh
      SUNDANCE_VERB_LOW(" create Peano Mesh 3D ... ");
      _dimension = 3;
      // here we create the Peano grid
      _peanoMesh = new SundancePeanoInterface3D( position , offset , res );
      _uniqueResolution = _peanoMesh->returnResolution(0);

      SUNDANCE_VERB_LOW(" Peano Mesh created 3D ... \n");
      _peanoMesh->plotVTK("Peano3D");
      SUNDANCE_VERB_LOW(" After Plot 3D ... \n");
}

PeanoMesh3D::~PeanoMesh3D() {
	//delete _peanoMesh;
}


int PeanoMesh3D::numCells(int dim) const  {
	//printf("PeanoMesh3D::numCells(int dim):%d   dim:%d \n",_peanoMesh->numCells(dim),dim);
	return _peanoMesh->numCells(dim);
}

Point PeanoMesh3D::nodePosition(int i) const {
	//SUNDANCE_VERB_HIGH("nodePosition(int i)");
	//printf("PeanoMesh3D::nodePosition(int i)   i:%d \n", i);
	double* coords;
	coords = _peanoMesh->nodePositionView(i);
	// set the values
	PeanoMesh3D::returnPoint[0] = coords[0];
	PeanoMesh3D::returnPoint[1] = coords[1];
	PeanoMesh3D::returnPoint[2] = coords[2];
    return PeanoMesh3D::returnPoint;
}

const double* PeanoMesh3D::nodePositionView(int i) const {
	//printf("PeanoMesh3D::nodePositionView(int i)   i:%d \n", i);
	//SUNDANCE_VERB_HIGH("nodePosition(int i)");
	nodePosition(i);
	return &(PeanoMesh3D::returnPoint[0]);
}

void PeanoMesh3D::getJacobians(int cellDim, const Array<int>& cellLID,
                          CellJacobianBatch& jBatch) const
{
	  //printf("cellDim:%d  _uniqueResolution:%f ",cellDim, _uniqueResolution);
	  SUNDANCE_VERB_HIGH("getJacobians()");
	  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	  int nCells = cellLID.size();
 	  int tmp_index , tmp;
 	  int tmp_index1 , tmp_index2;
 	  Point pnt(0.0,0.0,0.0);
 	  Point pnt1(0.0,0.0,0.0);
 	  Point pnt2(0.0,0.0,0.0);
	  jBatch.resize(cellLID.size(), spatialDim(), cellDim);
	  if (cellDim < spatialDim()) // they need the Jacobian of a lower dinemsional element
	  {
		  //printf("PeanoMesh3D::getJacobians() cellDim < spatialDim() \n");
		   for (int i=0; i<nCells; i++)
		    {
		      //printf("PeanoMesh3D::getJacobian() cellDim < spatialDim() cellDim:%d , ret:%f \n",cellDim , _uniqueResolution);
		      double* detJ = jBatch.detJ(i);
		      switch(cellDim)
		      {
		        case 0: *detJ = 1.0;
		          break;
		        case 1:
			    	 tmp_index = this->facetLID(cellDim,  cellLID[i] , 0 , 0 , tmp );
			    	 tmp_index1= this->facetLID(cellDim,  cellLID[i] , 0 , 1 , tmp );
			    	 pnt = nodePosition(tmp_index);
			    	 pnt1 = nodePosition(tmp_index1);
			    	 pnt = pnt1 - pnt;
		             *detJ = sqrt(pnt * pnt); // the length of the edge
		        break;
		        case 2:{
			    	 tmp_index = this->facetLID(cellDim,  cellLID[i] , 0 , 0 , tmp );
			    	 tmp_index1= this->facetLID(cellDim,  cellLID[i] , 0 , 1 , tmp );
			    	 tmp_index2= this->facetLID(cellDim,  cellLID[i] , 0 , 2 , tmp );
			    	 pnt = nodePosition(tmp_index);
			    	 pnt1 = nodePosition(tmp_index1);
			    	 pnt2 = nodePosition(tmp_index2);
			    	 Point directedArea = cross( pnt1 - pnt , pnt2 - pnt );
		             *detJ = sqrt(directedArea * directedArea); // the are of the face
		        break;}
		        default:
		          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
		            "cellDim=" << cellDim << " in PeanoMesh3D::getJacobians()");
		      }
		    }
	  }else{ // they request the complete Jacoby matrix for this bunch of elements
		    //Array<double> J(cellDim*cellDim);
		    SUNDANCE_VERB_HIGH("cellDim == spatialDim()");
		    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
		    {
			  //printf("PeanoMesh3D::getJacobian() cellDim == spatialDim() cellDim:%d , ret:%f \n",cellDim , _uniqueResolution);
		      double* J = jBatch.jVals(i);
		      switch(cellDim)
		      {
		        case 3:
		          J[0] = _peanoMesh->returnResolution(0);
		          J[1] = 0.0; J[2] = 0.0; J[3] = 0.0;
		          J[4] = _peanoMesh->returnResolution(1);
		          J[5] = 0.0; J[6] = 0.0; J[7] = 0.0;
		          J[8] = _peanoMesh->returnResolution(2); // the Jacobi of the tet
		        break;
		        default:
		          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
		            "cellDim=" << cellDim
		            << " in PeanoMesh3D::getJacobians()");
		      }
		    }
	  }
}

void PeanoMesh3D::getCellDiameters(int cellDim, const Array<int>& cellLID,
                              Array<double>& cellDiameters) const {
	 TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	 SUNDANCE_VERB_HIGH("getCellDiameters()");
	  cellDiameters.resize(cellLID.size());

 	  int tmp_index , tmp;
 	  int tmp_index1;
 	  Point pnt(0.0,0.0,0.0);
 	  Point pnt1(0.0,0.0,0.0);

	  if (cellDim < spatialDim())
	  {
		//printf("PeanoMesh3D::getCellDiameters(), cellDim < spatialDim() \n ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
	        case 0:
	          cellDiameters[i] = 1.0;
	          break;
	        case 1:
		      tmp_index = this->facetLID(cellDim,  cellLID[i] , 0 , 0 , tmp );
		      tmp_index1= this->facetLID(cellDim,  cellLID[i] , 0 , 1 , tmp );
		      pnt = nodePosition(tmp_index);
		      pnt1 = nodePosition(tmp_index1);
		      pnt = pnt1 - pnt;
		      cellDiameters[i] = sqrt(pnt * pnt); // the length of the edge
	        break;
	        case 2:
		      tmp_index = this->facetLID(cellDim,  cellLID[i] , 0 , 0 , tmp );
		      tmp_index1= this->facetLID(cellDim,  cellLID[i] , 0 , 3 , tmp );
		      pnt = nodePosition(tmp_index);
		      pnt1 = nodePosition(tmp_index1);
		      pnt = pnt1 - pnt;
		      cellDiameters[i] = sqrt(pnt * pnt); // the length of the edge
	        break;
	        default:
	          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
	            "cellDim=" << cellDim << " in PeanoMesh3D::getCellDiameters()");
	      }
	    }
	  }
	  else
	  {
		//printf("PeanoMesh3D::getCellDiameters(), cellDim == spatialDim() \n ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
			case 3:
	          cellDiameters[i] = (_peanoMesh->returnResolution(0) + _peanoMesh->returnResolution(1)
	        		            + _peanoMesh->returnResolution(2))/3.0;
	        break;
	        default:
	          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
	            "cellDim=" << cellDim
	            << " in PeanoMesh3D::getCellDiameters()");
	      }
	    }
	  }
}

void PeanoMesh3D::pushForward(int cellDim, const Array<int>& cellLID,
                         const Array<Point>& refQuadPts,
                         Array<Point>& physQuadPts) const {

	  //printf("PeanoMesh3D::pushForward cellDim:%d\n",cellDim);
	  TEUCHOS_TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
	    "cellDim=" << cellDim
	    << " is not in expected range [0, " << spatialDim()
	    << "]");

	  int nQuad = refQuadPts.size();
	  double start_point[3] , end_point[3];
	  Array<double> J(cellDim*cellDim);

	  if (physQuadPts.size() > 0) physQuadPts.resize(0);
	  physQuadPts.reserve(cellLID.size() * refQuadPts.size());
	  for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	  {
	    int lid = cellLID[i];
	    _peanoMesh->pushForward( cellDim, lid ,start_point , end_point );
	    //std::cout << "_peanoMesh->pushForward: " << cellDim << "," << lid << std::endl;
	    //std::cout << "Start Point: "<< start_point[0] << " , " << start_point[1] << " , " << start_point[2] << std::endl;
	    //std::cout << "End Point: "<< end_point[0] <<" , " << end_point[1] << " , "<< end_point[2] << std::endl;
   	    Point pnt( start_point[0] , start_point[1] , start_point[2] );
   	    Point pnt1( end_point[0] , end_point[1] , end_point[2]);
	    switch(cellDim)
	    {
	      case 0: // integrate one point
	         physQuadPts.append(pnt);
	        break;
	      case 1:{ // integrate on one line
	         for (int q=0; q<nQuad; q++) {
	           physQuadPts.append(pnt + refQuadPts[q][0]*(pnt1 - pnt));
	         }
	      break;}
	      case 2:{
	         for (int q=0; q<nQuad; q++) {
	        	 if (fabs(pnt[0] - pnt1[0]) < 1e-8)
	          	      physQuadPts.append( pnt + Point(pnt[0] - pnt1[0],
	           	        		                      refQuadPts[q][0]*_peanoMesh->returnResolution(1),
	           	        		                      refQuadPts[q][1]*_peanoMesh->returnResolution(2)));
	        	 else
	        		 if (fabs(pnt[1] - pnt1[1]) < 1e-8)
		          	      physQuadPts.append( pnt + Point(refQuadPts[q][0]*_peanoMesh->returnResolution(0),
		          	    		                          pnt[1] - pnt1[1],
		           	        		                      refQuadPts[q][1]*_peanoMesh->returnResolution(2)));
	        		 else
		          	      physQuadPts.append( pnt + Point(refQuadPts[q][0]*_peanoMesh->returnResolution(0),
		          	    		                          refQuadPts[q][1]*_peanoMesh->returnResolution(1),
		          	    		                          pnt[2] - pnt1[2]));
	         }
	      break;}
	      case 3:{
	          for (int q=0; q<nQuad; q++) {
	          	      physQuadPts.append( pnt
	  	           	     + Point(refQuadPts[q][0]*_peanoMesh->returnResolution(0),
	  	           	        		refQuadPts[q][1]*_peanoMesh->returnResolution(1),
	  	           	        	    refQuadPts[q][2]*_peanoMesh->returnResolution(2)));
	          }
	      break;}
	      default:
	        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "impossible switch value "
	          "in PeanoMesh3D::getJacobians()");
	    }
	  }
}

int PeanoMesh3D::ownerProcID(int cellDim, int cellLID) const  {
	 SUNDANCE_VERB_HIGH("ownerProcID()"); return 0; }


int PeanoMesh3D::numFacets(int cellDim, int cellLID,
                      int facetDim) const  {
	SUNDANCE_VERB_HIGH("numFacets()");
    return _peanoMesh->numFacets(cellDim, cellLID, facetDim);
}


int PeanoMesh3D::facetLID(int cellDim, int cellLID,
                     int facetDim, int facetIndex,
                     int& facetOrientation) const  {
    int LID;
    LID = _peanoMesh->facetLID( cellDim,cellLID, facetDim, facetIndex, facetOrientation);
  	//printf("PeanoMesh3D::facetLID  cellDim: %d , cellLID: %d , facetDim %d , facetIndex:%d  %d\n" , cellDim , cellLID , facetDim , facetIndex , LID );
	return LID;
}

void PeanoMesh3D::getFacetLIDs(int cellDim,
                          const Array<int>& cellLID,
                          int facetDim,
                          Array<int>& facetLID,
                          Array<int>& facetSign) const {
	  SUNDANCE_VERB_HIGH("getFacetLIDs()");
	  //printf("PeanoMesh3D::getFacetLIDs()  cellDim:%d  cellLID.size():%d  facetDim:%d\n" , cellDim, (int)cellLID.size() , facetDim);
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

const int* PeanoMesh3D::elemZeroFacetView(int cellLID) const {
	return _peanoMesh->elemZeroFacetView(cellLID);
}

int PeanoMesh3D::numMaxCofacets(int cellDim, int cellLID) const  {
	  //SUNDANCE_VERB_HIGH("numMaxCofacets()");
      int coFacetCounter;
      coFacetCounter = _peanoMesh->numMaxCofacets( cellDim, cellLID);
	  //printf("numMaxCofacets:  cellDim:%d cellLID:%d ret:%d\n",cellDim, cellLID, coFacetCounter);
	  return coFacetCounter;
}

int PeanoMesh3D::maxCofacetLID(int cellDim, int cellLID,
                       int cofacetIndex,
                       int& facetIndex) const  {
	  int rtn;
	  rtn = _peanoMesh->maxCofacetLID(cellDim, cellLID, cofacetIndex, facetIndex);
	  //printf("maxCofacetLID() cellDim:%d,  cellLID:%d, cofacetIndex:%d , rtn:%d , facetIndex:%d\n",
		//	  cellDim,  cellLID, cofacetIndex , rtn , facetIndex);
	  return rtn;
}

void PeanoMesh3D::getMaxCofacetLIDs(const Array<int>& cellLIDs,
  MaximalCofacetBatch& cofacets) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error," PeanoMesh3D::getMaxCofacetLIDs() not implemented yet");
	//TODO: Implement this, uses only in ExodusWriter::writeMesh
}


void PeanoMesh3D::getCofacets(int cellDim, int cellLID,
                 int cofacetDim, Array<int>& cofacetLIDs) const {
    int tmpVect[12] , nrCofacets;
    _peanoMesh->getCofacets( cellDim, cellLID, cofacetDim, &tmpVect[0], nrCofacets);
    cofacetLIDs.resize(nrCofacets);
    for (int ii = 0 ; ii < nrCofacets ; ii++ ) cofacetLIDs[ii] = tmpVect[ii];
}


int PeanoMesh3D::mapGIDToLID(int cellDim, int globalIndex) const  {
	SUNDANCE_VERB_HIGH("mapGIDToLID()");
	// in the serial implementation GID = LID
	// in the parallel version this should be done differently
	return globalIndex;
}

bool PeanoMesh3D::hasGID(int cellDim, int globalIndex) const {
	SUNDANCE_VERB_HIGH("hasGID()");
	// since currently we have a serial implementation , this is always true
	// in the parallel version this function has to be implemented differetly
	return true;
}

int PeanoMesh3D::mapLIDToGID(int cellDim, int localIndex) const  {
	SUNDANCE_VERB_HIGH("mapLIDToGID()");
	// at the current stage we have only serial implementation,
	// parallel implementation will(should) come soon
	return localIndex;
}

CellType PeanoMesh3D::cellType(int cellDim) const  {
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

int PeanoMesh3D::label(int cellDim, int cellLID) const {
   return _peanoMesh->label( cellDim, cellLID);
}

void PeanoMesh3D::getLabels(int cellDim, const Array<int>& cellLID,
		Array<int>& labels) const {
    int tmpIndex;
    SUNDANCE_VERB_HIGH("getLabels()");
    // resize the array
	labels.resize(cellLID.size());

    for (tmpIndex = 0 ; tmpIndex < (int)cellLID.size() ; tmpIndex++){
    	labels[tmpIndex] = _peanoMesh->label( cellDim, cellLID[tmpIndex]);
    }
}

Set<int> PeanoMesh3D::getAllLabelsForDimension(int cellDim) const {
	  Set<int>                 rtn;
	  int                      tmpIndex;
	  SUNDANCE_VERB_HIGH("getAllLabelsForDimension()");

	  for (tmpIndex = 0 ; tmpIndex < _peanoMesh->numCells(cellDim) ; tmpIndex++){
		  rtn.put( _peanoMesh->label( cellDim, tmpIndex) );
	  }
	  return rtn;
}

void PeanoMesh3D::getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const {
    int                      tmpIndex , tmpLabel;
	SUNDANCE_VERB_HIGH("getLIDsForLabel()");
    for (tmpIndex = 0 ; tmpIndex < _peanoMesh->numCells(cellDim) ; tmpIndex++){
    	tmpLabel = this->label( cellDim , tmpIndex);
    	if (tmpLabel == label) cellLIDs.append( tmpIndex );
    }
}

void PeanoMesh3D::setLabel(int cellDim, int cellLID, int label) {
	_peanoMesh->setLabel(cellDim, cellLID, label);
}


void PeanoMesh3D::assignIntermediateCellGIDs(int cellDim) {
	SUNDANCE_VERB_HIGH("assignIntermediateCellGIDs()");
	// in this method we could do synchronization between processors, not usede now
}


bool PeanoMesh3D::hasIntermediateGIDs(int dim) const {
	SUNDANCE_VERB_HIGH("hasIntermediateGIDs()");
	return true; // true means they have been synchronized ... not used now
}

#endif
