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
 * SundanceCellCurvePredicate.cpp
 *
 *  Created on: Feb 19, 2010
 *      Author: benk
 */

#include "SundanceCellCurvePredicate.hpp"
#include "SundanceStdMathOps.hpp"

using namespace Sundance;

#define SIGN(X) ((X>0.0)?1:-1)

bool CellCurvePredicate::lessThan(const CellPredicateBase* other) const
{
  const CellCurvePredicate* S = dynamic_cast<const CellCurvePredicate*>(other);

  TEUCHOS_TEST_FOR_EXCEPTION( S== 0,
                     std::logic_error,
                     "argument " << other->toXML()
                     << " to CellCurvePredicate::lessThan() should be "
                     "a CellCurvePredicate pointer.");
  // This comparison is IMPORTANT, to determine RQCs
  return OrderedPair<ParametrizedCurve, int>(curve_, (int)filterMode_)
      < OrderedPair<ParametrizedCurve, int>(S->curve_, (int)S->filterMode_);
}

void CellCurvePredicate::testBatch(const Array<int>& cellLID,
                                        Array<int>& results) const
{
  results.resize(cellLID.size());
  double eps = 1e-14;

  if (cellDim()==0)
    {
	  switch (filterMode_){
	  case Outside_Curve:
	      for (int i=0; i<cellLID.size(); i++)
	    	  if ( curve_.curveEquation( mesh().nodePosition(cellLID[i])) >= 0.0 )
	              results[i] = true;
	    	  else
	    		  results[i] = false;
		  break;
	  case Inside_Curve:
	      for (int i=0; i<cellLID.size(); i++)
	    	  if ( curve_.curveEquation( mesh().nodePosition(cellLID[i])) <= 0.0 )
	              results[i] = true;
	    	  else
	    		  results[i] = false;
		  break;
	  case On_Curve:
	      for (int i=0; i<cellLID.size(); i++)
	    	  if ( fabs(curve_.curveEquation( mesh().nodePosition(cellLID[i]))) < eps )
	              results[i] = true;
	    	  else
	    		  results[i] = false;
		  break;
	  }
    }
  else
    {

	  switch (mesh().spatialDim()) {
	  case 2:{
		  // for 2D we test the intersection , in order to determine if the cell is cut or not
    	  // check for intersection of the edges, rework this
    	  // - we must have al least two DIFFERENT intersection point in order to have a clear intersection
	      Array<int> facetLIDs;
	      Array<int> facetSigns;
	      Array<double> equSum( cellLID.size() , 0.0 );
	      int nf = mesh().numFacets(cellDim(), cellLID[0], 0);
	      int spaceDim = mesh().spatialDim() ;
	      mesh().getFacetLIDs(cellDim(), cellLID, 0, facetLIDs, facetSigns);
	      for (int c=0; c<cellLID.size(); c++)
	        {
	          results[c] = true;
	          if (filterMode_ == On_Curve) results[c] = false;
	          int curve_sign = 0;
	          for (int f=0; f<nf; f++)
	            {
	        	  int fLID = facetLIDs[c*nf + f];
	        	  double curveEqu = curve_.curveEquation( mesh().nodePosition(fLID) );
	        	  //SUNDANCE_MSG3( 4 , " node: " << mesh().nodePosition(fLID) << " curveEq:" << curveEqu);
	        	  switch (filterMode_){
	        	  case Outside_Curve:{
	        		  if ( curveEqu <= 0.0 ) {
	        		     results[c] = false;
	        		  }
	        		  equSum[c] = equSum[c] + curveEqu;
	        		  break;}
	        	  case Inside_Curve: {
	        		  if ( curveEqu >= 0.0 ) {
	        		     results[c] = false;
	        		  }
	        		  equSum[c] = equSum[c] + curveEqu;
	        		  break; }
	        	  case On_Curve:
	        		  if (f == 0){
	        			  curve_sign = SIGN(curveEqu);
	        			  equSum[c] = equSum[c] + curveEqu;
	        		  } else {
	        			  if ( curve_sign != SIGN(curveEqu) ){
	             		     results[c] = true;
	        			  }
	             		  equSum[c] = equSum[c] + curveEqu;
	        		  }
	        		  break;
	        	  }
	            } // from for loop over points
	        } // loop over cells

	        if ( (spaceDim == 2) && (cellDim() == 2) )
	        {
	        	int nrEdge = mesh().numFacets(cellDim(), cellLID[0], 1 ) , nrPoints = 0 ;
	        	Array<int> edgeLIDs;
	        	Array<int> edgeLID(1);
	        	Array<int> cLID(1);
	        	Array<int> edgeSigns;
	        	Array<int> pointsOfEdgeLIDs;
	        	Array<bool> hasIntersectionPoints( cellLID.size() , false );
	        	Array<Point> intPoints;
	        	// for each cell look get all edge and test for intersection
	        	for (int c=0; c<cellLID.size(); c++)
	        	{
	        		Point p0 , start, end;
	        		int nrIntPoint = 0;
	        		cLID[0] = cellLID[c];
		        	mesh().getFacetLIDs(cellDim(), cLID, 1 , edgeLIDs, edgeSigns);
		        	// for each edge test intersection
	        		for (int edgeI = 0 ; edgeI < nrEdge ; edgeI++ ){
	        			edgeLID[0] = edgeLIDs[edgeI];
	        			mesh().getFacetLIDs(1, edgeLID, 0 , pointsOfEdgeLIDs, edgeSigns);
	        			start = mesh().nodePosition(pointsOfEdgeLIDs[0]);
	        			end = mesh().nodePosition(pointsOfEdgeLIDs[1]);
	        			// once we have the start and end point of the edge then test for intersection points
	        			curve_.returnIntersectPoints( start , end , nrPoints , intPoints);
	        			for (int p = 0 ;p < nrPoints ; p++){
	        				// test if we have more than two intersection points
	        				if (nrIntPoint == 0){
	        					p0 = intPoints[p]; nrIntPoint++;
	        				} else {
	        					double dist = ::sqrt( (p0 - intPoints[p])*(p0 - intPoints[p]) );
	        					if (dist > eps){
	        						// here we found two different intersection points
	        						hasIntersectionPoints[c] = true;
	        						continue;
	        					}
	        				}
	        			}
	        			// if we already found intersection points then stop
	        			if (hasIntersectionPoints[c]) continue;
	        		}
	        	} // for loop over cells

	        	// for each cell cell overwrite the result if there are intersection points
	        	int count = 0;
	        	switch (filterMode_){
	        	case Outside_Curve:{
	        	  for (int c=0; c<cellLID.size(); c++){
	        		  //SUNDANCE_MSG3( 4 , "filterMode_ :" << filterMode_ << " , BEFORE result:" << results[c] );
	        		  if (hasIntersectionPoints[c]) results[c] = false;
	        		  // if no intersection point then the cell must be complete in or out, this we can test with the sum
	        		  else results[c] = (equSum[c] >= 0);
	        		  //SUNDANCE_MSG3( 4 , "filterMode_ :" << filterMode_ << " , result:" << results[c] << " , equSum[c]:" << equSum[c]
	        		  //                   << " , hasIntersectionPoints[c]:" << hasIntersectionPoints[c]);
	        	  } }
	        	  break;
	        	case Inside_Curve:{
	        	  for (int c=0; c<cellLID.size(); c++){
	        		  //SUNDANCE_MSG3( 4 , "filterMode_ :" << filterMode_ << " , BEFORE result:" << results[c] );
	        		  if (hasIntersectionPoints[c]) results[c] = false;
	        		  // if no intersection point then the cell must be complete in or out, this we can test with the sum
	        		  else results[c] = (equSum[c] <= 0);
	        		  //SUNDANCE_MSG3( 4 , "filterMode_ :" << filterMode_ << " , result:" << results[c] << " , equSum[c]:" << equSum[c]
	        		  //                   << " , hasIntersectionPoints[c]:" << hasIntersectionPoints[c]);
	        	  } }
	        	  break;
	        	case On_Curve:{
	        	  for (int c=0; c<cellLID.size(); c++){
	        		  //SUNDANCE_MSG3( 4 , "filterMode_ :" << filterMode_ << " , BEFORE result:" << results[c] );
	        		  if (hasIntersectionPoints[c]) results[c] = true;
	        		  else results[c] = false;
	        		  if (results[c]) {count++;}
	        		  //SUNDANCE_MSG3( 4 , "filterMode_ :" << filterMode_ << " cellLID["<<c<<"]=" << cellLID[c] <<
	        			//	  " , result:" << results[c] << " , equSum[c]:" << equSum[c]
	        		    //      << " , hasIntersectionPoints[c]:" << hasIntersectionPoints[c] << " cellNr=" << count);
	        	  } }
	        	  break;
	        	}
        	}
	  } break;

	  // for 3D or 1D use the
	  default:{
	      Array<int> facetLIDs;
	      Array<int> facetSigns;
	      int nf = mesh().numFacets(cellDim(), cellLID[0], 0);
	      mesh().getFacetLIDs(cellDim(), cellLID, 0, facetLIDs, facetSigns);
	      for (int c=0; c<cellLID.size(); c++)
	        {
	          results[c] = true;
	          if (filterMode_ == On_Curve) results[c] = false;
	          int curve_sign = 0;
	          for (int f=0; f<nf; f++)
	            {
	        	  int fLID = facetLIDs[c*nf + f];
	        	  switch (filterMode_){
	        	  case Outside_Curve:
	        		  if ( curve_.curveEquation( mesh().nodePosition(fLID) ) <= 0.0 ) {
	        		     results[c] = false;
	        		     continue;
	        		  }
	        		  break;
	        	  case Inside_Curve:
	        		  if ( curve_.curveEquation( mesh().nodePosition(fLID) ) >= 0.0 ) {
	        		     results[c] = false;
	        		     continue;
	        		  }
	        		  break;
	        	  case On_Curve:
	        		  if (f == 0){
	        			  curve_sign = SIGN(curve_.curveEquation( mesh().nodePosition(fLID)));
	        		  } else {
	        			  if ( curve_sign != SIGN(curve_.curveEquation( mesh().nodePosition(fLID))) ){
	             		     results[c] = true;
	             		     continue;
	        			  }
	        		  }
	        		  break;
	        	  }
	            } // from for loop
	          //SUNDANCE_MSG3( 4 , "filterMode_ :" << filterMode_ << " cellLID=" << cellLID[c] << " , result:" << results[c]);
	        } // loop over cells
	  } break;
	  }


    }
}

XMLObject CellCurvePredicate::toXML() const
{
  XMLObject rtn("SundanceCellCurvePredicate");
  return rtn;
}
