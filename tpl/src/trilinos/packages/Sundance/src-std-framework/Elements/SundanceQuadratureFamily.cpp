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


#include "SundanceQuadratureFamily.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Teuchos;
using Playa::Handle;
using Playa::Handleable;
using std::ios_base;
using std::setw;
using std::endl;
using std::setprecision;

int QuadratureFamily::getNumPoints( const CellType & cellType ) const
{
  const QuadratureFamilyBase* q 
    = dynamic_cast<const QuadratureFamilyBase*>(ptr().get());
  return q->getNumPoints( cellType );
}

void QuadratureFamily::getPoints(const CellType& cellType, 
                                 Array<Point>& quadPoints,
                                 Array<double>& quadWeights) const 
{
  const QuadratureFamilyBase* q 
    = dynamic_cast<const QuadratureFamilyBase*>(ptr().get());
  
  TEUCHOS_TEST_FOR_EXCEPTION(q==0, std::logic_error, 
                     "QuadratureFamilyStub pointer" << toXML().toString() 
                     << " could not be cast to a QuadratureFamilyBase ptr");
                     
  q->getPoints(cellType, quadPoints, quadWeights);
}

void QuadratureFamily::getAdaptedWeights(const CellType& cellType ,
		                         int cellDim,
	                             int celLID ,
	            	             int facetIndex ,
                                 const Mesh& mesh ,
                                 const ParametrizedCurve& globalCurve ,
                                 Array<Point>& quadPoints ,
                                 Array<double>& quadWeights ,
                                 bool& isCut) const
{

  const QuadratureFamilyBase* q
   = dynamic_cast<const QuadratureFamilyBase*>(ptr().get());

  TEUCHOS_TEST_FOR_EXCEPTION(q==0, std::logic_error,
                     "QuadratureFamilyStub pointer" << toXML().toString()
                      << " could not be cast to a QuadratureFamilyBase ptr");

  q->getAdaptedWeights(cellType, cellDim , celLID , facetIndex ,
  	  mesh , globalCurve , quadPoints, quadWeights , isCut);

}

XMLObject QuadratureFamily::toXML() const 
{
  return ptr()->toXML();
}

int QuadratureFamily::order() const 
{
  return ptr()->order();
}



void QuadratureFamily::getFacetPoints(const CellType& cellType, 
                                      int facetDim,
                                      int facetIndex,
                                      Array<Point>& quadPoints,
                                      Array<double>& quadWeights) const
{
  switch(cellType)
    {
    case LineCell:
      getLineFacetQuad(facetDim, facetIndex, quadPoints, quadWeights);
      break;
    case TriangleCell:
      getTriangleFacetQuad(facetDim, facetIndex, quadPoints, quadWeights);
      break;
    case QuadCell:
      getQuadFacetQuad(facetDim, facetIndex, quadPoints, quadWeights);
      break;
    case TetCell:
      getTetFacetQuad(facetDim, facetIndex, quadPoints, quadWeights);
      break;
    case BrickCell:
      getBrickFacetQuad(facetDim, facetIndex, quadPoints, quadWeights);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                         "getFacetPoints() not implemented for cell type "
                         << cellType);
    }
}


void QuadratureFamily::getLineFacetQuad(int facetDim,
                                        int facetIndex,
                                        Array<Point>& quadPoints,
                                        Array<double>& quadWeights) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(facetDim > 0, std::runtime_error,
                     "Invalid facet dimension " << facetDim 
                     << " in getLineFacetQuad()");
  TEUCHOS_TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 1, std::runtime_error,
                     "Invalid facet index " << facetIndex  
                     << " in getLineFacetQuad()");

  quadPoints.resize(1);
  quadWeights.resize(1);
  quadWeights[0] = 1.0;  
  quadPoints[0] = Point((double) facetIndex);
}




void QuadratureFamily::getTriangleFacetQuad(int facetDim,
                                            int facetIndex,
                                            Array<Point>& quadPoints,
                                            Array<double>& quadWeights) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(facetDim > 1, std::runtime_error,
                     "Invalid facet dimension " << facetDim 
                     << " in getTriangleFacetQuad()");
  TEUCHOS_TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 2, std::runtime_error,
                     "Invalid facet index " << facetIndex  
                     << " in getTriangleFacetQuad()");
  if (facetDim==1)
    {
      Array<Point> facetPts;
      Array<double> facetWts;
      getPoints(LineCell, facetPts, facetWts);
      quadPoints.resize(facetPts.size());
      quadWeights.resize(facetWts.size());
      for (int i=0; i<facetPts.size(); i++)
        {
          if (facetIndex==2)
            {
              quadPoints[i] = Point(facetPts[i][0], 0.0);
              quadWeights[i] = facetWts[i];
            }
          else if (facetIndex==0)
            {
              quadPoints[i] = Point(1.0-facetPts[i][0], facetPts[i][0]);
              quadWeights[i] = facetWts[i];
            }
          else
            {
              quadPoints[i] = Point(0.0, facetPts[i][0]);
              quadWeights[i] = facetWts[i];
            }
        }
    }
  else
    {
      quadPoints.resize(1);
      quadWeights.resize(1);
      quadWeights[0] = 1.0;  
      if (facetIndex==0) quadPoints[0] = Point(0.0, 0.0);
      else if (facetIndex==1) quadPoints[0] = Point(1.0, 0.0);
      else quadPoints[0] = Point(0.0, 1.0);
    }
}

void QuadratureFamily::getQuadFacetQuad(int facetDim,
                           int facetIndex,
                           Array<Point>& quadPoints,
                           Array<double>& quadWeights) const {

	  TEUCHOS_TEST_FOR_EXCEPTION(facetDim > 1, std::runtime_error,
	                     "Invalid facet dimension " << facetDim
	                     << " in getQuadFacetQuad()");

	  TEUCHOS_TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 3, std::runtime_error,
	                     "Invalid facet index " << facetIndex
	                     << " in getQuadFacetQuad()");
	  if (facetDim==1)
	    {
	      Array<Point> facetPts;
	      Array<double> facetWts;
	      getPoints(LineCell, facetPts, facetWts);
	      quadPoints.resize(facetPts.size());
	      quadWeights.resize(facetWts.size());
	      for (int i=0; i<facetPts.size(); i++)
	        {
	          if (facetIndex==0) // the 4 edges of the Quad
	            {               // Numbering is important, if Edge Numbering in QuadCell changes, change this as well
	              quadPoints[i] = Point(facetPts[i][0], 0.0);
	              quadWeights[i] = facetWts[i];
	            }
	          else if (facetIndex==1)
	            {
	              quadPoints[i] = Point(0.0 , facetPts[i][0]);
	              quadWeights[i] = facetWts[i];
	            }
	          else if (facetIndex==2)
	          	{
	              quadPoints[i] = Point(1.0, facetPts[i][0]);
	              quadWeights[i] = facetWts[i];
	          	}
	          else
	            {
	          	  quadPoints[i] = Point(facetPts[i][0], 1.0);
	          	  quadWeights[i] = facetWts[i];
	            }
	        }
	    }
	  else
	    {
	      quadPoints.resize(1);
	      quadWeights.resize(1);
	      quadWeights[0] = 1.0;
	      if (facetIndex==0) quadPoints[0] = Point(0.0, 0.0);
	      else if (facetIndex==1) quadPoints[0] = Point(1.0, 0.0);
	      else if (facetIndex==2) quadPoints[0] = Point(0.0, 1.0);
	      else quadPoints[0] = Point(1.0, 1.0);
	    }
}

void QuadratureFamily::getTetFacetQuad(int facetDim,
                                       int facetIndex,
                                       Array<Point>& quadPoints,
                                       Array<double>& quadWeights) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(facetDim > 2, std::runtime_error,
                     "Invalid facet dimension " << facetDim 
                     << " in getTetFacetQuad()");
  TEUCHOS_TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 4, std::runtime_error,
                     "Invalid facet index " << facetIndex  
                     << " in getTetFacetQuad()");
  if (facetDim==2)
    {
      Array<Point> facetPts;
      Array<double> facetWts;
      getPoints(TriangleCell, facetPts, facetWts);
      quadPoints.resize(facetPts.size());
      quadWeights.resize(facetWts.size());
      for (int i=0; i<facetPts.size(); i++)
        {
          double s = facetPts[i][0];
          double t = facetPts[i][1];
          double x,y,z;
          if (facetIndex==2) 
            {
              x = s;
              y = 0.0;
              z = t;
            }
          else if (facetIndex==0) 
            {
              x = 1.0-s-t;
              y = s;
              z = t;
            }
          else if (facetIndex==1) 
            {
              x = 0.0;
              y = t;
              z = s;
            }
          else
            {
              x = s;
              y = 1.0-s-t;
              z = 0.0;
            }
          quadPoints[i] = Point(x, y, z);
	  if (facetIndex != 1)
	    quadWeights[i] = facetWts[i];
	  else 
	    quadWeights[i] = facetWts[i] * sqrt(3.0);

        }
    }
  else if (facetDim==0)
    {
      quadPoints.resize(1);
      quadWeights.resize(1);
      quadWeights[0] = 1.0;  
      if (facetIndex==0) quadPoints[0] = Point(0.0, 0.0, 0.0);
      else if (facetIndex==1) quadPoints[0] = Point(1.0, 0.0, 0.0);
      else if (facetIndex==2) quadPoints[0] = Point(0.0, 1.0, 0.0);
      else quadPoints[0] = Point(0.0, 0.0, 1.0);
    }
  TEUCHOS_TEST_FOR_EXCEPT(facetDim==1);
}


void QuadratureFamily::getBrickFacetQuad(int facetDim,
                                       int facetIndex,
                                       Array<Point>& quadPoints,
                                       Array<double>& quadWeights) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(facetDim > 2, std::runtime_error,
                     "Invalid facet dimension " << facetDim
                     << " in getBrickFacetQuad()");
  if (facetDim==2)
    {
	  TEUCHOS_TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 5, std::runtime_error,
	                      "Invalid facet index " << facetIndex
	                      << " in getBrickFacetQuad()");
      Array<Point> facetPts;
      Array<double> facetWts;
      getPoints(QuadCell, facetPts, facetWts);
      quadPoints.resize(facetPts.size());
      quadWeights.resize(facetWts.size());
      for (int i=0; i<facetPts.size(); i++)
        {
          if (facetIndex==0) // the 6 faces of the Brick
            {// Numbering is important, if Edge Numbering in BrickCell changes, change this as well
              quadPoints[i] = Point(facetPts[i][0], facetPts[i][1] , 0.0);
              quadWeights[i] = facetWts[i];
            }
          else if (facetIndex==1)
            {
              quadPoints[i] = Point(facetPts[i][0], 0.0 , facetPts[i][1] );
              quadWeights[i] = facetWts[i];
            }
          else if (facetIndex==2)
          	{
              quadPoints[i] = Point( 0.0 , facetPts[i][0] , facetPts[i][1] );
              quadWeights[i] = facetWts[i];
          	}
          else if (facetIndex==3)
          	{
              quadPoints[i] = Point( 1.0 , facetPts[i][0] , facetPts[i][1] );
              quadWeights[i] = facetWts[i];
          	}
          else if (facetIndex==4)
          	{
              quadPoints[i] = Point( facetPts[i][0] , 1.0 , facetPts[i][1] );
              quadWeights[i] = facetWts[i];
          	}
          else
            {
          	  quadPoints[i] = Point( facetPts[i][0] , facetPts[i][1] , 1.0 );
          	  quadWeights[i] = facetWts[i];
            }

        }
    }
  else if (facetDim==1)
     {
	  TEUCHOS_TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 11, std::runtime_error,
	                      "Invalid facet index " << facetIndex
	                      << " in getBrickFacetQuad()");
      Array<Point> facetPts;
      Array<double> facetWts;
      getPoints(LineCell, facetPts, facetWts);
      quadPoints.resize(facetPts.size());
      quadWeights.resize(facetWts.size());
      for (int i=0; i<facetPts.size(); i++)
        {
          if (facetIndex==0) // the 6 faces of the Brick
            {// Numbering is important, if Edge Numbering in BrickCell changes, change this as well
              quadPoints[i] = Point(facetPts[i][0], 0.0 , 0.0);  quadWeights[i] = facetWts[i];
            }
          else if (facetIndex==1)
          { quadPoints[i] = Point( 0.0 , facetPts[i][0], 0.0 ); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==2)
          { quadPoints[i] = Point( 0.0 , 0.0 , facetPts[i][0]); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==3)
          { quadPoints[i] = Point( 1.0 , facetPts[i][0] , 0.0 ); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==4)
          { quadPoints[i] = Point( 1.0 , 0.0 , facetPts[i][0] ); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==5)
          { quadPoints[i] = Point( facetPts[i][0] , 1.0 , 0.0 ); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==6)
          { quadPoints[i] = Point( 0.0 , 1.0 , facetPts[i][0]); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==7)
          { quadPoints[i] = Point( 1.0 , 1.0 , facetPts[i][0]); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==8)
          { quadPoints[i] = Point( facetPts[i][0] , 0.0 , 1.0); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==9)
          { quadPoints[i] = Point( 0.0 , facetPts[i][0] , 1.0); quadWeights[i] = facetWts[i]; }
          else if (facetIndex==10)
          { quadPoints[i] = Point( 1.0 , facetPts[i][0] , 1.0); quadWeights[i] = facetWts[i]; }
          else
          { quadPoints[i] = Point( facetPts[i][0] , 1.0 , 1.0); quadWeights[i] = facetWts[i]; }
        }
     }
  else if (facetDim==0)
    {
      quadPoints.resize(1);
      quadWeights.resize(1);
      quadWeights[0] = 1.0;
      if (facetIndex==0) quadPoints[0] = Point(0.0, 0.0 , 0.0);
      else if (facetIndex==1) quadPoints[0] = Point(1.0, 0.0, 0.0);
      else if (facetIndex==2) quadPoints[0] = Point(0.0, 1.0, 0.0);
      else if (facetIndex==3) quadPoints[0] = Point(1.0, 1.0, 0.0);
      else if (facetIndex==4) quadPoints[0] = Point(0.0, 0.0, 1.0);
      else if (facetIndex==5) quadPoints[0] = Point(1.0, 0.0, 1.0);
      else if (facetIndex==6) quadPoints[0] = Point(0.0, 1.0, 1.0);
      else quadPoints[0] = Point(1.0, 1.0, 1.0);
    }
}


namespace Sundance
{
void printQuad(std::ostream& os, 
  const Array<Point>& pts, const Array<double>& wgts)
{
  Tabs tab(0);
  static Array<string> names = tuple<string>("x", "y", "z");
  TEUCHOS_TEST_FOR_EXCEPT(pts.size() != wgts.size());
  
  TEUCHOS_TEST_FOR_EXCEPT(pts.size() < 1);

  int dim = pts[0].dim();

  int prec=10;
  int w = prec+5;
  ios_base::fmtflags oldFlags = os.flags();
  os.setf(ios_base::internal);    
  os.setf(ios_base::showpoint);

  os << tab << setw(5) << "i" << setw(w) << "w";;
  
  for (int i=0; i<pts[0].dim(); i++)
  {
    os << setw(w) << names[i] ;
  }
  os << std::endl;
  os.setf(ios_base::right);    
  for (int i=0; i<pts.size(); i++)
  {
    os << tab << setw(5) << i << setw(w) << setprecision(prec) << wgts[i];
    for (int d=0; d<dim; d++)
    {
      os << setw(w) << setprecision(prec) << pts[i][d];
    }
    os << std::endl;
  }
  os.flags(oldFlags);
}
}
