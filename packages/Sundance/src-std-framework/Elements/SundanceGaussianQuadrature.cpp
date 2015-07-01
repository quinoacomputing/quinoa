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

#include "SundanceGaussianQuadrature.hpp"
#include "SundanceGauss1D.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceQuadQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"
#include "SundanceBrickQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;

GaussianQuadrature::GaussianQuadrature(int order)
  : QuadratureFamilyBase(order)
{
  
}

XMLObject GaussianQuadrature::toXML() const 
{
  XMLObject rtn("GaussianQuadrature");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}



void GaussianQuadrature::getLineRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const 
{
  int p = order() + 1;
  p = p + (p%2);
  int n = p/2;
			
  quadPoints.resize(n);
  quadWeights.resize(n);
			
  Gauss1D q1(n, 0.0, 1.0);
			
  for (int i=0; i<n; i++)
    {
      quadWeights[i] = q1.weights()[i];
      quadPoints[i] = Point(q1.nodes()[i]);
    }
}

void GaussianQuadrature::getTriangleRule(Array<Point>& quadPoints,
                                          Array<double>& quadWeights) const 
{
  Array<double> x;
  Array<double> y;
  Array<double> w;
			
  TriangleQuadrature::getPoints(order(), w, x, y);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = 0.5*w[i];
      quadPoints[i] = Point(x[i], y[i]);
    }  
}

void GaussianQuadrature::getQuadRule(Array<Point>& quadPoints,
                                          Array<double>& quadWeights) const
{
  Array<double> x;
  Array<double> y;
  Array<double> w;

  QuadQuadrature::getPoints(order(), w, x, y);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = w[i];
      quadPoints[i] = Point(x[i], y[i]);
    }
}

void GaussianQuadrature::getTetRule(Array<Point>& quadPoints,
                                    Array<double>& quadWeights) const 
{
  Array<double> x;
  Array<double> y;
  Array<double> z;
  Array<double> w;
			
  TetQuadrature::getPoints(order(), w, x, y, z);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = w[i]/6.0;
      quadPoints[i] = Point(x[i], y[i], z[i]);
    }  
}

void GaussianQuadrature::getBrickRule(Array<Point>& quadPoints,
                                    Array<double>& quadWeights) const
{
  Array<double> x;
  Array<double> y;
  Array<double> z;
  Array<double> w;

  BrickQuadrature::getPoints(order(), w, x, y, z);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = w[i];
      quadPoints[i] = Point(x[i], y[i], z[i]);
    }
}

void GaussianQuadrature::getAdaptedWeights(const CellType& cellType ,
									       int cellDim,
	                                       int celLID ,
	                	                   int facetIndex ,
	                                       const Mesh& mesh ,
	                                       const ParametrizedCurve& globalCurve ,
	                                       Array<Point>& quadPoints ,
	                                       Array<double>& quadWeights ,
	                                       bool &isCut) const {

	isCut = true; // todo: this not always might be true to save computations we might test if this is the case
	if (mesh.IsSpecialWeightValid() && mesh.hasSpecialWeight(cellDim, celLID))
	{
		mesh.getSpecialWeight(cellDim, celLID, quadWeights);
		//SUNDANCE_MSG3(verb, tabs << "GaussianQuadrature::getAdaptedWeights Found cached weights for cell LID " << cellLID)
		return;
	}

	// if we have no quad points then get them
	if (quadPoints.size() <= 1) getPoints(cellType,  quadPoints, quadWeights);

	// we say that the cell is always cut by the curve so these weights will be used
	isCut = true;
	Array<Point> tmpPoints = quadPoints;

    Array<int> celLIDs(1);
    celLIDs[0] = celLID;
    // transform the ref points to real coordinates
	mesh.pushForward( cellDim, celLIDs, quadPoints, tmpPoints );

	quadWeights.resize(tmpPoints.size());
    // simple weight calculation
	for (int i=0 ; i < tmpPoints.size() ; i++){
		quadWeights[i] = quadWeights[i] * globalCurve.integrationParameter(tmpPoints[i]);
	}

	// store the weights in the mesh
	mesh.setSpecialWeight(cellDim, celLID, quadWeights);
/*
	Array<Point> vertices;
    Array<int>  nodeLIDs;
    Array<int>  orient;
    int testcount = 1;

    // Get coordinates
    mesh.getFacetArray(cellDim, celLID, 0, nodeLIDs, orient);

    vertices.resize(nodeLIDs.size());
    for (int i=0; i<nodeLIDs.size(); i++) {
    	vertices[i] = mesh.nodePosition(nodeLIDs[i]);
    	SUNDANCE_MSG5(6, "Points [" << testcount << "/" << i << "] " << vertices[i]);
    	testcount++;
    }*/
}
