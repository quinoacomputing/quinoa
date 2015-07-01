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

#include "SundancePolygonQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;

int PolygonQuadrature::nrMaxLinePerCell_ = 5;

PolygonQuadrature::PolygonQuadrature( const QuadratureFamily& quad )
  : QuadratureFamilyBase(quad.order()) , quad_(quad)
{
  
}

XMLObject PolygonQuadrature::toXML() const
{
  XMLObject rtn("PolygonQuadrature");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}



void PolygonQuadrature::getLineRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const 
{
	Array<Point> quadPoints_tmp;
	Array<double> quadWeights_tmp;
	quad_.getPoints( LineCell , quadPoints_tmp , quadWeights_tmp );

	// the nr. of points per line segments
	int nrPointPerLine = quadPoints_tmp.size() , ind = 0;

	// resize the point arrays and the weight arrays
	quadPoints.resize( nrMaxLinePerCell_*nrPointPerLine );
	quadWeights.resize( nrMaxLinePerCell_*nrPointPerLine );

	// each line segment
	for (int nrl = 0 ; nrl < nrMaxLinePerCell_ ; nrl++ ){
		// loop over each quadrature point
		for (int q = 0 ; q < nrPointPerLine ; q++ ){
			// copy the points and the weights several times
			quadPoints[ind] = quadPoints_tmp[q];
			quadWeights[ind] = quadWeights_tmp[q]/((double)nrMaxLinePerCell_);
			ind++;
		}
	}
}
