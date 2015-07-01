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

#include "SundanceClosedNewtonCotes.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;

ClosedNewtonCotes::ClosedNewtonCotes(int order)
  : QuadratureFamilyBase(order)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(order==2 || order==3), std::runtime_error, "order " 
    << order << " not supported by ClosedNewtonCotes");
}

XMLObject ClosedNewtonCotes::toXML() const 
{
  XMLObject rtn("ClosedNewtonCotes");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}


void ClosedNewtonCotes::getLineRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const 
{
  Array<double> x;
  if (order()==2)
  {
    x = tuple(0.0, 1.0);
    quadWeights = tuple(0.5, 0.5);
  }
  else 
  {
    x = tuple(0.0, 0.5, 1.0);
    quadWeights = tuple(1.0/6.0, 2.0/3.0, 1.0/6.0);
  }
  quadPoints.resize(x.size());
    
  for (int i=0; i<x.length(); i++)
  {
    quadPoints[i] = Point(x[i]);
  }
}

void ClosedNewtonCotes::getTriangleRule(Array<Point>& quadPoints,
                                          Array<double>& quadWeights) const 
{
  if (order()==2)
  {
    quadPoints=tuple(Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0));
    quadWeights=tuple(1.0/3.0, 1.0/3.0, 1.0/3.0);
  }
  else
  {
    quadPoints = tuple(Point(0.0, 0.0), Point(0.5, 0.0), Point(1.0, 0.0),
      Point(0.5, 0.5), Point(0.0, 1.0), Point(0.0, 0.5),
      Point(1.0/3.0, 1.0/3.0));
    quadWeights = tuple(3.0/60.0, 8.0/60.0, 3.0/60.0,
      8.0/60.0, 3.0/60.0, 8.0/60.0,
      27.0/60.0);
  }
	for (int i=0; i<quadWeights.size(); i++)
		{
			quadWeights[i] *=  0.5;
		}
}

