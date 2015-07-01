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

#include "SundanceQuadratureFamilyBase.hpp"


using namespace Sundance;
using namespace Teuchos;

void QuadratureFamilyBase::getPoints(const CellType& cellType, 
                                     Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const
{
  switch(cellType)
    {
    case PointCell:
      quadPoints = tuple(Point());
      quadWeights = tuple(1.0);
      break;
    case LineCell:
      getLineRule(quadPoints, quadWeights);
      break;
    case TriangleCell:
      getTriangleRule(quadPoints, quadWeights);
      break;
    case QuadCell:
      getQuadRule(quadPoints, quadWeights);
      break;
    case TetCell:
      getTetRule(quadPoints, quadWeights);
      break;
    case BrickCell:
      getBrickRule(quadPoints, quadWeights);
      break;
    default:
#ifndef TRILINOS_7
      SUNDANCE_ERROR("cell type " << cellType << " not handled in "
                     "QuadratureFamilyBase::getPoints()");
#else
      SUNDANCE_ERROR7("cell type " << cellType << " not handled in "
                     "QuadratureFamilyBase::getPoints()");
#endif
    }
}

void QuadratureFamilyBase::getAdaptedWeights(const CellType& cellType ,
									 int cellDim,
	                                 int celLID ,
	                	             int facetIndex ,
	                                 const Mesh& mesh ,
	                                 const ParametrizedCurve& globalCurve ,
	                                 Array<Point>& quadPoints ,
	                                 Array<double>& quadWeights ,
	                                 bool &isCut) const {
#ifndef TRILINOS_7
  SUNDANCE_ERROR("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for " << toXML());
#else
  SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getLineRule(Array<Point>& /* quadPoints */,
                                       Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Line rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Line rule not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getTriangleRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Triangle rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Triangle rule not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getQuadRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Quad cell rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Quad cell rule not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getTetRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Tet cell rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Tet cell rule not available for " << toXML());
#endif
}

void QuadratureFamilyBase::getBrickRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
#ifndef TRILINOS_7
  SUNDANCE_ERROR("Brick rule not available for " << toXML());
#else
  SUNDANCE_ERROR7("Brick rule not available for " << toXML());
#endif
}
