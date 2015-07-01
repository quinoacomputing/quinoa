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

#include "SundanceFeketeQuadratureType.hpp"
#include "SundanceFeketeQuadrature.hpp"
#include "SundanceGaussLobatto1D.hpp"
#include "SundanceFeketeTriangleQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;

FeketeQuadratureType::FeketeQuadratureType() :
	QuadratureTypeBase()
{

}

XMLObject FeketeQuadratureType::toXML() const
{
	XMLObject rtn("FeketeQuadratureType");
	return rtn;
}

bool FeketeQuadratureType::supportsCellType(const CellType& cellType) const
{
	switch (cellType)
	{
	case PointCell:
	case LineCell:
	case TriangleCell:
	case QuadCell:
	case BrickCell:
		return true;
	default:
		return false;
	}
}

int FeketeQuadratureType::maxOrder(const CellType& cellType) const
{
	switch (cellType)
	{
	case TriangleCell:
		return FeketeTriangleQuadrature::maxOrder();
	default:
		return -1;
	}
}

bool FeketeQuadratureType::hasLimitedOrder(const CellType& cellType) const
{
	switch (cellType)
	{
	case TriangleCell:
		return true;
	default:
		return false;
	}
}

bool FeketeQuadratureType::supports(const CellType& cellType, int order) const
{
	if (order <= 0)
		return false;

	switch (cellType)
	{
	case PointCell:
		return true;
	case LineCell:
		return true;
	case QuadCell:
		return true;
	case BrickCell:
		return true;
	case TriangleCell:
		return FeketeTriangleQuadrature::supportsOrder(order);
	default:
		return false;
	}
}

QuadratureFamily FeketeQuadratureType::createQuadFamily(int order) const
{
	return new FeketeQuadrature(order);
}
