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

#include "SundanceEdgeLocalizedBasis.hpp"
#include "SundanceADReal.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


EdgeLocalizedBasis::EdgeLocalizedBasis()
{}

bool EdgeLocalizedBasis::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case TriangleCell:
    case TetCell:
    case QuadCell:
    case BrickCell:
      switch(cellType)
      {
        case LineCell:
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}

void EdgeLocalizedBasis::print(std::ostream& os) const 
{
  os << "EdgeLocalizedBasis()";
}

int EdgeLocalizedBasis::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(cellType)
  {
    case LineCell:
      return 1;
    default:
      return 0;
  }
}

void EdgeLocalizedBasis::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  typedef Array<int> Aint;
  switch(cellType)
  {
    case LineCell:
      dofs.resize(2);
      dofs[0] = Array<Array<int> >();
      dofs[1] = tuple<Aint>(tuple<int>(0));
      return;
    case TriangleCell:
      dofs.resize(3);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Aint>(tuple(0), tuple(1), tuple(2));
      dofs[2] = tuple(Array<int>());
      return;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
        << cellType << " not implemented in EdgeLocalizedBasis basis");
  }
}



void EdgeLocalizedBasis::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    std::runtime_error,
    "cannot evaluate spatial derivative " << sds << " on EdgeLocalizedBasis basis");
  const MultiIndex& deriv = sds.mi();
  typedef Array<double> Adouble;
  result.resize(1);
  result[0].resize(pts.length());

  int dim = dimension(cellType);
  
  if (dim==0)
  {
    result[0] = tuple<Adouble>(tuple(1.0));
  }
  else if (dim==1)
  {
    for (int i=0; i<pts.length(); i++)
    {
      evalOnLine(pts[i], deriv, result[0][i]);
    }
  }
  else if (dim==2)
  {
    for (int i=0; i<pts.length(); i++)
    {
      evalOnTriangle(pts[i], deriv, result[0][i]);
    }
  }
  else if (dim==3)
  {
    for (int i=0; i<pts.length(); i++)
    {
      evalOnTet(pts[i], deriv, result[0][i]);
    }
  }
}

/* ---------- evaluation on different cell types -------------- */


void EdgeLocalizedBasis::evalOnLine(const Point& pt, 
													const MultiIndex& deriv,
													Array<double>& result) const
{
	ADReal one(1.0, 1);
	result.resize(1);
	Array<ADReal> tmp(result.length());

  tmp[0] = one;

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}

void EdgeLocalizedBasis::evalOnTriangle(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  ADReal x = ADReal(pt[0], 0, 2);
	ADReal y = ADReal(pt[1], 1, 2);
	ADReal one(1.0, 2);
	ADReal zero(0.0, 2);

  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
    << y.value());

  result.resize(3);
  tmp.resize(3);

  bool onEdge2 = std::fabs(pt[1]) < 1.0e-14;
  bool onEdge0 = std::fabs(1.0-pt[0]-pt[1]) < 1.0e-14;
  bool onEdge1 = std::fabs(pt[0]) < 1.0e-14;
  
  TEUCHOS_TEST_FOR_EXCEPTION(!(onEdge0 || onEdge1 || onEdge2),
    std::runtime_error,
    "EdgeLocalizedBasis should not be evaluated at points not on edges");
  
  TEUCHOS_TEST_FOR_EXCEPTION((onEdge0 && onEdge1) || (onEdge1 && onEdge2)
    || (onEdge2 && onEdge0), std::runtime_error,
    "Ambiguous edge in EdgeLocalizedBasis::evalOnTriangle()");

  if (onEdge0)
  {
    tmp[0] = one;
    tmp[1] = zero;
    tmp[2] = zero;
  }
  if (onEdge1)
  {
    tmp[0] = zero;
    tmp[1] = one;
    tmp[2] = zero;
  }
  if (onEdge2)
  {
    tmp[0] = zero;
    tmp[1] = zero;
    tmp[2] = one;
  }


	for (int i=0; i<tmp.length(); i++)
  {
    SUNDANCE_OUT(this->verb() > 3,
      "tmp[" << i << "]=" << tmp[i].value() 
      << " grad=" << tmp[i].gradient());
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
  
}




void EdgeLocalizedBasis::evalOnTet(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "EdgeLocalizedBasis::evalOnTet not implemented");
}
