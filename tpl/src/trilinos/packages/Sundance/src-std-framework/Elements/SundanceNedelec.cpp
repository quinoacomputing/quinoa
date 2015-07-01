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

#include "SundanceNedelec.hpp"
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


Nedelec::Nedelec(int spatialDim)
  : HCurlVectorBasis(spatialDim)
{}

bool Nedelec::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case TriangleCell:
      switch(cellType)
      {
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    case TetCell:
      // tets not yet implemented
      return false;
      switch(cellType)
      {
        case TetCell:
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}

void Nedelec::print(std::ostream& os) const 
{
  os << "Nedelec()";
}

int Nedelec::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(maximalCellType != TriangleCell);
  switch(cellType)
    {
    case PointCell:
      return 0;
    case LineCell:
      return 1;
    case TriangleCell:
      return 0;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
                         << cellType << " not implemented in Nedelec basis");
      return -1; // -Wall
    }
}

void Nedelec::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  switch(cellType)
    {
    case PointCell:
      dofs.resize(1);
      dofs[0] = tuple(Array<int>());
      return;
    case LineCell:
      dofs.resize(2);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0));
      return;
    case TriangleCell:
      dofs.resize(3);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0), tuple(1), tuple(2));
      dofs[2] = tuple(Array<int>());
      return;
    case TetCell:
      dofs.resize(4);
      dofs[0] = tuple(Array<int>());
      dofs[1] = tuple<Array<int> >(tuple(0), tuple(1), tuple(2), tuple(3));
      dofs[2] = tuple(Array<int>());
      dofs[3] = tuple(Array<int>());
      return;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
                         << cellType << " not implemented in Nedelec basis");
    }
}



void Nedelec::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& deriv,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "evaluation of Nedelec elements not yet supported");
}


