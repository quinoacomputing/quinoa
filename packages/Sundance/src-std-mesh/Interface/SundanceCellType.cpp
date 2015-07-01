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

#include "SundanceCellType.hpp"
#include "PlayaExceptions.hpp"

using namespace Sundance;
using namespace Sundance;

namespace Sundance
{

  std::string toString(const CellType& cellType)
  {
    switch(cellType)
      {
      case NullCell:
        return "NullCell";
      case PointCell:
        return "PointCell";
      case LineCell:
        return "LineCell";
      case TriangleCell:
        return "TriangleCell";
      case QuadCell:
        return "QuadCell";
      case TetCell:
        return "TetCell";
      case BrickCell:
        return "BrickCell";
      case PrismCell:
        return "PrismCell";
      }
    return "NullCell"; // -Wall
  }

  int dimension(const CellType& cellType)
  {
    switch(cellType)
      {
      case NullCell:
        return -1;
      case PointCell:
        return 0;
      case LineCell:
        return 1;
      case TriangleCell:
      case QuadCell:
        return 2;
      case TetCell:
      case BrickCell:
      case PrismCell:
        return 3;
      }
    return -1; // -Wall
  }

  int numFacets(const CellType& cellType, int facetDim)
  {
    int d = dimension(cellType);
    if (facetDim == d) return 1;

    TEUCHOS_TEST_FOR_EXCEPTION(facetDim > d, std::runtime_error,
                       "invalid facet dim " << facetDim << " for cell "
                       << toString(cellType));

    switch(cellType)
      {
      case NullCell:
      case PointCell:
        return -1;
      case LineCell:
        return 2;
      case TriangleCell:
        return 3;
      case TetCell:
        if (facetDim==0 || facetDim==2) return 4;
        return 6;
      case QuadCell:
        return 4;
      case BrickCell:
        if (facetDim==0) return 8;
        else if (facetDim==1) return 12;
        else return 6;
      case PrismCell:
        if (facetDim==0) return 6;
        else if (facetDim==1) return 9;
        else return 5;
      }
    return -1; // -Wall
  }


  CellType facetType(const CellType& cellType, int facetDim, int facetIndex)
  {
    int d = dimension(cellType);
    if (facetDim == d) return cellType;
    TEUCHOS_TEST_FOR_EXCEPTION(facetDim > d, std::runtime_error,
                       "invalid facet dim " << facetDim << " for cell "
                       << toString(cellType));

    if (facetDim==0) return PointCell;
    if (facetDim==1) return LineCell;

    switch(cellType)
      {
      case NullCell:
      case PointCell:
      case LineCell:
      case TriangleCell:
      case QuadCell:
        return NullCell;

      case TetCell:
        return TriangleCell;
      case BrickCell:
        return QuadCell;
      case PrismCell:
        if (facetIndex==0 || facetIndex==4) return TriangleCell;
        return QuadCell;
      }
    return NullCell;
  }
  
}
