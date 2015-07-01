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

#ifndef SUNDANCE_VECTORBASISCOMPONENT_H
#define SUNDANCE_VECTORBASISCOMPONENT_H

#include "SundanceDefs.hpp"
#include "SundanceBasisFamily.hpp"

namespace Sundance {

using namespace Teuchos;

/** 
 * This class is for the representation of a single component
 * of a vector-valued basis family. 
 */
class VectorBasisComponent : public BasisFamilyBase
{
public:
  /** */
  VectorBasisComponent(const BasisFamily& master, int direction);

  /** */
  bool lessThan(const BasisFamilyBase* other) const ;

  /** */
  int order() const {return master_.order();}

  /** */
  int dim() const 
    {return master_.dim();}

  /** */
  bool isCovariantBasis() const 
    {return master_.isCovariantBasis();}

  /** */
  bool isContravariantBasis() const 
    {return master_.isContravariantBasis();}

  /** */
  int direction() const {return direction_;}

  /** */
  bool supportsCellTypePair(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const
    {
      return master_.ptr()->supportsCellTypePair(maximalCellType, 
        cellType);
    }


  /** */
  void getReferenceDOFs(
    const CellType& maximalCellType,
    const CellType& cellType,
    Array<Array<Array<int> > >& dofs
    ) const 
    {
      master_.ptr()->getReferenceDOFs(maximalCellType, 
        cellType, dofs);
    }


  /** */
  int nReferenceDOFsWithFacets(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const
    {
      return master_.ptr()->nReferenceDOFsWithFacets(maximalCellType, 
        cellType);
    }

  /** */
  int nReferenceDOFsWithoutFacets(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const
    {
      return master_.ptr()->nReferenceDOFsWithoutFacets(maximalCellType, 
        cellType);
    }

  /** */
  void refEval(
    const CellType& maximalCellType,
    const CellType& cellType,
    const Array<Point>& pts,
    const MultiIndex& deriv,
    Array<Array<Array<double> > >& result
    ) const
    {
      master_.ptr()->refEval(maximalCellType, cellType,
        pts, deriv, result);
    }


private:
  BasisFamily master_;
  int direction_;
};


} // namespace Sundance


#endif
