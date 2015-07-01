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

#include "SundanceFieldBase.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceMesh.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

using Sundance::CellSet;
using Sundance::CellIterator;
using Sundance::CellFilter;
using Sundance::MaximalCellFilter;


void FieldBase::getDataBatch(
  int cellDim, 
  const Array<int>& cellID,
  const Array<int>& funcElem,
  Array<double>& batch) const
{
  /* This is a silly default implementation */
  
  int nFunc = funcElem.size();
  batch.resize(cellID.size() * funcElem.size());

  for (int c=0; c<cellID.size(); c++)
  {
    for (int f=0; f<funcElem.size(); f++)
    {
      batch[c*nFunc + f] = getData(cellDim, cellID[c], f);
    }
  }
}

const CellFilter& FieldBase::domain() const 
{
  static CellFilter dum = new MaximalCellFilter();
  return dum;
}


namespace Sundance
{

using Sundance::CellSet;

CellSet connectedNodeSet(const CellFilter& f, const Mesh& mesh)
{
  CellSet cells = f.getCells(mesh);
  int dim = cells.dimension();
  if (dim==0) return cells;


  Array<int> cellLID;

  for (CellIterator i=cells.begin(); i!=cells.end(); i++)
  {
    cellLID.append(*i);
  }

  Array<int> nodes;
  Array<int> fo;

  mesh.getFacetLIDs(dim, cellLID, 0, nodes, fo);

  Set<int> nodeSet;

  for (int i=0; i<nodes.size(); i++)
  {
    nodeSet.put(nodes[i]);
  }
  
  return CellSet(mesh, 0, PointCell, nodeSet);
}


RCP<Array<int> > cellSetToLIDArray(const CellSet& cs)
{
  RCP<Array<int> > cellLID = rcp(new Array<int>());

  for (CellIterator i=cs.begin(); i!=cs.end(); i++)
  {
    cellLID->append(*i);
  }
  
  return cellLID;
}

}
