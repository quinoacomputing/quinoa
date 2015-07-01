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

#include "SundanceInhomogeneousEdgeLocalizedDOFMap.hpp"

#include "SundanceEdgeLocalizedBasis.hpp"
#include "SundanceOrderedTuple.hpp"

#include <numeric>
#include <algorithm>

namespace Sundance {

using Teuchos::null;

InhomogeneousEdgeLocalizedDOFMap::InhomogeneousEdgeLocalizedDOFMap(const Mesh& mesh,
    const Array<Map<Set<int>, CellFilter> >& funcSetToDomainMap,
    int setupVerb) :
  DOFMapBase(mesh, setupVerb),
  funcDomains_(),
  edgeDofs_()
{
  SUNDANCE_MSG1(setupVerb, "in InhomogeneousEdgeLocalizedDOFMap ctor");

  TEUCHOS_TEST_FOR_EXCEPTION(mesh.comm().getNProc() != 1,
    std::runtime_error,
    "distributed inhomogeneous edge localized DOF maps not yet supported");

  SUNDANCE_MSG4(setupVerb, "func set to domain map " << funcSetToDomainMap);

  TEUCHOS_ASSERT(funcSetToDomainMap.size() > 1);
  TEUCHOS_ASSERT(funcSetToDomainMap[0].empty());

  TEUCHOS_ASSERT(funcSetToDomainMap.size() <= this->meshDimension() + 1);

  Map<int, Array<Array<CellFilter> > > funcFilters;
  Map<CellFilter, Array<int> > filterEdges;
  for (int dim = 1; dim <= this->meshDimension(); ++dim) {
    const Map<Set<int>, CellFilter> &funcsOnDom = funcSetToDomainMap[dim];
    for (Map<Set<int>, CellFilter>::const_iterator it = funcsOnDom.begin(),
      it_end = funcsOnDom.end();
      it != it_end;
      ++it)
    {
      const CellFilter &subdomain = it->second;
      filterEdges[subdomain] = getEdgeLIDs(subdomain);

      const Set<int> &funcs = it->first;
      for (Set<int>::const_iterator jt = funcs.begin(),
          jt_end = funcs.end();
          jt != jt_end;
          ++jt)
      {
        Array<Array<CellFilter> > &filters = funcFilters[*jt];
        filters.resize(this->meshDimension() + 1);
        filters[dim].append(subdomain);
      }
    }
  }

  SUNDANCE_MSG4(setupVerb, "subdomains to edges " << filterEdges);
  SUNDANCE_MSG4(setupVerb, "func ids to subdomains " << funcFilters);

  const int functionCount = funcFilters.empty() ? 0 : funcFilters.rbegin()->first + 1;
  funcDomains_.resize(functionCount);

  for (Map<int, Array<Array<CellFilter> > >::const_iterator funcIt = funcFilters.begin(),
      funcItEnd = funcFilters.end();
      funcIt != funcItEnd;
      ++funcIt)
  {
    const Array<CellFilter> &maxDimFilters = (funcIt->second)[this->meshDimension()];
    funcDomains_[funcIt->first] = std::accumulate(maxDimFilters.begin(), maxDimFilters.end(), CellFilter());
  }

  SUNDANCE_MSG4(setupVerb, "func ids to max dim subdomains " << funcDomains_);

  typedef OrderedPair<int, int> DofId; // DofId = (EdgeId, FuncId)
  Array<DofId> dofSet;
  for (Map<int, Array<Array<CellFilter> > >::const_iterator funcIt = funcFilters.begin(),
      funcItEnd = funcFilters.end();
      funcIt !=funcItEnd;
      ++funcIt)
  {
    const int funcId = funcIt->first;
    Set<int> dofBearingEdges;
    for (int dim = 1; dim <= this->meshDimension(); ++dim)
    {
      const Array<CellFilter> &dimFilters = (funcIt->second)[dim];
      for (Array<CellFilter>::const_iterator filtIt  = dimFilters.begin(),
          filtItEnd = dimFilters.end();
          filtIt != filtItEnd;
          ++filtIt)
      {
        const Array<int> &edges = filterEdges[*filtIt];
        dofBearingEdges.insert(edges.begin(), edges.end());
      }
    }
    for (Set<int>::const_iterator edgeIt = dofBearingEdges.begin(),
         edgeItEnd = dofBearingEdges.end();
         edgeIt != edgeItEnd;
         ++edgeIt)
    {
      dofSet.push_back(DofId(*edgeIt, funcId));
    }
  }

  const int dofCount = dofSet.size();

  SUNDANCE_MSG1(setupVerb, "number of degrees of freedom = " << dofCount);
  SUNDANCE_MSG4(setupVerb, "set of degrees of freedom " << dofSet);

  const int edgeCount = mesh.numCells(1);
  edgeDofs_.resize(edgeCount, Array<int>(functionCount, -1));

  std::sort(dofSet.begin(), dofSet.end()); // Sort by edge id, then by function id
  int dofRank = 0;
  for (Array<DofId>::const_iterator dofIt = dofSet.begin(),
      dofItEnd = dofSet.end();
      dofIt != dofItEnd;
      ++dofIt)
  {
    const int edgeId = dofIt->first();
    const int funcId = dofIt->second();
    edgeDofs_[edgeId][funcId] = dofRank++;
  }

  SUNDANCE_MSG4(setupVerb, "degrees of freedom on edges " << edgeDofs_);

  this->setLowestLocalDOF(0);
  this->setNumLocalDOFs(dofCount);
  this->setTotalNumDOFs(dofCount);
}

RCP<const MapStructure>
InhomogeneousEdgeLocalizedDOFMap::getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLIDs,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verb) const
{
  TEUCHOS_ASSERT(requestedFuncSet.setDifference(*this->allowedFuncsOnCellBatch(cellDim, cellLIDs)).empty());

  int edgesPerCell;
  if (cellDim == 0)
  {
    // No dofs on nodes
    edgesPerCell = 0;
    const Array<int> empty;
    getDOFsForEdgeBatch(empty, requestedFuncSet, dofs, verb);
  }
  else if (cellDim == 1)
  {
    // Cells are actually edges
    edgesPerCell = 1;
    getDOFsForEdgeBatch(cellLIDs, requestedFuncSet, dofs, verb);
  }
  else
  {
    // Edges are facets of the cells
    const int edgeDim = 1;
    edgesPerCell = this->mesh().numFacets(cellDim, 0, edgeDim);
    Array<int> edgeLIDs;
    {
      Array<int> dummyOrientations;
      this->mesh().getFacetLIDs(cellDim, cellLIDs, edgeDim, edgeLIDs, dummyOrientations);
      TEUCHOS_ASSERT(edgeLIDs.size() == edgesPerCell * cellLIDs.size());
    }

    getDOFsForEdgeBatch(edgeLIDs, requestedFuncSet, dofs, verb);
  }

  nNodes.resize(1);
  nNodes[0] = edgesPerCell;

  const Array<int> requestedFuncArray(requestedFuncSet.begin(), requestedFuncSet.end());
  return rcp(new MapStructure(requestedFuncSet.size(), rcp(new EdgeLocalizedBasis()), tuple(requestedFuncArray)));
}

void
InhomogeneousEdgeLocalizedDOFMap::getDOFsForEdgeBatch(const Array<int>& edgeLIDs,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    int verb) const
{
  dofs.resize(1); // One basis topology only
  Array<int> &dofChunk = dofs[0];
  dofChunk.resize(edgeLIDs.size() * requestedFuncSet.size());

  Array<int>::iterator entryIt = dofChunk.begin();
  for (Array<int>::const_iterator edgeIt = edgeLIDs.begin(),
      edgeItEnd = edgeLIDs.end();
      edgeIt != edgeItEnd;
      ++edgeIt)
  {
    for (Set<int>::const_iterator funcIt = requestedFuncSet.begin(),
        funcItEnd = requestedFuncSet.end();
        funcIt != funcItEnd;
        ++funcIt)
    {
      (*entryIt++) = edgeDofs_[*edgeIt][*funcIt];
    }
  }
}

RCP<const Set<int> >
InhomogeneousEdgeLocalizedDOFMap::allowedFuncsOnCellBatch(int cellDim,
    const Array<int>& cellLIDs) const
{
  if (cellDim == 0)
  {
    // No functions allowed on nodes
    return rcp(new Set<int>());
  }

  if (cellDim == this->meshDimension())
  {
    // For compatibility and as a special case,
    // allow all functions on cells of maximum dimension
    return this->allFuncIDs();
  }

  const int edgeDim = 1;
  if (cellDim == edgeDim)
  {
    // Cells are actually edges
    return this->allowedFuncsOnEdgeBatch(cellLIDs);
  }

  // Edges are facets of the cells
  Array<int> edgeLIDs;
  Array<int> dummyOrientations;
  this->mesh().getFacetLIDs(cellDim, cellLIDs, edgeDim, edgeLIDs, dummyOrientations);
  return this->allowedFuncsOnEdgeBatch(edgeLIDs);
}

RCP<Set<int> >
InhomogeneousEdgeLocalizedDOFMap::allowedFuncsOnEdgeBatch(const Array<int> &edgeLIDs) const
{
  const RCP<Set<int> > result = this->allFuncIDs();

  for (Array<int>::const_iterator edgeIt = edgeLIDs.begin(),
      edgeItEnd = edgeLIDs.end();
      edgeIt != edgeItEnd;
      ++edgeIt)
  {
    const Array<int> &dofs = edgeDofs_[*edgeIt];

    for (Set<int>::iterator it = result->begin(); it != result->end(); /* Nothing */)
    {
      if (dofs[*it] < 0)
      {
        result->erase(it++);
      }
      else
      {
        ++it;
      }
    }

    if (result->empty())
    {
      break;
    }
  }

  return result;
}

RCP<Set<int> >
InhomogeneousEdgeLocalizedDOFMap::allFuncIDs() const
{
  const RCP<Set<int> > result(new Set<int>());
  for (int i = 0; i < funcDomains_.size(); ++i)
  {
    result->insert(result->end(), i);
  }
  return result;
}

void
InhomogeneousEdgeLocalizedDOFMap::print(std::ostream& os) const
{
  os << "InhomogeneousEdgeLocalizedDOFMap\n";
  os << "Edges in mesh = " << edgeDofs_.size() << "\n";
  os << "Functions = " << funcDomains_.size() << "\n";

  for (int dim = 1; dim <= this->meshDimension(); ++dim)
  {
    os << "Dofs on cells of dimension " << dim << ":\n";
    const int cellCount = this->mesh().numCells(dim);
    for (int iCell = 0; iCell < cellCount; ++iCell)
    {
      os << "  Cell LID " << iCell << "\n";
      const RCP<const Set<int> > funcs = this->allowedFuncsOnCellBatch(dim, tuple(iCell));
      for (Set<int>::const_iterator it = funcs->begin(),
          it_end = funcs->end();
          it != it_end;
          ++it)
      {
        const int funcId = *it;
        os << "    Function " << funcId << " : {";
        Array<int> dofs;
        this->getDOFsForCell(dim, iCell, funcId, dofs);
        for (Array<int>::const_iterator jt = dofs.begin(),
            jt_end = dofs.end();
            jt != jt_end;
            ++jt)
        {
          if (jt != dofs.begin())
          {
            os << ", ";
          }
          os << *jt;
        }
        os << "}\n";
      }
    }
  }
}

int
InhomogeneousEdgeLocalizedDOFMap::meshDimension() const
{
  return this->mesh().spatialDim();
}

Array<int>
InhomogeneousEdgeLocalizedDOFMap::getEdgeLIDs(const CellFilter &filter) const
{
  const int cellDim = filter.dimension(this->mesh());

  const int nodeDim = 0;
  if (cellDim == nodeDim)
  {
    // Ignore isolated nodes
    return Array<int>();
  }

  const CellSet cells = filter.getCells(this->mesh());
  const Array<int> cellLIDs(cells.begin(), cells.end());

  const int edgeDim = 1;
  if (cellDim == edgeDim)
  {
    // Cells are actually edges
    return cellLIDs;
  }

  // Edges are cofacets of the cells
  Array<int> edgeLIDs;
  Array<int> dummyOrientations;
  this->mesh().getFacetLIDs(cellDim, cellLIDs, edgeDim, edgeLIDs, dummyOrientations);

  return edgeLIDs;
}

} // namespace Sundance
