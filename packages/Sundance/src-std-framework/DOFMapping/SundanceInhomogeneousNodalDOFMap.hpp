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

#ifndef SUNDANCE_INHOMOGENEOUSNODALDOFMAP_H
#define SUNDANCE_INHOMOGENEOUSNODALDOFMAP_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceObjectWithVerbosity.hpp"


namespace Sundance
{
using namespace Teuchos;

/** 
 * 
 */
class InhomogeneousNodalDOFMap : public DOFMapBase
{
public:
  /** */
  InhomogeneousNodalDOFMap(const Mesh& mesh, 
    const Array<Map<Set<int>, CellFilter> >& funcSetToDomainMap, 
    int setupVerb);
      
  /** */
  virtual ~InhomogeneousNodalDOFMap(){;}

  /** */
  RCP<const MapStructure> 
  getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLID,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verb) const ;

protected:
  /** */
  void getFunctionDofs(int cellDim,
    const Array<int>& cellLID,
    const Array<int>& facetLID,
    const Array<int>& funcs,
    Array<Array<int> >& dofs) const ;

public:
  /** */
  RCP<const Set<int> >
  allowedFuncsOnCellBatch(int cellDim,
    const Array<int>& cellLID) const ;

  /** */
  const Array<CellFilter>& funcDomains() const {return funcDomains_;}

  /** */
  virtual void print(std::ostream& os) const ;


protected:

  /** */
  Array<int> dofsOnCell(int cellDim, int cellLID, const Set<int>& reqFuncs) const ;
                              

  void init();

  void computeOffsets(int localCount)  ;

  void shareRemoteDOFs(const Array<Array<int> >& remoteNodes);

  void assignNode(int fLID,
    int funcComboIndex,
    int dofOffset,
    int nFuncs,
    Array<Array<int> >& remoteNodes,
    Array<int>& hasProcessedCell,
    int& nextDOF) ;

  int dim_;
  RCP<BasisDOFTopologyBase> basis_;
  int nTotalFuncs_;
  Array<CellFilter> funcDomains_;

  Array<Array<int> > nodeDofs_;
  Array<Array<int> > elemDofs_;
  Array<int> nodeToFuncSetIndexMap_;
  Array<int> elemToFuncSetIndexMap_;
  Array<Set<int> > elemFuncSets_;
  Array<Set<int> > nodalFuncSets_;
  Array<int> nodeToOffsetMap_;
  Array<int> elemToOffsetMap_;

  Array<Array<int> > funcIndexWithinNodeFuncSet_;

  Array<RCP<const MapStructure> > elemStructure_;
  Array<RCP<const MapStructure> > nodeStructure_;
};
}


#endif
