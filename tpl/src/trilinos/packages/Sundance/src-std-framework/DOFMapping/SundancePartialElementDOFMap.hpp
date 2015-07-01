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

#ifndef SUNDANCE_PARTIALELEMENTDOFMAP_H
#define SUNDANCE_PARTIALELEMENTDOFMAP_H

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
 * PartialElementDOFMap is a DOF map specialized to the case of element-based
 * DOFs on a subset of cells in the domain. All elements must have the same
 * set of functions. 
 */
class PartialElementDOFMap : public DOFMapBase
{
public:
  /** */
  PartialElementDOFMap(const Mesh& mesh, 
    const CellFilter& subdomain,
    int nFuncs,
    int setupVerb);
      
  /** */
  virtual ~PartialElementDOFMap(){;}

  /** */
  RCP<const MapStructure> 
  getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLID,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verbosity) const ;

  /** */
  RCP<const Set<int> >
  allowedFuncsOnCellBatch(int cellDim,
    const Array<int>& cellLID) const ;

  /** */
  const Array<CellFilter>& funcDomains() const {return funcDomains_;}

  /** */
  virtual void print(std::ostream& os) const ;


protected:

  void init();

  void computeOffsets(int localCount)  ;

  void shareRemoteDOFs(const Array<Array<int> >& remoteElems);

  int dim_;
  int nFuncs_;
  int nElems_;
  CellFilter subdomain_;
  Array<CellFilter> funcDomains_;
  Array<int> elemDofs_;
  RCP<MapStructure> structure_;
  RCP<const Set<int> > allFuncs_;
};

}


#endif
