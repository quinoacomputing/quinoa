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

#ifndef SUNDANCE_HOMOGENEOUSDOFMAP_H
#define SUNDANCE_HOMOGENEOUSDOFMAP_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceDOFMapBase.hpp"


namespace Sundance
{
using namespace Teuchos;

/** 
 * A HomogeneousDOFMap is a DOF map for the special (and common)
 * case in which every function has the same basis and is defined
 * on every cell in the mesh. 
 */
class HomogeneousDOFMap : public DOFMapBase
{
public:
  /** */
  HomogeneousDOFMap(const Mesh& mesh, 
    const BasisFamily& basis,
    int numFuncs);
                        
  /** */
  HomogeneousDOFMap(const Mesh& mesh, 
    const BasisFamily& basis,
    const Array<CellFilter>& subregions,
    int numFuncs);

  /** */
  virtual ~HomogeneousDOFMap(){;}


     

      

  /** */
  virtual void getDOFsForCellBatch(int cellDim, const Array<int>& cellLID,
    Array<int>& dofs,  
    Array<Array<int> >& funcIDs,
    Array<int>& nNodes) const ;


  /** */
  virtual void print(std::ostream& os) const ;

private:

  /** */
  void allocate(const Mesh& mesh, 
    const BasisFamily& basis,
    int numFuncs);
      
  /** */
  void buildMaximalDofTable() const ;

  /** */
  bool hasBeenAssigned(int cellDim, int cellLID) const 
    {return dofs_[cellDim][cellLID][0] != uninitializedVal();}

  /** */
  void initMap();

  /** */
  void setDOFs(int cellDim, int cellLID, 
    int& nextDOF, bool isRemote=false);

  /** */
  void shareDOFs(int cellDim,
    const Array<Array<int> >& outgoingCellRequests);

  /** */
  void computeOffsets(int dim, int localCount);

  /** */
  const Array<int>& funcIDList() const {return funcIDOnCellSet(0);}

  static int uninitializedVal() {return -1;}

  int dim_;

  Array<Array<Array<int> > > dofs_;

  mutable Array<int> maximalDofs_;

  mutable bool haveMaximalDofs_;

  Array<Array<Array<Array<int> > > > localNodePtrs_;

  Array<int> nNodesPerCell_;

  Array<int> totalNNodesPerCell_;

  Array<Array<int> > numFacets_;

  Array<Array<int> > originalFacetOrientation_;

  bool basisIsContinuous_;

      
};
}


                  


#endif
