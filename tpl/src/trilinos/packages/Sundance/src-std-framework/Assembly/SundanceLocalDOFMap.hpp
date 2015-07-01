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



#ifndef SUNDANCE_LOCALDOFMAP_H
#define SUNDANCE_LOCALDOFMAP_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"


namespace Sundance
{

/** 
 * LocalDOFMap bundles several compact tables used for fast lookup of
 * local DOFs. 
 */
class LocalDOFMap
{
public:
  /** */
  LocalDOFMap(int numBlocks, int verb);

  /** */
  void setCells(int cellDim, int maxCellDim,
    const RCP<const Array<int> >& cellLID);

  /** */
  int nCells() const ;

  /** */
  bool isUsed(int b) const {return isUsed_[b];}

  /** */
  bool isUnused(int b) const {return !isUsed(b);}

  /** */
  bool isUnused() const ;

  /** */
  void markAsUnused() ;

  /** */
  bool hasCells() const {return hasCells_;}

  /** */
  const RCP<const Array<int> >& cellLIDs() const {return cellLID_;}

  /** */
  void markAsUsed(int b) {isUsed_[b] = true;}

  /** */
  int numBlocks() const {return mapStruct_->size();}

  /** */
  const Array<int>& nLocalNodesPerChunk(int b) const 
    {return (*nLocalNodesPerChunk_)[b];}

  /** */
  const RCP<const MapStructure>& mapStruct(int b) const 
    {return (*mapStruct_)[b];}

  /** */
  const Array<Array<int> >& localDOFs(int b) const 
    {return (*localDOFs_)[b];}

  /** */
  std::ostream& print(std::ostream& os) const ;

  /** */
  void fillBlock(int b, const RCP<DOFMapBase>& globalMap,
    const Array<Set<int> >& requiredFunc);

  /** */
  void setVerb(int v) const {verb_ = v;}

private:


  /** */
  Array<int>& nLocalNodesPerChunk(int b)
    {return (*nLocalNodesPerChunk_)[b];}

  /** */
  RCP<const MapStructure>& mapStruct(int b) 
    {return (*mapStruct_)[b];}

  /** */
  Array<Array<int> >& localDOFs(int b) 
    {return (*localDOFs_)[b];}

  /** */
  void verifyValidBlock(int b) const ;


  mutable int verb_;
  Array<int> isUsed_;
  bool hasCells_;
  RCP<Array<Array<int> > > nLocalNodesPerChunk_;
  RCP<Array<RCP<const MapStructure> > > mapStruct_;  
  RCP<Array<Array<Array<int> > > > localDOFs_;
  RCP<const Array<int> > cellLID_;
  int activeCellDim_;
  int maxCellDim_;
};


/** \relates LocalDOFMap */
inline std::ostream& operator<<(std::ostream& os,
  const LocalDOFMap& m)
{
  return m.print(os);
}

}



#endif
