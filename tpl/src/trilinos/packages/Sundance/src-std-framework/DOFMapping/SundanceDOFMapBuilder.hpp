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

#ifndef SUNDANCE_DOFMAPBUILDER_H
#define SUNDANCE_DOFMAPBUILDER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceFunctionSupportResolver.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCFMeshPair.hpp"
#include "SundanceMap.hpp"
#include "SundanceObjectWithVerbosity.hpp"

namespace Sundance
{


/** 
 * 
 */
class DOFMapBuilder 
{
public:
  /** */
  DOFMapBuilder(int setupVerb);
  /** */
  DOFMapBuilder(const Mesh& mesh, const RCP<FunctionSupportResolver>& fsr, 
    bool findBCCols, int setupVerb);

  /** */
  const Array<RCP<DOFMapBase> >& rowMap() const {return rowMap_;}

  /** */
  const Array<RCP<DOFMapBase> >& colMap() const {return colMap_;}

  /** */
  const Array<RCP<Array<int> > >& isBCRow() const {return isBCRow_;}

  /** */
  const Array<RCP<Array<int> > >& isBCCol() const {return isBCCol_;}


  /** */
  const Array<RCP<std::set<int> > >& remoteBCCols() const 
    {return remoteBCCols_;}

  Array<Array<RCP<BasisDOFTopologyBase> > > testBasisTopologyArray() const ;

  Array<Array<RCP<BasisDOFTopologyBase> > > unkBasisTopologyArray() const ;

  Array<Array<Set<CellFilter> > > testCellFilters() const ;

  Array<Array<Set<CellFilter> > > unkCellFilters() const ;

  const Mesh& mesh() const {return mesh_;}



  RCP<DOFMapBase> makeMap(const Mesh& mesh,
    const Array<RCP<BasisDOFTopologyBase> >& basis,
    const Array<Set<CellFilter> >& filters) ;

  bool hasOmnipresentNodalMap(const Array<RCP<BasisDOFTopologyBase> >& basis,
    const Mesh& mesh,
    const Array<Set<CellFilter> >& filters) const ;

  bool hasCommonDomain(const Array<Set<CellFilter> >& filters) const ;

  bool hasNodalBasis(const Array<RCP<BasisDOFTopologyBase> >& basis) const ;
  
  bool hasEdgeLocalizedBasis(const Array<RCP<BasisDOFTopologyBase> >& basis) const ;

  bool hasCellBasis(const Array<RCP<BasisDOFTopologyBase> >& basis) const ;

  bool filtersAreZeroDimensional(const Mesh& mesh,
    const Array<Set<CellFilter> >& filters) const ;

  bool allFuncsAreOmnipresent(const Mesh& mesh,
    const Array<Set<CellFilter> >& filters) const ;

  bool isWholeDomain(const Mesh& mesh,
    int maxFilterDim,
    const Set<CellFilter>& filters) const ;

  CellFilter getMaxCellFilter(const Array<Set<CellFilter> >& filters) const ;

  static bool& allowNodalMap() {static bool rtn=true; return rtn;}

  /** */
  void extractUnkSetsFromFSR(const FunctionSupportResolver& fsr,
    Array<Set<int> >& funcSets,
    Array<CellFilter>& regions) const ;

  /** */
  void extractVarSetsFromFSR(const FunctionSupportResolver& fsr,
    Array<Set<int> >& funcSets,
    Array<CellFilter>& regions) const ;

  /** */
  const RCP<FunctionSupportResolver>& fsr() const {return fsr_;}

  /** */
  Sundance::Map<Set<int>, Set<CellFilter> > 
  buildFuncSetToCFSetMap(const Array<Set<int> >& funcSets,
    const Array<CellFilter>& regions,
    const Mesh& mesh) const ;
        
  void getSubdomainUnkFuncMatches(const FunctionSupportResolver& fsr,
    Array<Sundance::Map<CellFilter, Set<int> > >& fmap) const ;
        
  void getSubdomainVarFuncMatches(const FunctionSupportResolver& fsr,
    Array<Sundance::Map<CellFilter, Set<int> > >& fmap) const ;

  Array<Sundance::Map<Set<int>, CellFilter> > 
  funcDomains(const Mesh& mesh,
    const Sundance::Map<CellFilter, Set<int> >& fmap,
    Sundance::Map<CellFilter, Sundance::Map<Set<int>, CellSet> >& inputToChildrenMap) const ;

  Sundance::Map<CellFilter, Set<int> > domainToFuncSetMap(const Array<Set<CellFilter> >& filters) const ;

private:

  Set<CellFilter> reduceCellFilters(const Mesh& mesh,
    const Set<CellFilter>& inputSet) const ;

  bool hasUnks() const ;

  bool unksAreOmnipresent() const ;

  bool testsAreOmnipresent() const ;

  bool regionIsMaximal(int r) const ;

  bool isSymmetric(int block) const ;

  void markBCRows(int block) ;

  void markBCCols(int block) ;

  const MPIComm& comm() const {return mesh().comm();}

  void init(bool findBCCols);

  int verb_;

  Mesh mesh_;

  RCP<FunctionSupportResolver> fsr_;

  Array<RCP<DOFMapBase> > rowMap_;

  Array<RCP<DOFMapBase> > colMap_;

  Array<RCP<Array<int> > > isBCRow_;

  Array<RCP<Array<int> > > isBCCol_;

  Array<RCP<std::set<int> > > remoteBCCols_;

};

/** \relates DOFMapBuilder */
Array<Array<BasisFamily> > testBasisArray(const RCP<FunctionSupportResolver>& fsr) ;

/** \relates DOFMapBuilder */
Array<Array<BasisFamily> > unkBasisArray(const RCP<FunctionSupportResolver>& fsr) ;

}


#endif
