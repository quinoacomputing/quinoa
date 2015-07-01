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

#include "SundanceDOFMapBuilder.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceEdgeLocalizedBasis.hpp"
#include "SundanceMixedDOFMap.hpp"
#include "SundanceMixedDOFMapHN.hpp"
#include "SundanceNodalDOFMap.hpp"
#include "SundanceSubmaximalNodalDOFMap.hpp"
#include "SundanceNodalDOFMapHN.hpp"
#include "SundancePartialElementDOFMap.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceInhomogeneousNodalDOFMap.hpp"
#include "SundanceInhomogeneousEdgeLocalizedDOFMap.hpp"
#include "SundanceInhomogeneousDOFMapHN.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCFMeshPair.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"


using namespace Sundance;
using namespace Teuchos;
using Playa::MPIComm;
using Playa::MPIDataType;
using Playa::MPIOp;

static Time& DOFBuilderCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("DOF map building"); 
  return *rtn;
}

static Time& cellFilterReductionTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("cell filter reduction"); 
  return *rtn;
}

static Time& findFuncDomainTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("finding func domains"); 
  return *rtn;
}

DOFMapBuilder::DOFMapBuilder(const Mesh& mesh, 
  const RCP<FunctionSupportResolver>& fsr, bool findBCCols,
  int setupVerb)
  : verb_(setupVerb),
    mesh_(mesh),
    fsr_(fsr),
    rowMap_(),
    colMap_(),
    isBCRow_(),
    isBCCol_(),
    remoteBCCols_()
{
  init(findBCCols);
}

DOFMapBuilder::DOFMapBuilder(int setupVerb)
  : verb_(setupVerb),
    mesh_(),
    fsr_(),
    rowMap_(),
    colMap_(),
    isBCRow_(),
    isBCCol_(),
    remoteBCCols_()
{}

RCP<DOFMapBase> DOFMapBuilder::makeMap(const Mesh& mesh,
  const Array<RCP<BasisDOFTopologyBase> >& basis,
  const Array<Set<CellFilter> >& filters) 
{
  verb_=0;
  TimeMonitor timer(DOFBuilderCtorTimer());
  SUNDANCE_MSG1(verb_, "in DOFMapBuilder::makeMap()");
  for (int i=0; i<basis.size(); i++)
  {
    SUNDANCE_MSG2(verb_, "i=" << i 
      << " basis=" << basis[i]->description()
      << " filters=" << filters[i]);
  }

  RCP<DOFMapBase> rtn;

  if (allowNodalMap() && hasOmnipresentNodalMap(basis, mesh, filters))
  {
    SUNDANCE_MSG2(verb_, "creating omnipresent nodal map");
    CellFilter maxCells = getMaxCellFilter(filters);
    // if the mesh allows hanging nodes then create different DOF Map
    if (mesh.allowsHangingHodes()){
      rtn = rcp(new NodalDOFMapHN(mesh, basis.size(), maxCells, verb_));
    }
    else {
      rtn = rcp(new NodalDOFMap(mesh, basis.size(), maxCells, verb_));
    }
  }
  else if (hasNodalBasis(basis) && filtersAreZeroDimensional(mesh, filters))
  {
    SUNDANCE_MSG2(verb_, "creating submaximal nodal map");
    TEUCHOS_TEST_FOR_EXCEPT(filters.size() != 1);
    TEUCHOS_TEST_FOR_EXCEPT(filters[0].size() != 1);
    rtn = rcp(new SubmaximalNodalDOFMap(mesh, *filters[0].begin(), basis.size(), verb_));
  }
  else if (hasCellBasis(basis) && hasCommonDomain(filters))
  {
    TEUCHOS_TEST_FOR_EXCEPTION(filters[0].size() != 1, std::runtime_error,
      "only a single domain expected in construction of an element "
      "DOF map");
    rtn = rcp(new PartialElementDOFMap(mesh, *filters[0].begin(), basis.size(), verb_));
  }
  else if (allFuncsAreOmnipresent(mesh, filters))
  {
    SUNDANCE_MSG2(verb_, "creating omnipresent mixed map");
    CellFilter maxCells = getMaxCellFilter(filters);
    if (mesh.allowsHangingHodes()){
    	rtn = rcp(new MixedDOFMapHN(mesh, basis, maxCells, verb_));
    }else
    {
    	rtn = rcp(new MixedDOFMap(mesh, basis, maxCells, verb_));
    }
  }
  else if ( hasNodalBasis(basis) && (!mesh.allowsHangingHodes()) )
  {
    SUNDANCE_MSG2(verb_, "creating inhomogeneous nodal map");
    Sundance::Map<CellFilter, Set<int> > fmap = domainToFuncSetMap(filters);
    Sundance::Map<CellFilter, Sundance::Map<Set<int>, CellSet> > inputChildren;

    Array<Sundance::Map<Set<int>, CellFilter> > disjoint 
      = DOFMapBuilder::funcDomains(mesh, fmap, inputChildren);

    rtn = rcp(new InhomogeneousNodalDOFMap(mesh, disjoint, verb_));
  }
  else if ( hasEdgeLocalizedBasis(basis) && (!mesh.allowsHangingHodes()) )
  {
    SUNDANCE_MSG2(verb_, "creating inhomogeneous edge-localized map");

    Sundance::Map<CellFilter, Set<int> > fmap = domainToFuncSetMap(filters);
    
    Sundance::Map<CellFilter, Sundance::Map<Set<int>, CellSet> > inputChildren;
    Array<Sundance::Map<Set<int>, CellFilter> > disjoint 
      = DOFMapBuilder::funcDomains(mesh, fmap, inputChildren);

    rtn = rcp(new InhomogeneousEdgeLocalizedDOFMap(mesh, disjoint, verb_));
  }
  else
  {
	SUNDANCE_MSG2(verb_, "creating inhomogeneous HN map (the last possibility)");
	Sundance::Map<CellFilter, Set<int> > fmap = domainToFuncSetMap(filters);
	Sundance::Map<CellFilter, Sundance::Map<Set<int>, CellSet> > inputChildren;

	Array<Sundance::Map<Set<int>, CellFilter> > disjoint
	      = DOFMapBuilder::funcDomains(mesh, fmap, inputChildren);
    // the last option inhomogeneous mixed map, which supports hanging nodes
	rtn = rcp(new InhomogeneousDOFMapHN(mesh, basis , disjoint, verb_));
  }

  if (verb_ > 0)
  {
    Out::os() << "done building DOF map" << std::endl;
    if (verb_ > 1) Out::os() << "num DOFs" << rtn->numDOFs() << std::endl;
    if (verb_ > 1) Out::os() << "num local DOFs" 
                             << rtn->numLocalDOFs() << std::endl;
    if (verb_ > 4) rtn->print(Out::os());
  }
  return rtn;
}


Sundance::Map<CellFilter, Set<int> > DOFMapBuilder::domainToFuncSetMap(const Array<Set<CellFilter> >& filters) const 
{
  SUNDANCE_MSG2(verb_, "in DOFMapBuilder::domainToFuncSetMap()");
  Map<CellFilter, Set<int> > rtn;
  for (int i=0; i<filters.size(); i++)
  {
    const Set<CellFilter>& s = filters[i];
    for (Set<CellFilter>::const_iterator j=s.begin(); j!=s.end(); j++)
    {
      const CellFilter& cf = *j;
      if (rtn.containsKey(cf)) 
      {
        rtn[cf].put(i);
      }
      else
      {
              
        rtn.put(cf, makeSet((int) i));
      }
    }
  }
  for (Map<CellFilter, Set<int> >::const_iterator 
         i=rtn.begin(); i!=rtn.end(); i++)
  {
    SUNDANCE_MSG2(verb_, "subdomain=" << i->first << ", functions="
      << i->second);
  }
  return rtn;
}


void DOFMapBuilder
::getSubdomainUnkFuncMatches(const FunctionSupportResolver& fsr,
  Array<Sundance::Map<CellFilter, Set<int> > >& fmap) const 
{
  fmap.resize(fsr.numUnkBlocks());
  
  for (int r=0; r<fsr.numRegions(); r++)
  {
    CellFilter subreg = fsr.region(r);
    Set<int> funcs = fsr.unksOnRegion(r).setUnion(fsr.bcUnksOnRegion(r));
    for (Set<int>::const_iterator i=funcs.begin(); i!=funcs.end(); i++)
    {
      int block = fsr.blockForUnkID(*i);
      if (fmap[block].containsKey(subreg))
      {
        fmap[block][subreg].put(*i);
      }
      else
      {
        fmap[block].put(subreg, makeSet(*i));
      }
    }
  }
}

void DOFMapBuilder
::getSubdomainVarFuncMatches(const FunctionSupportResolver& fsr,
  Array<Sundance::Map<CellFilter, Set<int> > >& fmap) const 
{
  fmap.resize(fsr.numVarBlocks());
  
  for (int r=0; r<fsr.numRegions(); r++)
  {
    CellFilter subreg = fsr.region(r);
    Set<int> funcs = fsr.varsOnRegion(r).setUnion(fsr.bcVarsOnRegion(r));
    for (Set<int>::const_iterator i=funcs.begin(); i!=funcs.end(); i++)
    {
      int block = fsr.blockForVarID(*i);
      if (fmap[block].containsKey(subreg))
      {
        fmap[block][subreg].put(*i);
      }
      else
      {
        fmap[block].put(subreg, makeSet(*i));
      }
    }
  }
}

Array<Sundance::Map<Set<int>, CellFilter> > DOFMapBuilder
::funcDomains(const Mesh& mesh,
  const Sundance::Map<CellFilter, Set<int> >& fmap,
  Sundance::Map<CellFilter, Sundance::Map<Set<int>, CellSet> >& inputToChildrenMap) const 
{
  TimeMonitor timer(findFuncDomainTimer());
  Array<Array<CellFilter> > filters(mesh.spatialDim()+1);
  Array<Array<Set<int> > > funcs(mesh.spatialDim()+1);

  for (Sundance::Map<CellFilter, Set<int> >::const_iterator 
         i=fmap.begin(); i!=fmap.end(); i++)
  {
    int d = i->first.dimension(mesh);
    filters[d].append(i->first);
    funcs[d].append(i->second);
  }
  Array<Array<CFMeshPair> > tmp(mesh.spatialDim()+1);
  for (int d=0; d<tmp.size(); d++)
  {
    if (filters[d].size() != 0)
      tmp[d] = findDisjointFilters(filters[d], funcs[d], mesh);
  }

  for (int d=0; d<tmp.size(); d++)
  {
    for (int r=0; r<tmp[d].size(); r++)
    {
      for (int p=0; p<filters[d].size(); p++)
      {
        if (tmp[d][r].filter().isSubsetOf(filters[d][p], mesh)) 
        {
          if (inputToChildrenMap.containsKey(filters[d][p]))
          {
            Sundance::Map<Set<int>, CellSet>& m 
              = inputToChildrenMap[filters[d][p]];
            if (m.containsKey(tmp[d][r].funcs()))
            {
              m.put(tmp[d][r].funcs(), m[tmp[d][r].funcs()].setUnion(tmp[d][r].cellSet())); 
            }
            else
            {
              m.put(tmp[d][r].funcs(), tmp[d][r].cellSet()); 
            }
          }
          else
          {
            Sundance::Map<Set<int>, CellSet> m;
            m.put(tmp[d][r].funcs(), tmp[d][r].cellSet());
            inputToChildrenMap.put(filters[d][p], m);
          }
        }
      }
    }
  }

  Array<Sundance::Map<Set<int>, CellFilter> > rtn(mesh.spatialDim()+1);
  for (int d=0; d<tmp.size(); d++)
  {
    if (tmp[d].size() == 0) continue;
    for (int i=0; i<tmp[d].size(); i++)
    {
      const Set<int>& f = tmp[d][i].funcs();
      const CellFilter& cf = tmp[d][i].filter();
      if (rtn[d].containsKey(f))
      {
        rtn[d].put(f, rtn[d][f] + cf);
      }
      else
      {
        rtn[d].put(f, cf);
      }
    }
  }

  return rtn;
}


void DOFMapBuilder::init(bool findBCCols)
{
  SUNDANCE_MSG1(verb_, "in DOFMapBuilder::init()");
  SUNDANCE_MSG2(verb_, "num var blocks=" << fsr_->numVarBlocks());
  SUNDANCE_MSG2(verb_, "num unk blocks=" << fsr_->numUnkBlocks());

  rowMap_.resize(fsr_->numVarBlocks());
  colMap_.resize(fsr_->numUnkBlocks());
  isBCRow_.resize(fsr_->numVarBlocks());
  isBCCol_.resize(fsr_->numUnkBlocks());

  Array<Array<RCP<BasisDOFTopologyBase> > > testBasis = testBasisTopologyArray();
  Array<Array<Set<CellFilter> > > testRegions = testCellFilters();

  for (int br=0; br<fsr_->numVarBlocks(); br++)
  {
    SUNDANCE_MSG2(verb_, "making map for block row=" << br);
    rowMap_[br] = makeMap(mesh_, testBasis[br], testRegions[br]);
    SUNDANCE_MSG2(verb_, "marking BC rows for block row=" << br);
    markBCRows(br);
  }      


  Array<Array<RCP<BasisDOFTopologyBase> > > unkBasis = unkBasisTopologyArray();
  Array<Array<Set<CellFilter> > > unkRegions = unkCellFilters();

  for (int bc=0; bc<fsr_->numUnkBlocks(); bc++)
  {
    if (isSymmetric(bc))
    {
      colMap_[bc] = rowMap_[bc];
    }
    else
    {
      SUNDANCE_MSG2(verb_, "making map for block col=" << bc);
      colMap_[bc] = makeMap(mesh_, unkBasis[bc], unkRegions[bc]);
    }
    SUNDANCE_MSG2(verb_, "marking BC cols for block col=" << bc);
    if (findBCCols) markBCCols(bc);
  }
}

void DOFMapBuilder::extractUnkSetsFromFSR(const FunctionSupportResolver& fsr,
  Array<Set<int> >& funcSets,
  Array<CellFilter>& regions) const
{
  funcSets.resize(fsr.numRegions());
  regions.resize(fsr.numRegions());
  for (int r=0; r<fsr.numRegions(); r++)
  {
    regions[r] = fsr.region(r);
    funcSets[r] = fsr.unksOnRegion(r).setUnion(fsr.bcUnksOnRegion(r));
  }
}

void DOFMapBuilder::extractVarSetsFromFSR(const FunctionSupportResolver& fsr,
  Array<Set<int> >& funcSets,
  Array<CellFilter>& regions) const
{
  funcSets.resize(fsr.numRegions());
  regions.resize(fsr.numRegions());
  for (int r=0; r<fsr.numRegions(); r++)
  {
    regions[r] = fsr.region(r);
    funcSets[r] = fsr.varsOnRegion(r).setUnion(fsr.bcVarsOnRegion(r));
  }
}

Sundance::Map<Set<int>, Set<CellFilter> > 
DOFMapBuilder::buildFuncSetToCFSetMap(const Array<Set<int> >& funcSets,
  const Array<CellFilter>& regions,
  const Mesh& mesh) const 
{
  Sundance::Map<Set<int>, Set<CellFilter> > tmp;
  
  for (int r=0; r<regions.size(); r++)
  {
    const CellFilter& reg = regions[r];
    if (!tmp.containsKey(funcSets[r]))
    {
      tmp.put(funcSets[r], Set<CellFilter>());
    }
    tmp[funcSets[r]].put(reg);
  }
  
  /* eliminate overlap between cell filters */
  
  Sundance::Map<Set<int>, Set<CellFilter> > rtn=tmp;
  /*
    for (Sundance::Map<Set<int>, Set<CellFilter> >::const_iterator 
    i=tmp.begin(); i!=tmp.end(); i++)
    {
    rtn.put(i->first, reduceCellFilters(mesh, i->second));
    }
  */
  return rtn;
}

bool DOFMapBuilder::hasOmnipresentNodalMap(const Array<RCP<BasisDOFTopologyBase> >& basis,
  const Mesh& mesh,
  const Array<Set<CellFilter> >& filters) const 
{
  return hasNodalBasis(basis) && allFuncsAreOmnipresent(mesh, filters);
}
                                           

                                           
bool DOFMapBuilder::hasCommonDomain(const Array<Set<CellFilter> >& filters) const
{
  Set<CellFilter> first = filters[0];
  for (int i=1; i<filters.size(); i++) 
  {
    if (! (filters[i] == first) ) return false;
  }
  return true;
}                           

bool DOFMapBuilder::hasNodalBasis(const Array<RCP<BasisDOFTopologyBase> >& basis) const
{
  for (int i=0; i<basis.size(); i++)
  {
    const Lagrange* lagr 
      = dynamic_cast<const Lagrange*>(basis[i].get());
    if (lagr==0 || lagr->order()!=1) return false;
  }
  return true;
}

bool DOFMapBuilder::hasEdgeLocalizedBasis(const Array<RCP<BasisDOFTopologyBase> >& basis) const
{
  for (int i=0; i<basis.size(); i++)
  {
    const RCP<const EdgeLocalizedBasis> downcasted = rcp_dynamic_cast<EdgeLocalizedBasis>(basis[i]);
    if (is_null(downcasted)) return false;
  }
  return true;
}

bool DOFMapBuilder::hasCellBasis(const Array<RCP<BasisDOFTopologyBase> >& basis) const
{
  for (int i=0; i<basis.size(); i++)
  {
    const Lagrange* lagr 
      = dynamic_cast<const Lagrange*>(basis[i].get());
    if (lagr==0 || lagr->order()!=0) return false;
  }
  return true;
}

bool DOFMapBuilder::filtersAreZeroDimensional(const Mesh& mesh, 
  const Array<Set<CellFilter> >& filters) const
{
  for (int i=0; i<filters.size(); i++)
  {
    for (Set<CellFilter>::const_iterator 
           j=filters[i].begin(); j!=filters[i].end(); j++)
    {
      if (j->dimension(mesh) != 0) return false;
    }
  }
  return true;
}


bool DOFMapBuilder::allFuncsAreOmnipresent(const Mesh& mesh, 
  const Array<Set<CellFilter> >& filters) const
{
  int rtn = 0;

  int maxFilterDim = 0;
  Set<Set<CellFilter> > distinctSets;
  for (int i=0; i<filters.size(); i++)
  {
    for (Set<CellFilter>::const_iterator iter=filters[i].begin();
         iter != filters[i].end(); iter++)
    {
      int dim = iter->dimension(mesh);
      if (dim > maxFilterDim) maxFilterDim = dim;
    }
    distinctSets.put(filters[i]);
  }

  for (Set<Set<CellFilter> >::const_iterator 
         iter=distinctSets.begin(); iter != distinctSets.end(); iter++)
  {
    if (!isWholeDomain(mesh, maxFilterDim, *iter)) rtn = 1;
  }

  // make synchronization with the other processors
  int omniPresent = rtn;
  mesh.comm().allReduce((void*) &omniPresent, (void*) &rtn, 1,
    MPIDataType::intType(), MPIOp::sumOp());

  // it is true only when the summed value is zero
  return (rtn < 1);
}

bool DOFMapBuilder::isWholeDomain(const Mesh& mesh, 
  int maxFilterDim,
  const Set<CellFilter>& filters) const
{
  CellFilter allMax;
  if (maxFilterDim==mesh.spatialDim()) allMax = new MaximalCellFilter();
  else allMax = new DimensionalCellFilter(maxFilterDim);

  CellSet remainder = allMax.getCells(mesh);

  for (Set<CellFilter>::const_iterator 
         i=filters.begin(); i!=filters.end(); i++)
  {
    const CellFilter& cf = *i;
    if (maxFilterDim==mesh.spatialDim())
    {
      if (0 != dynamic_cast<const MaximalCellFilter*>(cf.ptr().get()))
      {
        return true;
      }
    }
    else
    {
      const DimensionalCellFilter* dcf 
        = dynamic_cast<const DimensionalCellFilter*>(cf.ptr().get());
      if (0 != dcf && dcf->dimension(mesh) == maxFilterDim) 
      {
        return true;
      }
    }
    if (cf.dimension(mesh) != maxFilterDim) continue;
    CellSet cells = cf.getCells(mesh);
    SUNDANCE_MSG2(verb_, "found " << cells.numCells() << " in CF " << cf);
    remainder = remainder.setDifference(cells);
    if (remainder.begin() == remainder.end()) return true;
  }

  SUNDANCE_MSG2(verb_, "num remaining cells: " << remainder.numCells());

  return false;
}



CellFilter DOFMapBuilder::getMaxCellFilter(const Array<Set<CellFilter> >& filters) const
{
  for (int i=0; i<filters.size(); i++)
  {
    const Set<CellFilter>& cfs = filters[i];
    if (cfs.size() != 1) continue;
    const CellFilter& cf = *cfs.begin();
    if (0 != dynamic_cast<const MaximalCellFilter*>(cf.ptr().get()))
      return cf;
  }
//  TEUCHOS_TEST_FOR_EXCEPT(true);
  return new MaximalCellFilter();
}



Array<Array<Set<CellFilter> > > DOFMapBuilder::testCellFilters() const
{
  Array<Array<Set<CellFilter> > > rtn(fsr_->numVarBlocks());

  for (int b=0; b<fsr_->numVarBlocks(); b++)
  {
    for (int i=0; i<fsr_->numVarIDs(b); i++) 
    {
      int testID = fsr_->unreducedVarID(b, i);
      Set<CellFilter> s;
      const Set<OrderedHandle<CellFilterStub> >& cfs 
        = fsr_->regionsForTestFunc(testID);
      for (Set<OrderedHandle<CellFilterStub> >::const_iterator 
             j=cfs.begin(); j!=cfs.end(); j++)
      {
        RCP<CellFilterBase> cfb 
          = rcp_dynamic_cast<CellFilterBase>(j->ptr());
        TEUCHOS_TEST_FOR_EXCEPT(cfb.get()==0);
        CellFilter cf = j->ptr();
        s.put(cf);
      }
      //         Set<CellFilter> reducedS = reduceCellFilters(mesh(), s);
//          rtn[b].append(reducedS);
      rtn[b].append(s);
    }
  }
  return rtn;
}

Array<Array<Set<CellFilter> > > DOFMapBuilder::unkCellFilters() const
{
  Array<Array<Set<CellFilter> > > rtn(fsr_->numUnkBlocks());

  for (int b=0; b<fsr_->numUnkBlocks(); b++)
  {
    for (int i=0; i<fsr_->numUnkIDs(b); i++) 
    {
      int unkID = fsr_->unreducedUnkID(b, i);
      Set<CellFilter> s;
      const Set<OrderedHandle<CellFilterStub> >& cfs 
        = fsr_->regionsForUnkFunc(unkID);
      for (Set<OrderedHandle<CellFilterStub> >::const_iterator 
             j=cfs.begin(); j!=cfs.end(); j++)
      {
        RCP<CellFilterBase> cfb 
          = rcp_dynamic_cast<CellFilterBase>(j->ptr());
        TEUCHOS_TEST_FOR_EXCEPT(cfb.get()==0);
        CellFilter cf = j->ptr();
        s.put(cf);
      }
//          Set<CellFilter> reducedS = reduceCellFilters(mesh(), s);
//          rtn[b].append(reducedS);
      rtn[b].append(s);
    }
  }
  return rtn;
}

Array<Array<RCP<BasisDOFTopologyBase> > > DOFMapBuilder::testBasisTopologyArray() const 
{
  Array<Array<RCP<BasisDOFTopologyBase> > > rtn(fsr_->numVarBlocks());
  for (int b=0; b<fsr_->numVarBlocks(); b++)
  {
    for (int i=0; i<fsr_->numVarIDs(b); i++) 
    {
      rtn[b].append(BasisFamily::getBasisTopology(fsr_->varFuncData(b, i)));
    }
  }
  return rtn;
}

Array<Array<RCP<BasisDOFTopologyBase> > > DOFMapBuilder::unkBasisTopologyArray() const 
{
  Array<Array<RCP<BasisDOFTopologyBase> > > rtn(fsr_->numUnkBlocks());
  for (int b=0; b<fsr_->numUnkBlocks(); b++)
  {
    for (int i=0; i<fsr_->numUnkIDs(b); i++) 
    {
      rtn[b].append(BasisFamily::getBasisTopology(fsr_->unkFuncData(b, i)));
    }
  }
  return rtn;
}


Set<CellFilter> DOFMapBuilder
::reduceCellFilters(const Mesh& mesh, 
  const Set<CellFilter>& inputSet) const
{
  TimeMonitor timer(cellFilterReductionTimer());
  Set<CellFilter> rtn;
  /* If the input set explicitly contains all maximal cells, we're done */
  CellFilter m = new MaximalCellFilter();
  if (inputSet.contains(m))
  {
    rtn.put(m);
    return rtn;
  }

  /* Next, see if combining the maximal-dimension filters in the
   * input set gives us all maximal cells. */
  CellFilter myMaxFilters;
  for (Set<CellFilter>::const_iterator 
         i=inputSet.begin(); i!=inputSet.end(); i++)
  {
    CellFilter f = *i;
    if (f.dimension(mesh) != mesh.spatialDim()) continue;
    myMaxFilters = myMaxFilters + f;
  }

  CellSet myMax = myMaxFilters.getCells(mesh);
//  if (myMax.size() == mesh.numCells(mesh.spatialDim()))
//  {
//    rtn.put(m);
//    return rtn;
//  }

  CellSet allMax = m.getCells(mesh);
  CellSet diff = allMax.setDifference(myMax);
  /* if the difference between the collected max cell set and the known
   * set of all max cells is empty, then we're done */
  if (diff.begin() == diff.end())
  {
    rtn.put(m);
    return rtn;
  }
  
  /* Otherwise, we return the remaining max cell filters, and possibly
   * some lower-dimensional filters to be identified below. */
  if (myMax.begin() != myMax.end()) rtn.put(myMaxFilters);

  /* At this point, we can eliminate as redundant any lower-dimensional
   * cell filters all of whose cells are facets of our subset of
   * maximal cell filters. Any other lower-dim cell filters must be 
   * appended to our list. */
  for (Set<CellFilter>::const_iterator 
         i=inputSet.begin(); i!=inputSet.end(); i++)
  {
    CellFilter f = *i;
    if (f.dimension(mesh) == mesh.spatialDim()) continue;
    CellSet s = f.getCells(mesh);
    if (s.areFacetsOf(myMax)) continue;
    /* if we're here, then we have a lower-dimensional cell filter
     * whose cells are not facets of cells in our maximal cell filters.
     * These will need to be processed separately in assigning DOFs.
     */
    rtn.put(f);
  }
  return rtn;
}


bool DOFMapBuilder::isSymmetric(int b) const 
{
  if ((int)fsr_->numVarBlocks() < b || (int)fsr_->numUnkBlocks() < b) return false;

  if (fsr_->numVarIDs(b) != fsr_->numUnkIDs(b)) return false;

  for (int i=0; i<fsr_->numVarIDs(b); i++) 
  {
    BasisFamily basis1 = BasisFamily::getBasis(fsr_->varFuncData(b,i));
    BasisFamily basis2 = BasisFamily::getBasis(fsr_->unkFuncData(b,i));
    if (!(basis1 == basis2)) return false;
  }
  return true;
}

bool DOFMapBuilder::regionIsMaximal(int r) const 
{
  const CellFilterStub* reg = fsr_->region(r).get();
  return (dynamic_cast<const MaximalCellFilter*>(reg) != 0);
}

void DOFMapBuilder::markBCRows(int block)
{
  isBCRow_[block] = rcp(new Array<int>(rowMap_[block]->numLocalDOFs()));
  int ndof = rowMap_[block]->numLocalDOFs();
  Array<int>& isBC = *isBCRow_[block];
  for (int i=0; i<ndof; i++) isBC[i] = false;

  RCP<Array<int> > cellLID = rcp(new Array<int>());
  Array<RCP<Array<int> > > cellBatches;
  const RCP<DOFMapBase>& rowMap = rowMap_[block];

  for (int r=0; r<fsr_->numRegions(); r++)
  {
    /* find the cells in this region */
    CellFilter region = fsr_->region(r);

    if (!fsr_->isBCRegion(r)) continue;

    int dim = region.dimension(mesh_);
    CellSet cells = region.getCells(mesh_);
    cellLID->resize(0);
    for (CellIterator c=cells.begin(); c != cells.end(); c++)
    {
      cellLID->append(*c);
    }
    if (cellLID->size() == 0) continue;
      
    /* find the functions that appear in BCs on this region */
    const Set<int>& allBcFuncs = fsr_->bcVarsOnRegion(r);

    Set<int> bcFuncs;
    for (Set<int>::const_iterator 
           i=allBcFuncs.begin(); i != allBcFuncs.end(); i++)
    {
      if (block == fsr_->blockForVarID(*i)) 
      {
        bcFuncs.put(fsr_->reducedVarID(*i));
      }
    }
    if (bcFuncs.size()==0) continue;
    Array<int> bcFuncID = bcFuncs.elements();

    Array<Array<int> > dofs;
    Array<int> nNodes;

    RCP<const MapStructure> s 
      = rowMap->getDOFsForCellBatch(dim, *cellLID, bcFuncs, dofs, nNodes,0);
    int offset = rowMap->lowestLocalDOF();
    int high = offset + rowMap->numLocalDOFs();
      
    for (int c=0; c<cellLID->size(); c++)
    {
      for (int b=0; b< s->numBasisChunks(); b++)
      {
        int nFuncs = s->numFuncs(b);
        for (int n=0; n<nNodes[b]; n++)
        {
          for (int f=0; f<bcFuncID.size(); f++)
          {
            int chunk = s->chunkForFuncID(bcFuncID[f]);
            if (chunk != b) continue;
            int funcOffset = s->indexForFuncID(bcFuncID[f]);
            int dof = dofs[b][(c*nFuncs + funcOffset)*nNodes[b]+n];
            if (dof < offset || dof >= high) continue;
            (*isBCRow_[block])[dof-offset]=true;
          }
        }
      }
    }
  }
}
        


void DOFMapBuilder::markBCCols(int block)
{
  isBCCol_[block] = rcp(new Array<int>(colMap_[block]->numLocalDOFs()));
  int ndof = colMap_[block]->numLocalDOFs();
  Array<int>& isBC = *isBCCol_[block];
  for (int i=0; i<ndof; i++) isBC[i] = false;

  RCP<Array<int> > cellLID = rcp(new Array<int>());
  Array<RCP<Array<int> > > cellBatches;
  const RCP<DOFMapBase>& colMap = colMap_[block];

  for (int r=0; r<fsr_->numRegions(); r++)
  {
    /* find the cells in this region */
    CellFilter region = fsr_->region(r);

    if (!fsr_->isBCRegion(r)) continue;

    int dim = region.dimension(mesh_);
    CellSet cells = region.getCells(mesh_);
    cellLID->resize(0);
    for (CellIterator c=cells.begin(); c != cells.end(); c++)
    {
      cellLID->append(*c);
    }
    if (cellLID->size() == 0) continue;
      
    /* find the functions that appear in BCs on this region */
    const Set<int>& allBcFuncs = fsr_->bcUnksOnRegion(r);

    Set<int> bcFuncs;
    for (Set<int>::const_iterator 
           i=allBcFuncs.begin(); i != allBcFuncs.end(); i++)
    {
      if (block == fsr_->blockForUnkID(*i)) 
      {
        bcFuncs.put(fsr_->reducedUnkID(*i));
      }
    }
    if (bcFuncs.size()==0) continue;
    Array<int> bcFuncID = bcFuncs.elements();

    Array<Array<int> > dofs;
    Array<int> nNodes;

    RCP<const MapStructure> s 
      = colMap->getDOFsForCellBatch(dim, *cellLID, bcFuncs, dofs, nNodes,0);
    int offset = colMap->lowestLocalDOF();
    int high = offset + colMap->numLocalDOFs();
      
    for (int c=0; c<cellLID->size(); c++)
    {
      for (int b=0; b< s->numBasisChunks(); b++)
      {
        int nFuncs = s->numFuncs(b);
        for (int n=0; n<nNodes[b]; n++)
        {
          for (int f=0; f<bcFuncID.size(); f++)
          {
            int chunk = s->chunkForFuncID(bcFuncID[f]);
            if (chunk != b) continue;
            int funcOffset = s->indexForFuncID(bcFuncID[f]);
            int dof = dofs[b][(c*nFuncs + funcOffset)*nNodes[b]+n];
            if (dof < offset || dof >= high) 
            {
              remoteBCCols_[block]->insert(dof);
            }
            else
            {
              (*isBCCol_[block])[dof-offset]=true;
            }
          }
        }
      }
    }
  }
}



namespace Sundance
{

Array<Array<BasisFamily> > testBasisArray(const RCP<FunctionSupportResolver>& fsr) 
{
  Array<Array<BasisFamily> > rtn(fsr->numVarBlocks());
  for (int b=0; b<fsr->numVarBlocks(); b++)
  {
    for (int i=0; i<fsr->numVarIDs(b); i++) 
    {
      rtn[b].append(BasisFamily::getBasis(fsr->varFuncData(b, i)));
    }
  }
  return rtn;
}


Array<Array<BasisFamily> > unkBasisArray(const RCP<FunctionSupportResolver>& fsr) 
{
  Array<Array<BasisFamily> > rtn(fsr->numUnkBlocks());
  for (int b=0; b<fsr->numUnkBlocks(); b++)
  {
    for (int i=0; i<fsr->numUnkIDs(b); i++) 
    {
      rtn[b].append(BasisFamily::getBasis(fsr->unkFuncData(b, i)));
    }
  }
  return rtn;
}


}        
