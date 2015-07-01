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

#include "SundanceExodusWriter.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceVertexSort.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_XMLObject.hpp"

#ifdef HAVE_SUNDANCE_EXODUS 
#include "exodusII.h"
#endif


using namespace Sundance;
using namespace Teuchos;
using namespace std;


void ExodusWriter::write() const 
{
//  Out::os() << "in ExodusWriter::write()" << endl;
  std::string exoFile = filename();
  std::string parFile = filename();
  if (nProc() > 1) 
  {
    exoFile = exoFile + "-" + Teuchos::toString(nProc()) + "-" + Teuchos::toString(myRank());
    parFile = parFile + "-" + Teuchos::toString(nProc()) + "-" + Teuchos::toString(myRank());
  }
  exoFile = exoFile + ".exo";
  parFile = parFile + ".pxo";

  if (nProc() > 1) writeParallelInfo(parFile);
#ifdef HAVE_SUNDANCE_EXODUS
  int ws = 8;
  int exoid = ex_create(exoFile.c_str(), EX_CLOBBER, &ws, &ws);

  TEUCHOS_TEST_FOR_EXCEPTION(exoid < 0, std::runtime_error, "failure to create file "
    << filename());


  Array<CellFilter> nsFilters;
  Array<int> omniNodalFuncs;
  Array<RCP<Array<int> > > funcsForNodeset;
  Array<RCP<Array<int> > > nodesForNodeset;
  Array<int> nsID;
  Array<int> nsNodesPerSet;
  Array<int> nsNodePtr;
  RCP<Array<int> > allNodes=rcp(new Array<int>());

  Array<CellFilter> blockFilters;
  Array<int> omniElemFuncs;
  Array<RCP<Array<int> > > funcsForBlock;
  Array<RCP<Array<int> > > elemsForBlock;
  Array<int> blockID;
  Array<int> nElemsPerBlock;
  Array<int> blockElemPtr;
  RCP<Array<int> > allElems=rcp(new Array<int>());

  
  findNodeSets(nsFilters, omniNodalFuncs, funcsForNodeset,
    nodesForNodeset, nsID, nsNodesPerSet, nsNodePtr, allNodes);

  findBlocks(blockFilters, omniElemFuncs, funcsForBlock,
    elemsForBlock, blockID, nElemsPerBlock, blockElemPtr, allElems);
  


  writeMesh(exoid, nsFilters, nsID, nsNodesPerSet, nsNodePtr, allNodes );
  
  writeFields(exoid, nsFilters, omniNodalFuncs, omniElemFuncs, 
    funcsForNodeset,
    nodesForNodeset, nsID);


  ex_close(exoid);
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Exodus not enabled");
#endif
}


void ExodusWriter::offset(Array<int>& x) const
{
  for (int i=0; i<x.size(); i++) x[i]++;
}


void ExodusWriter::writeMesh(int exoid, 
  const Array<CellFilter>& nodesetFilters,
  const Array<int>& nsID,
  const Array<int>& nNodesPerSet,
  const Array<int>& nsNodePtr,
  const RCP<Array<int> >& allNodes) const
{
#ifdef HAVE_SUNDANCE_EXODUS
//  Out::os() << "in ExodusWriter::writeMesh()" << endl;
  int ierr = 0;

  int dim = mesh().spatialDim();
  int nElems = mesh().numCells(dim);
  
  Array<int> ssLabels = mesh().getAllLabelsForDimension(dim-1).elements();
  int numSS = 0;

  for (int ss=0; ss<ssLabels.size(); ss++) 
  {
    if (ssLabels[ss] != 0) numSS++;
  }
  
  int numNS = nsID.size();

  
  /* initialize the output file */
  int nElemBlocks = mesh().numLabels(dim);

  ierr = ex_put_init(
    exoid, 
    filename().c_str(), 
    dim,
    mesh().numCells(0), 
    nElems,
    nElemBlocks, 
    numNS,
    numSS
    );
  TEUCHOS_TEST_FOR_EXCEPT(ierr<0);


  /* write the vertices */
  int nPts = mesh().numCells(0);

  Array<double> x(nPts);
  Array<double> y(nPts);
  Array<double> z(nPts * (dim > 2));

  for (int n=0; n<nPts; n++)
  {
    const Point& P = mesh().nodePosition(n);
    x[n] = P[0];
    y[n] = P[1];
    if (dim==3) z[n] = P[2];
  }

  if (dim==2)
  {
    ierr = ex_put_coord(exoid, &(x[0]), &(y[0]), (void*) 0);
  }
  else
  {
    ierr = ex_put_coord(exoid, &(x[0]), &(y[0]), &(z[0]));
  }

  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);

  if (dim==2)
  {
    Array<std::string> cn;
    cn.append("x");
    cn.append("y");
    Array<const char*> pp;
    getCharpp(cn, pp);
    ierr = ex_put_coord_names(exoid,(char**) &(pp[0]));
  }
  else
  {
    Array<std::string> cn;
    cn.append("x");
    cn.append("y");
    cn.append("z");
    Array<const char*> pp;
    getCharpp(cn, pp);
    ierr = ex_put_coord_names(exoid, (char**)&(pp[0]));
  }

  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);

  /* write the element blocks */

  Array<int> blockLabels = mesh().getAllLabelsForDimension(dim).elements();
  int nodesPerElem = dim+1;
  std::string eType = elemType(mesh().cellType(dim));
  
  bool minBlockIsZero = false;
  for (int b=0; b<blockLabels.size(); b++)
  {
    if (blockLabels[b] < 1) {minBlockIsZero = true; break;}
  }
  

  for (int b=0; b<blockLabels.size(); b++)
  {
//    Out::os() << "writing element block" << b << endl;
    int numBlockAttr = 0;
    Array<int> blockElemLIDs;
    Array<int> nodeLIDs;
    Array<int> orient;
    mesh().getLIDsForLabel(dim, blockLabels[b], blockElemLIDs);
    int numElemsThisBlock = blockElemLIDs.size();
//    Out::os() << "writing block #" << blockLabels[b] << " with " << numElemsThisBlock << " elems"  
//              << endl;
    mesh().getFacetLIDs(dim, blockElemLIDs, 0, nodeLIDs, orient);
    offset(nodeLIDs);
    ierr = ex_put_elem_block(
      exoid, blockLabels[b]+minBlockIsZero, eType.c_str(), 
      numElemsThisBlock, nodesPerElem, numBlockAttr
      );

    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
//    Out::os() << "block: " << blockLabels[b]+minBlockIsZero << endl;
    ierr = ex_put_elem_conn(exoid, blockLabels[b]+minBlockIsZero, 
      &(nodeLIDs[0]));
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  }
//  Out::os() << "done all element blocks" << endl;

  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  
  /* write the side sets */
  
//  Out::os() << "writing side sets "<< endl;  
  for (int ss=0; ss<ssLabels.size(); ss++)
  {
//    Out::os() << "writing side set=" << ss << " of " << ssLabels.size() << endl;
    if (ssLabels[ss]==0) continue;
    Array<int> sideLIDs;
    RCP<Array<int> > elemLIDs = rcp(new Array<int>());
    RCP<Array<int> > facets = rcp(new Array<int>());
    MaximalCofacetBatch maxCofacetBatch;

//    Out::os() << "getting LIDs "<< endl;  
    mesh().getLIDsForLabel(dim-1, ssLabels[ss], sideLIDs);
    sort(&(sideLIDs[0]), &(sideLIDs[sideLIDs.size()-1]));
//    Out::os() << "LIDs= "<< sideLIDs << endl;  
//    Out::os() << "getting cofacets "<< endl;  
    mesh().getMaxCofacetLIDs(sideLIDs, maxCofacetBatch);
    //Out::os() << "getting specified cofacets "<< endl;  
    maxCofacetBatch.getSpecifiedCofacets(0, elemLIDs, facets);

    int numSides = sideLIDs.size();
    int numDists = 0;

    for (int i=0; i<elemLIDs->size(); i++)
    {
      (*facets)[i] = ufcFacetIndexToExFacetIndex(dim,  (*facets)[i]);
    }

    offset(sideLIDs);
    offset(*elemLIDs);
    offset(*facets);


    ierr = ex_put_side_set_param(exoid, ssLabels[ss], numSides, numDists);
    ierr = ex_put_side_set(exoid, ssLabels[ss], &((*elemLIDs)[0]), &((*facets)[0]));
  }
//  Out::os() << "done all side sets" << endl;
  
  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);


//  Out::os() << "writing node sets "<< endl;  
  if (nsID.size() > 0)
  {
    /* write the node sets */

    Array<int> nsDistPerSet(nsID.size(), 0);
    Array<int> nsDistPtr(nsID.size(), 0);
    Array<int> emptyDist(1, 0);

    offset(*allNodes);
  
    ierr = ex_put_concat_node_sets( exoid,
      (int*) &(nsID[0]),
      (int*) &(nNodesPerSet[0]),
      &(nsDistPerSet[0]),
      (int*) &(nsNodePtr[0]),
      &(nsDistPtr[0]),
      &((*allNodes)[0]),
      &(emptyDist[0]));

    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  }
//  Out::os() << "done all node sets" << endl;
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Exodus not enabled");
#endif
}


void ExodusWriter::writeFields(int exoid, 
  const Array<CellFilter>& nodesetFilters,
  const Array<int>& omnipresentNodalFuncs,
  const Array<int>& omnipresentElemFuncs,
  const Array<RCP<Array<int> > >& funcsForNodeset,
  const Array<RCP<Array<int> > >& nodesForNodeset,
  const Array<int>& nsID) const 
{

#ifdef HAVE_SUNDANCE_EXODUS
  Tabs tab0(0);
  Tabs tab1;
  int verb=3;

  
  int nNodalFuncs = omnipresentNodalFuncs().size();
  int nElemFuncs = omnipresentElemFuncs().size();

  int nNodesetFuncs = pointScalarFields().size() - nNodalFuncs;

  int nNodesets = funcsForNodeset.size();
  
  PLAYA_ROOT_MSG1(verb, tab0 << "ExodusWriter::writeFields()");
  PLAYA_ROOT_MSG2(verb, tab1 << "nNodalFuncs = " << nNodalFuncs);
  PLAYA_ROOT_MSG2(verb, tab1 << "nElemFuncs = " << nElemFuncs);
  PLAYA_ROOT_MSG2(verb, tab1 << "nNodesetFuncs = " << nNodesetFuncs);
  PLAYA_ROOT_MSG2(verb, tab1 << "nNodesets = " << nNodesets);


  Set<int> nsFuncSet;
  Map<int, Array<int> > funcToNSMap;

  for (int i=0; i<nNodesets; i++)
  {
    const Array<int>& f = *(funcsForNodeset[i]);
    for (int j=0; j<f.size(); j++)
    {
      nsFuncSet.put(f[j]);
      if (funcToNSMap.containsKey(f[j]))
      {
        funcToNSMap[f[j]].append(i);
      }
      else
      {
        funcToNSMap.put(f[j], tuple(i));
      }
    }
  }
  Array<int> nsFuncs = nsFuncSet.elements();
  TEUCHOS_TEST_FOR_EXCEPT(nsFuncs.size() != nNodesetFuncs);

  Map<int, int > funcIDToNSFuncIndex;
  for (int i=0; i<nNodesetFuncs; i++) funcIDToNSFuncIndex.put(nsFuncs[i],i);

  Array<Array<int> > nsFuncNodesets(nsFuncs.size());
  for (int i=0; i<nNodesetFuncs; i++)
  {
    nsFuncNodesets[i] = funcToNSMap.get(nsFuncs[i]);
  }

  Array<int> nodesetFuncTruthTable(nNodesetFuncs * nNodesets, 0);
  for (int i=0; i<nNodesetFuncs; i++)
  {
    for (int j=0; j<nsFuncNodesets[i].size(); j++)
    {
      int ns = nsFuncNodesets[i][j];
      nodesetFuncTruthTable[ns*nNodesetFuncs + i] = 1;
    }

    nsFuncNodesets[i] = funcToNSMap.get(nsFuncs[i]);
  }

  

  int ierr;

  if (nNodalFuncs > 0)
  {
    ierr = ex_put_var_param(exoid, "N", nNodalFuncs);
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  }

  if (nElemFuncs > 0)
  {
    ierr = ex_put_var_param(exoid, "E", nElemFuncs);
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  }


  if (nNodesets > 0)
  {
    ierr = ex_put_var_param(exoid, "M", nNodesetFuncs);
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    
    ierr = ex_put_nset_var_tab(exoid, nNodesets, 
      nNodesetFuncs, &(nodesetFuncTruthTable[0]));
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    
    Array<std::string> nsFuncNames(nNodesetFuncs);
    Array<const char*> nsNameP;
    
    for (int i=0; i<nNodesetFuncs; i++)
    {
      nsFuncNames[i] = pointScalarNames()[nsFuncs[i]];
    }
    getCharpp(nsFuncNames, nsNameP);  
    
    ierr = ex_put_var_names(exoid, "M", nNodesetFuncs, (char**)&(nsNameP[0]));
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  }


  Array<std::string> nodalFuncNames(nNodalFuncs);
  Array<std::string> elemFuncNames(nElemFuncs);
  Array<const char*> nNameP;
  Array<const char*> eNameP;
  
  for (int i=0; i<nNodalFuncs; i++)
  {
    nodalFuncNames[i] = pointScalarNames()[omnipresentNodalFuncs[i]];
  }
  getCharpp(nodalFuncNames, nNameP);  
  
  for (int i=0; i<nElemFuncs; i++)
  {
    elemFuncNames[i] = cellScalarNames()[omnipresentElemFuncs[i]];
  }
  getCharpp(elemFuncNames, eNameP);  
  
  if (nNodalFuncs > 0)
  {
    ierr = ex_put_var_names(exoid, "N", nNodalFuncs, (char**)&(nNameP[0]));
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    
    Array<double> funcVals;
    Array<int> nodeID(mesh().numCells(0));
    for (int i=0; i<mesh().numCells(0); i++) nodeID[i]=i;
    
    for (int i=0; i<nNodalFuncs; i++)
    {
      int f = omnipresentNodalFuncs[i];
      pointScalarFields()[f]->getDataBatch(0, nodeID, tuple(f), funcVals);
      int t = 1;
      int numNodes = funcVals.size();
      ierr = ex_put_nodal_var(exoid, t, i+1, numNodes, &(funcVals[0]));
      TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    }
    
    for (int i=0; i<nNodesetFuncs; i++)
    {
      const Array<int>& ns = nsFuncNodesets[i];
      int fid = nsFuncs[i];
      
      for (int s=0; s<ns.size(); s++)
      {
        const Array<int>& nodes = *(nodesForNodeset[ns[s]]);
        pointScalarFields()[fid]->getDataBatch(0, nodes, tuple(fid), funcVals);
        int t = 1;
        int numNodes = funcVals.size();
        int id = nsID[ns[s]];
        ierr = ex_put_nset_var(exoid, t, i+1, id, numNodes, &(funcVals[0]));
        TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
      }
    }
  }

  if (nElemFuncs > 0)
  {
    ierr = ex_put_var_names(exoid, "E", nElemFuncs, (char**)&(eNameP[0]));
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    
    Array<double> funcVals;
    int dim = mesh().spatialDim();
    Array<int> elemID(mesh().numCells(dim));
    for (int i=0; i<mesh().numCells(dim); i++) elemID[i]=i;
    
    for (int i=0; i<nElemFuncs; i++)
    {
      int f = omnipresentElemFuncs[i];
      cellScalarFields()[f]->getDataBatch(dim, elemID, tuple(f), funcVals);
      int t = 1;
      int numElems = funcVals.size();
      ierr = ex_put_elem_var(exoid, t, i+1, 1, numElems, &(funcVals[0]));
      TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    }
  }


#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Exodus not enabled");
#endif
  
}



std::string ExodusWriter::elemType(const CellType& type) const
{
  switch(type)
  {
    case TriangleCell:
      return "TRIANGLE";
    case TetCell:
      return "TETRA";
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "cell type=" << type << " cannot be used as a "
        "maximal-dimension cell in exodus");
  }
  return "NULL"; //-Wall
}



void ExodusWriter::writeParallelInfo(const std::string& parfile) const 
{
  std::ofstream os(parfile.c_str());

  int dim = mesh().spatialDim();
  int nCells = mesh().numCells(dim);
  int nPts = mesh().numCells(0);

  os << myRank() << " " << nProc() << std::endl;

  os << nPts << std::endl;
  for (int i=0; i<nPts; i++)
  {
    os << i << " " << mesh().mapLIDToGID(0,i) 
       << " " << mesh().ownerProcID(0,i) << std::endl;
  }

  os << nCells << std::endl;
  for (int c=0; c<nCells; c++)
  {
    os << c << " " << mesh().mapLIDToGID(dim,c) 
       << " " << mesh().ownerProcID(dim,c) << std::endl;
  }

  for (int i=0; i<comments().length(); i++)
  {
    os << "# " << comments()[i] << std::endl;
  }
}





void ExodusWriter::findNodeSets(
  Array<CellFilter>& nodesetFilters,
  Array<int>& omnipresentFuncs,
  Array<RCP<Array<int> > >& funcsForNodeset,
  Array<RCP<Array<int> > >& nodesForNodeset,
  Array<int>& nsID,
  Array<int>& nNodesPerSet,
  Array<int>& nsNodePtr,
  RCP<Array<int> > allNodes
  ) const 
{
//  Out::os() << "in ExodusWriter::findNodeSets()" << endl;
  int verb = 2;

  const Array<RCP<FieldBase> >& f = pointScalarFields();
  CellFilter maximal = new MaximalCellFilter();

  nNodesPerSet.resize(0);
  nsNodePtr.resize(0);
  nsID.resize(0);

  Map<CellFilter, RCP<Array<int> > > tmp;

  for (int i=0; i<f.size(); i++)
  {
    const CellFilter& cf = f[i]->domain(); 
    if (cf==maximal) 
    {
      SUNDANCE_MSG2(verb, "function #" << i << " is defined on all nodes");
      omnipresentFuncs.append(i);
      continue;
    }
    if (!tmp.containsKey(cf))
    {
      RCP<Array<int> > a = rcp(new Array<int>());
      tmp.put(cf, a);
    }
    SUNDANCE_MSG2(verb, "function #" << i << " is defined on CF " << cf);
    tmp[cf]->append(i);
  }

  int nodesetID=1;
  int nodePtr=0;
  nodesetFilters.resize(0);
  funcsForNodeset.resize(0);
  nodesForNodeset.resize(0);

  for (Map<CellFilter, RCP<Array<int> > >::const_iterator
         i=tmp.begin(); i!=tmp.end(); i++)
  {
    const CellFilter& cf = i->first;
    nodesetFilters.append(cf);
    funcsForNodeset.append(i->second);
    RCP<Array<int> > cells 
      = cellSetToLIDArray(connectedNodeSet(cf, mesh()));
    nodesForNodeset.append(cells);
    int nn = cells->size();
    nNodesPerSet.append(nn);
    nsID.append(nodesetID++);
    nsNodePtr.append(nodePtr);
    nodePtr += nn;
  }

  SUNDANCE_MSG2(verb, "node set IDs = " << nsID);
  SUNDANCE_MSG2(verb, "num nodes = " << nNodesPerSet);
  SUNDANCE_MSG2(verb, "node set pointers = " << nsNodePtr);


  int numNodes = nodePtr;
  allNodes->resize(numNodes);

  int k=0;
  for (int i=0; i<nsID.size(); i++)
  {
    SUNDANCE_MSG4(verb, "node set " << i << " funcs = " 
      << *funcsForNodeset[i]);
    SUNDANCE_MSG4(verb, "node set " << i 
      << " nodes = " << *nodesForNodeset[i]);
    const Array<int>& myCells = *(nodesForNodeset[i]);
    for (int c=0; c<myCells.size(); c++)
    {
      (*allNodes)[k++] = myCells[c];
    }
  }

  SUNDANCE_MSG4(verb, "all nodes = " << *allNodes);
}



void ExodusWriter::findBlocks(
  Array<CellFilter>& blockFilters,
  Array<int>& omnipresentFuncs,
  Array<RCP<Array<int> > >& funcsForBlock,
  Array<RCP<Array<int> > >& elemsForBlock,
  Array<int>& blockIDs,
  Array<int>& nElemsPerBlock,
  Array<int>& elemBlockPtr,
  RCP<Array<int> > allElems
  ) const 
{
//  Out::os() << "in ExodusWriter::findBlocks()" << endl;
  int verb=0;
  const Array<RCP<FieldBase> >& f = cellScalarFields();
  CellFilter maximal = new MaximalCellFilter();

  nElemsPerBlock.resize(0);
  elemBlockPtr.resize(0);
  blockIDs.resize(0);

  Map<CellFilter, RCP<Array<int> > > tmp;

  for (int i=0; i<f.size(); i++)
  {
    const CellFilter& cf = f[i]->domain(); 
    if (cf==maximal) 
    {
      SUNDANCE_MSG2(verb, "function #" << i << " is defined on all nodes");
      omnipresentFuncs.append(i);
      continue;
    }
    if (!tmp.containsKey(cf))
    {
      RCP<Array<int> > a = rcp(new Array<int>());
      tmp.put(cf, a);
    }
    SUNDANCE_MSG2(verb, "function #" << i << " is defined on CF " << cf);
    tmp[cf]->append(i);
  }

  int blockID=1;
  int blockPtr=0;
  blockFilters.resize(0);
  funcsForBlock.resize(0);
  elemsForBlock.resize(0);

  for (Map<CellFilter, RCP<Array<int> > >::const_iterator
         i=tmp.begin(); i!=tmp.end(); i++)
  {
    const CellFilter& cf = i->first;
    blockFilters.append(cf);
    funcsForBlock.append(i->second);
    RCP<Array<int> > cells 
      = cellSetToLIDArray(cf.getCells(mesh()));
    elemsForBlock.append(cells);
    int nn = cells->size();
    nElemsPerBlock.append(nn);
    blockIDs.append(blockID++);
    elemBlockPtr.append(blockPtr);
    blockPtr += nn;
  }

  SUNDANCE_MSG2(verb, "block IDs = " << blockIDs);
  SUNDANCE_MSG2(verb, "num elems = " << nElemsPerBlock);
  SUNDANCE_MSG2(verb, "block pointers = " << elemBlockPtr);


  int numElems = blockPtr;
  allElems->resize(numElems);

  int k=0;
  for (int i=0; i<blockIDs.size(); i++)
  {
    SUNDANCE_MSG2(verb, "block " << i << " funcs = " 
      << *funcsForBlock[i]);
    SUNDANCE_MSG2(verb, "block " << i 
      << " elems = " << *elemsForBlock[i]);
    const Array<int>& myCells = *(elemsForBlock[i]);
    for (int c=0; c<myCells.size(); c++)
    {
      (*allElems)[k++] = myCells[c];
    }
  }

  SUNDANCE_MSG2(verb, "all elems = " << *allElems);
}

void ExodusWriter::getCharpp(const Array<std::string>& s, Array<const char*>& p) const
{
  p.resize(s.size());
  for (int i=0; i<p.size(); i++) p[i] = s[i].c_str();
}


