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


#include "SundanceSerialPartitionerBase.hpp"
#include "SundanceMap.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;

using std::max;
using std::min;
using std::endl;
using std::cout;
using std::cerr;



Set<int> SerialPartitionerBase::arrayToSet(const Array<int>& a) const
{
  Set<int> rtn;
  for (int i=0; i<a.size(); i++) 
  {
    rtn.put(a[i]);
  }
  return rtn;
}

int SerialPartitionerBase::max(const Set<int>& s) const
{
  return *(s.rbegin());
}

void SerialPartitionerBase::getNeighbors(const Mesh& mesh, 
  Array<Array<int> >& neighbors, int& nEdges) const
{
  int dim = mesh.spatialDim();
  int nElems = mesh.numCells(dim);
  int nVerts = mesh.numCells(0);

  neighbors.resize(nElems);
  nEdges = 0;

  elemVerts_.resize(nElems);
  elemEdgewiseNbors_.resize(nElems);
  vertElems_.resize(nVerts);

  
  for (int c=0; c<nElems; c++)
  {
    Array<int> facetLID;
    Array<int> facetDir;
    /* Find the vertices associated with this element */
    mesh.getFacetArray(dim, c, 0, facetLID, facetDir);
    elemVerts_[c] = arrayToSet(facetLID);

    /* Find the neighboring elements that share an edge. Needed for Rivara refinement */
    facetLID.resize(0);
    facetDir.resize(0);
    mesh.getFacetArray(dim, c, 1, facetLID, facetDir);
    Set<int> adjElems;
    for (int f=0; f<facetLID.size(); f++)
    {
      Array<int> maxCofs;
      mesh.getCofacets(1, facetLID[f], dim, maxCofs);
      for (int cf=0; cf<maxCofs.size(); cf++)
      {
        if (maxCofs[cf] != c) adjElems.put(maxCofs[cf]);
      }
    }
    elemEdgewiseNbors_[c] = adjElems;
  }    

  for (int v=0; v<nVerts; v++)
  {
    Array<int> cofacetLID;
    mesh.getCofacets(0, v, dim, cofacetLID);
    vertElems_[v] = arrayToSet(cofacetLID);
  }    


  for (int i=0; i<nElems; i++)
  {
    Set<int> allNbors;
    const Set<int>& ev = elemVerts_[i];
    for (Set<int>::const_iterator v=ev.begin(); v!=ev.end(); v++)
    {
      allNbors = allNbors.setUnion(vertElems_[*v]);
    }
    allNbors.erase(i);
    
    Array<int> fullNbors;
    for (Set<int>::const_iterator j=allNbors.begin(); j!=allNbors.end(); j++)
    {
      Set<int> numCommonNodes = elemVerts_[i].intersection(elemVerts_[*j]);
      if ((int) numCommonNodes.size() == dim) fullNbors.append(*j);
    }
    nEdges += fullNbors.size();
    neighbors[i] = fullNbors;
  }

  nEdges = nEdges/2;
}


void SerialPartitionerBase::getNodeAssignments(int nProc, 
  const Array<int>& elemAssignments,
  Array<int>& nodeAssignments,
  Array<int>& nodeOwnerElems,
  Array<int>& nodesPerProc) const
{
  nodesPerProc.resize(nProc);
  nodeAssignments.resize(vertElems_.size());
  nodeOwnerElems.resize(vertElems_.size());

  for (int i=0; i<nodesPerProc.size(); i++) nodesPerProc[i] = 0;

  for (int v=0; v<vertElems_.size(); v++)
  {
    /* assign the vertex to the highest-numbered cofacet */
    int ownerElem = max(vertElems_[v]);
    nodeOwnerElems[v] = ownerElem;
    nodeAssignments[v] = elemAssignments[ownerElem];
    nodesPerProc[nodeAssignments[v]]++;
  }
}

void SerialPartitionerBase::getElemsPerProc(int nProc, 
  const Array<int>& elemAssignments,
  Array<int>& elemsPerProc) const 
{
  elemsPerProc.resize(nProc);
  for (int i=0; i<elemsPerProc.size(); i++) elemsPerProc[i] = 0;

  for (int e=0; e<elemAssignments.size(); e++)
  {
    elemsPerProc[elemAssignments[e]]++;
  }
}

void SerialPartitionerBase::remapEntities(const Array<int>& assignments, 
  int nProc,
  Array<int>& entityMap) const 
{
  Array<Array<int> > procBuckets(nProc);
  entityMap.resize(assignments.size());

  for (int i=0; i<assignments.size(); i++)
  {
    procBuckets[assignments[i]].append(i);
  }

  int g = 0;
  
  for (int p=0; p<nProc; p++)
  {
    for (int i=0; i<procBuckets[p].size(); i++)
    {
      int lid = procBuckets[p][i];
      entityMap[lid] = g++;
    }
  }
}


void SerialPartitionerBase::getOffProcData(int p, 
  const Array<int>& elemAssignments,
  const Array<int>& nodeAssignments,
  Set<int>& offProcNodes,
  Set<int>& offProcElems) const
{
  Out::root() << "ignore ghosts = " << ignoreGhosts_ << endl;
  /* first pass: find off-proc nodes and elems required by on-proc elems */
  for (int e=0; e<elemAssignments.size(); e++)
  {
    if (elemAssignments[e] != p) continue;
    for (Set<int>::const_iterator v=elemVerts_[e].begin(); v!=elemVerts_[e].end(); v++)
    {
      if (nodeAssignments[*v]!=p) offProcNodes.put(*v);
    }
    if (!ignoreGhosts_)
    {
      for (Set<int>::const_iterator 
             n=elemEdgewiseNbors_[e].begin(); n!=elemEdgewiseNbors_[e].end(); n++)
      {
        if (elemAssignments[*n]!=p) offProcElems.put(*n);
      }
    }
  }

  if (!ignoreGhosts_)
  {
    /* second pass: find off-proc elems required by the on-proc nodes */
    for (int v=0; v<nodeAssignments.size(); v++)
    {
      if (nodeAssignments[v] != p) continue;
      const Set<int>& v2e = vertElems_[v];
      for (Set<int>::const_iterator e=v2e.begin(); e!=v2e.end(); e++)
      {
        if (elemAssignments[*e] != p) offProcElems.put(*e);
      }
    }

    /* third pass: find the additional nodes required by the off-proc
     * elems found in the previous step */
    for (Set<int>::const_iterator e=offProcElems.begin(); e!=offProcElems.end(); e++)
    {
      for (Set<int>::const_iterator v=elemVerts_[*e].begin(); v!=elemVerts_[*e].end(); v++)
      {
        if (nodeAssignments[*v]!=p) offProcNodes.put(*v);
      }
    }
  }  
}


Array<Mesh> SerialPartitionerBase::makeMeshParts(const Mesh& mesh, int np,
  Array<Sundance::Map<int, int> >& oldElemLIDToNewLIDMaps,
  Array<Sundance::Map<int, int> >& oldVertLIDToNewLIDMaps) const 
{
  int dim = mesh.spatialDim();

  Array<int> elemAssignments;
  Array<int> nodeAssignments;
  Array<int> nodeOwnerElems;
  Array<int> nodesPerProc;
  
  getAssignments(mesh, np, elemAssignments);

  getNodeAssignments(np, elemAssignments, nodeAssignments, nodeOwnerElems,
    nodesPerProc);

  Array<Array<int> > offProcNodes(np);
  Array<Array<int> > offProcElems(np);

  Array<Set<int> > offProcNodeSet(np);
  Array<Set<int> > offProcElemSet(np);

  Array<Array<int> > procNodes(np);
  Array<Array<int> > procElems(np);

  for (int p=0; p<np; p++)
  {
    getOffProcData(p, elemAssignments, nodeAssignments, 
      offProcNodeSet[p], offProcElemSet[p]);
    offProcNodes[p] = offProcNodeSet[p].elements();
    offProcElems[p] = offProcElemSet[p].elements();
    procElems[p].reserve(elemAssignments.size()/np);
    procNodes[p].reserve(nodeAssignments.size()/np);
  }

  Array<int> remappedElems;
  Array<int> remappedNodes;

  remapEntities(elemAssignments, np, remappedElems);
  remapEntities(nodeAssignments, np, remappedNodes);

  for (int e=0; e<elemAssignments.size(); e++)
  {
    procElems[elemAssignments[e]].append(e);
  }

  for (int n=0; n<nodeAssignments.size(); n++)
  {
    procNodes[nodeAssignments[n]].append(n);
  }

  /* Now we make NP submeshes */
  Array<Mesh> rtn(np);
  oldVertLIDToNewLIDMaps.resize(np);
  oldElemLIDToNewLIDMaps.resize(np);
  for (int p=0; p<np; p++)
  {
    Sundance::Map<int, int>& oldVertLIDToNewLIDMap 
      = oldVertLIDToNewLIDMaps[p];
    Sundance::Map<int, int>& oldElemLIDToNewLIDMap 
      = oldElemLIDToNewLIDMaps[p];
    MeshType type = new BasicSimplicialMeshType();
    rtn[p] = type.createEmptyMesh(mesh.spatialDim(), MPIComm::world());

    Set<int> unusedVertGID;

    /* add the on-processor nodes */
    for (int n=0; n<procNodes[p].size(); n++)
    {
      int oldLID = procNodes[p][n];
      int newGID = remappedNodes[oldLID];
      unusedVertGID.put(newGID);
      int newLID = rtn[p].addVertex(newGID, mesh.nodePosition(oldLID), p, 0);
      oldVertLIDToNewLIDMap.put(oldLID, newLID);
    }

    /* add the off-processor nodes */
    for (int n=0; n<offProcNodes[p].size(); n++)
    {
      int oldLID = offProcNodes[p][n];
      int nodeOwnerProc = nodeAssignments[oldLID];
      int newGID = remappedNodes[oldLID];
      unusedVertGID.put(newGID);
      int newLID = rtn[p].addVertex(newGID, mesh.nodePosition(oldLID), nodeOwnerProc, 0);
      oldVertLIDToNewLIDMap.put(oldLID, newLID);
    }
    
    /* add the on-processor elements */
    for (int e=0; e<procElems[p].size(); e++)
    {
      int oldLID = procElems[p][e];
      int newGID = remappedElems[oldLID];
      Array<int> vertGIDs;
      Array<int> orientations;
      mesh.getFacetArray(dim, oldLID, 0, vertGIDs, orientations);
      for (int v=0; v<vertGIDs.size(); v++)
      {
        vertGIDs[v] = remappedNodes[vertGIDs[v]];
        if (unusedVertGID.contains(vertGIDs[v])) unusedVertGID.erase(newGID);
      }
      int newLID = rtn[p].addElement(newGID, vertGIDs, p, 1);
      oldElemLIDToNewLIDMap.put(oldLID, newLID);
    }
    
    /* add the off-processor elements */
    for (int e=0; e<offProcElems[p].size(); e++)
    {
      int oldLID = offProcElems[p][e];
      int newGID = remappedElems[oldLID];
      int elemOwnerProc = elemAssignments[oldLID];
      Array<int> vertGIDs;
      Array<int> orientations;
      mesh.getFacetArray(dim, oldLID, 0, vertGIDs, orientations);
      for (int v=0; v<vertGIDs.size(); v++)
      {
        vertGIDs[v] = remappedNodes[vertGIDs[v]];
        if (unusedVertGID.contains(vertGIDs[v])) unusedVertGID.erase(newGID);
      }
      int newLID = rtn[p].addElement(newGID, vertGIDs, elemOwnerProc, 1);
      oldElemLIDToNewLIDMap.put(oldLID, newLID);

//      TEUCHOS_TEST_FOR_EXCEPTION(unusedVertGID.size() != 0, std::logic_error,
//        "unused vertices=" << unusedVertGID);
    }

    /* Now, propagate the labels of any intermediate-dimension cells
     * to the submesh */
    for (int d=1; d<dim; d++)
    {
      Set<int> labels = mesh.getAllLabelsForDimension(d);
      for (Set<int>::const_iterator i=labels.begin(); i!=labels.end(); i++)
      {
        Array<int> labeledCells;
        int label = *i;
        if (label == 0) continue;
        mesh.getLIDsForLabel(d, label, labeledCells);
        for (int c=0; c<labeledCells.size(); c++)
        {
          int lid = labeledCells[c];
          Array<int> cofacets;
          mesh.getCofacets(d, lid, dim, cofacets);
          for (int n=0; n<cofacets.size(); n++)
          {
            int cofLID = cofacets[n];
            if (elemAssignments[cofLID]==p 
              || offProcElemSet[p].contains(cofLID))
            {
              /* at this point we need to find the facet index of the side
               * relative to the new element. */
              
              /* find vertices of old cell */
              Array<int> oldVerts;
              Array<int> newVerts;
              Array<int> orientation;
              mesh.getFacetArray(d, lid, 0, oldVerts, orientation);
              for (int v=0; v<oldVerts.size(); v++)
              {
                newVerts.append(remappedNodes[oldVerts[v]]);
              }
              /* Put the vertices in a set. This will let us compare to the
               * vertex sets in the new submesh. */
              Set<int> newVertSet = arrayToSet(newVerts);
              
              /* Find the cofacet in the new submesh */
              int newElemLID = oldElemLIDToNewLIDMap.get(cofLID);

              /* Get the d-facets of the element in the new submesh */
              Array<int> submeshFacets;
              rtn[p].getFacetArray(dim, newElemLID, d, 
                submeshFacets, orientation);
              int facetIndex = -1;
              for (int df=0; df<submeshFacets.size(); df++)
              {
                /* Get the vertices of this d-facet */
                int facetLID = submeshFacets[df];
                Array<int> verts;
                rtn[p].getFacetArray(d, facetLID, 0, verts, orientation);
                Array<int> vertGID(verts.size());
                for (int v=0; v<verts.size(); v++)
                {
                  vertGID[v] = rtn[p].mapLIDToGID(0, verts[v]);
                }
                Set<int> subVertSet = arrayToSet(vertGID);
                if (subVertSet==newVertSet)
                {
                  facetIndex = df;
                  break;
                }
              }
              TEUCHOS_TEST_FOR_EXCEPTION(facetIndex==-1, std::logic_error,
                "couldn't match new " << d << "-cell in submesh to old " << d
                << "cell. This should never happen");

              /* OK, now we have the d-cell's facet index relative to one
               * of its cofacets existing on the new submesh. We now
               * find the LID of the d-cell so we can set its label */
              int o;  // dummy orientation variable; not needed here
              int newFacetLID = rtn[p].facetLID(dim, newElemLID, d, 
                facetIndex, o);
              /* Set the label, finally! */
              rtn[p].setLabel(d, newFacetLID, label);
              break; /* no need to continue the loop over cofacets once 
                      * we've set the label */
            }
            else
            {
            }
          }
        }
      }
    }
  }

  return rtn;
}


