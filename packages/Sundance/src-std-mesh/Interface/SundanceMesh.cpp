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

#include "SundanceMesh.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceOut.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "PlayaTabs.hpp"
#include "PlayaMPIContainerComm.hpp"

using namespace Teuchos;
using namespace Sundance;
using Playa::Handle;
using Playa::Handleable;
using Playa::MPIComm;
using Playa::MPIContainerComm;


IncrementallyCreatableMesh* Mesh::creatableMesh()
{
  IncrementallyCreatableMesh* mci 
    = dynamic_cast<IncrementallyCreatableMesh*>(ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(mci==0, std::runtime_error, 
                     "Mesh::creatableMesh() could not convert mesh to "
                     "a type deriving from IncrementallyCreatableMesh");

  return mci;
}


void Mesh::dump(const std::string& filename) const
{
  FieldWriter w = new VerboseFieldWriter(filename);
  w.addMesh(*this);
  w.write();
}

bool Mesh::checkConsistency(const std::string& filename) const 
{
  std::string f = filename;

  if (comm().getNProc() > 1)
    {
      f = f + "." + Teuchos::toString(comm().getRank());
    }
  std::ofstream os(f.c_str());
  return checkConsistency(os);
}

bool Mesh::checkConsistency(std::ostream& os) const 
{
  /* eliminate the trivial serial case */
  if (comm().getNProc()==1) 
    {
      os << "Mesh is serial, thus it is consistent across processors" << std::endl;
      return true;
    }

  
  if (spatialDim() >= 2) 
    const_cast<Mesh*>(this)->assignIntermediateCellGIDs(1);
  if (spatialDim() >= 3)  
    const_cast<Mesh*>(this)->assignIntermediateCellGIDs(2);
  bool isOK = checkVertexConsistency(os);
  for (int d=1; d<=spatialDim(); d++)
    {
      isOK = checkCellConsistency(os, d) && isOK;
    }

  return isOK;
}


bool Mesh::checkCellConsistency(std::ostream& os, int dim) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(dim==0, std::runtime_error,
                     "Mesh::checkCellConsistency called for d=0. "
                     "checkVertexConsistency() should be called instead");

  int myRank = comm().getRank();
  int nProcs = comm().getNProc();
  Array<int> nFacets(dim);

  bool elemsAreOK = true;

  /* compute the amount of data associated with each element. The
   * first two entries are the element's GID and owner. The subsequent
   * entires are the GIDs for the facets of all dimensions */
  int dataSize = 2;
  for (int d=0; d<dim; d++) 
    {
      nFacets[d] = numFacets(dim, 0, d);
      dataSize += 2*nFacets[d];
      if (d > 0) dataSize += nFacets[d]*numFacets(d, 0, 0);
    }

  int nCells = numCells(dim);
  Array<int> elemData(dataSize*nCells);

  for (int c=0; c<nCells; c++)
    {
      elemData[dataSize*c] = mapLIDToGID(dim, c);
      elemData[dataSize*c + 1] = ownerProcID(dim, c);
      int base = dataSize*c + 2;
      int k=0;
      for (int d=0; d<dim; d++)
        {
          for (int f=0; f<nFacets[d]; f++)
            {
              int orientation;
              int lid = facetLID(dim, c, d, f, orientation);
              int facetGID = mapLIDToGID(d, lid);
              int facetOwner = ownerProcID(d, lid);
              elemData[base + k++] = facetGID;
              elemData[base + k++] = facetOwner;
              if (d > 0)
                {
                  int nNodes = numFacets(d, lid, 0);
                  for (int v=0; v<nNodes; v++)
                    {
                      int nodeLID = facetLID(d, lid, 0, v, orientation);
                      int nodeGID = mapLIDToGID(0, nodeLID);
                      elemData[base + k++] = nodeGID;
                    }
                }          
            }
        }
    }

  /* do an all-to-all. AllGather would be better, but we have a clean
   * all-to-all implemented */
  Array<Array<int> > outData(nProcs);
  for (int p=0; p<nProcs; p++)
    {
      outData[p] = elemData;
    }

  Array<Array<int> > inData(nProcs);
  
  MPIContainerComm<int>::allToAll(outData, inData, comm()); 
  
  for (int p=0; p<nProcs; p++)
    {
      Tabs tab1;
      if (p==myRank) continue;

      os << tab1 << "p=" << myRank << " testing " << dim 
           << "-cells from remote p=" << p << std::endl;
      
      const Array<int>& remoteData = inData[p];
      int nRemote = remoteData.size()/dataSize;
      
      for (int c=0; c<nRemote; c++)
        {
          Tabs tab2;
          int cellGID = remoteData[dataSize*c];
          int cellOwner = remoteData[dataSize*c+1];
          int cellLID = -1;
          bool elemIsOK = checkRemoteEntity(os, p, dim, cellGID, 
                                            cellOwner, false, cellLID);
          if (cellLID >= 0 && elemIsOK) 
            {
              int base = dataSize*c + 2;
              int k = 0;
              for (int d=0; d<dim; d++)
                {
                  Tabs tab3;
                  int dir;
                  os << tab3 << "checking " << d << "-facets" << std::endl;
                  /* The facets may be stored in permuted order on the
                   * different processors. We can get a common ordering
                   * by sorting them. We define a STL map from facet GID to
                   * facet owner, which serves the dual purpose of sorting
                   * the facets by GID and preserving the pairing of.
                   * GIDs and owners.
                   */
                  Sundance::Map<int, int> remoteFacetToOwnerMap;
                  Sundance::Map<int, int> localFacetToOwnerMap;
                  Sundance::Map<int, Set<int> > remoteFacetToNodesMap;
                  Sundance::Map<int, Set<int> > localFacetToNodesMap;

                  for (int f=0; f<nFacets[d]; f++)
                    {

                      int lfLID = facetLID(dim, cellLID, d, f, dir); 
                      int lfGID = mapLIDToGID(d, lfLID);
                      int lfOwner = ownerProcID(d, lfLID);
                      localFacetToOwnerMap.put(lfGID, lfOwner);
                      int rfGID = remoteData[base + k++];
                      int rfOwner = remoteData[base + k++];
                      remoteFacetToOwnerMap.put(rfGID, rfOwner);
                      if (d > 0)
                        {
                          int numNodes = numFacets(d, 0, 0);
                          Sundance::Set<int> localNodes;
                          Sundance::Set<int> remoteNodes;
                           for (int v=0; v<numNodes; v++)
                            {
                              int nodeLID = facetLID(d, lfLID, 0, v, dir);
                              int nodeGID = mapLIDToGID(0, nodeLID);
                              localNodes.put(nodeGID);
                              int remoteNodeGID = remoteData[base + k++];
                              remoteNodes.put(remoteNodeGID);
                            }
                           localFacetToNodesMap.put(lfGID, localNodes);
                           remoteFacetToNodesMap.put(rfGID, remoteNodes);
                        }
                    }
                  /* Now that the facets are in sorted order, we can 
                   * compare them */
                  for (Sundance::Map<int,int>::const_iterator 
                         rf=remoteFacetToOwnerMap.begin(), 
                         lf=localFacetToOwnerMap.begin();
                       rf != remoteFacetToOwnerMap.end();
                       rf++, lf++)
                    {
                      int lfGID = lf->first;
                      int lfOwner = lf->second;
                      int rfGID = lf->first;
                      int rfOwner = lf->second;
                      Tabs tab4;
                      os << tab4 << " local facet GID=" << lfGID
                         << " remote GID=" << rfGID 
                         << " local Own=" << lfOwner 
                         << " remote Own=" << rfOwner << std::endl;
                      elemIsOK = testIdentity(os, rfGID, lfGID, 
                                              "facet GIDs") && elemIsOK;
                      elemIsOK = testIdentity(os, rfOwner, 
                                              lfOwner, 
                                              "facet owners") && elemIsOK;
                      
                      if (d > 0)
                        {
                          Tabs tab5;
                          /* compare nodes */
                          const Set<int>& localNodes = 
                            localFacetToNodesMap.get(lfGID);
                          const Set<int>& remoteNodes = 
                            remoteFacetToNodesMap.get(rfGID);
                          elemIsOK = testIdentity(os, localNodes.elements(),
                                                  remoteNodes.elements(), 
                                                  "facet nodes") && elemIsOK;
                          os << tab5 << "facet node GIDs: local=" 
                             << localNodes << " remote=" 
                             << remoteNodes << std::endl;
                        }
                    }
                }
            }

          elemsAreOK = elemsAreOK && elemIsOK;
        }
    }

  return elemsAreOK;
} 

bool Mesh::testIdentity(std::ostream& os, int a, int b, const std::string& msg) const
{
  Tabs tab;
  if (a != b)
    {
      os << tab << "CONFLICT in " << msg << std::endl;
      return false;
    }
  return true;
}

bool Mesh::testIdentity(std::ostream& os, 
                        const Array<int>& a,
                        const Array<int>& b, const std::string& msg) const
{
  Tabs tab;
  bool ok = true;
  for (int i=0; i<a.size(); i++)
    {
      if (a[i] != b[i]) ok = false;
    }
  if (!ok)
    {
      os << tab << "CONFLICT in " << msg << std::endl;
    }
  return ok;
}

bool Mesh::checkRemoteEntity(std::ostream& os, int p, int dim, int gid, int owner,
                             bool mustExist, int& lid) const
{
  Tabs tab;
  bool isOK = true;

  int myRank = comm().getRank();
  lid = -1;
  
  if (hasGID(dim, gid)) 
    {
      lid = mapGIDToLID(dim, gid); 
      os << tab << "p=" << myRank << " got " 
           << dim << "-cell GID=" << gid << " from proc=" << p 
           << ", is LID=" << lid << " locally" << std::endl;
    }
  else
    {
      os << tab << "p=" << myRank << " got " << dim << "-cell GID=" << gid 
           << " from proc=" << p << ", doesn't exist locally" << std::endl;
      if (mustExist) isOK = false;
    }

  if (lid >= 0)
    {
      int localOwner = ownerProcID(dim, lid);
      os << tab << dim << "-cell GID=" << gid 
           << " is locally LID=" << lid << " and owned by "
           << localOwner << std::endl;
      if (localOwner != owner)
        {
          os << tab << "OWNERSHIP CONFLICT: local p=" << myRank
               << " thinks " << dim << "-cell GID=" << gid << " is owned by "
               << localOwner << " but remote proc=" << p
               << " says it's owned by " << owner << std::endl;
          isOK = false;
        }
    }

  return isOK;

}



bool Mesh::checkVertexConsistency(std::ostream& os) const 
{
  int dim = spatialDim();
  int myRank = comm().getRank();
  int nProcs = comm().getNProc();
  std::string pHdr = "p=" + Teuchos::toString(myRank);

  Out::println(pHdr + " testing consistency of vertices");

  int dataSize = 2;
  Array<int> vertData(dataSize*numCells(0));
  Array<double> vertPositions(numCells(0)*dim);
  for (int v=0; v<numCells(0); v++)
    {
      vertData[dataSize*v] = mapLIDToGID(0, v);
      vertData[dataSize*v + 1] = ownerProcID(0, v);
      Point x = nodePosition(v);
      for (int j=0; j<dim; j++)
        { 
          vertPositions[dim*v + j] = x[j];
        }
    }

  /* do an all-to-all. AllGather would be better, but we have a clean
   * all-to-all implemented */
  Array<Array<int> > outData(nProcs);
  Array<Array<double> > outPos(nProcs);
  for (int p=0; p<nProcs; p++)
    {
      outData[p] = vertData;
      outPos[p] = vertPositions;
    }

  Array<Array<int> > inData(nProcs);
  Array<Array<double> > inPos(nProcs);

  MPIContainerComm<int>::allToAll(outData,
                                  inData,
                                  comm()); 

  MPIContainerComm<double>::allToAll(outPos, inPos, comm()); 

  double tol = 1.0e-15;

  bool ok = true;
  bool allVertsOK = true;
  

  for (int p=0; p<nProcs; p++)
    {
      Tabs tab1;
      if (p==myRank) continue;

      os << tab1 << "p=" << myRank << " testing vertices from remote p="
           << p << std::endl;

      int nVerts = inData[p].size()/dataSize;

      for (int v=0; v<nVerts; v++)
        {
          Tabs tab2;
          int vertGID = inData[p][v*dataSize];
          int vertOwner = inData[p][v*dataSize+1];
          int vertLID = -1;
          bool vertIsOK = checkRemoteEntity(os, p, 0, vertGID, vertOwner,
                                            false, vertLID);

          if (vertLID >= 0)
            {
              Point localX = nodePosition(vertLID);
              /* copy pt to get size right */
              Point remoteX = localX;
              for (int j=0; j<dim; j++) remoteX[j] = inPos[p][dim*v + j];
              double dx = remoteX.distance(localX);
              if (dx > tol) 
                {
                  os << tab2 << "POSITION CONFLICT: local p=" << myRank
                       << " thinks GID=" << vertGID << " is at xLocal="
                       << localX << " but remote proc=" << p
                       << " says it's at xRemote" << remoteX << std::endl;
                  os << tab2 << "distance = " << dx << std::endl;
                  vertIsOK = false;
                }
            }
          allVertsOK = vertIsOK && allVertsOK;
          if (vertIsOK) 
            {
              os << tab2 << "p=" << myRank 
                   << " says vertex GID=" << vertGID << " is OK" << std::endl;
            }
        }
    }

  if (allVertsOK)
    {
      os << "p=" << myRank << " found all vertex data is OK" << std::endl;
    }
  else
    {
      os << "p=" << myRank << " detected conflicts in vertex data" << std::endl;
    }

  ok = allVertsOK && ok;


  
  return ok;
}

