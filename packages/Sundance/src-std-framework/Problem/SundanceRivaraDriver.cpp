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


#include "SundanceRivaraDriver.hpp"
#include "SundanceMesh.hpp"
#include "SundanceRivaraMesh.hpp"
#include "SundanceExprFieldWrapper.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Sundance::Rivara;
using Sundance::ExprFieldWrapper;
using std::endl;

void sort(const Array<int>& in, Array<int>& rtn)
{
  rtn.resize(in.size());
  const int* pIn = &(in[0]);
  int* pRtn = &(rtn[0]);
  for (int i=0; i<in.size(); i++) pRtn[i] = pIn[i];
  
  for (int j=1; j<in.size(); j++)
  {
    int key = pRtn[j];
    int i=j-1;
    while (i>=0 && pRtn[i]>key)
    {
      pRtn[i+1]=pRtn[i];
      i--;
    }
    pRtn[i+1]=key;
  }
  //std::sort(rtn.begin(), rtn.end());
}

static Time& refTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("mesh refinement"); 
  return *rtn;
}

static Time& m2rTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("mesh to rivara copy"); 
  return *rtn;
}

static Time& r2mTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("rivara to mesh copy"); 
  return *rtn;
}

static Time& volSetTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("refinement stack building"); 
  return *rtn;
}


Mesh RefinementTransformation::apply(const Mesh& inputMesh) const 
{
  TimeMonitor timer(refTimer());

  int dim = inputMesh.spatialDim();
  MPIComm comm = inputMesh.comm();
  int numElems = inputMesh.numCells(dim);

  RCP<RivaraMesh> rivMesh = rcp(new RivaraMesh(dim, comm));

  Array<int> lidMap;

  meshToRivara(inputMesh, lidMap,rivMesh);
  
  ExprFieldWrapper expr(errExpr_);
  TEUCHOS_TEST_FOR_EXCEPTION(expr.isPointData(), std::runtime_error,
    "Expected cell-based discrete function for area specification");

  {
    TimeMonitor timer(volSetTimer());
    numRefined_ = 0;
    for (int c=0; c<numElems; c++)
    {
      double err = expr.getData(dim,c,0);
      int rivLID = lidMap[c];
      Element* e = rivMesh->element(rivLID).get();
      double vol = e->volume();
      double newVol = vol * std::pow(reqErr_/(err+1.0e-12), 0.5*dim);
///      Out::os() << "c=" << c << " refine by " << newVol/vol << std::endl;
      if (newVol >= vol) continue;
      rivMesh->refinementSet().push(e);
      rivMesh->refinementAreas().push(newVol);
      numRefined_ ++;
    }
  }
  Out::os() << "refining n=" << numRefined_ << " cells" << std::endl;
  rivMesh->refine();


  Mesh outputMesh = rivaraToMesh(rivMesh, comm);

  return outputMesh;
}


void RefinementTransformation::meshToRivara(
  const Mesh& mesh, 
  Array<int>& lidMap,
  RCP<RivaraMesh>& rivMesh) const 
{
  TimeMonitor timer(m2rTimer());
  int dim = mesh.spatialDim();
  int numNodes = mesh.numCells(0);
  int numElems = mesh.numCells(dim);

  for (int n=0; n<numNodes; n++)
  {
    int gid = n;
    int label = mesh.label(0,n);
    int ownerProcID = mesh.ownerProcID(0,n);
    Point x = mesh.nodePosition(n);
    rivMesh->addVertex(gid, x, ownerProcID, label);
  }

  lidMap.resize(numElems);
  for (int e=0; e<numElems; e++)
  {
    int gid = e;
    int label = mesh.label(dim,e);
    int ownerProcID = mesh.ownerProcID(dim,e);
    Array<int> verts;
    Array<int> fo;
    mesh.getFacetArray(dim, e, 0, verts, fo);
    int elemLID = rivMesh->addElement(gid, verts, ownerProcID, label);
    lidMap[e] = elemLID;
    /* label edges or faces */
    if (dim==2)
    {
      Array<int> edgeLIDs;
      mesh.getFacetArray(dim, e, 1, edgeLIDs, fo);
      for (int f=0; f<3; f++)
      {
        int edgeLabel = mesh.label(1,edgeLIDs[f]);
        rivMesh->element(elemLID)->edge(f)->setLabel(edgeLabel);
      }
    }
    else if (dim==3)
    {
      Array<int> faceLIDs;
      mesh.getFacetArray(dim, e, 2, faceLIDs, fo);
      for (int f=0; f<4; f++)
      {
        int faceLabel = mesh.label(2,faceLIDs[f]);
        rivMesh->element(elemLID)->face(f)->setLabel(faceLabel);
      }
    }
  }
}


Mesh RefinementTransformation::rivaraToMesh(const RCP<RivaraMesh>& rivMesh,
  const MPIComm& comm) const 
{
  TimeMonitor timer(r2mTimer());
  int dim = rivMesh->spatialDim();
  int numNodes = rivMesh->numNodes();

  Mesh mesh = meshType_.createEmptyMesh(dim, comm);

  for (int n=0; n<numNodes; n++)
  {
    const RCP<Node>& node = rivMesh->node(n);
    const Point& x = node->pt();
    int gid = node->globalIndex();
    int ownerProcID = node->ownerProc();
    int label = node->label();
    mesh.addVertex(gid, x, ownerProcID, label);
  }


  ElementIterator iter(rivMesh.get());

  int gid=0;

  Array<int> verts(dim+1);
  Array<int> fo;
  Array<int> edgeLIDs;
  Array<int> faceLIDs;
      
  while (iter.hasMoreElements())
  {
    Tabs tab;
    const Element* e = iter.getNextElement();
    int ownerProcID = e->ownerProc();
    int label = e->label();
    const Array<RCP<Node> >& nodes = e->nodes();

    for (int i=0; i<nodes.size(); i++)
    {
      verts[i] = nodes[i]->globalIndex();
    }
    int lid = mesh.addElement(gid, verts, ownerProcID, label);
    gid++;

    Array<int> sortedVerts(verts.size());
    sort(verts, sortedVerts);

    Out::os() << tab << "elem LID=" << lid << " verts=" << sortedVerts << endl; 
    /* label edges or faces */
    if (dim==2)
    {
      if (e->hasNoEdgeLabels()) continue;
      mesh.getFacetArray(dim, lid, 1, edgeLIDs, fo);
      for (int f=0; f<3; f++)
      {
        int edgeLabel = e->edge(f)->label();
        if (edgeLabel < 0) continue;
        mesh.setLabel(1, edgeLIDs[f], edgeLabel);
      }
    }
    else if (dim==3)
    {
      Tabs tab1;
      if (e->hasNoFaceLabels()) continue;
      mesh.getFacetArray(dim, lid, 2, faceLIDs, fo);
      for (int f=0; f<4; f++)
      {
        /* the faces have been renumbered by the mesh */
        int faceLabel = e->face(f)->label();
        if (faceLabel <= 0) continue;
        const FaceNodes& nodes = e->face(f)->nodes();
        Array<int> nodeLIDs = nodes.nodes().elements();
        Out::os() << tab1 << "face nodes = " << nodeLIDs << endl;
        int faceNum = -1;
        if (nodeLIDs[0]==sortedVerts[0] && nodeLIDs[1]==sortedVerts[1]
          && nodeLIDs[2]==sortedVerts[2])
        {
          faceNum = 3;
        }
        else if (nodeLIDs[0]==sortedVerts[0] && nodeLIDs[1]==sortedVerts[1]
          && nodeLIDs[2]==sortedVerts[3])
        {
          faceNum = 2;
        }
        else if (nodeLIDs[0]==sortedVerts[0] && nodeLIDs[1]==sortedVerts[2]
          && nodeLIDs[2]==sortedVerts[3])
        {
          faceNum = 1;
        }
        else if (nodeLIDs[0]==sortedVerts[1] && nodeLIDs[1]==sortedVerts[2]
          && nodeLIDs[2]==sortedVerts[3])
        {
          faceNum = 0;
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPT(true);
        }
        Out::os() << tab1 << "faceLID=" << faceLIDs[faceNum] << " UFC face num=" 
                  << faceNum 
                  << " label = " << faceLabel 
                  << " parent = " << lid << std::endl;
        mesh.setLabel(2, faceLIDs[faceNum], faceLabel);
      }
    }
  }

  return mesh;
  
}

  
