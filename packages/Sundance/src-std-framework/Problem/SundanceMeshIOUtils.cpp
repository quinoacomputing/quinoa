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


#include "SundanceMeshIOUtils.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
  using std::ifstream;
  using std::ofstream;


Expr readNodalFields(const MeshSource& mesher, const Mesh& mesh,
  const VectorType<double>& vecType, int verb)
{
  Tabs tab0(0);
  PLAYA_ROOT_MSG1(verb, tab0 << "begin readNodalFields()");
  Tabs tab1;

  /* Read the attributes */
  RCP<Array<Array<double> > > elemAttr;
  RCP<Array<Array<double> > > vertAttr;

  PLAYA_ROOT_MSG2(verb, tab1 << "reading attributes");
  mesher.getAttributes(vertAttr, elemAttr);
  int nAttrs = vertAttr->size();
  PLAYA_ROOT_MSG2(verb, tab1 << "found " << nAttrs << " attributes");
  const Array<Array<double> >& funcVals = *vertAttr;
  
  /* create an empty (zero-valued) discrete function */
  Array<BasisFamily> bas(nAttrs);
  for (int i=0; i<bas.size(); i++) bas[i] = new Lagrange(1);
  PLAYA_ROOT_MSG2(verb, tab1 << "forming discrete space");
  DiscreteSpace discSpace(mesh, bas, vecType);
  PLAYA_ROOT_MSG2(verb, tab1 << "forming discrete function");
  Expr u0 = new DiscreteFunction(discSpace, 0.0, "u0");
  
  /* get from the discrete function a pointer to the underlying
   * vector and a pointer to the node-to-dof map. */
  Vector<double> vec = DiscreteFunction::discFunc(u0)->getVector();
  const RCP<DOFMapBase>& dofMap 
    = DiscreteFunction::discFunc(u0)->map();
  
  /* run through the data, putting each entry into its correct place in
   * the vector as indexed by the dof number, NOT the node number. */
  PLAYA_ROOT_MSG2(verb, tab1 << "filling discrete function");
  Array<int> dofs(1);
  for (int i=0; i<mesh.numCells(0); i++)
  {
    Tabs tab2;
    PLAYA_ROOT_MSG3(verb, tab2 << "node i=" << i);
    for (int f=0; f<nAttrs; f++)
    {
      Tabs tab3;
      /* look up the dof for the f-th function on this node */
      dofMap->getDOFsForCell(0, i, f, dofs);
      int dof = dofs[0];
      PLAYA_ROOT_MSG3(verb, tab3 << "f=" << f << " dof=" << dof << " val=" 
        << funcVals[f][i]);
      loadable(vec)->setElement(dof, funcVals[f][i]);
    }
  }
  
  /* Reset the vector */
  DiscreteFunction::discFunc(u0)->setVector(vec);

  PLAYA_ROOT_MSG1(verb, tab0 << "done readNodalFields()");
  return u0;
}





Expr readSerialGridField(const std::string& gridFile, 
  double ax, double bx,
  double ay, double by,
  const VectorType<double>& vecType,
  const MeshType& meshType,
  Mesh& mesh)
{
  ifstream is(gridFile.c_str());
  int nNodes ;
  int nx;
  int ny;
  int nAttrs;
  is >> nx >> ny >> nAttrs;
  nNodes = nx*ny;

  Array<Array<double> > funcVals(nAttrs);
  for (int i=0; i<nAttrs; i++) funcVals[i].resize(nNodes);
  for (int i=0; i<nNodes; i++) 
  {
    for (int f=0; f<nAttrs; f++) 
    {
      is >> funcVals[f][i];
    }
  }
  
  
  MeshSource mesher = new PartitionedRectangleMesher(ax, bx, nx-1, 1,
    ay, by, ny-1, 1,
    meshType);

  mesh = mesher.getMesh();

  
  
  /* create an empty (zero-valued) discrete function */
  Array<BasisFamily> bas(nAttrs);
  for (int i=0; i<bas.size(); i++) bas[i] = new Lagrange(1);
  DiscreteSpace discSpace(mesh, bas, vecType);
  Expr u0 = new DiscreteFunction(discSpace, 0.0, "u0");

  
  /* get from the discrete function a pointer to the underlying
   * vector and a pointer to the node-to-dof map. */
  Vector<double> vec = DiscreteFunction::discFunc(u0)->getVector();
  const RCP<DOFMapBase>& dofMap 
    = DiscreteFunction::discFunc(u0)->map();
  
  /* run through the data, putting each entry into its correct place in
   * the vector as indexed by the dof number, NOT the node number. */
  Array<int> dofs(1);
  for (int i=0; i<mesh.numCells(0); i++)
  {
    for (int f=0; f<nAttrs; f++)
    {
      /* look up the dof for the f-th function on this node */
      dofMap->getDOFsForCell(0, i, f, dofs);
      int dof = dofs[0];
      loadable(vec)->setElement(dof, funcVals[f][i]);
    }
  }
  
  /* Reset the vector */
  DiscreteFunction::discFunc(u0)->setVector(vec);

  return u0;
}

double readbackTester(const std::string& infile, const MPIComm& comm, int verb) 
{
  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();
  
  /* Read a mesh */
  Out::root() << "starting to read mesh " << std::endl;
  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new ExodusMeshReader(infile, meshType, verb, comm);
  Mesh mesh = mesher.getMesh();

  Out::root() << "done reading mesh " << std::endl;
  int dim = mesh.spatialDim();
  
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);
  Expr z = new CoordExpr(2);
  

  /* Create a discrete function on the mesh */
  DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
  Out::root() << "done making discfunc " << std::endl;
  Expr u;
  if (dim==3)
  {
    L2Projector proj(discSpace, x*sin(x)*y*y+exp(z)*cos(y));
    u = proj.project();
  }
  else if (dim==2)
  {
    L2Projector proj(discSpace, x*sin(x)+x*x*cos(y));
    u = proj.project();
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPT(dim != 2 && dim != 3);
  }
  Out::root() << "done projecting function " << std::endl;
  /* Compute some functional using the mesh */
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(2);
  Expr F = Integral(interior, u*u, quad);
  double fVal = evaluateIntegral(mesh, F);
  Out::root() << "done evaluating functional on original mesh" << std::endl;
  
  /* write the solution */
  FieldWriter w = new ExodusWriter("./readbackTest-out");
  w.addMesh(mesh);
  w.addField("u", new ExprFieldWrapper(u));
  w.write();

  /* Read back the solution */
  MeshSource mesher2 = new ExodusMeshReader("./readbackTest-out", meshType, verb, comm);
  Mesh mesh2 = mesher2.getMesh();
  Expr u2 = readNodalFields(mesher2, mesh2, vecType, verb);
  Out::root() << "done readback " << std::endl;
  
  /* Compute the functional using the second mesh */
  Expr F2 = Integral(interior, u2*u2, quad);
  double fVal2 = evaluateIntegral(mesh2, F2);
    
  Out::root() << "functional on original mesh = " << fVal << std::endl;
  Out::root() << "functional on readback mesh = " << fVal2 << std::endl;

  double diff = std::fabs(fVal - fVal2);
  Out::root() << "diff = " << diff << std::endl;

  return diff;
}

   

void invertMap(const Map<int,int>& in, Map<int,int>& out)
{
  typedef Map<int,int>::const_iterator iter;

  for (iter i=in.begin(); i!=in.end(); i++)
  {
    out.put(i->second, i->first);
  }
}
 
void serialPartition(
  const RCP<SerialPartitionerBase>& part,
  int numProc,
  const MeshSource& mesher, 
  const std::string& outfile)
{
  /* This should only be run on a single-process communicator. If run in 
   * parallel, this function should be called only by the "MPI_COMM_SELF" 
   * communicator of a single processor. */
  TEUCHOS_TEST_FOR_EXCEPTION(mesher.comm().getNProc() > 1, std::runtime_error,
    "serialPartition() should only be called from a "
    "single-process communicator");

  /* Make the mesh */
  Mesh mesh = mesher.getMesh();
  int dim = mesh.spatialDim();

  /* */
  FieldWriter origWriter = new ExodusWriter("orig-writback");
  origWriter.addMesh(mesh);
  origWriter.write();

  /* Set up maps between old (unpartitioned) and new (partitioned) element
   * and vertex indices. These will be needed to transfer fields, if any,
   * to the partitioned mesh */
  Array<Sundance::Map<int, int> > oldElemLIDToNewLIDMaps;
  Array<Sundance::Map<int, int> > oldVertLIDToNewLIDMaps;

  /* Carry out the partitioning */
  Array<Mesh> submesh = part->makeMeshParts(mesh, numProc,
    oldElemLIDToNewLIDMaps,
    oldVertLIDToNewLIDMaps
    );

  /* Read the fields, if any, from the old mesh */
  RCP<Array<Array<double> > > oldElemData;
  RCP<Array<Array<double> > > oldNodeData;
  mesher.getAttributes(oldNodeData, oldElemData);

  /* Find the process IDs that view each vertex */
  Array<Set<int> > vertViewers(mesh.numCells(0));
  for (int p=0; p<numProc; p++)
  {
    for (int v=0; v<submesh[p].numCells(0); v++)
    {
      int gid = submesh[p].mapLIDToGID(0,v);
      vertViewers[gid].put(p);
    }
  }

  /* Now write the submeshes using the specified writer */
  for (int p=0; p<numProc; p++)
  {
    Out::os() << "writing part=" << p << " of " << numProc << std::endl; 
    FieldWriter writer = new ExodusWriter(outfile);
    /* Set the mesh for writing */
    writer.impersonateParallelProc(numProc, p);
    writer.addMesh(submesh[p]);

    /* Prepare to map the field data to the new submeshes */
    Map<int, int> newElemLIDToOldLIDMap;
    Map<int, int> newVertLIDToOldLIDMap;
    Out::os() << "preparing maps to submeshes" << std::endl; 
    invertMap(oldElemLIDToNewLIDMaps[p], newElemLIDToOldLIDMap);
    invertMap(oldVertLIDToNewLIDMaps[p], newVertLIDToOldLIDMap);

    /* Map the element-based field data */
    Out::os() << "mapping element data to submeshes" << std::endl; 
    RCP<Array<double> > elemData;
    int nElemFuncs = oldElemData->size();
    int nElems = submesh[p].numCells(dim);
    for (int f=0; f<nElemFuncs; f++)
    {
      elemData = rcp(new Array<double>(nElems));
      for (int lid=0; lid<nElems; lid++)
      {
        int oldLID = newElemLIDToOldLIDMap.get(lid);
        (*elemData)[lid] = (*oldElemData)[f][oldLID];
      }
      /* Write the element-based data */
      writer.addField("uElem["+ Teuchos::toString(f)+"]", 
        new CellLIDMappedFieldWrapper(dim, 1, elemData));
    }

    /* Map the node-based field data */
    Out::os() << "mapping node data to submeshes" << std::endl; 
    RCP<Array<double> > nodeData;
    int nNodeFuncs = oldNodeData->size();
    int nNodes = submesh[p].numCells(0);
    for (int f=0; f<nNodeFuncs; f++)
    {
      nodeData = rcp(new Array<double>(nNodes));
      for (int lid=0; lid<nNodes; lid++)
      {
        int oldLID = newVertLIDToOldLIDMap.get(lid);
        (*nodeData)[lid] = (*oldNodeData)[f][oldLID];
      }
      /* Write the node-based data */
      writer.addField("uNode["+ Teuchos::toString(f)+"]", 
        new CellLIDMappedFieldWrapper(0, 1, nodeData));
    }
        
    /* Complete the write of the p-th segment */
    Out::os() << "doing write()" << std::endl; 
    writer.write();
    Out::os() << "done part=" << p << " of " << numProc << std::endl; 
  }

  /* Write the vertex view files */
  for (int p=0; p<numProc; p++)
  {
    Out::os() << "writing part=" << p << " of " << numProc << std::endl; 

    string vvFile = outfile + "-" + Teuchos::toString(numProc)
      + "-" + Teuchos::toString(p) + ".pvv";
    ofstream of(vvFile.c_str());
    of << submesh[p].numCells(0) << endl;
    for (int v=0; v<submesh[p].numCells(0); v++)
    {
      int gid = submesh[p].mapLIDToGID(0, v);
      const Set<int>& viewers = vertViewers[gid];
      of << v << " " << gid << " " << viewers.size();
      for (Set<int>::const_iterator 
             i=viewers.begin(); i!=viewers.end(); i++)
      {
        of << " " << *i;
      }
      of << endl;
    }
  }
}

  
  


}
