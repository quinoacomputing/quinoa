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


#include "SundanceExodusMeshReader.hpp"
#include "SundanceVertexSort.hpp"
#include "SundanceOut.hpp"
#include "PlayaExceptions.hpp"
#include "SundancePathUtils.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifdef HAVE_SUNDANCE_EXODUS
#include "exodusII.h"
#endif 

using namespace Teuchos;
using namespace Sundance;
using namespace std;

static Time& getExoTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("raw exodus reader functions"); 
  return *rtn;
}


static Time& getExoFillTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("exodus reader fillMesh"); 
  return *rtn;
}


ExodusMeshReader::ExodusMeshReader(const std::string& fname,
  const MeshType& meshType,
  int verbosity,
  const MPIComm& comm)
  : MeshReaderBase(fname, meshType, verbosity, comm), 
    exoFilename_(fname),
    parFilename_(fname)
{
  Tabs tab0(0);

  if (nProc() > 1)
    {
      std::string suffix =  "-" + Teuchos::toString(nProc()) 
        + "-" + Teuchos::toString(myRank());
      exoFilename_ = exoFilename_ + suffix;
      parFilename_ = parFilename_ + suffix;
    }
  exoFilename_ = exoFilename_ + ".exo";
  parFilename_ = parFilename_ + ".pxo";
  
  PLAYA_MSG1(verb(), "exodus filename = " << exoFilename_);
  
  if (nProc() > 1)
  {
    SUNDANCE_OUT(this->verb() > 1,
      "par filename = " << parFilename_);
  }
}





Mesh ExodusMeshReader::fillMesh() const 
{
  TimeMonitor fillTimer(getExoFillTimer());
  Tabs tab0(0);
  Tabs tab1;
  Mesh mesh;
#ifndef HAVE_SUNDANCE_EXODUS
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
    "ExodusMeshReader called for a build without ExodusII");
#else
  PLAYA_ROOT_MSG1(verb(), tab0 << "in ExodusMeshReader::fillMesh()");

  int CPU_word_size = 8;
  int IO_word_size = 0;
  float version;

  Array<int> ptGID;
  Array<int> ptOwner;
  Array<int> elemGID;
  Array<int> elemOwner;

  readParallelInfo(ptGID, ptOwner, elemGID, elemOwner);

  if (verb() > 2) ex_opts(EX_DEBUG | EX_VERBOSE);

  PLAYA_MSG3(verb(), tab1 << "opening file");
  std::string resolvedName = searchForFile(exoFilename_);
  int exoID = ex_open(resolvedName.c_str(), EX_READ, 
    &CPU_word_size, &IO_word_size, &version);

  TEUCHOS_TEST_FOR_EXCEPTION(exoID < 0, std::runtime_error, "ExodusMeshReader unable to "
    "open file: " << exoFilename_);

  PLAYA_MSG3(verb(), tab1 << "exoID=" << exoID);

  TEUCHOS_TEST_FOR_EXCEPT(IO_word_size != 8 || CPU_word_size != 8);

  char title[MAX_LINE_LENGTH+1];

  int dim = 0;
  int numNodes = 0;
  int numElems = 0;
  int numElemBlocks = 0;
  int numNodeSets = 0;
  int numSideSets = 0;

  int ierr = ex_get_init(exoID, title, &dim, &numNodes, &numElems,
    &numElemBlocks, &numNodeSets, &numSideSets);

  TEUCHOS_TEST_FOR_EXCEPTION(numNodes <= 0, std::runtime_error, "invalid numNodes=" 
    << numNodes);
  TEUCHOS_TEST_FOR_EXCEPTION(numElems <= 0, std::runtime_error, "invalid numElems=" 
    << numElems);

  PLAYA_MSG2(verb(), tab1 << "numNodes=" << numNodes);
  PLAYA_MSG2(verb(), tab1 << "numElems=" << numElems);

  /* */
  if (nProc()==1)
  {
    ptGID.resize(numNodes);
    ptOwner.resize(numNodes);
    for (int i=0; i<numNodes; i++)
    {
      ptGID[i] = i;
      ptOwner[i] = 0;
    }
  }
  else
  {
    /* If we're running in parallel, we'd better have consistent numbers
     * of points in the .exo and .pex files. */
    TEUCHOS_TEST_FOR_EXCEPTION((int)ptGID.size() != numNodes, std::runtime_error,
      "ExodusMeshReader::getMesh() found inconsistent "
      "numbers of points in .exo file and .pex files. Exo "
      "file " << exoFilename_ << " had nPoints=" 
      << numNodes << " but .pex file " 
      << parFilename_ << " had nPoints=" << ptGID.size());
  }

  /* now we can build the mesh */
  mesh = createMesh(dim);



  /* Read the points */
  PLAYA_MSG2(verb(), tab1 << "reading vertices");
  Array<double> x(numNodes);
  Array<double> y(numNodes);
  Array<double> z(numNodes * (dim > 2));

  if (dim == 2)
  {
    ierr = ex_get_coord(exoID, &(x[0]), &(y[0]), (void*) 0);
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    TEUCHOS_TEST_FOR_EXCEPT(ierr > 0);
  }
  else if (dim==3)
  {
    ierr = ex_get_coord(exoID, (void*) &(x[0]), (void*) &(y[0]), (void*) &(z[0]));
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  }
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(dim < 2 || dim > 3, std::runtime_error, 
      "invalid dimension=" << dim << " in ExodusMeshReader");
  }


  /* add the points to the mesh */
  PLAYA_MSG2(verb(), tab1 << "adding vertices to mesh");
  for (int n=0; n<numNodes; n++)
  {
    Point p;
    if (dim==2)
    {
      p = Point(x[n], y[n]);
    }
    else
    {
      p = Point(x[n], y[n], z[n]);
    }
    mesh.addVertex(ptGID[n], p, ptOwner[n], 0);
  }


  /* Set up the element numbering */
  if (nProc()==1)
  {
    elemGID.resize(numElems);
    elemOwner.resize(numElems);
    for (int i=0; i<numElems; i++)
    {
      elemGID[i] = i;
      elemOwner[i] = 0;
    }
  }
  else
  {
    /* If we're running in parallel, we'd better have consistent numbers
     * of elements in the .exo and .pex files. */
    TEUCHOS_TEST_FOR_EXCEPTION((int)elemGID.size() != numElems, std::runtime_error,
      "ExodusMeshReader::readElems() found inconsistent "
      "numbers of elements in .exo file and .pex files. Exodus "
      "file " << exoFilename_ << " had nElems=" 
      << numElems << " but .pex file " 
      << parFilename_ << " had nElems=" << elemGID.size());
  }

  /* Read the elements for each block */

  PLAYA_MSG2(verb(), tab1 << "reading elements");
  Array<int> blockIDs(numElemBlocks);
  if (numElemBlocks > 0)
  {
    TimeMonitor timer(getExoTimer());
    ierr = ex_get_elem_blk_ids(exoID, &(blockIDs[0]));
  }
  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  PLAYA_MSG2(verb(), tab1 << "file contains block IDs " << blockIDs);
  int count = 0;
  Array<int> permKey;
  permKey.resize(numElems);
//  Sundance::Map<int,Array<int> > elemVerts;
  Array<int> blockIsSimplicial(numElemBlocks);
  bool allBlocksAreSimplicial = true;

  for (int b=0; b<numElemBlocks; b++)
  {
    Tabs tab2;
    char elemType[MAX_LINE_LENGTH+1];
    int elsInBlock=0;
    int nodesPerEl=0;
    int numAttrs=0;
    int bid = blockIDs[b];
    PLAYA_MSG2(verb(), tab2 << "reading elems for block ID=" << bid);
    {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_elem_block(exoID, bid, elemType, &elsInBlock,
        &nodesPerEl, &numAttrs);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
      "ex_get_elem_block returned error code ierr=" << ierr);

    bool blockIsSimplicial = true;
    if (nodesPerEl != dim+1) 
    {
      blockIsSimplicial=false;
      allBlocksAreSimplicial=false;
    }


    /* get element connectivity */
    PLAYA_MSG2(verb(), tab2 << "reading connectivity for block ID=" << bid);
    Array<int> connect(elsInBlock * nodesPerEl);

    {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_elem_conn(exoID, bid, &(connect[0]));
    }
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    int n=0;
    Array<int> orderedVerts(nodesPerEl);
    Array<int> exVerts(nodesPerEl);

    for (int e=0; e<elsInBlock; e++, n+=nodesPerEl, count++)
    {

      for (int v=0; v<nodesPerEl; v++)
      {
        orderedVerts[v] = ptGID[connect[n+v]-1];
      }
      exVerts = orderedVerts;
      int key = -1;
      if (blockIsSimplicial) 
      {
        vertexSort(orderedVerts, &key);
//        elemVerts.put(elemGID[count], orderedVerts); 
      }

      int elemLID 
        = mesh.addElement(elemGID[count], orderedVerts, 
          elemOwner[count], bid);
      if (blockIsSimplicial) 
      {
//        elemVerts.put(elemLID, orderedVerts); 
        permKey[elemLID]=key;
      }
    }
  }

  /* The topology of the mesh is now fully defined and we can call
   * topological functions */
  mesh.freezeTopology();
 

  /* Read the node sets */
  Array<int> nsIDs(numNodeSets);
  if (numNodeSets > 0)
  {
    TimeMonitor timer(getExoTimer());
    ierr = ex_get_node_set_ids(exoID, &(nsIDs[0]));
  }
  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  for (int ns=0; ns<numNodeSets; ns++)
  {
    int nNodes;
    int nDist;
    int nsID = nsIDs[ns];
    {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_node_set_param(exoID, nsID, &nNodes, &nDist);
    }
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    Array<int> nodes(nNodes);
    {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_node_set(exoID, nsID, &(nodes[0]));
    }
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    for (int n=0; n<nNodes; n++)
    {
      mesh.setLabel(0, nodes[n]-1, nsID);
    }
  }

    
  /* Read the side sets */
  Array<int> ssIDs(numSideSets);
  if (numSideSets > 0)
  {
    TimeMonitor timer(getExoTimer());
    ierr = ex_get_side_set_ids(exoID, &(ssIDs[0]));
  }
  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  for (int ss=0; ss<numSideSets; ss++)
  {
    int nSides;
    int nDist;
    int ssID = ssIDs[ss];
    {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_side_set_param(exoID, ssID, &nSides, &nDist);
    }
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    Array<int> sides(nSides);
    Array<int> elems(nSides);
    {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_side_set(exoID, ssID, &(elems[0]), &(sides[0]));
    }
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
    for (int n=0; n<nSides; n++)
    {
      int elemID = elems[n]-1;
      int facetNum = sides[n]-1;
      int fsign;
      if (allBlocksAreSimplicial)
      {
        int key = permKey[elemID];
        facetNum = exFacetIndexToUFCFacetIndex(dim, key, facetNum);
      }
      int sideLID = mesh.facetLID(dim, elemID, dim-1, 
        facetNum, fsign);
      mesh.setLabel(dim-1, sideLID, ssID);
    }
  }


  /* Read the nodal attributes */
  int nNodalVars = 0;
  {
    TimeMonitor timer(getExoTimer());
    ierr = ex_get_var_param(exoID, "N", &nNodalVars);
  }
  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);

  Array<Array<double> >& funcVals = *nodeAttributes();
  funcVals.resize(nNodalVars);

  for (int i=0; i<nNodalVars; i++)
  {
    int t = 1;
    funcVals[i].resize(mesh.numCells(0));
    {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_nodal_var(exoID, t, i+1, mesh.numCells(0), &(funcVals[i][0]));
    }
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  }

  /* Read the element attributes */
  int nElemVars = 0;
  {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_var_param(exoID, "E", &nElemVars);
  }
  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);

  Array<Array<double> >& eFuncVals = *elemAttributes();
  eFuncVals.resize(nElemVars);

  for (int i=0; i<nElemVars; i++)
  {
    int t = 1;
    eFuncVals[i].resize(mesh.numCells(mesh.spatialDim()));
    {
      TimeMonitor timer(getExoTimer());
      ierr = ex_get_elem_var(exoID, t, i+1, 1, mesh.numCells(mesh.spatialDim()), &(eFuncVals[i][0]));
    }
    TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);
  }

  ierr = ex_close(exoID);
  TEUCHOS_TEST_FOR_EXCEPT(ierr < 0);

#endif
	return mesh;
}




void ExodusMeshReader::readParallelInfo(Array<int>& ptGID, 
  Array<int>& ptOwner,
  Array<int>& cellGID, 
  Array<int>& cellOwner) const
{
  int nPoints;
  int nElems;
  std::string line;
  Array<string> tokens;
  try
  {
    ptGID.resize(0);
    ptOwner.resize(0);
    cellGID.resize(0);
    cellOwner.resize(0);

    /* if we're running in parallel, read the info on processor 
     * distribution */
    if (nProc() > 1)
    {
      RCP<std::ifstream> parStream 
        = openFile(parFilename_, "parallel info");
     
      /* read the number of processors and the processor rank in 
       * the file. These must be consistent with the current number of
       * processors and the current rank */
      getNextLine(*parStream, line, tokens, '#');
      
      TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 2, std::runtime_error,
        "ExodusMeshReader::getMesh() expects 2 entries "
        "on the first line of .par file. In " 
        << parFilename_ << " I found \n[" << line << "]\n");

      int np = atoi(tokens[1]);
      int pid = atoi(tokens[0]);

      /* check consistency with the current number of
       * processors and the current rank */
      
      TEUCHOS_TEST_FOR_EXCEPTION(np != nProc(), std::runtime_error,
        "ExodusMeshReader::getMesh() found "
        "a mismatch between the current number of processors="
        << nProc() << 
        "and the number of processors=" << np
        << "in the file " << parFilename_);

      TEUCHOS_TEST_FOR_EXCEPTION(pid != myRank(), std::runtime_error,
        "ExodusMeshReader::getMesh() found "
        "a mismatch between the current processor rank="
        << myRank() << "and the processor rank="
        << pid << " in the file " << parFilename_);

      /* read the number of points */
      getNextLine(*parStream, line, tokens, '#');

      TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 1, std::runtime_error,
        "ExodusMeshReader::getMesh() requires 1 entry "
        "on the second line of .pxo file. Found line \n[" 
        << line << "]\n in file " << parFilename_);
      
      nPoints = StrUtils::atoi(tokens[0]);

      /* read the global ID and the owner PID for each point */
      ptGID.resize(nPoints);
      ptOwner.resize(nPoints);

      for (int i=0; i<nPoints; i++)
      {
        getNextLine(*parStream, line, tokens, '#');

        TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 3, std::runtime_error,
          "ExodusMeshReader::getMesh() requires 3 "
          "entries on each line of the point section in "
          "the .pxo file. Found line \n[" << line
          << "]\n in file " << parFilename_);

        ptGID[i] = StrUtils::atoi(tokens[1]);
        ptOwner[i] = StrUtils::atoi(tokens[2]);
      }


      /* Read the number of elements */

      getNextLine(*parStream, line, tokens, '#');

      TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 1, std::runtime_error,
        "ExodusMeshReader::getMesh() requires 1 entry "
        "on the cell count line of .pxo file. Found line \n[" 
        << line << "]\n in file " << parFilename_);

      nElems = StrUtils::atoi(tokens[0]);

      SUNDANCE_OUT(this->verb() > 1,
        "read nElems = " << nElems);


      /* read the global ID and the owner PID for each element */

      cellGID.resize(nElems);
      cellOwner.resize(nElems);
      for (int i=0; i<nElems; i++)
      {
        getNextLine(*parStream, line, tokens, '#');

        TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 3, std::runtime_error,
          "ExodusMeshReader::getMesh() requires 3 "
          "entries on each line of the element section in "
          "the .pxo file. Found line \n[" << line
          << "]\n in file " << parFilename_);

        cellGID[i] = StrUtils::atoi(tokens[1]);
        cellOwner[i] = StrUtils::atoi(tokens[2]);
      }
    }

    nPoints = ptGID.length();
    nElems = cellGID.length();
  }
  catch(std::exception& e)
  {
    SUNDANCE_TRACE(e);
  }
}
