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


#include "SundanceTriangleMeshReader.hpp"
#include "SundanceOut.hpp"
#include "SundanceGeomUtils.hpp"
#include "PlayaExceptions.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


TriangleMeshReader::TriangleMeshReader(const std::string& fname,
				       const MeshType& meshType,
				       int verbosity,
				       const MPIComm& comm)
  : MeshReaderBase(fname, meshType, verbosity, comm),
    nodeFilename_(filename()),
    elemFilename_(filename()),
    parFilename_(filename()),
    sideFilename_(filename()),
    edgeFilename_(filename()),
    offset_(0)
{
  

  if (nProc() > 1)
    {
      std::string suffix =  "." + Teuchos::toString(nProc()) 
        + "." + Teuchos::toString(myRank());
      nodeFilename_ = nodeFilename_ + suffix;
      parFilename_ = parFilename_ + suffix;
      elemFilename_ = elemFilename_ + suffix;
      sideFilename_ = sideFilename_ + suffix;
      edgeFilename_ = edgeFilename_ + suffix;
    }
  nodeFilename_ = nodeFilename_ + ".node";
  elemFilename_ = elemFilename_ + ".ele";
  parFilename_ = parFilename_ + ".par";
  sideFilename_ = sideFilename_ + ".side";
  edgeFilename_ = edgeFilename_ + ".edge";
  
  //  verbosity() = 5;
  SUNDANCE_OUT(this->verb() > 1,
               "node filename = " << nodeFilename_);
  
  SUNDANCE_OUT(this->verb() > 1,
               "elem filename = " << elemFilename_);
  
}
TriangleMeshReader::TriangleMeshReader(const ParameterList& params)
  : MeshReaderBase(params),
    nodeFilename_(filename()),
    elemFilename_(filename()),
    parFilename_(filename()),
    sideFilename_(filename()),
    edgeFilename_(filename()),
    offset_(0)
{
  

  if (nProc() > 1)
    {
      std::string suffix =  "." + Teuchos::toString(nProc()) 
        + "." + Teuchos::toString(myRank());
      nodeFilename_ = nodeFilename_ + suffix;
      parFilename_ = parFilename_ + suffix;
      elemFilename_ = elemFilename_ + suffix;
      sideFilename_ = sideFilename_ + suffix;
      edgeFilename_ = edgeFilename_ + suffix;
    }
  nodeFilename_ = nodeFilename_ + ".node";
  elemFilename_ = elemFilename_ + ".ele";
  sideFilename_ = sideFilename_ + ".side";
  edgeFilename_ = edgeFilename_ + ".edge";
  parFilename_ = parFilename_ + ".par";
  
  SUNDANCE_OUT(this->verb() > 1,
               "node filename = " << nodeFilename_);
  
  SUNDANCE_OUT(this->verb() > 1,
               "elem filename = " << elemFilename_);
  
}


Mesh TriangleMeshReader::fillMesh() const 
{
  Mesh mesh;

  Array<int> ptGID;
  Array<int> ptOwner;
  Array<int> cellGID;
  Array<int> cellOwner;

  readParallelInfo(ptGID, ptOwner, cellGID, cellOwner);

  mesh = readNodes(ptGID, ptOwner);

  readElems(mesh, ptGID, cellGID, cellOwner);

  mesh.freezeTopology();

  readSides(mesh);

  readEdges(mesh);

  return mesh;
}

void TriangleMeshReader::readParallelInfo(Array<int>& ptGID, 
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
				     "TriangleMeshReader::getMesh() expects 2 entries "
				     "on the first line of .par file. In " 
				     << parFilename_ << " I found \n[" << line << "]\n");

          int np = atoi(tokens[1]);
          int pid = atoi(tokens[0]);

          /* check consistency with the current number of
           * processors and the current rank */
      
          TEUCHOS_TEST_FOR_EXCEPTION(np != nProc(), std::runtime_error,
				     "TriangleMeshReader::getMesh() found "
				     "a mismatch between the current number of processors="
				     << nProc() << 
				     "and the number of processors=" << np
				     << "in the file " << parFilename_);

          TEUCHOS_TEST_FOR_EXCEPTION(pid != myRank(), std::runtime_error,
				     "TriangleMeshReader::getMesh() found "
				     "a mismatch between the current processor rank="
				     << myRank() << "and the processor rank="
				     << pid << " in the file " << parFilename_);

          /* read the number of points */
          getNextLine(*parStream, line, tokens, '#');

          TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 1, std::runtime_error,
				     "TriangleMeshReader::getMesh() requires 1 entry "
				     "on the second line of .par file. Found line \n[" 
				     << line << "]\n in file " << parFilename_);
      
          nPoints = StrUtils::atoi(tokens[0]);

          /* read the global ID and the owner PID for each point */
          ptGID.resize(nPoints);
          ptOwner.resize(nPoints);

          for (int i=0; i<nPoints; i++)
            {
              getNextLine(*parStream, line, tokens, '#');

              TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 3, std::runtime_error,
					 "TriangleMeshReader::getMesh() requires 3 "
					 "entries on each line of the point section in "
					 "the .par file. Found line \n[" << line
					 << "]\n in file " << parFilename_);

              ptGID[i] = StrUtils::atoi(tokens[1]);
              ptOwner[i] = StrUtils::atoi(tokens[2]);
            }


          /* Read the number of elements */

          getNextLine(*parStream, line, tokens, '#');

          TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 1, std::runtime_error,
				     "TriangleMeshReader::getMesh() requires 1 entry "
				     "on the cell count line of .par file. Found line \n[" 
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
					 "TriangleMeshReader::getMesh() requires 3 "
					 "entries on each line of the element section in "
					 "the .par file. Found line \n[" << line
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

Mesh TriangleMeshReader::readNodes(Array<int>& ptGID,
                                   Array<int>& ptOwner) const 
{
  Mesh mesh;
  std::string line;
  Array<string> tokens;
  int nPoints = -1;

  /* Open the node file so we can read in the nodes */
	
  RCP<std::ifstream> nodeStream = openFile(nodeFilename_, "node info");
	
  /* read the header line */
  getNextLine(*nodeStream, line, tokens, '#');
  TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 4, std::runtime_error,
			     "TriangleMeshReader::getMesh() requires 4 "
			     "entries on the header line in "
			     "the .node file. Found line \n[" << line
			     << "]\n in file " << nodeFilename_);
  string headerLine = line;
  SUNDANCE_OUT(this->verb() > 2,
               "read point header " << line);
  
  
  if (nProc()==1)
    {
      nPoints = atoi(tokens[0]);
      ptGID.resize(nPoints);
      ptOwner.resize(nPoints);
      for (int i=0; i<nPoints; i++)
        {
          ptGID[i] = i;
          ptOwner[i] = 0;
        }
    }
  else
    {
      /* If we're running in parallel, we'd better have consistent numbers
       * of points in the .node and .par file. */
      nPoints = ptGID.length();
      TEUCHOS_TEST_FOR_EXCEPTION(atoi(tokens[0]) != nPoints, std::runtime_error,
				 "TriangleMeshReader::getMesh() found inconsistent "
				 "numbers of points in .node file and par file. Node "
				 "file " << nodeFilename_ << " had nPoints=" 
				 << atoi(tokens[0]) << " but .par file " 
				 << parFilename_ << " had nPoints=" << nPoints);
    }

  SUNDANCE_OUT(this->verb() > 3,
               "expecting to read " << nPoints << " points");
  
  int dimension = atoi(tokens[1]);
  int nAttributes = atoi(tokens[2]);
  int nBdryMarkers = atoi(tokens[3]);

  /* now that we know the dimension, we can build the mesh object */
  mesh = createMesh(dimension);

  /* size point-related arrays */
  Array<int> ptIndices(nPoints);
  Array<bool> usedPoint(nPoints);
  nodeAttributes()->resize(nAttributes);
  for (int i=0; i<nAttributes; i++)
    {
      (*nodeAttributes())[i].resize(nPoints);
    }
  offset_=-1;
  bool first = true;

  /* read all the points */
  for (int count=0; count<nPoints; count++)
    {
      getNextLine(*nodeStream, line, tokens, '#');
      
      TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() 
				 != (1 + dimension + nAttributes + nBdryMarkers),
				 std::runtime_error,
				 "TriangleMeshReader::getMesh() found bad node input "
				 "line. Expected " 
				 << (1 + dimension + nAttributes + nBdryMarkers)
				 << " entries but found line \n[" << 
				 line << "]\n in file " << nodeFilename_);
      /* Triangle files can use either 0-offset or 1-offset numbering. We'll
       * inspect the first node line to decide which numbering to use. */
      if (first)
        {
          offset_ = atoi(tokens[0]);
          TEUCHOS_TEST_FOR_EXCEPTION(offset_ < 0 || offset_ > 1, std::runtime_error,
				     "TriangleMeshReader::getMesh() expected "
				     "either 0-offset or 1-offset numbering. Found an "
				     "initial offset of " << offset_ << " in line \n["
				     << line << "]\n of file " << nodeFilename_);
          first = false;
        }
      
      /* now we can add the point to the mesh */
      double x = atof(tokens[1]); 
      double y = atof(tokens[2]);
      double z = 0.0;
      Point pt;
      int ptLabel = 0;

      if (dimension==3)
	{
	  z = atof(tokens[3]);
          pt = Point(x,y,z);
	}
      else 
	{
          pt = Point(x,y);
	}

      ptIndices[count] 
        = mesh.addVertex(ptGID[count], pt, ptOwner[count], ptLabel);

      for (int i=0; i<nAttributes; i++)
	{
	  (*nodeAttributes())[i][count] = atof(tokens[dimension+1+i]);
	}
    }
  return mesh;
}

void TriangleMeshReader::readElems(Mesh& mesh,
                                   const Array<int>& ptGID,
                                   Array<int>& elemGID,
                                   Array<int>& elemOwner) const 
{
  try
    {
      std::string line;  
      Array<string> tokens;
      /* Open the element file */
	
      RCP<std::ifstream> elemStream = openFile(elemFilename_, "element info");

      getNextLine(*elemStream, line, tokens, '#');

      TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 3, std::runtime_error,
				 "TriangleMeshReader::getMesh() requires 3 "
				 "entries on the header line in "
				 "the .ele file. Found line \n[" << line
				 << "]\n in file " << elemFilename_);
                   
      int nElems = -1;

      if (nProc()==1)
        {
          nElems = atoi(tokens[0]);
          elemGID.resize(nElems);
          elemOwner.resize(nElems);
          for (int i=0; i<nElems; i++)
            {
              elemGID[i] = i;
              elemOwner[i] = 0;
            }
        }
      else
        {
          /* If we're running in parallel, we'd better have consistent numbers
           * of points in the .node and .par file. */
          nElems = elemGID.length();
          TEUCHOS_TEST_FOR_EXCEPTION(atoi(tokens[0]) != nElems, std::runtime_error,
				     "TriangleMeshReader::readElems() found inconsistent "
				     "numbers of elements in .ele file and par file. Elem "
				     "file " << elemFilename_ << " had nElems=" 
				     << atoi(tokens[0]) << " but .par file " 
				     << parFilename_ << " had nElems=" << nElems);
        }

      int ptsPerElem = atoi(tokens[1]);

      TEUCHOS_TEST_FOR_EXCEPTION(ptsPerElem != mesh.spatialDim()+1, std::runtime_error,
				 "TriangleMeshReader::readElems() found inconsistency "
				 "between number of points per element=" << ptsPerElem 
				 << " and dimension=" << mesh.spatialDim() << ". Number of pts "
				 "per element should be dimension + 1");

      int nAttributes = atoi(tokens[2]);
      elemAttributes()->resize(nElems);

      int dim = mesh.spatialDim();
      Array<int> nodes(dim+1);

      for (int count=0; count<nElems; count++)
        {
          getNextLine(*elemStream, line, tokens, '#');
      
          TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() 
				     != (1 + ptsPerElem + nAttributes),
				     std::runtime_error,
				     "TriangleMeshReader::readElems() found bad elem "
				     "input line. Expected " 
				     << (1 + ptsPerElem + nAttributes)
				     << " entries but found line \n[" << 
				     line << "]\n in file " << elemFilename_);

          for (int d=0; d<=dim; d++)
            {
              nodes[d] = ptGID[atoi(tokens[d+1])-offset_];
            }

          int elemLabel = 0;
          try
            {
              mesh.addElement(elemGID[count], nodes, elemOwner[count], elemLabel);
            }
          catch(std::exception& ex1)
            {
              SUNDANCE_TRACE(ex1);
            }
			
          (*elemAttributes())[count].resize(nAttributes);
          for (int i=0; i<nAttributes; i++)
            {
              (*elemAttributes())[count][i] = atof(tokens[1+ptsPerElem+i]);
            }
        }
    }
  catch(std::exception& ex)
    {
      SUNDANCE_TRACE(ex);
    }
}


void TriangleMeshReader::readSides(Mesh& mesh) const 
{
  try
    {
      std::string line;  
      Array<string> tokens;
      /* Open the side file */
      RCP<std::ifstream> sideStream;
      bool fileOK = false;

      try
        {
          sideStream = openFile(sideFilename_, "side info");
          fileOK = true;
        }
      catch(std::exception& e) {;}

      /* Not all meshes will have sides files.
       * If the sides file doesn't exist, return. */
      if (!fileOK)
        {
          SUNDANCE_VERB_LOW("side file [" << sideFilename_ << "] not found");
          return;
        }

      getNextLine(*sideStream, line, tokens, '#');

      TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 1, std::runtime_error,
				 "TriangleMeshReader::readSides() requires 1 "
				 "entry on the header line in "
				 "the .side file. Found line \n[" << line
				 << "]\n in file " << sideFilename_);

      int nSides = atoi(tokens[0]);

      int elemDim = mesh.spatialDim();
      int sideDim = elemDim - 1;

      for (int i=0; i<nSides; i++)
        {
          getNextLine(*sideStream, line, tokens, '#');
      
          TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 4,
				     std::runtime_error,
				     "TriangleMeshReader::readSides() found bad side "
				     "input line. Expected 4 entries but found line \n[" 
				     << line << "]\n in file " << sideFilename_);

          int elemGID = atoi(tokens[1]);
          int elemFacet = atoi(tokens[2]);
          TEUCHOS_TEST_FOR_EXCEPTION(!mesh.hasGID(elemDim, elemGID), std::runtime_error,
				     "element GID " << elemGID << " not found");
          int elemLID = mesh.mapGIDToLID(elemDim, elemGID);
          int o=0; // dummy orientation variable; not needed here
          int sideLID = mesh.facetLID(elemDim, elemLID, sideDim, 
				      elemFacet, o);
          
          int sideLabel = atoi(tokens[3]);

          mesh.setLabel(sideDim, sideLID, sideLabel);
        }
    }
  catch(std::exception& ex)
    {
      SUNDANCE_TRACE(ex);
    }
}



void TriangleMeshReader::readEdges(Mesh& mesh) const 
{
  try
    {
      std::string line;  
      Array<string> tokens;
      /* Open the edge file */
      RCP<std::ifstream> edgeStream;
      bool fileOK = false;

      try
        {
          edgeStream = openFile(edgeFilename_, "edge info");
          fileOK = true;
        }
      catch(std::exception& e) {;}

      /* Not all meshes will have edges files.
       * If the edges file doesn't exist, return. */
      if (!fileOK)
        {
          SUNDANCE_VERB_LOW("edge file [" << edgeFilename_ << "] not found");
          return;
        }

      getNextLine(*edgeStream, line, tokens, '#');

      TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 2, std::runtime_error,
				 "TriangleMeshReader::readEdges() requires 2 "
				 "entries on the header line in "
				 "the .edge file. Found line \n[" << line
				 << "]\n in file " << edgeFilename_);

      int nEdges = atoi(tokens[0]);

      for (int i=0; i<nEdges; i++)
        {
          getNextLine(*edgeStream, line, tokens, '#');
      
          TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 4,
				     std::runtime_error,
				     "TriangleMeshReader::readEdges() found bad edge "
				     "input line. Expected 4 entries but found line \n[" 
				     << line << "]\n in file " << edgeFilename_);

          int v1 = atoi(tokens[1])-offset_;
          int v2 = atoi(tokens[2])-offset_;
          int label = atoi(tokens[3]);
	  if (label==0) continue;

	  int edgeLID = lookupEdgeLIDFromVerts(mesh, v1, v2);

          TEUCHOS_TEST_FOR_EXCEPTION(edgeLID < 0, std::runtime_error,
				     "edge(" << v1 << ", " << v2 <<") not found");
          mesh.setLabel(1, edgeLID, label);
        }
    }
  catch(std::exception& ex)
    {
      SUNDANCE_TRACE(ex);
    }
}
