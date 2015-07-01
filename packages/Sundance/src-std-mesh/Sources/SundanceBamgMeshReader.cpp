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


#include "SundanceBamgMeshReader.hpp"
#include "SundanceOut.hpp"
#include "PlayaExceptions.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Teuchos;
using namespace Sundance;


BamgMeshReader::BamgMeshReader(const std::string& fname,
                               const MeshType& meshType, bool bbAttr,
  int verbosity,
                               const MPIComm& comm)
  : MeshReaderBase(fname, meshType, verbosity, comm),
    nodeFilename_(filename()),
    elemFilename_(filename()),
    parFilename_(filename()),
    meshFilename_(filename()), //new
    bbFilename_(filename()), //new
    bbAttr_(bbAttr) //new
  
{
  
  //keep, but assume nProc() = 1
  if (nProc() > 1)
    {
      nodeFilename_ = nodeFilename_ + Teuchos::toString(myRank());
      parFilename_ = parFilename_ + Teuchos::toString(myRank());
      elemFilename_ = elemFilename_ + Teuchos::toString(myRank());
    }
  //

  nodeFilename_ = nodeFilename_ + ".node";
  elemFilename_ = elemFilename_ + ".ele";
  parFilename_ = parFilename_ + ".par";
  meshFilename_ = meshFilename_ + ".mesh"; //new
  bbFilename_ = bbFilename_ + ".bb"; //new
  
  SUNDANCE_OUT(this->verb() > 1,
               "node filename = " << nodeFilename_);
  
  SUNDANCE_OUT(this->verb() > 1,
               "elem filename = " << elemFilename_);

  SUNDANCE_OUT(this->verb() > 1,
               "mesh filename = " << meshFilename_); //new

  SUNDANCE_OUT(this->verb() > 1,
               "bb filename = " << bbFilename_); //new
  
}
BamgMeshReader::BamgMeshReader(const ParameterList& params)
  : MeshReaderBase(params),
    nodeFilename_(filename()),
    elemFilename_(filename()),
    parFilename_(filename()),
    meshFilename_(filename()) //new
{
  
  //keep, but assume nProc() = 1
  if (nProc() > 1)
    {
      nodeFilename_ = nodeFilename_ + Teuchos::toString(myRank());
      parFilename_ = parFilename_ + Teuchos::toString(myRank());
      elemFilename_ = elemFilename_ + Teuchos::toString(myRank());
    }
  //

  nodeFilename_ = nodeFilename_ + ".node";
  elemFilename_ = elemFilename_ + ".ele";
  parFilename_ = parFilename_ + ".par";
  meshFilename_ = meshFilename_ + ".mesh"; //new
  
  SUNDANCE_OUT(this->verb() > 1,
               "node filename = " << nodeFilename_);
  
  SUNDANCE_OUT(this->verb() > 1,
               "elem filename = " << elemFilename_);

  SUNDANCE_OUT(this->verb() > 1,
               "mesh filename = " << meshFilename_); //new

  SUNDANCE_OUT(this->verb() > 1,
               "bb filename = " << bbFilename_); //new
  
}


Mesh BamgMeshReader::fillMesh() const 
{
  Mesh mesh;

  Array<int> ptGID;
  Array<int> ptOwner;
  Array<int> cellGID;
  Array<int> cellOwner;

  readParallelInfo(ptGID, ptOwner, cellGID, cellOwner);

  //mesh = readNodes(ptGID, ptOwner); //replaced by readMesh

  //readElems(mesh, ptGID, cellGID, cellOwner); //readMesh reads nodes+elems

  mesh = readMesh(ptGID, ptOwner); //new -- replaces readNodes, readElems

  //  mesh.assignGlobalIndices();

	return mesh;
}

void BamgMeshReader::readParallelInfo(Array<int>& ptGID, 
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

      //keep, but assume nProc() = 1
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


//Mesh TriangleMeshReader::readNodes(Array<int>& ptGID, 
//                                   Array<int>& ptOwner) const 
Mesh BamgMeshReader::readMesh(Array<int>& ptGID, 
                              Array<int>& ptOwner) const 
{
  Mesh mesh;
  std::string line;
  Array<string> tokens;
  int nPoints=0;

  /* Open the mesh file so we can read in the nodes and elements */
	
  ///////////////////////////////////////////////////////////////////
    // skip this part and plug in the corresponding Bamg reader code //
    ///////////////////////////////////////////////////////////////////

    /*

    RCP<std::ifstream> nodeStream = openFile(nodeFilename_, "node info");
    // read the header line //
    getNextLine(*nodeStream, line, tokens, '#');
    TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 4, std::runtime_error,
    "TriangleMeshReader::getMesh() requires 4 "
    "entries on the header line in "
    "the .node file. Found line \n[" << line
    << "]\n in file " << nodeFilename_);
    std::string headerLine = line;
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
    // If we're running in parallel, we'd better have consistent numbers
    // of points in the .node and .par file. //
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

    // now that we know the dimension, we can build the mesh object //
    mesh = createMesh(dimension);

    // size point-related arrays //
    Array<int> ptIndices(nPoints);
    Array<bool> usedPoint(nPoints);
    nodeAttributes()->resize(nPoints);

    int offset=0;
    bool first = true;

    // read all the points //
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
    // Triangle files can use either 0-offset or 1-offset numbering. We'll
    // inspect the first node line to decide which numbering to use. //
    if (first)
    {
    offset = atoi(tokens[0]);
    TEUCHOS_TEST_FOR_EXCEPTION(offset < 0 || offset > 1, std::runtime_error,
    "TriangleMeshReader::getMesh() expected "
    "either 0-offset or 1-offset numbering. Found an "
    "initial offset of " << offset << " in line \n["
    << line << "]\n of file " << nodeFilename_);
    first = false;
    }
    } //added this '}' because the for loop continued in TriangleMeshReader
  
    */

    /////////// here's the Bamg reader for reading nodes /////////////
    std::cerr << "starting to read meshFile" << std::endl;

    RCP<std::ifstream> meshStream = openFile(meshFilename_, 
                                                "node & elem info");
    //Array<string> liners = StrUtils::readFile(meshStream, '#');
    Array<string> liners = StrUtils::readFile(*meshStream, '#');
    
    std::cerr << "defined liners" << std::endl;
    //nodeStream.close(); //old

    SUNDANCE_OUT(this->verb() > 3, 
                 "reading nodes from " + meshFilename_);

    // extract dimension, nodes, elements (triangles) from lines
    int dimension=0;
    int dimensionIndex = 0;
    int lineIndex = 0;
    int verticesIndex = 0;
    int edgesIndex = 0;
    int triangleIndex = 0;
    int nCells=0;
    std::cerr << "DimensionIndex = " << dimensionIndex << std::endl;

    //int linerssize = liners.size();
    int linerssize = liners.length();
    std::cerr << "number of lines in liners = " << linerssize << std::endl;

    for (int i = lineIndex; i < linerssize; i++)
      {
        tokens = StrUtils::stringTokenizer(liners[i]);
        std::cerr << "test: i = " << i << "tokens[0] = " << tokens[0] << std::endl;
        if (i < 20)
          {
            std::cerr << "i = " << i << ": number of tokens = " << tokens.length() 
                 << std::endl;
            for(int j = 0; j < tokens.length();j++)
              {
                std::cerr << "     token " << j << "  = " << tokens[j] << std::endl;
                std::cerr << "     length of token[" << j << "] = " 
                     << tokens[j].length() << std::endl;
              }
          }
        if (tokens[0] == "Dimension") 
          {
            dimensionIndex = i;
            std::cerr << "token[0] = Dimension" << std::endl;
            break;
          }
      }
    std::cerr << "DimensionIndex = " << dimensionIndex << std::endl;
    if (dimensionIndex > 0)
      {
        tokens = StrUtils::stringTokenizer(liners[dimensionIndex + 1]);
        dimension = atoi(tokens[0]);
        std::cerr << "dimension = " << dimension << std::endl;
        lineIndex = dimensionIndex + 2;
      }

    for (int i = lineIndex; i < liners.size(); i++)
      {
        tokens = StrUtils::stringTokenizer(liners[i]);
        if (tokens[0] == "Vertices") verticesIndex = i;
        if (tokens[0] == "Vertices") break;
      }
    if (verticesIndex > 0)
      {
        tokens = StrUtils::stringTokenizer(liners[verticesIndex + 1]);
        nPoints = atoi(tokens[0]);
        std::cerr << "nPoints = " << nPoints << std::endl;
        ptGID.resize(nPoints);
        ptOwner.resize(nPoints);
        for (int i=0; i<nPoints; i++)
          {
            ptGID[i] = i;
            ptOwner[i] = 0;
          }
        lineIndex = verticesIndex + 2; 
        //assume node data starts 2 lines below "Vertices"
      }
	
    Array<int> ptIndices(nPoints);
    Array<int> rawIndices(nPoints);

    Array<double> velVector(2*nPoints);
    int numbbAttr = 0; //will reset if bbAttr_ = true
    
    
    ///// read the velocity data from bb file if bbAttr_ = true ////////////
      if (bbAttr_)
        {
          RCP<std::ifstream> bbStream = openFile(bbFilename_, 
                                                    "velocity info");
          Array<string> bbliners = StrUtils::readFile(*bbStream, '#');

          //extract dimension, # solutions, # vertices, solution type 
          //from first line

          int bbdimension;
          int bbnumSolns;
          int bbnumPoints=0;
          int bbsolnType;
          int bbdimensionIndex = 0;
          int bblineIndex = 0;
          //          int bbverticesIndex = 0;
          std::cerr << "bbDimensionIndex = " << bbdimensionIndex << std::endl;

          int bblinerssize = bbliners.size();
          std::cerr << "number of lines in bbliners = " << bblinerssize << std::endl;

          for (int i = bblineIndex; i < bblinerssize; i++)
            {
              Array<string> bbtokens = StrUtils::stringTokenizer(bbliners[i]);
              if(bbtokens.length() > 0) //read first nonblank line
                {
                  if(bbtokens.length() != 4)
                    {
                      std::cerr << 
                        "warning: .bb header line should have 4 tokens, read " 
                           << bbtokens.length() << std::endl;
                    }				   
                  std::cerr << "bbline " << i << ": read bbtokens " << bbtokens[0] 
                       << " " << bbtokens[1] << " " << bbtokens[2] << " " 
                       << bbtokens[3] << std::endl;
                  bbdimensionIndex = i;
                  break;
                }
            }
          bblineIndex = bbdimensionIndex + 1;
          std::cerr << "bblineIndex = " << bblineIndex << std::endl;
          if (bblineIndex > 0)
            {
              Array<string> bbtokens = 
                StrUtils::stringTokenizer(bbliners[bbdimensionIndex]);
              bbdimension = StrUtils::atoi(bbtokens[0]);
              std::cerr << "bbdimension = " << bbdimension << std::endl;
              if (bbdimension != 2) 
                std::cerr << "Error! bbdimension should be 2" << std::endl;
              bbnumSolns = StrUtils::atoi(bbtokens[1]);
              numbbAttr = bbnumSolns; //sets the number of attributes
              std::cerr << "number of solutions per vertex = " << bbnumSolns << std::endl;
              //if (bbnumSolns != nAttributes) 
              //  std::cerr << "Error: bbnumSolns not equal to nAttributes!" << std::endl;
              bbnumPoints = StrUtils::atoi(bbtokens[2]);
              std::cerr << "number of vertices = " << bbnumPoints << std::endl;
              if (bbnumPoints != nPoints) 
                std::cerr << "Error!! number of bb points != nPoints" << std::endl;
              bbsolnType = StrUtils::atoi(bbtokens[3]);
              std::cerr << "bbsolution type = " << bbsolnType  << std::endl;
            }

          //assume solution data starts immediately below header 
          //and has no blank lines
          for (int i = bblineIndex; i < bbnumPoints + bblineIndex; i++) 
            {
              Array<string> bbtokens = StrUtils::stringTokenizer(bbliners[i]);
              if (bbtokens.length() == 0) 
                std::cerr << "warning: encountered a blank soln line" << std::endl;
              int ii = i - bblineIndex;
              //pointAttributes_[ii].resize(nAttributes);
              velVector[ii] = StrUtils::atof(bbtokens[0]);
              velVector[ii+bbnumPoints] = StrUtils::atof(bbtokens[1]);
              //pointAttributes_[ii][0] = StrUtils::atof(bbtokens[0]);
              //pointAttributes_[ii][1] = StrUtils::atof(bbtokens[1]);
      
              if( i < 5 || i > bbnumPoints - 5) //check velocities
                {
                  //std::cerr << "vel0[" << ii << "] = " 
                  //     << pointAttributes_[ii][0] << std::endl;
                  //std::cerr << "vel1[" << ii << "] = " << pointAttributes_[ii][1] 
                  //     << std::endl;
                  std::cerr << "vel0[" << ii << "] = " << velVector[ii] << std::endl;
                  std::cerr << "vel1[" << ii << "] = " << velVector[ii+bbnumPoints] 
                       << std::endl;
                }
	    
            }
        }
      //////////// got velocity vector ///////////////////////
  

      int nAttributes = 0; // value from Triangle .node file
      if (bbAttr_) nAttributes = numbbAttr; // expect two velocities per node
      std::cerr << "nAttributes = " << nAttributes << std::endl;
      //      int nBdryMarkers = 1; // unused - commented out - KL
      //value from Triangle .node file; consistent with Bamg .mesh file

      //Mesh mesh(dimension); /old
      mesh = createMesh(dimension);
      int count=0;
      int offset=0; 
      //initialization--will later set = 1 since Bamg #'s count from 1, not 0
      Array<bool> usedPoint(nPoints);

      //pointAttributes_.resize(nPoints); //old
      nodeAttributes()->resize(nPoints);

      //bool first = true;// unused - commented out - KL

      /* now we can add the point to the mesh */
      for (int i = lineIndex; i < liners.size(); i++)  
        //proceed to read nodes, forget bdry markers
        {
          tokens = StrUtils::stringTokenizer(liners[i]);
          if (tokens.length() > 1)
            {
              double x = atof(tokens[0]); 
              double y = atof(tokens[1]);
              count = i - lineIndex;//assume no blank lines once node data begins
              rawIndices[count] = count;
              if ((i == lineIndex) || (i == lineIndex + 1) || 
                  (i > lineIndex + nPoints - 2))
                {
                  std::cerr << "i = " << i << ";  node = (" << x << "," << y << ")" 
                       << std::endl;
                  std::cerr << "count = " << count << std::endl;
                }

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
	  
              (*nodeAttributes())[count].resize(nAttributes);
              for (int i=0; i<nAttributes; i++)
                {
                  //(*nodeAttributes())[count][i] = atof(tokens[dimension+1+i]);
                  if(i == 0) (*nodeAttributes())[count][i] = velVector[count];
                  if(i == 1) (*nodeAttributes())[count][i] = 
                    velVector[count + nPoints];
                }
            }

          if (tokens[0] == "Edges") 
            {
              edgesIndex = i;
              break;
            }
        }
      std::cerr << "edgesIndex = " << edgesIndex << std::endl;
      std::cerr << "count = " << count << "; npoints - 1 = " << nPoints - 1 << std::endl;
  
      if (count != nPoints - 1) std::cout << "error: # of nodes != # of vertices" 
                                     << std::endl;
      else std::cerr << "successfully read node data" << std::endl;
      lineIndex = edgesIndex + 1;

      // done reading nodes; now read connectivity data.

      std::cerr << "lineIndex = " << lineIndex << "; triangleIndex = " 
           << triangleIndex << std::endl;
      for (int i = lineIndex; i < liners.size(); i++)
        {
          tokens = StrUtils::stringTokenizer(liners[i]);
          if (tokens[0] == "Triangles") {triangleIndex = i; break;}
          if (tokens[0] == "CrackedEdges") {triangleIndex = i+2; break;} 
          //skip over CrackedEdges info 
        }
      std::cerr << "lineIndex = " << lineIndex << "; triangleIndex = " 
           << triangleIndex << std::endl;
      std::cerr << "tokens[0] = " << tokens[0] << std::endl;
    
      lineIndex = triangleIndex + 1;

      Array<int> elemGID;
      Array<int> elemOwner;

      if (triangleIndex > 0)
        {
          tokens = StrUtils::stringTokenizer(liners[lineIndex]);
          std::cerr << "lineIndex = " << lineIndex << "; tokens[0] = " << tokens[0] 
               << std::endl;
          nCells = atoi(tokens[0]);
          elemGID.resize(nCells);
          elemOwner.resize(nCells);
          std::cerr << "nCells = " << nCells << std::endl;
          for (int i=0; i<nCells; i++)
            {
              elemGID[i] = i;
              elemOwner[i] = 0;
            }
          lineIndex = lineIndex + 1;
        }


      // continue to read element data from mesh file//
      /*
        void TriangleMeshReader::readElems(Mesh& mesh,
        const Array<int>& ptGID,
        Array<int>& elemGID,
        Array<int>& elemOwner) const 
      */

      //string line; //already defined
      //Array<string> tokens; //already defined

      ///////////////////////////////////////////////////////////////////
        // skip this part and plug in the corresponding Bamg reader code //
        ///////////////////////////////////////////////////////////////////

        /*
   
        // Open the element file //
        RCP<std::ifstream> elemStream = openFile(elemFilename_, "element info");

        getNextLine(*elemStream, line, tokens, '#');
     
        TEUCHOS_TEST_FOR_EXCEPTION(tokens.length() != 3, std::runtime_error,
        "TriangleMeshReader::getMesh() requires 3 "
        "entries on the header line in "
        "the .ele file. Found line \n[" << line
        << "]\n in file " << elemFilename_);
                   
        int nElems = 0;
        int offset = 1;
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
        // If we're running in parallel, we'd better have consistent numbers
        // of points in the .node and .par file. //
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
        nodes[d] = ptGID[atoi(tokens[d+1])-offset];
        }

        int elemLabel = 0;
        mesh.addElement(elemGID[count], nodes, elemOwner[count], elemLabel);
       
        (*elemAttributes())[count].resize(nAttributes);
        for (int i=0; i<nAttributes; i++)
        {
        (*elemAttributes())[count][i] = atof(tokens[1+ptsPerElem+i]);
        }
        }
        */

        //////// here's the Bamg code for reading triangles (elements) ////////
        std::cerr << "# of triangles = " << nCells 
             << "; starting to read triangle data" << std::endl;
        SUNDANCE_OUT(this->verb() > 3, 
                     "done reading nodes, ready to read elements from " + meshFilename_);

        int nElems = nCells;
        int ptsPerElem = dimension + 1;
        elemGID.resize(nElems);
        elemOwner.resize(nElems);
        offset = 1; //assume Bamg node #'s start with 1, not 0

        for (int i=0; i<nElems; i++)
          {
            elemGID[i] = i;
            elemOwner[i] = 0;
          }
        nAttributes = 0;
        //cellAttributes_.resize(nCells); //old
        elemAttributes()->resize(nElems);
        count = 0;

        int dim = mesh.spatialDim(); //should equal dimension
        Array<int> nodes(dim+1); 
        if (dim != dimension) std::cerr << "ERROR: dim = " << dim << "!= dimension = "
                                   << dimension << std::endl;
      
        std::cerr << "lineIndex = " << lineIndex << std::endl;
        std::cerr << "size of liners = " << liners.size() << std::endl;
        std::cerr << "lineIndex+nCells-1 = " << lineIndex+nCells-1 << std::endl;
        for (int i = lineIndex; i < liners.size();i++) 
          //proceed to read elements, forget bdry markers
          {
            tokens = StrUtils::stringTokenizer(liners[i]);
            if (tokens.length() > 1)
              {
                if (i < lineIndex + 5) std::cerr << "i = " << i << " first triangles: " 
                                            << tokens[0] << " " << tokens[1] 
                                            << " " << tokens[2] << std::endl;
                if (i == lineIndex+nCells-1) std::cerr << "i = " << i 
                                                  << " last triangle: " << tokens[0] 
                                                  << " " << tokens[1] << " " 
                                                  << tokens[2] << std::endl;

                for (int d=0; d<=dim; d++)
                  {
                    //nodes[d] = ptGID[atoi(tokens[d+1])-offset]; 
                    //a Triangle .ele file reads the element no. first, 
                    //so must offset d by 1
                    nodes[d] = ptGID[atoi(tokens[d])-offset]; 
                    //no need to offset d since a Bamg .mesh file doesn't list 
                    //node numbers; offset=1 since Bamg #'s start from 1, not 0
                  }

                count = i - lineIndex; 
                int elemLabel = 0;
                mesh.addElement(elemGID[count], nodes, elemOwner[count], elemLabel);
       
                (*elemAttributes())[count].resize(nAttributes);
                for (int i=0; i<nAttributes; i++)
                  {
                    //(*elemAttributes())[count][i] = atof(tokens[1+ptsPerElem+i]);
                    //offset by 1 since a Triangle .ele file reads the 
                    //element no. first
                    (*elemAttributes())[count][i] = atof(tokens[ptsPerElem+i]);
                    //no need to offset by 1 since a Bamg .mesh file doesn't 
                    //list node numbers;
                  }
              }
            if (tokens[0] == "SubDomainFromMesh") break; //we're done
          }
  
        std::cerr << "# of cells read = " << count + 1 << std::endl;
        if (count != nCells - 1) std::cerr << "error: # of cells read != nCells" << std::endl;

        // done reading elements

        return mesh;
}

///////////////////////////////////////////////////////////////////////////
//ignore this -- we read velocities from bb file within readMesh method
/*

Array<Array<double> > BamgMeshReader::getVelocityField(const std::string& bbFilename) const
// read .bb file for 2-D velocity field & return a ListExpr for the field //
{
RCP<std::ifstream> bbStream = openFile(bbFilename_, "velocity info");
Array<string> liners = StrUtils::readFile(*bbStream, '#');

//extract dimension, # solutions, # vertices, solution type from first line

int dimension;
int numSolns;
int numPoints;
int solnType;
int dimensionIndex = 0;
int lineIndex = 0;
int verticesIndex = 0;
std::cerr << "DimensionIndex = " << dimensionIndex << std::endl;

int linerssize = liners.size();
std::cerr << "number of lines in liners = " << linerssize << std::endl;

for (int i = lineIndex; i < linerssize; i++)
{
Array<string> tokens = StrUtils::stringTokenizer(liners[i]);
//replaced 'PlayaArray' with 'Array'
if(tokens.length() > 0) //read first nonblank line
{
if(tokens.length() != 4)
{
std::cerr << "warning: .bb header line should have 4 tokens, read " 
<< tokens.length() << std::endl;
}				   
std::cerr << "line " << i << ": read tokens " << tokens[0] << " " 
<< tokens[1] << " " << tokens[2] << " " << tokens[3] << std::endl;
dimensionIndex = i;
break;
}
}
lineIndex = dimensionIndex + 1;
std::cerr << "lineIndex = " << lineIndex << std::endl;
if (lineIndex > 0)
{
Array<string> tokens = StrUtils::stringTokenizer(liners[dimensionIndex]);
//replaced 'PlayaArray' with 'Array'
dimension = StrUtils::atoi(tokens[0]);
std::cerr << "dimension = " << dimension << std::endl;
if (dimension != 2) std::cerr << "Error! dimension should be 2" << std::endl;
numSolns = StrUtils::atoi(tokens[1]);
std::cerr << "number of solutions per vertex = " << numSolns << std::endl;
//if (numSolns != nAttributes) 
//  std::cerr << "Error: numSolns not equal to nAttributes!" << std::endl;
numPoints = StrUtils::atoi(tokens[2]);
std::cerr << "number of vetices = " << numPoints << std::endl;
solnType = StrUtils::atoi(tokens[3]);
std::cerr << "solution type = " << solnType  << std::endl;
}

//BasisFamily lagr = new Lagrange(1);
//VectorType<double> petra = new EpetraVectorType();
//Array<BasisFamily> lagrVelBasis;
//for (int n = 0; n < dimension; n++) lagrVelBasis.append(lagr);
//DiscreteSpace dStateSpace(mesh, lagrVelBasis, petra);

//Array<Expr> vel0(numPoints);
//Array<Expr> vel1(numPoints);

Array<double> vel(2*numPoints); //concatenated velocity vectors

//assume solution data starts immediately below header and has no blank lines
for (int i = lineIndex; i < numPoints + lineIndex; i++) 
{
Array<string> tokens = StrUtils::stringTokenizer(liners[i]);
if (tokens.length() == 0) 
std::cerr << "warning: encountered a blank soln line" << std::endl;
int ii = i - lineIndex;
//pointAttributes_[ii].resize(nAttributes);
vel[ii] = StrUtils::atof(tokens[0]);
vel[ii+numPoints] = StrUtils::atof(tokens[1]);
//pointAttributes_[ii][0] = StrUtils::atof(tokens[0]);
//pointAttributes_[ii][1] = StrUtils::atof(tokens[1]);
      
if( i < 5 || i > numPoints - 5) //check to see if we got velocities
{
std::cerr << "vel0[" << ii << "] = " << pointAttributes_[ii][0] << std::endl;
std::cerr << "vel1[" << ii << "] = " << pointAttributes_[ii][1] << std::endl;
}
 
}

return vel;
}

*/
