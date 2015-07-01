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


#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceExodusMeshReader.hpp"
#include "SundanceExodusWriter.hpp"
#include "SundanceTriangleWriter.hpp"
#include "SundanceFileIOChacoPartitioner.hpp"
#include "SundanceVTKWriter.hpp"
#include "Teuchos_CommandLineProcessor.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

using namespace std;

static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

int main(int argc, char** argv)
{
  
  try
  {
    GlobalMPISession session(&argc, &argv);

    CommandLineProcessor clp;
    TimeMonitor t(totalTimer());

    std::string infile="wheel";
    std::string outfile=infile;
    int numProc = 4;
    bool help=false;
    clp.setOption("i", &infile, "Input mesh filename");
    clp.setOption("o", &outfile, "Output mesh filename");
    clp.setOption("np", &numProc, "Number of partitions");
    clp.setOption("h", "nohelp", &help, "Help");
      
    clp.throwExceptions(false);

    CommandLineProcessor::EParseCommandLineReturn rtn 
      = clp.parse(argc, (char**) argv);

    TEUCHOS_TEST_FOR_EXCEPTION(rtn != CommandLineProcessor::PARSE_SUCCESSFUL,
      std::runtime_error,
      "Command-line parsing failed");

    if (help)
    {
      cout << "Usage: partitionExo --i=inputFilename --o=outputFilename "
        "--np=numberOfPartitions" << std::endl;
      cout << "Do not include .exo suffix on filenames" << std::endl;
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPT(infile.length()==0);
      TEUCHOS_TEST_FOR_EXCEPT(outfile.length()==0);
      TEUCHOS_TEST_FOR_EXCEPT(numProc<=1);


      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusMeshReader(infile, meshType);

      Mesh mesh = mesher.getMesh();

      RCP<SerialPartitionerBase> part 
        = rcp(new FileIOChacoPartitioner("part"));

      Array<Mesh> submesh = part->makeMeshParts(mesh, numProc);


      for (int p=0; p<numProc; p++)
      {
        FieldWriter w = new ExodusWriter(outfile);
        w.impersonateParallelProc(numProc, p);
        w.addMesh(submesh[p]);
        w.write();
      }
    }
    TimeMonitor::summarize();
  }
	catch(std::exception& e)
  {
    std::cerr << "Detected exception: " << e.what() << std::endl;
  }
}

