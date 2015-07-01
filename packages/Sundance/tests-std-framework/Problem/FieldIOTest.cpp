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



#include "Sundance.hpp"
#include "SundanceMeshIOUtils.hpp"




#if defined(HAVE_SUNDANCE_EXODUS)

int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);


    MPIComm world = MPIComm::world();
#ifdef HAVE_MPI
    MPIComm self = MPIComm(MPI_COMM_SELF);
#else 
    MPIComm self = world;
#endif

    std::string infile = "tetPrism";
//    std::string infile = "./prism";
    int numProc = world.getNProc();
    

#ifdef HAVE_SUNDANCE_CHACO
    Out::os() << "have chaco!" << std::endl;
    /* If in parallel, partition the mesh */
    if (world.getNProc() > 1)
    {
      /* Partition on the root only */
      if (world.getRank()==0)
      {
        Out::os() << "root processor is working on the partitioning" << std::endl;
        MeshType meshType = new BasicSimplicialMeshType();
        MeshSource mesher = new ExodusMeshReader(infile, meshType, self);

        RCP<SerialPartitionerBase> part 
          = rcp(new FileIOChacoPartitioner("part"));
        
        serialPartition(part, numProc, mesher, infile);
      }
      /* everyone else waits for the root processor to finish partitioning */
      MPIComm::world().synchronize();
    }

    /* Now do a readback test on the mesh */
    MPIComm::world().synchronize();
    Out::root() << "starting readback" << endl;
    MPIComm::world().synchronize();
    double err = readbackTester(infile, world,4);
#else
    Out::root() << "starting readback" << endl;
    double err = 0.0;
    if (world.getRank()==0)
    {
      Out::os() << "partitioning..." << std::endl;
      err = readbackTester(infile, self, 2);
      Out::os() << "done partitioning..." << std::endl;
    }
    else
    {
      Out::os() << "waiting..." << std::endl;
    }
    Out::os() << "synching..." << std::endl;
    MPIComm::world().synchronize();
    Out::os() << "sharing erro..." << std::endl;
    MPIComm::world().bcast(&err, 1, MPIDataType::intType(), 0);
#endif
    
    /* Use a very tight tolerance because the readback should be
     * essentially exact */
    Sundance::passFailTest(err, 1.0e-15);
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}
    

#else // don't have exodus


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy FieldIOTest PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif

    
    

    
    
