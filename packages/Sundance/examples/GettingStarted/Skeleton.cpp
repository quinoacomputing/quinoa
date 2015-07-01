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

int main(int argc, char** argv)
{
  try
  {
    /* Initialization */
    Sundance::init(&argc, &argv);

    /* ---- BEGIN CODE BODY --- */

    /* The main simulation code goes here. In this example, all we do
     * is to print some information about the processor ranks. */

    MPIComm comm = MPIComm::world();

    /* Print a header from the root processor only. Although this executes on
     * all processors, anything written to the output stream Out::root() 
     * is ignored on all non-root processors (rank != 0). 
     * After writing, synchronize to keep this message from getting jumbled
     * together with the subsequent messages. 
     */
    Out::root() << "Example: getting started" << endl;
    comm.synchronize();

    /* Every processor now speaks up and identifies itself */
    int myRank = comm.getRank();
    int nProc = comm.getNProc();
    Out::os() << "Processor " << myRank 
              << " of " << nProc << " checking in" << endl;

    /* ---- END CODE BODY --- */

    /* Test success or failure. Most examples you'll see will do this 
     * as part of the Trilinos regression testing system. 
     * If you write a simulation code that won't become part of Trilinos,
     * you often can bypass this step.
     *
     * Here the test is a trival one: every processor's rank must be
     * smaller than the total number of processors. If this fails,
     * your MPI installation is probably broken!
     * */
    Sundance::passFailTest(myRank < nProc);
  }
	catch(std::exception& e) /* exception handling */
  {
    cerr << "exception!" << endl;
    Sundance::handleException(e);
  }
  /* Finalization */
  Sundance::finalize();

  return Sundance::testStatus();
}
