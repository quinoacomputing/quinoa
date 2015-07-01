/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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

#include "PlayaEpetraVectorType.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaOut.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaVectorOpsImpl.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaMPIComm.hpp"


/* ------------------------------------------------------------------------
 *
 * This example shows the construction and use of simple vectors
 * 
 * ------------------------------------------------------------------------ */

using namespace Playa;

int main(int argc, char *argv[])
{
  int rtn = 0;

  try
  {
    /* Initialize MPI */
    GlobalMPISession session(&argc, &argv);
    Tabs::showDepth() = false;


    /* The VectorType object will be used when we create vector spaces, 
     * specifying what type of low-level linear algebra implementation
     * will be used. */
    VectorType<double> vecType = new EpetraVectorType();

    /* Construct a vector space  */
    int n = 4;
    VectorSpace<double> vs = vecType.createEvenlyPartitionedSpace(MPIComm::world(), n);

    /* */
    Rand::setLocalSeed(MPIComm::world(), 12345);

    /* Make some vectors */
    Vector<double> x = vs.createMember();
    x.randomize();

    Vector<double> y = vs.createMember();
    y.randomize();

    Vector<double> z = vs.createMember();
    z.randomize();

    Vector<double> u = vs.createMember();
    u.randomize();

    Vector<double> v = vs.createMember();
    v.randomize();

    
    double err = 1.0e10;
    double tol = 1.0e-10;
    int fails = 0;

    /* Print a vector */
    Out::root() << "x = " << x.description() << endl;
    Out::os() << x << endl;
    /* Make sure we've not created a zero-norm random vector */
    Out::root() << "||x||=" << norm2(x) << endl;

    


    /* Test the deep copying of a vector */
    Vector<double> xCopy = x.copy();
    err = norm2(x-xCopy);
    Out::root() << "||copy error|| = " << err << endl;
    if (!(err < tol)) fails++;

    /* Do the computation w/o making a temporary */
    err = norm2Dist(x, xCopy); 
    Out::root() << "||copy error|| = " << err << endl;
    if (!(err < tol)) fails++;

    /* Test scalar multiplication */
    x *= 2.0;
 
    err = ::fabs( (x-xCopy).norm2() - xCopy.norm2() );
    Out::root() << "|norm(2x-x) - norm(x)| = " << err << endl;
    if (!(err < tol)) fails++;

    /* Test scalar division */
    x /= 2.0;
 
    err = norm2Dist(x, xCopy);
    Out::root() << "norm(2x/2-x)= " << err << endl;
    if (!(err < tol)) fails++;

    
     
    if (fails == 0)
    {
      Out::root() << "test PASSED" << endl;
    }
    else
    {
      Out::root() << "test FAILED" << endl;
    }
  }
  catch(std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << endl;
    rtn = -1;
  }
  return rtn;
}
