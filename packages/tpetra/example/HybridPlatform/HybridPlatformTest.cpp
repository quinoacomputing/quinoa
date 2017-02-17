/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_HybridPlatform.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_FileInputSource.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLParameterListReader.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_MatrixIO.hpp>

std::string fnMatrix("bcsstk17.rsa");
bool testPassed;

template <class Node, class Scalar, class Ordinal>
Scalar power_method(const Teuchos::RCP<const Tpetra::Operator<Scalar,Ordinal,Ordinal,Node> > &A, size_t niters, typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance, bool verbose) {
  using Teuchos::RCP;
  typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node> Vector;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  // create three vectors
  RCP<Vector> z = Tpetra::createVector<Scalar>(A->getRangeMap()),
              q = Tpetra::createVector<Scalar>(A->getRangeMap()),
              r = Tpetra::createVector<Scalar>(A->getRangeMap());
  // Fill z with random numbers
  z->randomize();
  // Variables needed for iteration
  const Scalar ONE  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar lambda = static_cast<Scalar>(0.0);
  Magnitude normz, residual = static_cast<Magnitude>(0.0);
  // power iteration
  for (size_t iter = 0; iter < niters; ++iter) {
    normz = z->norm2();                            // Compute 2-norm of z
    q->scale(ONE/normz, *z);                       // Set q = z / normz
    A->apply(*q, *z);                              // Compute z = A*q
    lambda = q->dot(*z);                           // Approximate maximum eigenvalue: lamba = dot(q,z)
    if ( iter % 100 == 0 || iter + 1 == niters ) {
      r->update(ONE, *z, -lambda, *q, ZERO);       // Compute A*q - lambda*q
      residual = Teuchos::ScalarTraits<Scalar>::magnitude(r->norm2() / lambda);
      if (verbose) {
        std::cout << "Iter = " << iter << "  Lambda = " << lambda 
                  << "  Residual of A*q - lambda*q = " 
                  << residual << std::endl;
      }
    } 
    if (residual < tolerance) {
      break;
    }
  }
  return lambda;
}

template <class Node>
class runTest {
  public:
  static void run(Teuchos::ParameterList &myMachPL, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) {
    using std::cout; 
    using std::endl;
    cout << "Running test with Node==" << Teuchos::typeName(*node) << " on rank " << comm->getRank() << "/" << comm->getSize() << endl;
    //
    // Get the data from the HB file and build the Map,Matrix
    //
    // we prefer float for this test, as it is more likely to exercise the GPU
#if   defined(HAVE_TPETRA_INST_FLOAT)
    typedef float  TestScalar;
#elif defined(HAVE_TPETRA_INST_DOUBLE)
    typedef double TestScalar;
#endif
    Teuchos::RCP< Tpetra::CrsMatrix<TestScalar,int,int,Node> > A;
    try {
      Tpetra::Utils::readHBMatrix(fnMatrix,comm,node,A);
    }
    catch (std::runtime_error &e) {
      if (comm->getRank() == 0) {
        cout << "Tpetra::Utils::readHBMatrix() threw exception: " << endl << e.what() << endl;
      }
      testPassed = false;      
      return;
    }
    (void)power_method<Node,TestScalar,int>(A,100,1e-4f,comm->getRank() == 0);
    testPassed = true;
  }
};

int main(int argc, char **argv) {
  using std::string;
  using std::cout;
  using std::endl;
  using Teuchos::FileInputSource;
  using Teuchos::XMLObject;
  using Teuchos::ParameterList;
  using Teuchos::XMLParameterListReader;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;

  Teuchos::GlobalMPISession mpisess(&argc,&argv,&cout);
  RCP<const Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));

  //
  // Get test parameters from command-line processor
  //  
  Teuchos::CommandLineProcessor cmdp(false,true);
  string fnMachine("mpionly.xml");
  cmdp.setOption("matrix-file",&fnMatrix,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("machine-file",&fnMachine,"Filename for XML machine description file.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // read machine file and initialize platform
  // 
  FileInputSource fileSrc(fnMachine);
  XMLObject machXML = fileSrc.getObject();
  XMLParameterListReader pl2xml;
  ParameterList machPL = pl2xml.toParameterList(machXML);
  Tpetra::HybridPlatform platform(comm,machPL);
  platform.runUserCode<runTest>();

  if (testPassed == false) {
    if (comm->getRank() == 0) {
      cout << "End Result: TEST FAILED" << endl;
      return -1;
    }
  }

  if (comm->getRank() == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;
}
