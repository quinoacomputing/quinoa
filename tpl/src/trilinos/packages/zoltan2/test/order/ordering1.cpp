// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Zoltan2_OrderingProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <MatrixMarket_Tpetra.hpp>

using Teuchos::RCP;
using namespace std;

/////////////////////////////////////////////////////////////////////////////
// Program to demonstrate use of Zoltan2 to order a TPetra matrix 
// (read from a MatrixMarket file or generated by Galeri::Xpetra).
// Usage:
//     a.out [--inputFile=filename] [--outputFile=outfile] [--verbose] 
//           [--x=#] [--y=#] [--z=#] [--matrix={Laplace1D,Laplace2D,Laplace3D}
// Karen Devine, 2011
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Eventually want to use Teuchos unit tests to vary z2TestLO and
// GO.  For now, we set them at compile time based on whether Tpetra
// is built with explicit instantiation on.  (in Zoltan2_TestHelpers.hpp)

typedef lno_t z2TestLO;
typedef gno_t z2TestGO;
typedef scalar_t Scalar;

typedef KokkosClassic::DefaultNode::DefaultNodeType Node;
typedef Tpetra::CrsMatrix<Scalar, z2TestLO, z2TestGO> SparseMatrix;
typedef Tpetra::Vector<Scalar, z2TestLO, z2TestGO> Vector;

typedef Zoltan2::XpetraCrsMatrixAdapter<SparseMatrix> SparseMatrixAdapter;

#define epsilon 0.00000001

int validatePerm(size_t n, z2TestLO *perm)
// returns 0 if permutation is valid
{
  std::vector<int> count(n);
  int status = 0;
  size_t i;

  for (i=0; i<n; i++)
    count[i]=0;

  for (i=0; i<n; i++){
    if ((perm[i]<0) || (perm[i]>=z2TestLO(n)))
      status = -1;
    else
      count[perm[i]]++;
  }

  // Each index should occur exactly once (count==1)
  for (i=0; i<n; i++){
    if (count[i] != 1){
      status = -2;
      break;
    }
  }

  return status;
}

size_t computeBandwidth(RCP<SparseMatrix> A, z2TestLO *perm)
// Returns the bandwidth of the (local) permuted matrix
// Note we need to use the inverse permutation, but assume
// the direct permutation is passed in.
{
  z2TestLO ii, i, j, k;
  z2TestLO *iperm = 0;
  ArrayView<const z2TestLO> indices;
  ArrayView<const Scalar> values;
  z2TestLO bw_left = 0;
  z2TestLO bw_right = 0;

  z2TestLO  n = A->getNodeNumRows();

  // Construct inverse perm
  if (perm){
    iperm = new z2TestLO [n];
    for (ii=0; ii<n; ii++) {
      iperm[perm[ii]] = ii;
    }
  }

  // Loop over rows of matrix
  for (ii=0; ii<n; ii++) {
    A->getLocalRowView (ii, indices, values);
    for (k=0; k< indices.size(); k++){
      if (perm){
        i = iperm[ii];
        j = iperm[indices[k]];
      } else {
        i = ii;
        j = indices[k];
      }
      if (j-i > bw_right)
        bw_right = j-i;
      if (i-j > bw_left)
        bw_left = i-j;
    }
  }

  if (iperm)
    delete [] iperm;

  // Total bandwidth is the sum of left and right + 1
  return (bw_left + bw_right + 1);
}

/////////////////////////////////////////////////////////////////////////////
int main(int narg, char** arg)
{
  std::string inputFile = "";            // Matrix Market file to read
  std::string outputFile = "";           // Output file to write
  bool verbose = false;                  // Verbosity of output
  int testReturn = 0;

  ////// Establish session.
  Teuchos::GlobalMPISession mpiSession(&narg, &arg, NULL);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  int me = comm->getRank();

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("inputFile", &inputFile,
                 "Name of a Matrix Market file in the data directory; "
                 "if not specified, a matrix will be generated by Galeri.");
  cmdp.setOption("outputFile", &outputFile,
                 "Name of file to write the permutation");
  cmdp.setOption("verbose", "quiet", &verbose,
                 "Print messages and results.");
  cout << "Starting everything" << endl;

  //////////////////////////////////
  // Even with cmdp option "true", I get errors for having these
  //   arguments on the command line.  (On redsky build)
  // KDDKDD Should just be warnings, right?  Code should still work with these
  // KDDKDD params in the create-a-matrix file.  Better to have them where
  // KDDKDD they are used.
  int xdim=10;
  int ydim=10;
  int zdim=10;
  std::string matrixType("Laplace3D");

  cmdp.setOption("x", &xdim,
                "number of gridpoints in X dimension for "
                "mesh used to generate matrix.");
  cmdp.setOption("y", &ydim,
                "number of gridpoints in Y dimension for "
                "mesh used to generate matrix.");
  cmdp.setOption("z", &zdim,              
                "number of gridpoints in Z dimension for "
                "mesh used to generate matrix.");
  cmdp.setOption("matrix", &matrixType,
                "Matrix type: Laplace1D, Laplace2D, or Laplace3D");

  //////////////////////////////////
  // Ordering options to test.
  //////////////////////////////////
  std::string orderMethod("rcm"); // TODO: Allow "RCM" as well
  cmdp.setOption("order_method", &orderMethod,
                "order_method: natural, random, rcm, sorted_degree");
  
  //////////////////////////////////
  cmdp.parse(narg, arg);


  RCP<UserInputForTests> uinput;

  if (inputFile != ""){ // Input file specified; read a matrix
    uinput = rcp(new UserInputForTests(testDataFilePath, inputFile, comm, true));
  }
  else                  // Let Galeri generate a matrix

    uinput = rcp(new UserInputForTests(xdim, ydim, zdim, matrixType, comm, true));

  RCP<SparseMatrix> origMatrix = uinput->getTpetraCrsMatrix();

  if (me == 0) 
    cout << "NumRows     = " << origMatrix->getGlobalNumRows() << endl
         << "NumNonzeros = " << origMatrix->getGlobalNumEntries() << endl
         << "NumProcs = " << comm->getSize() << endl;

  ////// Create a vector to use with the matrix.
  RCP<Vector> origVector, origProd;
  origProd   = Tpetra::createVector<Scalar,z2TestLO,z2TestGO>(
                                    origMatrix->getRangeMap());
  origVector = Tpetra::createVector<Scalar,z2TestLO,z2TestGO>(
                                    origMatrix->getDomainMap());
  origVector->randomize();

  ////// Specify problem parameters
  Teuchos::ParameterList params;
  params.set("order_method", orderMethod);

  ////// Create an input adapter for the Tpetra matrix.
  SparseMatrixAdapter adapter(origMatrix);

  ////// Create and solve ordering problem
  try
  {
  Zoltan2::OrderingProblem<SparseMatrixAdapter> problem(&adapter, &params);
  cout << "Going to solve" << endl;
  problem.solve();

  ////// Basic metric checking of the ordering solution
  size_t checkLength;
  z2TestGO *checkGIDs;
  z2TestLO *checkPerm;
  Zoltan2::OrderingSolution<z2TestGO, z2TestLO> *soln = problem.getSolution();

  cout << "Going to get results" << endl;
  // Check that the solution is really a permutation
  checkLength = soln->getPermutationSize();
  checkGIDs = soln->getGids(); // NOT computed; this should be done in ApplyOrderingSolution()
  checkPerm = soln->getPermutation();

  if (outputFile != "") {
    ofstream permFile;

    // Write permutation (0-based) to file
    permFile.open(outputFile.c_str());
    for (size_t i=0; i<checkLength; i++){
      permFile << " " << checkPerm[i] << endl;
    }
    permFile.close();

  }

  cout << "Going to validate the soln" << endl;
  // Verify that checkPerm is a permutation
  testReturn = validatePerm(checkLength, checkPerm);

  cout << "Going to compute the bandwidth" << endl;
  // Compute original bandwidth 
  cout << "Original Bandwidth: " << computeBandwidth(origMatrix, 0) << endl;
  // Compute permuted bandwidth
  cout << "Permuted Bandwidth: " << computeBandwidth(origMatrix, checkPerm) << endl;

  } catch (std::exception &e){
      if (comm->getSize() != 1)
      {
          std::cout << "Ordering does not support distributed matrices."
             << std::endl;
          std::cout << "PASS" << std::endl;
      }
      else
      {
          std::cout << "Exception caught in ordering" << std::endl;
          std::cout << e.what() << std::endl;
          std::cout << "FAIL" << std::endl;
      }
      return 0;
  }

  if (me == 0) {
    if (testReturn)
      std::cout << "Solution is not a permutation; FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

}

