// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
//  This test is for the block Krylov-Schur eigensolver
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziBlockKrylovSchur.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziSolverUtils.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

using namespace Teuchos;
using namespace Anasazi;

typedef double                              ScalarType;
typedef ScalarTraits<ScalarType>                   SCT;
typedef SCT::magnitudeType               MagnitudeType;
typedef Epetra_MultiVector                 MV;
typedef Epetra_Operator                    OP;
typedef MultiVecTraits<ScalarType,MV>     MVT;
typedef OperatorTraits<ScalarType,MV,OP>  OPT;

class get_out : public std::logic_error {
  public: get_out(const std::string &whatarg) : std::logic_error(whatarg) {}
};

void checks( RCP<BlockKrylovSchur<ScalarType,MV,OP> > solver, int blocksize, int numblocks, 
             RCP<Eigenproblem<ScalarType,MV,OP> > problem,
             RCP<MatOrthoManager<ScalarType,MV,OP> > ortho,
             SolverUtils<ScalarType,MV,OP> &msutils) {
  BlockKrylovSchurState<ScalarType,MV> state = solver->getState();

  // Remember that block Krylov-Schur needs to keep an extra vector for F, 
  // the next block of the factorization (AV=VH+FB^T).
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.V)-solver->getBlockSize() != solver->getMaxSubspaceDim(),get_out,"getMaxSubspaceDim() does not match allocated size for V");

  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize, get_out,"Solver block size does not match specified block size.");  

  TEUCHOS_TEST_FOR_EXCEPTION(&solver->getProblem() != problem.get(),get_out,"getProblem() did not return the submitted problem.");

  if (solver->getProblem().isHermitian()) {
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim() != blocksize*numblocks,get_out,
        "BlockKrylovSchur::getMaxSubspaceDim() should be blocksize times the number of blocks.");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim() != blocksize*numblocks+1,get_out,
        "BlockKrylovSchur::getMaxSubspaceDim() should be one vector more than the blocksize times the number of blocks.");
  }

  if (solver->isInitialized()) 
  {
    // check Ritz values
    solver->computeRitzValues();
    
    std::vector<Anasazi::Value<ScalarType> > ritzValues = solver->getRitzValues();
    
    // check Ritz residuals
    std::vector<MagnitudeType> ritzResids = solver->getRitzRes2Norms();
    
    // get Ritz index
    std::vector<int> ritzIndex = solver->getRitzIndex();
    
    // get Ritz vector
    RCP<const MV> ritzVectors = solver->getRitzVectors();
    
    // check Ritz vector
    const int numRitzVecs = MVT::GetNumberVecs( *ritzVectors );
    if (solver->getCurSubspaceDim() >= numRitzVecs ) {
      solver->computeRitzVectors();
    }

    if (solver->isRitzValsCurrent() && solver->isRitzVecsCurrent()) {
      
      // Compute Ritz residuals like R = OP*X - X*T 
      Teuchos::SerialDenseMatrix<int,ScalarType> T(numRitzVecs,numRitzVecs);
      Teuchos::RCP<MV> ritzResiduals = MVT::Clone( *ritzVectors, numRitzVecs );
      for (int i=0; i<T.numRows(); i++) T(i,i) = ritzValues[i].realpart;
      OPT::Apply( *(problem->getOperator()), *ritzVectors, *ritzResiduals );
      MVT::MvTimesMatAddMv(-1.0,*ritzVectors,T,1.0,*ritzResiduals);
      
      // Compute the norm of the Ritz residual vectors
      std::vector<MagnitudeType> ritzVecNrm( numRitzVecs );
      MVT::MvNorm( *ritzVectors, ritzVecNrm );
      MagnitudeType error;
      for (int i=0; i<numRitzVecs; i++) {
        error = Teuchos::ScalarTraits<MagnitudeType>::magnitude( ritzVecNrm[i] - 1.0 );
        TEUCHOS_TEST_FOR_EXCEPTION(error > 1e-14,get_out,"Ritz vectors are not normalized.");
      }      

      /* TO DO: Fix the iteration to compute residuals in the event of a non-euclidean normalization of the basis.

      std::vector<MagnitudeType> ritzResNrm( MVT::GetNumberVecs( *ritzResiduals ) );
      MVT::MvNorm( *ritzResiduals, ritzResNrm );
      for (int i=0; i<(int)ritzResNrm.size(); i++) {
      error = Teuchos::ScalarTraits<MagnitudeType>::magnitude( ritzResids[i] - ritzResNrm[i] );
      cout << error/ritzResNrm[i] << endl;
      TEUCHOS_TEST_FOR_EXCEPTION(error > 1e-13,get_out,"Ritz residuals from iteration do not compare to those computed.");
      }
      */
    }
  }
  else {
    // not initialized
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 0,get_out,"In unitialized state, getCurSubspaceDim() should be 0.");
  }
}

void testsolver( RCP<BasicEigenproblem<ScalarType,MV,OP> > problem,
                 RCP< OutputManager<ScalarType> > printer,
                 RCP< MatOrthoManager<ScalarType,MV,OP> > ortho,
                 RCP< SortManager<MagnitudeType> > sorter,
                 ParameterList &pls,bool invalid,
                 BlockKrylovSchurState<ScalarType,MV> initstate, bool invalidinit)
{
  // create a status tester
  RCP< StatusTest<ScalarType,MV,OP> > tester = rcp( new StatusTestMaxIters<ScalarType,MV,OP>(1) );

  // create the solver
  RCP< BlockKrylovSchur<ScalarType,MV,OP> > solver;
  try {
    solver = rcp( new BlockKrylovSchur<ScalarType,MV,OP>(problem,sorter,printer,tester,ortho,pls) );
    TEUCHOS_TEST_FOR_EXCEPTION(invalid, get_out, "Instantiating with invalid parameters failed to throw exception.")
  }
  catch (const std::invalid_argument &ia) {
    if (!invalid) {
      printer->stream(Warnings) << "Error thrown at instantiation: " << ia.what() << endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!invalid, get_out, "Instantiating with valid parameters unexpectadly threw exception.");

    // caught expected exception
    return;
  }

  const int  blocksize = pls.get<int>("Block Size");
  const int  numblocks = pls.get<int>("Num Blocks");
  const int  numritzvecs = pls.get<int>("Number of Ritz Vectors");

  SolverUtils<ScalarType,MV,OP> msutils;

  // solver should be uninitialized
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != false,get_out,"Solver should be un-initialized after instantiation.");  
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations after initialization should be zero after init.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  checks(solver,blocksize,numblocks,problem,ortho,msutils);

  // initialize solver and perform checks
  try {
    solver->initialize(initstate);
    TEUCHOS_TEST_FOR_EXCEPTION(invalidinit, get_out, "Initializing with invalid data failed to throw exception.")
  }
  catch (const std::invalid_argument &ia) {
    TEUCHOS_TEST_FOR_EXCEPTION(!invalidinit, get_out, "Initializing with valid data unexpectadly threw exception.");
    // caught expected exception
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");  
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*(solver->getRitzVectors())) != numritzvecs,get_out,"Number of Ritz vectors in storage is incorrect.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 0,get_out,"after init, getCurSubspaceDim() should be zero, only the kernel was generated.");
  checks(solver,blocksize,numblocks,problem,ortho,msutils);

  // call iterate(); solver should perform exactly one iteration and return; status test should be passed
  solver->iterate();
  TEUCHOS_TEST_FOR_EXCEPTION(tester->getStatus() != Passed,get_out,"Solver returned from iterate() but getStatus() not Passed.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");  
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 1,get_out,"Number of iterations should be zero.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*(solver->getRitzVectors())) != numritzvecs,get_out,"Number of Ritz vectors in storage is incorrect.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != blocksize,get_out,"after one step, getCurSubspaceDim() should be blocksize.");
  checks(solver,blocksize,numblocks,problem,ortho,msutils);

  // reset numiters, call iterate(); solver should perform exactly one iteration and return; status test should be passed
  solver->resetNumIters();
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero after resetNumIters().")
  solver->iterate();
  TEUCHOS_TEST_FOR_EXCEPTION(tester->getStatus() != Passed,get_out,"Solver returned from iterate() but getStatus() not Passed.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");  
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 1,get_out,"Number of iterations should be zero.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*(solver->getRitzVectors())) != numritzvecs,get_out,"Number of Ritz vectors in storage is incorrect.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 2*blocksize,get_out,"after two steps, getCurSubspaceDim() should be 2*blocksize.");
  checks(solver,blocksize,numblocks,problem,ortho,msutils);
}

int main(int argc, char *argv[]) 
{

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool testFailed;
  bool verbose = false;
  bool debug = false;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging output from iteration.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (debug) verbose = true;

  // create the output manager
  int verbosity = Anasazi::Errors;
  if (verbose) {
    verbosity += Anasazi::Warnings;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }
  RCP< OutputManager<ScalarType> > printer = 
    rcp( new BasicOutputManager<ScalarType>( verbosity ) );

  printer->stream(Debug) << Anasazi_Version() << endl;

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100+1;

  // Create problem
  RCP<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RCP<const Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RCP<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  //
  // Create the initial vectors
  const int nev = 5;
  RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), nev) );
  ivec->Random();
  //
  // Create eigenproblem: one standard and one generalized
  RCP<BasicEigenproblem<ScalarType,MV,OP> > probstd = rcp( new BasicEigenproblem<ScalarType, MV, OP>(K, ivec) );
  RCP<BasicEigenproblem<ScalarType,MV,OP> > probgen = rcp( new BasicEigenproblem<ScalarType, MV, OP>(K, M, ivec) );
  //
  // Inform the eigenproblem that the operator A is not symmetric (even though it is)
  probstd->setHermitian(false);
  probgen->setHermitian(false);
  //
  // Set the number of eigenvalues requested
  probstd->setNEV( nev );
  probgen->setNEV( nev );
  //
  // Inform the eigenproblem that you are finishing passing it information
  if ( probstd->setProblem() != true || probgen->setProblem() != true ) {
    printer->stream(Warnings) << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
      << "End Result: TEST FAILED" << endl;
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // create the orthogonalization managers: one standard and one M-based
  RCP< MatOrthoManager<ScalarType,MV,OP> > orthostd = rcp( new SVQBOrthoManager<ScalarType,MV,OP>() );
  RCP< MatOrthoManager<ScalarType,MV,OP> > orthogen = rcp( new SVQBOrthoManager<ScalarType,MV,OP>(M) );
  // create the sort manager
  RCP< SortManager<MagnitudeType> > sorter = rcp( new BasicSort<MagnitudeType>("LM") );
  // create the parameter list specifying blocksize > nev and full orthogonalization
  ParameterList pls;

  // begin testing 
  testFailed = false;

  try 
  {
    BlockKrylovSchurState<ScalarType,MV> istate;

    pls.set<int>("Block Size",nev);
    pls.set<int>("Num Blocks",3);
    pls.set<int>("Step Size", 2);
    pls.set<int>("Number of Ritz Vectors",nev);
    printer->stream(Warnings) << "Testing solver(nev,3) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Num Blocks",3);
    printer->stream(Warnings) << "Testing solver(nev,3) with generalized eigenproblem..." << endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",nev);
    pls.set<int>("Num Blocks",4);
    printer->stream(Warnings) << "Testing solver(nev,4) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Num Blocks",4);
    printer->stream(Warnings) << "Testing solver(nev,4) with generalized eigenproblem..." << endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",2*nev);
    pls.set<int>("Num Blocks",3);
    printer->stream(Warnings) << "Testing solver(2*nev,3) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Num Blocks",3);
    printer->stream(Warnings) << "Testing solver(2*nev,3) with generalized eigenproblem..." << endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",2*nev);
    pls.set<int>("Num Blocks",4);
    printer->stream(Warnings) << "Testing solver(2*nev,4) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Num Blocks",4);
    printer->stream(Warnings) << "Testing solver(2*nev,4) with generalized eigenproblem..." << endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",nev/2);
    pls.set<int>("Num Blocks",3);
    printer->stream(Warnings) << "Testing solver(nev/2,3) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Num Blocks",3);
    printer->stream(Warnings) << "Testing solver(nev/2,3) with generalized eigenproblem..." << endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",nev/2);
    pls.set<int>("Num Blocks",4);
    printer->stream(Warnings) << "Testing solver(nev/2,4) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Num Blocks",4);
    printer->stream(Warnings) << "Testing solver(nev/2,4) with generalized eigenproblem..." << endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    // try with an invalid blocksize
    pls.set<int>("Block Size",0);
    pls.set<int>("Num Blocks",4);
    printer->stream(Warnings) << "Testing solver(0,4) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with an invalid numblocks: invalid because it is less than the minimum (3) allowed
    pls.set<int>("Block Size",nev);
    pls.set<int>("Num Blocks",2);
    printer->stream(Warnings) << "Testing solver(4,2) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with a too-large subspace: Hermitian
    // subspace will be BlockSize*NumBlocks
    probstd->setHermitian(true);
    probstd->setProblem();
    pls.set<int>("Block Size",4);
    pls.set<int>("Num Blocks",100/4+1);
    printer->stream(Warnings) << "Testing solver(4,toomany,Hermitian) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with a too-large subspace: non-Hermitian, blocksize != 1
    // allocated subspace will be BlockSize*NumBlocks+1 >= 100
    probstd->setHermitian(false);
    probstd->setProblem();
    pls.set<int>("Block Size",4);
    pls.set<int>("Num Blocks",25);
    printer->stream(Warnings) << "Testing solver(4,toomany,non-Hermitian) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with a too-large subspace: non-Hermitian, blocksize == 1
    // allocated subspace will be BlockSize*NumBlocks+1 >= 100
    probstd->setHermitian(false);
    probstd->setProblem();
    pls.set<int>("Block Size",1);
    pls.set<int>("Num Blocks",100);
    printer->stream(Warnings) << "Testing solver(1,toomany,non-Hermitian) with standard eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with an unset problem
    // setHermitian will mark the problem as unset
    probstd->setHermitian(true);
    printer->stream(Warnings) << "Testing solver with unset eigenproblem..." << endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);
    // set problem: now hermitian
    probstd->setProblem();

    // try with a too-small initial basis
    printer->stream(Warnings) << "Initializing solver with too-small basis..." << endl;
    pls.set("Block Size",4);
    pls.set("Num Blocks",3);
    istate.V  = MVT::Clone(*ivec,3);
    istate.H = rcp( new SerialDenseMatrix<int,ScalarType>(3,3) );
    istate.curDim = 3;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,true);

    // try with a too-large initial basis
    // problem is Hermitian, maxSubSpaceSize is BlockSize*NumBlocks
    printer->stream(Warnings) << "Initializing solver with too-large basis..." << endl;
    pls.set("Block Size",4);
    pls.set("Num Blocks",3);
    istate.V  = MVT::Clone(*ivec,4*3+1);
    istate.H = rcp( new SerialDenseMatrix<int,ScalarType>(13,13) );
    istate.curDim = 13;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,true);

    // try with a inconsistent curDim and H
    printer->stream(Warnings) << "Initializing solver with inconsistent H and V..." << endl;
    pls.set("Block Size",4);
    pls.set("Num Blocks",3);
    istate.V  = MVT::Clone(*ivec,4);
    istate.H = rcp( new SerialDenseMatrix<int,ScalarType>(3,3) );
    istate.curDim = 4;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,true);

    // create a dummy status tester
    RCP< StatusTest<ScalarType,MV,OP> > dumtester = rcp( new StatusTestMaxIters<ScalarType,MV,OP>(1) );

    // try with a null problem
    printer->stream(Warnings) << "Testing solver with null eigenproblem..." << endl;
    try {
      RCP< BlockKrylovSchur<ScalarType,MV,OP> > solver 
        = rcp( new BlockKrylovSchur<ScalarType,MV,OP>(Teuchos::null,sorter,printer,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null sortman
    printer->stream(Warnings) << "Testing solver with null sort manager..." << endl;
    try {
      RCP< BlockKrylovSchur<ScalarType,MV,OP> > solver 
        = rcp( new BlockKrylovSchur<ScalarType,MV,OP>(probstd,Teuchos::null,printer,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a output man problem
    printer->stream(Warnings) << "Testing solver with null output manager..." << endl;
    try {
      RCP< BlockKrylovSchur<ScalarType,MV,OP> > solver 
        = rcp( new BlockKrylovSchur<ScalarType,MV,OP>(probstd,sorter,Teuchos::null,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null status test
    printer->stream(Warnings) << "Testing solver with null status test..." << endl;
    try {
      RCP< BlockKrylovSchur<ScalarType,MV,OP> > solver 
        = rcp( new BlockKrylovSchur<ScalarType,MV,OP>(probstd,sorter,printer,Teuchos::null,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null orthoman
    printer->stream(Warnings) << "Testing solver with null ortho manager..." << endl;
    try {
      RCP< BlockKrylovSchur<ScalarType,MV,OP> > solver 
        = rcp( new BlockKrylovSchur<ScalarType,MV,OP>(probstd,sorter,printer,dumtester,Teuchos::null,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

  }
  catch (const get_out &go) {
    printer->stream(Errors) << "Test failed: " << go.what() << endl;
    testFailed = true;
  }
  catch (const std::exception &e) {
    printer->stream(Errors) << "Caught unexpected exception: " << e.what() << endl;
    testFailed = true;
  }

  
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    printer->stream(Warnings) << endl << "End Result: TEST FAILED" << endl;
    return -1;
  }
  //
  // Default return value
  //
  printer->stream(Warnings) << endl << "End Result: TEST PASSED" << endl;
  return 0;

}
