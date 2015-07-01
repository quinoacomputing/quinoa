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

#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaDenseSerialMatrixFactory.hpp"
#include "PlayaSerialVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"  
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_BLAS.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaVectorImpl.hpp"
#endif

extern "C"
{
void dgesv_(int *n, int *nrhs, double *a, int* lda, 
  int *ipiv, double *b, int *ldb, int *info);

void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a,
  int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
  double* work, int* lwork, int* info );
}
using std::max;
using std::min;

using namespace Playa;
using namespace Teuchos;

using std::setw;

DenseSerialMatrix::DenseSerialMatrix(
  const VectorSpace<double>& domain,
  const VectorSpace<double>& range)
  : LinearOpWithSpaces<double>(domain, range),
    nRows_(range.dim()),
    nCols_(domain.dim()),
    data_(nRows_*nCols_)
{}


void DenseSerialMatrix::apply(
  Teuchos::ETransp transApplyType,
  const Vector<double>& in,
  Vector<double> out) const
{
  const SerialVector* rvIn = SerialVector::getConcrete(in);
  SerialVector* rvOut = SerialVector::getConcrete(out);

  Teuchos::BLAS<int, double> blas;
  int lda = numRows();
  blas.GEMV(transApplyType, numRows(), numCols(), 1.0, dataPtr(), 
    lda, rvIn->dataPtr(), 1, 1.0, rvOut->dataPtr(), 1);
}

void DenseSerialMatrix::addToRow(int globalRowIndex,
  int nElemsToInsert,
  const int* globalColumnIndices,
  const double* elementValues)
{
  int r = globalRowIndex;
  for (int k=0; k<nElemsToInsert; k++)
  {
    int c = globalColumnIndices[k];
    double x = elementValues[k];
    data_[r + c*numRows()] = x;
  }
}

void DenseSerialMatrix::zero()
{
  for (int i=0; i<data_.size(); i++) data_[i] = 0.0;
}


void DenseSerialMatrix::print(std::ostream& os) const
{
  if (numCols() <= 4)
  {
    for (int i=0; i<numRows(); i++)
    {
      for (int j=0; j<numCols(); j++)
      {
        os << setw(16) << data_[i+numRows()*j];
      }
      os << std::endl;
    }
  }
  else
  {
    for (int i=0; i<numRows(); i++)
    {
      for (int j=0; j<numCols(); j++)
      {
        os << setw(6) << i << setw(6) << j << setw(20) << data_[i+numRows()*j]
           << std::endl;
      }
    }
  }
}

void DenseSerialMatrix::setRow(int row, const Array<double>& rowVals)
{
  TEUCHOS_TEST_FOR_EXCEPT(rowVals.size() != numCols());
  TEUCHOS_TEST_FOR_EXCEPT(row < 0);
  TEUCHOS_TEST_FOR_EXCEPT(row >= numRows());

  for (int i=0; i<rowVals.size(); i++)
  {
    data_[row+numRows()*i] = rowVals[i];
  }
}


namespace Playa
{


SolverState<double> denseSolve(const LinearOperator<double>& A,
  const Vector<double>& b,
  Vector<double>& x)
{
  const DenseSerialMatrix* Aptr 
    = dynamic_cast<const DenseSerialMatrix*>(A.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(Aptr==0);
  /* make a working copy, because dgesv will overwrite the matrix */
  DenseSerialMatrix tmp = *Aptr;
  /* Allocate a vector for the solution */
  x = b.copy();
  SerialVector* xptr 
    = dynamic_cast<SerialVector*>(x.ptr().get());
  
  int N = Aptr->numRows();
  int nRHS = 1;
  int LDA = N;
  Array<int> iPiv(N);
  int LDB = N;
  int info = 0;
  dgesv_(&N, &nRHS, tmp.dataPtr(), &LDA, &(iPiv[0]), xptr->dataPtr(),
    &LDB, &info);

  if (info == 0)
  {
    return SolverState<double>(SolveConverged, "solve OK",
      0, 0.0);
  }
  else 
  {
    return SolverState<double>(SolveCrashed, "solve crashed with dgesv info="
      + Teuchos::toString(info),
      0, 0.0);
  }
}


void denseSVD(const LinearOperator<double>& A,
  LinearOperator<double>& U,  
  Vector<double>& Sigma,
  LinearOperator<double>& Vt)
{
  VectorSpace<double> mSpace = A.range();
  VectorSpace<double> nSpace = A.domain();

  const DenseSerialMatrix* Aptr 
    = dynamic_cast<const DenseSerialMatrix*>(A.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(Aptr==0);
  /* make a working copy, because dgesvd will overwrite the matrix */
  DenseSerialMatrix ATmp = *Aptr;

  int M = ATmp.numRows();
  int N = ATmp.numCols();
  int S = min(M, N);
  
  VectorSpace<double> sSpace;
  if (S==M) sSpace = mSpace;
  else sSpace = nSpace;

  Sigma = sSpace.createMember();
  SerialVector* sigPtr
    = dynamic_cast<SerialVector*>(Sigma.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(sigPtr==0);

  DenseSerialMatrixFactory umf(sSpace, mSpace);
  DenseSerialMatrixFactory vtmf(nSpace, sSpace);
  
  U = umf.createMatrix();
  Vt = vtmf.createMatrix();

  DenseSerialMatrix* UPtr 
    = dynamic_cast<DenseSerialMatrix*>(U.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(UPtr==0);

  DenseSerialMatrix* VtPtr 
    = dynamic_cast<DenseSerialMatrix*>(Vt.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(VtPtr==0);
  
  double* uData = UPtr->dataPtr();
  double* vtData = VtPtr->dataPtr();
  double* aData = ATmp.dataPtr();
  double* sData = sigPtr->dataPtr();

  char jobu = 'S';
  char jobvt = 'S';
 
  int LDA = M;
  int LDU = M;
  int LDVT = S;

  int LWORK = max(1, max(3*min(M,N)+max(M,N), 5*min(M,N)));
  Array<double> work(LWORK);
  
  int info = 0;

  dgesvd_(&jobu, &jobvt, &M, &N, aData, &LDA, sData, uData, &LDU, 
    vtData, &LDVT, &(work[0]), &LWORK, &info);

  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
    "dgesvd failed with error code info=" << info);

  
  
}

}
