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


/*
 * SundanceTransformationHN.cpp
 *
 *  Created on: Mar 16, 2010
 *      Author: benk
 */

#include "SundanceTransformationHN.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"

using std::ios_base;
using std::setw;
using std::endl;

using namespace Sundance;
using namespace Teuchos;

extern "C"
{
  int dgemm_(const char* transA, const char* transB,
	     const int* M, const int *N, const int* K,
	     const double* alpha,
	     const double* A, const int* ldA,
	     const double* B, const int* ldB,
	     const double* beta,
	     double* C, const int* ldC);
}


TransformationHN::TransformationHN(const HNDoFMapBase* dofMap ,
				   const int nrCol , const int nrRaw):
  dofMap_(dofMap) , nrRow_((nrRaw==0)?1:nrRaw) , nrCol_((nrCol==0)?1:nrCol) {
}

TransformationHN::~TransformationHN() {
}

void TransformationHN::preApply( const int funcID,
		         int cellDim ,
				 const CellJacobianBatch& JTrans,
				 const CellJacobianBatch& JVol,
				 const Array<int>& facetIndex,
				 const RCP<Array<int> >& cellLIDs,
				 RCP<Array<double> >& A
				 ) const 
{

  // if we do not operate on the maxCellType then do nothing
  //if (dofMap_->getSpacialMeshDim() != cellDim ) return;

  int nrRow = (nrRow_ == 0)? 1 : nrRow_;
  int nrCol = (nrCol_ == 0)? 1 : nrCol_;
  int                 cellLID , matrixSize , i;
  bool                doTransform = false;
  Array<double>       M;
  const Array<int>*   cellLIDAray = cellLIDs.get();
  double*             matrixAray = &((*A)[0]);
  Array<double>       tmpArray(nrRow*nrCol);

  for (i = 0 ; i < cellLIDs->size() ; i++){
    // get the Cell ID
    cellLID = (*cellLIDAray)[i];
    //get the Matrix transformation
    SUNDANCE_MSG2( verb() , "TransformationHN::preApply() cellLID:" << cellLID );
    if (cellDim == dofMap_->getSpacialMeshDim())
      dofMap_->getTrafoMatrixForCell( cellLID , funcID , matrixSize, doTransform, M );
    else
      dofMap_->getTrafoMatrixForFacet( cellDim , cellLID , facetIndex[i] , funcID , matrixSize, doTransform, M );

    if (doTransform)
    {
      for (int ii = 0 ; ii < nrRow*nrCol ; ii++) tmpArray[ii] = matrixAray[ i*nrRow*nrCol + ii];

      //SUNDANCE_MSG2( verb() , " TransformationHN::preApply() cellLID: " << cellLID << "  doTransform:" <<
      //		 doTransform << "  nrRow:" << nrRow << " nrCol:" << nrCol);
      //SUNDANCE_MSG2( verb() ," TransformationHN::preApply() Matrix:" << M);
 	  //SUNDANCE_MSG2( verb() , " TransformationHN::postApply() nrRow:" << nrRow << " nrCol:" << nrCol);
      // Make transformation with M^T from left
      multiplyFromLeftWithTransp( M , &(matrixAray[ i*nrRow*nrCol ]) , &(tmpArray[0]));
    }
  }
}

/** */
void TransformationHN::postApply( const int funcID,
		          int cellDim ,
				  const CellJacobianBatch& JTrans,
				  const CellJacobianBatch& JVol,
				  const Array<int>& facetIndex,
				  const RCP<Array<int> >& cellLIDs,
				  RCP<Array<double> >& A
				  ) const 
{

  // if we do not operate on the maxCellType then do nothing
  //if (dofMap_->getSpacialMeshDim() != cellDim ) return;

  int nrRow = (nrRow_ == 0)? 1 : nrRow_;
  int nrCol = (nrCol_ == 0)? 1 : nrCol_;
  int                 cellLID , matrixSize;
  bool                doTransform = false;
  Array<double>       M;
  const Array<int>*   cellLIDAray = cellLIDs.get();
  double*             matrixAray = &((*A)[0]);
  Array<double>       tmpArray(nrRow*nrCol);

  for (int i = 0 ; i < cellLIDs->size() ; i++){
    // get the Cell ID
    cellLID = (*cellLIDAray)[i];
    //get the Matrix transformation
    //SUNDANCE_MSG2( verb() ,"TransformationHN::postApply() cellLID:" << cellLID);
    if (cellDim == dofMap_->getSpacialMeshDim())
      dofMap_->getTrafoMatrixForCell( cellLID , funcID , matrixSize, doTransform, M );
    else
      dofMap_->getTrafoMatrixForFacet( cellDim , cellLID , facetIndex[i] , funcID , matrixSize, doTransform, M );

    if (doTransform)
    {
      for (int ii = 0 ; ii < nrRow*nrCol ; ii++) tmpArray[ii] = matrixAray[ i*nrRow*nrCol + ii];

      //SUNDANCE_MSG2( verb() , " TransformationHN::postApply() cellLID: " << cellLID << "  doTransform:" <<
      //		 doTransform << "  nrRow:" << nrRow << " nrCol:" << nrCol);
      //SUNDANCE_MSG2( verb() , " TransformationHN::postApply() Matrix:" << M );
	  // SUNDANCE_MSG2( verb() , " TransformationHN::postApply() nrRow:" << nrRow << " nrCol:" << nrCol);

      // Make transformation with M from right
      multiplyFromRight( &(matrixAray[ i*nrRow*nrCol ]) , M , &(tmpArray[0]));
    }
  }
}

void TransformationHN::preapplyTranspose( const int cellDim,
					  const int funcID,
					  const Array<int>& cellLIDs,
					  const Array<int>& facetIndex,
					  Array<double>& A
					  ) const 
{

  // if we do not operate on the maxCellType then do nothing
  //if (dofMap_->getSpacialMeshDim() != cellDim ) return;

  // Todo: in this function we can be sure that we we apply this to a vector
  int                 cellLID , matrixSize;
  bool                doTransform;
  Array<double>       M;
  double*             matrixAray = &(A[0]);
  // we can calculate the row size here
  int nrRow = A.size() / cellLIDs.size();
  int nrCol = 1;
  Array<double>       tmpArray(nrRow*nrCol);


  for (int i = 0 ; i < cellLIDs.size() ; i++){
    // get the Cell ID
    cellLID = cellLIDs[i];
    //get the Matrix transformation
    if (cellDim == dofMap_->getSpacialMeshDim())
      dofMap_->getTrafoMatrixForCell( cellLID , funcID , matrixSize, doTransform, M );
    else
      dofMap_->getTrafoMatrixForFacet( cellDim , cellLID , facetIndex[i] , funcID , matrixSize, doTransform, M );

    if (doTransform)
    {
      SUNDANCE_MSG2( verb() ," TransformationHN::postApply() A.size(): " << A.size() << "  cellLIDs.size():" <<
		     cellLIDs.size());

      for (int ii = 0 ; ii < nrRow*nrCol ; ii++) tmpArray[ii] = matrixAray[ i*nrRow*nrCol + ii];

      SUNDANCE_MSG2( verb() , " TransformationHN::postApply() cellLID: " << cellLID << "  doTransform:" <<
		     doTransform << "  nrRow:" << nrRow << " nrCol:" << nrCol );
      SUNDANCE_MSG2( verb() , " TransformationHN::postApply() Matrix:" << M );

      SUNDANCE_MSG2( verb() , " TransformationHN BEFORE tmpArray:" << tmpArray );
      // Make transformation with M from right
      multiplyFromLeft( M , &(matrixAray[ i*nrRow*nrCol] ) , &(tmpArray[0]) , nrRow , nrCol );

      for (int ii = 0 ; ii < nrRow*nrCol ; ii++) tmpArray[ii] = matrixAray[ i*nrRow*nrCol + ii];

      SUNDANCE_MSG2( verb() , " TransformationHN AFTER tmpArray:" << tmpArray );
    }
  }

}

void TransformationHN::multiplyFromLeftWithTransp(Array<double>& M , double* A_end ,const double* A_copy) const {
  const int nrRow = (nrRow_ == 0)? 1 : nrRow_;
  const int nrCol = (nrCol_ == 0)? 1 : nrCol_;

  const double* A = &(M[0]);
  const double* B =  A_copy;
  const double beta = 0.0;
  const double alpha = 1.0;
/*
  int ii,jj,kk;
  SUNDANCE_MSG2( 3 , " nrRow = " << nrRow << " nrCol = " << nrCol);
  SUNDANCE_MSG2( 3 , " A = ");
  for (ii = 0 ; ii < nrRow ; ii++){
      for(jj = 0 ; jj < nrRow ; jj++){ std::cout << setw(12) <<  A[ii*nrRow + jj] << " ,"; }
	  std::cout << std::endl;
  }
  SUNDANCE_MSG2( 3 , " B = ");
  for (ii = 0 ; ii < nrRow ; ii++){
      for(jj = 0 ; jj < nrCol ; jj++){ std::cout << setw(12) << B[ii*nrCol + jj] << " ,"; }
	  std::cout << std::endl;
  }
*/
  // C := alpha*op( A )*op( B ) + beta*C,
  //with op( A )  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  //dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    ::dgemm_("N", "T", &nrCol, &nrRow , &nrRow,  &alpha, B ,
    &nrCol, A , &nrRow, &beta,
    A_end, &nrCol);

    // ----------------- plotting ---------------
/*
    SUNDANCE_MSG2( 3 , " nrRow = " << nrRow << " nrCol = " << nrCol);
    SUNDANCE_MSG2( 3 , " A_end = ");
    for (ii = 0 ; ii < nrRow ; ii++){
        for(jj = 0 ; jj < nrCol ; jj++){ std::cout << setw(12) << A_end[ii*nrCol + jj] << " ,"; }
  	  std::cout << std::endl;
    }*/
/*double sum = 0;
  // M is nrRow X nrRow   and A is nrCol X nrRow
  for (ii = 0 ; ii < nrRow ; ii++){
    for(jj = 0 ; jj < nrCol ; jj++){
      sum = 0.0;
      for(kk = 0 ; kk < nrRow ; kk++){
	sum += M[kk*nrRow + ii] * A_copy[kk*nrCol + jj];
      }
      A_end[ii*nrCol + jj] = sum;
    }
  }*/
    /*  // ----------------- plotting ---------------
  SUNDANCE_MSG2( 3 , " nrRow = " << nrRow << " nrCol = " << nrCol);
  SUNDANCE_MSG2( 3 , " A_end = ");
  for (ii = 0 ; ii < nrRow ; ii++){
      for(jj = 0 ; jj < nrCol ; jj++){ std::cout << setw(12) << A_end[ii*nrCol + jj] << " ,"; }
	  std::cout << std::endl;
  }*/

}

void TransformationHN::multiplyFromRight(double* A_end , Array<double>& M , const double* A_copy) const{
  int nrRow = (nrRow_ == 0)? 1 : nrRow_;
  int nrCol = (nrCol_ == 0)? 1 : nrCol_;

  const double* A = A_copy;
  const double* B = &(M[0]);
  const double beta = 0.0;
  const double alpha = 1.0;

  // C := alpha*op( A )*op( B ) + beta*C,
  //with op( A )  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  //dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    ::dgemm_("N", "N", &nrCol, &nrRow , &nrCol,  &alpha, B ,
    &nrCol, A , &nrCol, &beta,
    A_end, &nrCol);

/* double sum = 0;
  // M is nrCol X nrCol   and A is nrCol X nrRow
  for (ii = 0 ; ii < nrRow ; ii++){
    for(jj = 0 ; jj < nrCol ; jj++){
      sum = 0.0;
      for(kk = 0 ; kk < nrCol ; kk++){
	sum += A_copy[ ii*nrCol + kk ] * M[ kk*nrCol + jj];
      }
      A_end[ii*nrCol + jj] = sum;
    }
  } */
}

void TransformationHN::multiplyFromLeft(Array<double>& M , double* A_end , const double* A_copy
					, const int nrRow , const int nrCol ) const{
  const double* A = &(M[0]);
  const double* B =  A_copy;
  const double beta = 0.0;
  const double alpha = 1.0;

  // we can assume that this is applied only to vectors,
  // C := alpha*op( A )*op( B ) + beta*C,
  //with op( A )  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  //dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    ::dgemm_("N", "N", &nrCol, &nrRow , &nrRow,  &alpha, B ,
    &nrCol, A , &nrRow, &beta,
    A_end, &nrCol);

/* int ii,jj,kk;
  double sum = 0; // todo: A_copy is always a column
  // M is nrRow X nrRow   and A is nrCol X nrRow
  for (ii = 0 ; ii < nrRow ; ii++){
    for(jj = 0 ; jj < nrCol ; jj++){
      sum = 0.0;
      for(kk = 0 ; kk < nrRow ; kk++){
	sum += M[ii*nrRow + kk] * A_copy[kk*nrCol + jj];
      }
      A_end[ii*nrCol + jj] = sum;
    }
  }*/
}
