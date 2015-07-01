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
 * SundanceMatrixStore.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: benk
 */

#include "SundanceMatrixStore.hpp"

using namespace Sundance;
using namespace Teuchos;

MatrixStore::MatrixStore() {
	nrChunk_ = 0;
	matrixLength_.resize(0);
	matrixStore_.resize(0);
}


void MatrixStore::init(int chunkNr){
	// initialize data structure
	nrChunk_ = chunkNr;
	matrixLength_.resize(chunkNr);
	matrixStore_.resize(chunkNr);
}

int MatrixStore::addMatrix(int chunkIndex, Array<double>& M){
   double L2_Vect_norm2 , tmp;
   int foundIndex = -1;
   TEUCHOS_TEST_FOR_EXCEPTION( (chunkIndex >= nrChunk_) , std::runtime_error, "MatrixStore::addMatrix" );
   int nrMatrix = matrixStore_[chunkIndex].size();

   for (int ii = 0 ; ii < nrMatrix ; ii++){
	   L2_Vect_norm2 = 0.0;
	   TEUCHOS_TEST_FOR_EXCEPTION( M.size() >  matrixStore_[chunkIndex][ii].size() , std::runtime_error, "MatrixStore::addMatrix" );
	   // calculate the L2 difference between the 2 matrixes
	   for (int jj = 0 ; jj < M.size() ; jj++){
		   tmp = (M[jj] - matrixStore_[chunkIndex][ii][jj]);
		   L2_Vect_norm2 += tmp*tmp;
	   }
	   // test if the two matrixes are equal, if yes then we can return the matrix Index
       if (L2_Vect_norm2 < 1e-10){
    	   foundIndex = ii;
    	   return foundIndex;
       }
   }
   // add the matrix (because this could not be found) and return
   matrixStore_[chunkIndex].append(M);
   return nrMatrix;
}

void MatrixStore::getMatrix(int chunkIndex, int matrixIndex, Array<double>& transfMatrix) const{
	// return the corresponding matrix
	TEUCHOS_TEST_FOR_EXCEPTION( (chunkIndex >= nrChunk_) || (matrixStore_[chunkIndex].size() <= matrixIndex)
			, std::runtime_error, "MatrixStore::getMatrix , invalid: " << nrChunk_ << " , " << matrixIndex);
	transfMatrix = matrixStore_[chunkIndex][matrixIndex];
}
