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
 * SundanceMatrixStore.hpp
 *
 *  Created on: Jun 21, 2010
 *      Author: benk
 */

#ifndef SUNDANCEMATRIXSTORE_HPP_
#define SUNDANCEMATRIXSTORE_HPP_

#include "SundanceDefs.hpp"
#include "PlayaExceptions.hpp"

namespace Sundance {

using namespace Teuchos;

/** Class to store the transformation matrixes globally. <br>
 * Previously we used them one per cell (with hanging nodes) */
class MatrixStore {
public:

	/** Empty ctor*/
	MatrixStore();

	/** empty dtor */
	virtual ~MatrixStore() {;}

	/** initializes the */
	void init(int chunkNr);

	/** returns the matrix index corresponding to the chunckIndex
	 * @param chunkIndex
	 * @param M , the linearized matrix */
	int addMatrix(int chunkIndex, Array<double>& M);

	/** returns the matrix */
	void getMatrix(int chunkIndex, int matrixIndex, Array<double>& transfMatrix) const;


private:

	/** nr chunck*/
	int nrChunk_;

	/** the length of the matrix for each chunck*/
	Array<int> matrixLength_;

	/** [chunckIndex][matrixIndex] -> trafo Matrix*/
	Array< Array< Array<double> > > matrixStore_;

};

}

#endif /* SUNDANCEMATRIXSTORE_HPP_ */
