// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER


#include "mex.h"
#include <Intrepid_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n defGrad ..... MEX interface for function defGrad.\n\n"
                    "\t defGrad(outputFields,inputData)\n\n"
                    "\t<1-in/out> outputData = Output fields array (4D array of size [#spaceDim x #spaceDim x #cubPoints x #cells])\n"
                    "\t<2-in> 	   inputData = Input data array (4D array of size [#spaceDim x #spaceDim x #numFields x #cells])\n");

    // Check the number of input arguments
    if(nInput != 2)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the output field values array
    Teuchos::Array<int> oFields_dims;
    m2iGetArrayDims(oFields_dims, pInput[0]);
    // Get the dimensions of the input vector data values array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);

    FieldContainer<double> outputData(oFields_dims, oFields_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);

    // get sizes
    int invalRank = inputData.rank();
    int outvalRank = outputData.rank();
    int numCells = inputData.dimension(0);
    int numQPs = inputData.dimension(1);
    int numDims = inputData.dimension(2);

    // Compute DefGrad tensor from displacement gradient
    for(int cell = 0; cell < numCells; ++cell)
    {
        for(int qp = 0; qp < numQPs; ++qp)
        {
            for(int i = 0; i < numDims; ++i)
            {
                for(int j = 0; j < numDims; ++j)
                {
                    outputData(cell, qp, i, j) = inputData(cell, qp, i, j);
                }
                outputData(cell, qp, i, i) += 1.0;
            }
        }
    }
}
