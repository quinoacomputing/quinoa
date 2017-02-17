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
            ("\n tangentStiffnessElastostaticsTL ..... MEX interface for function tangentStiffnessElastostaticsTL.\n\n"
                    "\t tangentStiffnessElastostaticsTL(internalForces,stress,inputFields)\n\n"
                    "\t<1-in/out> tangentStiffness = Output data array (3D array of size [#numFields x #numFields x #cells])\n"
                    "\t<2-in> 	   inputData = Input data array (4D array of size [#rows x #numFields x #cubPoints x #cells])\n"
                    "\t<3-in> 	   inputFields = Input fields array (4D array of size [#rows x #numFields x #cubPoints x #cells])\n");

    // Check the number of input arguments
    if(nInput != 3)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the output data values array
    Teuchos::Array<int> oData_dims;
    m2iGetArrayDims(oData_dims, pInput[0]);
    // Get the dimensions of the input data values array
    Teuchos::Array<int> iDataBF_dims;
    m2iGetArrayDims(iDataBF_dims, pInput[1]);
    // Get the dimensions of the input field values array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* oData_raw = mxGetPr(pInput[0]);
    double* iDataBF_raw = mxGetPr(pInput[1]);
    double* iFields_raw = mxGetPr(pInput[2]);

    FieldContainer<double> tangentStiffness(oData_dims, oData_raw);
    FieldContainer<double> inputDataBF(iDataBF_dims, iDataBF_raw);
    FieldContainer<double> wGradBF(iFields_dims, iFields_raw);

    // get sizes
    int numCells = wGradBF.dimension(0);
    int numQPs = wGradBF.dimension(1);
    int numDofs = wGradBF.dimension(2);
    int matDim = wGradBF.dimension(3);

    // compute tangent stiffness matrix
    for(int cell = 0; cell < numCells; ++cell)
    {
        for(int qp = 0; qp < numQPs; ++qp)
        {
            for(int i = 0; i < numDofs; ++i)
            {
                for(int j = 0; j < numDofs; j++)
                {
                    for(int k = 0; k < matDim; k++)
                    {
                        tangentStiffness(cell, i, j) += wGradBF(cell, qp, i, k) * inputDataBF(cell, qp, j, k);
                    }
                }
            }
        }
    }

}
