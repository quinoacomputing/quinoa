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
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_tensorMultiplyDataField ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::tensorMultiplyDataField.\n\n"
                    "\tintrepid_tensorMultiplyDataField(outputFields,inputData,inputFields,transpose)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (3D array of size [spaceDim x #cubPoints x #cells])\n"
                    "\t<2-in>     inputData = Input data array (4D array of size [spaceDim x spaceDim x #cubPoints x #cells] or 3D array of size [spaceDim x #cubPoints x #cells] or 2D array of size [#cubPoints x #cells])\n"
                    "\t<3-in>     inputFields = Input fields array (5D array of size [spaceDim x spaceDim x #cubPoints x #fields x #cells] or 4D array of size [spaceDim x #cubPoints x #fields x #cells])\n"
                    "\t<4-in>     transpose = If 'T', use transposed tensor; if 'N', no transpose. Default: 'N'\n\n");

    // Check the number of input arguments
    if((nInput != 3) && (nInput != 4))
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the physical values array
    Teuchos::Array<int> oFields_dims;
    m2iGetArrayDims(oFields_dims, pInput[0]);
    // Get the dimensions of the input data array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);
    // Get the dimensions of the input fields array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);
    double* iFields_raw = mxGetPr(pInput[2]);

    char transpose = 'N';
    if(nInput == 3)
    {
        transpose = 'N';
    }
    else if(nInput == 4)
    {
        const std::string si_str(mxArrayToString(pInput[4]));
        if(si_str == "T")
        {
            transpose = 'T';
        }
    }

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);
    FieldContainer<double> inputFields(iFields_dims, iFields_raw);

    try
    {
        FunctionSpaceTools::tensorMultiplyDataField<double>(outputFields, inputData, inputFields, transpose);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_tensorMultiplyDataField.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
