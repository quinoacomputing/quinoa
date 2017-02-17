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
            ("\nintrepid_HGRADtransformGRAD ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::HGRADtransformGRAD.\n\n"
                    "\tintrepid_HGRADtransformGRAD(physGrads,jacInvs,refGrads,transpose='T')\n\n"
                    "\t<1-in/out> physGrads = Gradients in physical space (4D array of size [spaceDim x #refPoints x #basisFields #cells])\n"
                    "\t<2-in>     jacInvs = Cell jacobian inverses (3D array of size [spaceDim x spacedim #cells])\n"
                    "\t<3-in>     refGrads = Gradients on the reference cell (3D array of size [spaceDim x #refPoints x #basisFields])\n"
                    "\t<4-in/opt> transpose = Cell number: 'T' | 't' (default -> transpose) | 'N' | 'n' (do not transpose) (string)\n\n");

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

    // Get the dimensions of the physical gradients array
    Teuchos::Array<int> pGrad_dims;
    m2iGetArrayDims(pGrad_dims, pInput[0]);
    // Get the dimensions of the inverse Jacobian array
    Teuchos::Array<int> invJac_dims;
    m2iGetArrayDims(invJac_dims, pInput[1]);
    // Get the dimensions of the reference gradients array
    Teuchos::Array<int> rGrad_dims;
    m2iGetArrayDims(rGrad_dims, pInput[2]);

    // Get the (pointers to) data
    double* pGrad_raw = mxGetPr(pInput[0]);
    double* invJac_raw = mxGetPr(pInput[1]);
    double* rGrad_raw = mxGetPr(pInput[2]);

    char transpose = 'Z';
    if(nInput == 3)
    {
        transpose = 'T';
    }
    else if(nInput == 4)
    {
        const std::string t_str(mxArrayToString(pInput[3]));
        if((t_str == "N") || (t_str == "n"))
        {
            transpose = 'N';
        }
        else if((t_str == "T") || (t_str == "t"))
        {
            transpose = 'T';
        }
    }

    FieldContainer<double> physGrads(pGrad_dims, pGrad_raw);
    FieldContainer<double> invJacs(invJac_dims, invJac_raw);
    FieldContainer<double> refGrads(rGrad_dims, rGrad_raw);

    try
    {
        FunctionSpaceTools::HGRADtransformGRAD<double>(physGrads, invJacs, refGrads, transpose);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_HGRADtransformGRAD.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
