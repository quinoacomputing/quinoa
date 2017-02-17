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
#include "Intrepid_CellTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_setJacobian ..... MEX interface for the Intrepid (Trilinos) function Intrepid::CellTools::setJacobian.\n\n"
                    "\tintrepid_setJacobian(jacobian,points,cellWorkset,cellType,whichCell=-1)\n\n"
                    "\t<1-in/out> jacobian = Jacobians of reference-to-physical map (4D array of size [spaceDim x spaceDim x #evalPoints x #cells] or 3D array of size [spaceDim x spaceDim x #evalPoints])\n"
                    "\t<2-in>     points = Evaluation points (2D array of size [spaceDim x #evalPoints] or 3D array of size [spaceDim x #evalPoints x #cells])\n"
                    "\t<3-in>     cellWorkset = Cell nodes (3D array of size [spaceDim x #nodes x #cells])\n"
                    "\t<4-in>     cellType = 'Line' | 'Triangle' | 'Quadrilateral' | 'Tetrahedron' | 'Hexahedron' (string)\n"
                    "\t<5-in/opt> whichCell = Cell number (default: -1 means ALL cells) (int)\n\n");

    // Check the number of input arguments
    if((nInput != 4) && (nInput != 5))
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the Jacobian array
    Teuchos::Array<int> jacs_dims;
    m2iGetArrayDims(jacs_dims, pInput[0]);
    // Get the dimensions of the evaluation point array
    Teuchos::Array<int> points_dims;
    m2iGetArrayDims(points_dims, pInput[1]);
    // Get the dimensions of the cell workset
    Teuchos::Array<int> workset_dims;
    m2iGetArrayDims(workset_dims, pInput[2], true);

    // Get the (pointers to) data
    double* jacs_raw = mxGetPr(pInput[0]);
    double* points_raw = mxGetPr(pInput[1]);
    double* workset_raw = mxGetPr(pInput[2]);
    const std::string cell_type(mxArrayToString(pInput[3]));
    int whichCell = -1;
    if(nInput == 5)
    {
        whichCell = (int) mxGetScalar(pInput[4]);
    }

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type, descriptor);

    FieldContainer<double> cellJacs(jacs_dims, jacs_raw);
    FieldContainer<double> evalPoints(points_dims, points_raw);
    FieldContainer<double> cellWorkset(workset_dims, workset_raw);

    try
    {
        CellTools<double>::setJacobian(cellJacs, evalPoints, cellWorkset, *(cellTopo.get()), whichCell);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_setJacobian.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
