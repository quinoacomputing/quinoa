/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/


#include <iostream>
#include <iomanip>


// Intrepid includes
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_Basis.hpp>
#include <Intrepid_CubatureDirectLineGauss.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_CubatureTensor.hpp>


// Teuchos includes
#include "Teuchos_RCP.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

#include <KokkosArray_DeviceHost.hpp>
#include <KokkosArray_DeviceHost_MDArrayView.hpp>
#include <KokkosArray_DeviceHost_MultiVectorView.hpp>
#include <KokkosArray_DeviceHost_ValueView.hpp>

#include <KokkosArray_DeviceCuda.hpp>
#include <KokkosArray_DeviceCuda_MDArrayView.hpp>
#include <KokkosArray_DeviceCuda_MultiVectorView.hpp>
#include <KokkosArray_DeviceCuda_ValueView.hpp>
#include <KokkosArray_DeviceCuda_ParallelFor.hpp>
#include <KokkosArray_DeviceCuda_ParallelReduce.hpp>

#include <KokkosArray_DeviceCuda_macros.hpp>
#include <Jacobian.hpp>
#include <Transform.hpp>
#include <TransformValue.hpp>
#include <simpleFill.hpp>
#include <Multiply.hpp>
#include <Integrate.hpp>
#include <computeCellMeasure.hpp>
#include <Invert.hpp>
#include <Determinant.hpp>
#include <Poisson_Driver.hpp>
#include <KokkosArray_DeviceClear_macros.hpp>

namespace Test {


void poisson_cuda(int beg , int end)
{
	KokkosArray::DeviceCuda::initialize();
	Test::poisson_run< KokkosArray::DeviceCuda>(beg , end);
};

} // namespace Test


