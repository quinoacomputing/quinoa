#if 0

// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "ConstrainedOptPack_MatrixHessianSuperBasicInitDiagonal.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "Midynamic_cast_verbose.h"

namespace ConstrainedOptPack {

MatrixHessianSuperBasicInitDiagonal::MatrixHessianSuperBasicInitDiagonal()
  : B_RR_init_(NULL)
{}

void MatrixHessianSuperBasicInitDiagonal::initialize(
  size_type            n
  ,size_type           n_R
  ,const size_type     i_x_free[]
  ,const size_type     i_x_fixed[]
  ,const EBounds       bnd_fixed[]
  ,const B_RR_ptr_t&   B_RR_ptr
  ,const B_RX_ptr_t&   B_RX_ptr
  ,BLAS_Cpp::Transp    B_RX_trans
  ,const B_XX_ptr_t&   B_XX_ptr
  )
{
  using Teuchos::dyn_cast;

  // Validate the B_RR supports this interface
#ifdef _WINDOWS
  B_RR_init_ = &dynamic_cast<MatrixSymInitDiag&>(
    const_cast<MatrixSymWithOpFactorized&>(*B_RR_ptr)
    );
#else
  B_RR_init_ = &dyn_cast<MatrixSymInitDiag>(
    const_cast<MatrixSymWithOpFactorized&>(*B_RR_ptr)
    );
#endif

  MatrixHessianSuperBasic::initialize(
    n,n_R,i_x_free,i_x_fixed,bnd_fixed
    ,B_RR_ptr,B_RX_ptr,B_RX_trans,B_XX_ptr 
    );
}

// Overridden from MatrixSymInitDiag

void MatrixHessianSuperBasicInitDiagonal::init_identity(
  size_type n, value_type alpha )
{
  assert_initialized();
  B_RR_init_->init_identity(n,alpha);
  MatrixHessianSuperBasic::initialize(
    n,n,NULL,NULL,NULL
    ,this->B_RR_ptr()
    ,this->B_RX_ptr(),this->B_RX_trans()
    ,this->B_XX_ptr()
    );
}

void MatrixHessianSuperBasicInitDiagonal::init_diagonal(
  const DVectorSlice& diag )
{
  assert_initialized();
  B_RR_init_->init_diagonal(diag);
  MatrixHessianSuperBasic::initialize(
    diag.size(),diag.size(),NULL,NULL,NULL
    ,this->B_RR_ptr()
    ,this->B_RX_ptr(),this->B_RX_trans()
    ,this->B_XX_ptr()
    );
}

} // end namespace ConstrainedOptPack

#endif // 0
