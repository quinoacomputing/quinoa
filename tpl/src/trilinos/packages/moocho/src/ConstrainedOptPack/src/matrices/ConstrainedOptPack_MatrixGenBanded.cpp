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

#include <assert.h>
#include <sstream>

#include "ConstrainedOptPack_MatrixGenBanded.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_BLAS_Cpp.hpp"
#include "MiWorkspacePack.h"

namespace ConstrainedOptPack {

MatrixGenBanded::MatrixGenBanded(
  size_type                         m
  ,size_type                        n
  ,size_type                        kl
  ,size_type                        ku
  ,DMatrixSlice                   *MB
  ,const release_resource_ptr_t&    MB_release_resource_ptr
  )
{
  initialize(m,n,kl,ku,MB,MB_release_resource_ptr);
}

void MatrixGenBanded::initialize(
  size_type                         m
  ,size_type                        n
  ,size_type                        kl
  ,size_type                        ku
  ,DMatrixSlice                   *MB
  ,const release_resource_ptr_t&    MB_release_resource_ptr
  )
{
  // Validate input

  if( m == 0 ) {
    if( n != 0 )
      throw std::invalid_argument(
        "MatrixGenBanded::initialize(...): Error, "
        "n must be 0 if m == 0" );
    if( kl != 0 )
      throw std::invalid_argument(
        "MatrixGenBanded::initialize(...): Error, "
        "kl must be 0 if m == 0" );
    if( ku != 0 )
      throw std::invalid_argument(
        "MatrixGenBanded::initialize(...): Error, "
        "ku must be 0 if m == 0" );
    if( MB != NULL )
      throw std::invalid_argument(
        "MatrixGenBanded::initialize(...): Error, "
        "MB must be NULL if m == 0" );
    if( MB_release_resource_ptr.get() != NULL )
      throw std::invalid_argument(
        "MatrixGenBanded::initialize(...): Error, "
        "MB_release_resource_ptr.get() must be NULL if m == 0" );
  }
  else {
    if( kl + 1 > m )
      throw std::invalid_argument(
        "MatrixGenBanded::initialize(...): Error, "
        "kl + 1 can not be larger than m" );
    if( ku + 1 > n )
      throw std::invalid_argument(
        "MatrixGenBanded::initialize(...): Error, "
        "ku + 1 can not be larger than n" );
    if( MB == NULL )
      throw std::invalid_argument(
        "MatrixGenBanded::initialize(...): Error, "
        "MB must not be NULL if n > 0" );
  }

  // Set the members

  if( m == 0 ) {
    m_                        = 0;
    n_                        = 0;
    kl_                       = 0;
    ku_                       = 0;
    MB_.bind(DMatrixSlice());
    MB_release_resource_ptr_  = NULL;
  }
  else {
    // Set the members
    m_                        = m;
    n_                        = n;
    kl_                       = kl;
    ku_                       = ku;
    MB_.bind(*MB);
  }
}

// Overridden from MatrixOp

size_type MatrixGenBanded::rows() const
{
  return m_;
}

size_type MatrixGenBanded::cols() const
{
  return n_;
}

size_type MatrixGenBanded::nz() const
{
  return (ku_ + kl_ + 1) * n_ - ( (ku_+1) * (ku_+1) - (ku_+1) )/2  - ( (kl_+1) * (kl_+1) - (kl_+1) )/2; // Is correct?
}

std::ostream& MatrixGenBanded::output(std::ostream& out) const
{
  return MatrixOp::output(out); // ToDo: Implement specialized version later!
}

void MatrixGenBanded::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
  , const DVectorSlice& x, value_type b) const
{
  assert_initialized();
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->size(), n_, n_, BLAS_Cpp::no_trans, x.size() );
  BLAS_Cpp::gbmv(M_trans,m_,n_,kl_,ku_,a,MB_.col_ptr(1),MB_.max_rows(),x.raw_ptr(),x.stride()
           ,b,y->raw_ptr(),y->stride());
}

void MatrixGenBanded::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
  , const SpVectorSlice& x, value_type b) const
{
  assert_initialized();
  MatrixOp::Vp_StMtV(y,a,M_trans,x,b); // ToDo: Implement spacialized operation when needed!
}

void MatrixGenBanded::Vp_StPtMtV(
  DVectorSlice* y, value_type a
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp M_trans
  , const DVectorSlice& x, value_type b) const
{
  assert_initialized();
  MatrixOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,b); // ToDo: Implement spacialized operation when needed!
}

void MatrixGenBanded::Vp_StPtMtV(
  DVectorSlice* y, value_type a
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp M_trans
  , const SpVectorSlice& x, value_type b) const
{
  assert_initialized();
  MatrixOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,b); // ToDo: Implement spacialized operation when needed!
}

// Private member functions

void MatrixGenBanded::assert_initialized() const
{
  if( m_ == 0 )
    throw std::logic_error("MatrixGenBanded::assert_initialized(): Error, "
                 "not initialized!" );
}

} // end namespace ConstrainedOptPack

#endif // 0
