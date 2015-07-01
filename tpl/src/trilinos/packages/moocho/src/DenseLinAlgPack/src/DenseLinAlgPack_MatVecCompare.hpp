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

#ifndef MAT_VEC_COMPARE_H
#define MAT_VEC_COMPARE_H

#include <limits>
#if defined(_GNU_GXX)
#include <cmath>
#else
#include <math.h>
#endif

#include "DenseLinAlgPack_Types.hpp"
#include "TestingHelperPack_update_success.hpp"

namespace DenseLinAlgPack {

using TestingHelperPack::update_success;

/* * @name DVectorSlice and DMatrixSlice comparison functions.
  *
  * These functions compare the elements of two DVectorSlice or DMatrixSlice
  * objects.  If any of the corresponding elements does not obey
  * abs(ele1 - ele2) < sqrt(eps) then the functions return false, otherwise
  * they return true.  An exact test (bits) is not performed to allow for some round-off
  * error to occur and still equate.
  */
// @{

/** \brief . */
const value_type sqrt_eps
#if defined(_GNU_GXX)
  = std::sqrt(std::numeric_limits<value_type>::epsilon());
#elif defined(_CPQ_CXX)
  = ::sqrt(std::numeric_limits<value_type>::epsilon());
#else
  = ::sqrt(std::numeric_limits<value_type>::epsilon());
#endif

/** \brief . */
bool comp(const DVectorSlice& vs1, const DVectorSlice& vs2);

/** \brief . */
bool comp(const DVectorSlice& vs, value_type alpha);

/** \brief . */
bool comp(const DMatrixSlice& gms1, BLAS_Cpp::Transp trans1
  , const DMatrixSlice& gms2, BLAS_Cpp::Transp trans2);

/////
//bool comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2);

inline
/** \brief . */
bool comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2)
{
  return comp(gms1, BLAS_Cpp::no_trans, gms2, BLAS_Cpp::no_trans);
}

/** \brief . */
bool comp(const DMatrixSlice& gms1, value_type alpha);

/** \brief . */
bool comp(const DMatrixSliceTriEle& tri_gms1, const DMatrixSliceTriEle& tri_gms2);

/** \brief . */
bool comp(const DMatrixSliceTriEle& tri_gms1, value_type alpha);

/** \brief . */
bool comp_less(const DVectorSlice& vs, value_type alpha);

// @}

}	// end namespace DenseLinAlgPack

// ////////////////////////////////////
// Inline definitions

//inline
//bool DenseLinAlgPack::comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2)
//{
//	return comp(gms1, BLAS_Cpp::no_trans, gms2, BLAS_Cpp::no_trans);
//}


#endif	// MAT_VEC_COMPARE_H
