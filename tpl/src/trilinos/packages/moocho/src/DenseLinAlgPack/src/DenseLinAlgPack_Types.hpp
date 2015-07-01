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

#ifndef LINALGPACK_TYPES_H
#define LINALGPACK_TYPES_H

#include "DenseLinAlgPack_Options.hpp"
#include "RangePack_Range1D.hpp"
#include "BLAS_Cpp_Types.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Assert.hpp"

namespace MemMngPack {}

namespace DenseLinAlgPack {

using Teuchos::typeName;
using Teuchos::TypeNameTraits;

/* * @name {\bf DenseLinAlgPack Type Declarations}.
  *
  * These are forward declarations of the types used with the DenseLinAlgPack
  * package (namespace).  In addition the BLAS_Cpp enumerations
  * \Ref{Transp}, \Ref{Side}, \Ref{Uplo}, and \Ref{Diag} and there values
  * are avalible using the qualifier #BLAS_Cpp#.
  */

// @{

/** \brief . */
using RangePack::Range1D;
#ifdef _INTEL_CXX
using RangePack::full_range;
#endif
/** \brief . */
using BLAS_Cpp::rows;
/** \brief . */
using BLAS_Cpp::cols;
/** \brief . */
using BLAS_Cpp::trans_not;

/// Enumeration for returning the amount of overlap between two objects
enum EOverLap { NO_OVERLAP = 0, SOME_OVERLAP, SAME_MEM };	

/** \brief . */
class IVector;
/** \brief . */
template<class T>
class VectorTmpl;
/** \brief . */
template<class T>
class VectorSliceTmpl;
/** \brief . */
typedef VectorTmpl<value_type>                DVector;
/** \brief . */
typedef VectorSliceTmpl<value_type>           DVectorSlice;
/** \brief . */
typedef VectorTmpl<extended_value_type>       VectorExt;
/** \brief . */
typedef VectorSliceTmpl<extended_value_type>  VectorSliceExt;
/** \brief . */
class TransVectorSlice;
/** \brief . */
class DMatrix;
/** \brief . */
class DMatrixSlice;
/** \brief . */
class TransGenMatrixSlice;
/** \brief . */
class DMatrixSliceTriEle;
/** \brief . */
class DMatrixSliceTri;
/** \brief . */
class DMatrixSliceSym;

// @}

}  // namespace DenseLinAlgPack

#endif // LINALGPACK_TYPES_H
