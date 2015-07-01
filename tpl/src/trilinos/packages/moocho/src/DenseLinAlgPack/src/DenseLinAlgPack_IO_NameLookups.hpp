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
//

#ifndef LINALGPACK_IO_NAME_LOOKUPS_H
#define LINALGPACK_IO_NAME_LOOKUPS_H

#ifdef _WINDOWS

//  MS VC++ 5.0 is not performing function look for these templated operator
//  functions as it should so I will do it for the compiler.  These
//	inline functions are injected into the local namespace.  These should
//  be removed once a standards conforming compiler is available.

namespace {

// Define some boiler plate macros

#define OPEATOR_FUNCTION(OPERATOR,STREAM_TYPE,FORMAT_TYPE,OBJECT_TYPE)								\
  inline STREAM_TYPE & OPERATOR ( STREAM_TYPE & s													\
    , DenseLinAlgPack::LinAlgPackIO:: ## FORMAT_TYPE ## <DenseLinAlgPack:: ## OBJECT_TYPE ## >& bf)	\
  {																								\
    return DenseLinAlgPack:: ## OPERATOR ## (s,bf);													\
  }

#define INPUT_OPEATOR_FUNCTION(FORMAT_TYPE,OBJECT_TYPE)												\
  OPEATOR_FUNCTION( operator>> , std::istream , FORMAT_TYPE , OBJECT_TYPE )						\

#define OUTPUT_OPEATOR_FUNCTION(FORMAT_TYPE,OBJECT_TYPE)											\
  OPEATOR_FUNCTION( operator<< , std::ostream , FORMAT_TYPE , OBJECT_TYPE )						\


INPUT_OPEATOR_FUNCTION(		bound_format		,	DVector			)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	DVector			)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	DVector			)

INPUT_OPEATOR_FUNCTION(		bound_format		,	DVectorSlice		)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	DVectorSlice		)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	DVectorSlice		)

INPUT_OPEATOR_FUNCTION(		bound_format		,	DMatrix		)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	DMatrix		)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	DMatrix		)

INPUT_OPEATOR_FUNCTION(		bound_format		,	DMatrixSlice	)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	DMatrixSlice	)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	DMatrixSlice	)

#undef OPEATOR_FUNCTION
#undef INPUT_OPEATOR_FUNCTION
#undef OUTPUT_OPEATOR_FUNCTION

}	// end namespace

#endif

#endif // LINALGPACK_IO_NAME_LOOKUPS_H
