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
// Options for DenseLinAlgPack compilation
//

#ifndef LINALGPACK_OPTIONS_H
#define LINALGPACK_OPTIONS_H

#include "DenseLinAlgPack_extended_value_type.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_F77_wrappers.h"

#if !defined(LINALGPACK_NO_CHECKS)

/* * @name {\bf DenseLinAlgPack Options}.
  *
  * The header file DenseLinAlgPack_Options.hpp contains the defines for several macros that
  * determine how the library is built.  The user should comment out any
  * macros that her or she does not want to be defined.  The definition of
  * these macros cause the library code to assert the preconditions documented
  * for each of the member and non-member functions and throw the listed exceptions
  * if they are not satisfied.  Precondtions are supposed to be the 
  * responcibility of the client code so the user may only want to define
  * these macros during debugging for better program verification.
  * If the user checks all of the preconditions listed in this documentation for the calls
  * to all functions then the checks performed by the library are redundant.
  */
// @{

/** \brief . */
/* * If defined the library code checks to see if subscripts are in bounds for element access
  * an subregion indexing.  If the preconditions for the subscripting operations are
  * not satisfied then the listed exceptions will be thrown.
  */
#ifndef LINALGPACK_CHECK_RANGE
#define LINALGPACK_CHECK_RANGE 1
#endif

/** \brief . */
/* * If defined the library code checks to see if the sizes of rhs arguments in expressions are compatible.
  * The exception std::length_error will be thrown if rhs sizes are not compatible.
  */
#ifndef LINALGPACK_CHECK_RHS_SIZES
#define LINALGPACK_CHECK_RHS_SIZES 1
#endif

/** \brief . */
/* * If defined the library code checks to see if DVectorSlice and DMatrixSlice objects have valid constructions.
  * If they do not have valid constructions then an exception will be thrown.  The operation of these
  * checks may depend on the definition of the macro \Ref{LINALGPACK_CHECK_RANGE}.
  */
#ifndef LINALGPACK_CHECK_SLICE_SETUP
#define LINALGPACK_CHECK_SLICE_SETUP 1
#endif

#endif

namespace DenseLinAlgPack{

/// Typedef for the value type of elements that is used for the library.
typedef FortranTypes::f_dbl_prec		value_type;
/// Typedef for the index type of elements that are used by the library
typedef Teuchos::Ordinal index_type;
/// Typedef for the size type of elements that are used by the library
typedef	Teuchos::Ordinal size_type;

}

// @}

#endif // LINALGPACK_OPTIONS_H
