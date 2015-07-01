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
// Types for BLAS_Cpp namespace
//

#ifndef BLAS_CPP_TYPES_H
#define BLAS_CPP_TYPES_H

#include "Moocho_ConfigDefs.hpp"

namespace BLAS_Cpp {

/** \defgroup BLAS_Cpp_grp Basic C++/Fortran BLAS declarations/utilities.
 * \ingroup Misc_grp
 *
 * These are declarations that are useful in BLAS (Basic Linear Algebra
 * Subroutines) types of codes.  They provide some foundation for C++/Fortran
 * compatibility.
 */
//@{

/// Size type
typedef size_t size_type;

/** \defgroup BLAS_Cpp_grp_enums BLAS Enumerations */
//@{

/// SIDE
enum Side{
  left
  ,right
};
/// TRANS
enum Transp{
  no_trans     ///< Not transposed
  ,trans       ///< Transposed
  ,conj_trans  ///< Conjugate transpose
};
/// UPLO
enum Uplo {upper, lower};
/// Return the opposite of Uplo argument
inline Uplo operator!(Uplo uplo)
{	return uplo == upper ? lower : upper; }
/// DIAG
enum Diag {unit, nonunit};

//@}

/** \defgroup BLAS_Cpp_grp_helper Helper functions */
//@{

/// Return Transp given a bool
inline Transp bool_to_trans(bool return_trans) {
  return return_trans ? trans : no_trans;
}

/// Returns true if _trans == trans
inline bool trans_to_bool(Transp _trans) {
  return ( _trans == trans );
}

/// Return the opposite of the transpose argument
inline Transp operator!(Transp _trans)	// does not work with conj_trans
{	return _trans == no_trans ? trans : no_trans; }

/// Return the opposite of the transpose argument
inline Transp trans_not(Transp _trans)	// does not work with conj_trans
{	return _trans == no_trans ? trans : no_trans; }

/// Return the transpose of the transpose argument
inline Transp trans_trans(Transp _trans1, Transp _trans2) // does not work with conj_trans
{	return _trans1 == _trans2 ? no_trans : trans; }

/// Give a string name to Transp value
inline const char*  trans_to_string(Transp _trans)
{
  return _trans == no_trans ? "no_trans" : "trans";
}

/// Return rows of a possible transposed matrix
inline size_type rows(size_type rows, size_type cols
  , BLAS_Cpp::Transp _trans)
{
  return _trans == BLAS_Cpp::no_trans ? rows : cols;
}

/// Return columns of a possible transposed matrix
inline size_type cols(size_type rows, size_type cols
  , BLAS_Cpp::Transp _trans)
{
  return _trans == BLAS_Cpp::no_trans ? cols : rows;
}

//@}

//@}

}	// end namespace BLAS_Cpp

#endif // BLAS_CPP_TYPES_H
