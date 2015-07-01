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

#ifndef ASSERT_PRINT_NAN_INF_H
#define ASSERT_PRINT_NAN_INF_H

#include <stdexcept>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief . */
class NaNInfException : public std::runtime_error
{public: NaNInfException(const std::string& what_arg) : std::runtime_error(what_arg) {}};

/** \brief This function asserts if a value_type scalare is a NaN or Inf and optionally
 * prints out these entires.
 *
 * @param  val             [in] Value the check
 * @param  name            [in] Name of the scale variable for output purposes
 * @param  throw_excpt     [in] If true and is found to be a NaN or Inf
 *                         then a NaNInfException excetion is thrown after
 *                         any output.
 * @param  out             [in/out] If out==NULL then not output is produced.
 *                         If out!=NULL and val is not
 *                         NaN or Inf, then no output is produced.
 *                         If out!=NULL and val is
 *                         NaN or Inf then this will be printed before any
 *                         execption is thrown.
 *
 * @return Returns true if val is not NaN or Inf.  If val
 * is NaN or Inf then false will be returned unless a
 * excetion NaNInfException was thrown (throw_except==true).
 */
bool assert_print_nan_inf( const value_type& val, const char name[]
  , bool throw_excpt, std::ostream* out );

/** \brief This function asserts if a vector has any NaN or inf entries and optionally
 * prints out these entires.
 *
 * @param  v              [in]	Vector slice to check
 * @param  name           [in]	Name of the vector for output purposes
 * @param  throw_excpt    [in]	If true and an entry is found to be a NaN or Inf
 *                        then a NaNInfException excetion is thrown after
 *                        any output.
 * @param  out            [in/out]	If out==NULL then not output is produced.
 *                        If out!=NULL and none of the entries is
 *                        NaN or Inf, then no output is produced.
 *                        If out!=NULL then any entries that are
 *                        NaN or Inf will be printed before any
 *                        execption is thrown.
 *
 * @return Returns true none of the entries are NaN or Inf.  If one of the
 * entries is NaN or Inf then false will be returned unless an
 * excetion was thrown (throw_except==true).
 */
bool assert_print_nan_inf( const Vector& v, const char name[]
  , bool throw_excpt, std::ostream* out );

}	// end namespace AbstractLinAlgPack

#endif // ASSERT_PRINT_NAN_INF_H
