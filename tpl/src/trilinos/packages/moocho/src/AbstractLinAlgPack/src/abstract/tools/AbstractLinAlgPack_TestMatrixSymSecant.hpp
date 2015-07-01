/*
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
*/

// ///////////////////////////////////////////////////////////
// AbstractLinAlgPack_TestMatrixSymSecant.hpp

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Checks the secant condition <tt>B*s = y</tt>.
 *
 * Call this function after calling <tt>B.secant_update(s,y,...)</tt>.
 *
 * Returns \c true if the secant condition satsifies \c error_tol
 * and \c false otherwise.
 *
 * @param  B  [in] The matrix object we are testing.
 * @param  s  [in] First secant vector
 * @param  y  [in] Second secant vector
 * @param  warning_tol
 *            [in] Any relative error above \c warning_tol will
 *            be noted and possibly printed.
 * @param  error_tol
 *            [in] Any relative error above \c error_tol will
 *            cause the tests to stop and false to
 *            be returned.
 * @param  print_all_warnings
 *            [in] If true then all relative errors greater than
 *            warning_tol will be printed.
 * @param  out
 *            [in/out] Stream that output or error messages are sent
 *            to follow the tests.
 * @param  trase
 *            [in] If \c trase==true then tests will be trased and
 *            if \c trase==false then only error messages will be
 *            output.
 */
bool TestMatrixSymSecant(
  const MatrixOp        &B
  ,const Vector       &s
  ,const Vector       &y
  ,value_type               warning_tol
  ,value_type               error_tol
  ,bool                     print_all_warnings
  ,std::ostream             *out
  ,bool                     trase                 = true
  );

} // end namespace AbstractLinAlgPack
