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

#ifndef INITIALIZE_Q_R_Q_X_H
#define INITIALIZE_Q_R_Q_X_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Initialize <tt>GenPermMatrixSlice</tt> mapping matrices for <tt>Q_R</tt> and <tt>Q_X</tt>.
 *
 *
 * @param  n_R        [in] Number of free variables
 * @param  n_X        [in] Number of fixed variables
 * @param  i_x_free   [in] array (length n_R) of indices of free variables.
 *                    If n_R == 0 then i_x_free can be NULL.  It is allowed
 *                    for i_x_free == NULL in which case it is determined to
 *                    from i_x_fixed[] (if n_X > 0) and i_x_free is assumed
 *                    to be sorted in assending order.
 * @param  i_x_fixed  [in] array (length n_X) of indices of fixed variables.
 *                    If n_X == 0 then i_x_fixed can be NULL.
 * @param  test_setup [in] If true then i_x_free[] and i_x_fixed[] will be
 *                    validated and if not okay then an exception will be
 *                    thown.
 * @param  Q_R_row_i  [out] array (length n_R) of row indices for Q_R.
 *                    If n_R == 0 or i_x_free == NULL and it is known that
 *                    i_x_fixed[l] > n_R, for l = 0...n_X-1, then Q_R_row_i
 *                    can be NULL and will not be accessed.  If Q_R_row_i
 *                    != NULL then it will always be set.
 *                    This array will be sorted in assending order on output
 *                    if it is set.
 * @param  Q_R_col_j  [out] array (length n_R) of column indices for Q_R.
 *                    Q_R_col_j can be NULL when Q_R_row_i is NULL.
 *                    If this array turns out to be sorted then Q_R will
 *                    be set to Q_R->ordered_by() == BY_ROW_AND_COL
 * @param  Q_R        [out] GenPermMatixSlice object initialized with 
 *                    Q_R_row_i and Q_R_col_j.  If n_R == 0 then Q_R
 *                    will be initialized to (n_X by 0).  If it turns out
 *                    that i_x_free == NULL and Q_R has the identity matrix
 *                    as its leading nonzero matrix, then Q_R->is_identity()
 *                    will be true on output.  In any case Q_R->ordered_by()
 *                    will be BY_ROW or BY_ROW_AND_COL on output.
 * @param  Q_X_row_i  [out] array (length n_X) of row indices for Q_X.
 *                    If n_X == 0 then Q_X_row_i can be NULL and will not be accessed.
 * @param  Q_X_col_j  [out] array (length n_X) of column indices for Q_X
 *                    If n_X == 0 then Q_X_col_j can be NULL and will not be accessed.
 *                    If this array turns out to be sorted then Q_X will
 *                    be set to Q_X->ordered_by() == BY_ROW_AND_COL
 * @param  Q_X        [out] GenPermMatixSlice object initialized with 
 *                    Q_X_row_i and Q_X_col_j  If n_X == 0 then Q_X
 *                    will be initialized to (n_X by 0.)  On output Q_X->ordered_by()
 *                    will be BY_ROW or BY_ROW_AND_COL.
 */
void initialize_Q_R_Q_X(
  size_type            n_R
  ,size_type           n_X
  ,const size_type     i_x_free[]
  ,const size_type     i_x_fixed[]
  ,bool                test_setup
  ,size_type           Q_R_row_i[]
  ,size_type           Q_R_col_j[]
  ,GenPermMatrixSlice  *Q_R
  ,size_type           Q_X_row_i[]
  ,size_type           Q_X_col_j[]
  ,GenPermMatrixSlice  *Q_X
  );

}  // end namespace ConstrainedOptPack

#endif // INITIALIZE_Q_R_Q_X_H
