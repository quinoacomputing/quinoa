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

#ifndef DELETE_ROW_COL_H
#define DELETE_ROW_COL_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

/** \brief . */
/* * Delete a symmetric row and a column form a triangular matrix.
 *
 * If #M# is a lower triangular matrix then we partition it
 * as:
 \verbatim

   1 |\
     |  \
     |    \
     | M11  \
     |________\ _
  kd |_________|_|
     |         | |\
     |         | |  \
     |         | |    \
     |   M31   | | M33  \
   n |         | |        \
     ----------------------
     1         kd         n

 \endverbatim
 * In order to delete row #kd# and column #kd# the rectangular
 * matrix #M31# is moved up one row and the triangular matrix
 * #M33# is moved up one row and to the left one column.
 *
 * If #M# is an upper triangular matrix then we partition it
 * as:
 \verbatim

  1         kd      n
  -------------------- 1
  \        | |       |
    \  M11 | |  M13  |
      \    | |       |
        \  | |       |
          \|_|_______|
           |_|_______| kd
             \       |
               \ M33 |
                 \   |
                   \ | n
 
 \endverbatim
 *
 * In order to delete row #kd# and column #kd# the matrix
 * #M13# is moved one column to the left and the upper
 * triangular matrix #M33# is moved one row up and
 * on column to the left.
 *
 * Preconditions:<ul>
 * <li> #M != NULL#
 * <li> #M->rows() >= 1#
 * <li> #1 <= kd <= M->rows()#
 * </ul>
 */
void delete_row_col( size_type kd, DMatrixSliceTriEle* M );

} // end namespace DenseLinAlgPack

#endif  // DELETE_ROW_COL_H
