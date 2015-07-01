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
// Utilities for pivoting VectorSlices and GenMatrixSlices.

#ifndef PIVOTVECMAT_H
#define PIVOTVECMAT_H

#include <stdexcept>

#include "DenseLinAlgPack_Types.hpp"

/* * @name {\bf DVector / Matrix Permutations}.
  *
  * These are functions for pivoting the elements of a vector and the 
  * rows and/or columns of a rectandular matrix.
  *
  * For DVector and DVectorSlice pivot funcitions the #perm# argument
  * gives the mapping from the new sequence to the old sequence.
  * The statement #i_old = perm(i_new)# gives the index in the
  * old vector.  After the permutation is performed the postcondition:
  \verbatim

    perm_vs(i) == vs(perm(i)), for i = 1,2,...,vs.size()

  \endverbatim
  * is true.
  *
  * For the DMatrix permutation functions the #row_perm# argument
  * gives the row permutations and the #col_perm# argument gives the column
  * permutations respectively.
  *
  * In each of these functions the permuted object can be the same
  * as the unpermuted object.
  */

// @{

// @Include: DenseLinAlgPack_IVector.hpp

// /////////////////////////////////////////////////////////////////////////////////////////
// Public Permutation functions

namespace DenseLinAlgPack {

/** \brief . */
/* * Initialize a permutation to the identiy permutation.
  *
  * Preconditions: <ul>
  * <li> #perm.size() > 0# (throw std::length_error)
  * </ul>
  *
  * Postconditions: <ul>
  * <li> #perm(i) = i#, for i = 1, 2, ... m #perm.size()#
  * </ul>
  */
void identity_perm(IVector* perm);

/** \brief . */
/* * Find the inverse permutation (inv_perm) given a permutation vector (perm).
  *
  * Postconditions: <ul>
  * <li> #inv_perm.size() == perm.size()#
  * </ul>
  */
void inv_perm(const IVector& perm, IVector* inv_perm);

/// Permute a DVectorSlice in place
void perm_ele(const IVector& perm, DVectorSlice* vs);

/// Perform y = P*x
void perm_ele(const DVectorSlice& x, const IVector& perm, DVectorSlice* y);

/// Perform x = P'*y
void inv_perm_ele(const DVectorSlice& y, const IVector& perm, DVectorSlice* x);

/// Permute a GenMatrixSlices rows
void perm_rows(const IVector& row_perm, DMatrixSlice* gms);

/// Permute a GenMatrixSlices columns
void perm_cols(const IVector& col_perm, DMatrixSlice* gms);

/// Permute a GenMatrixSlices rows and columns
void perm_rows_cols(const IVector& row_perm, const IVector& col_perm, DMatrixSlice* gms);

#ifdef TEUCHOS_DEBUG
extern bool PermVecMat_print;
#endif

} // end namespace DenseLinAlgPack

// @}

#endif // PIVOTVECMAT_H
