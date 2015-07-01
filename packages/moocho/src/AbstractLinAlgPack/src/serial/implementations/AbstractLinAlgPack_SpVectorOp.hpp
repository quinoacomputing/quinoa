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

#ifndef SP_VECTOR_OP_H
#define SP_VECTOR_OP_H

#include "AbstractLinAlgPack_SparseVectorSliceOp.hpp"
#include "AbstractLinAlgPack_SparseElement.hpp"

namespace AbstractLinAlgPack {

/** \brief Add elements from a dense vector to a sparse vector.
 *
 * Here sv_lhs is not resized and only elements are added.
 * The purpose of this function is to add elements from
 * a dense vector to a sparse vector.
 *
 * Postconditions:\begin{itemize}
 * \item sv_lhs->nz() == vs_rhs->size() + sv_lhs_before->size()
 * \item [sv_lhs_before->is_sorted() || sv_lhs_before->nz() == 0] sv_lhs->is_sorted() == true
 * \end{itemize}
 */
void add_elements( SpVector* sv_lhs, value_type alpha, const DVectorSlice& vs_rhs
           , size_type offset = 0, bool add_zeros = true );

/** \brief Add elements from a sparse vector to another sparse vector.
 *
 * Here sv_lhs is not resized and only elements are added.
 * The purpose of this function is to add elements from
 * a sparse vector to a sparse vector.
 *
 * Postconditions:\begin{itemize}
 * \item sv_lhs->nz() == sv_rhs->nz() + sv_lhs_before->size()
 * \item [(sv_lhs_before->is_sorted() || sv_lhs_before->nz() == 0)
 *        && (sv_rhs.is_sorted() || sv_rhs.nz() == 0)] sv_lhs->is_sorted() == true
 * \end{itemize}
 */
void add_elements( SpVector* sv_lhs, value_type alpha, const SpVectorSlice& sv_rhs
           , size_type offset = 0, bool add_zeros = true );

inline
/** \brief Create a dense representation of a sparse vector.
 *
 * The primary use if the function is to create a DVectorSlice
 * object that only represents the nonzero values of the
 * sparse vector.  This could have several different uses
 * but one of the most significant examples is when you want
 * to discard the indices when sv_rhs->size() == sv_rhs->nz() and
 * sv_rhs->is_sorted() == true.
 */
DVectorSlice dense_view( SpVectorSlice& sv_rhs );

inline
/** \brief . */
const DVectorSlice dense_view( const SpVectorSlice& sv_rhs );

} // end namespace AbstractLinAlgPack

// /////////////////////////////////////////
// Inline function definitions

inline
DenseLinAlgPack::DVectorSlice
AbstractLinAlgPack::dense_view( SpVectorSlice& sv_rhs )
{
  return sv_rhs.nz()
    ? DVectorSlice( &sv_rhs.begin()->value(), sv_rhs.nz(), 2 )
    : DVectorSlice( NULL, 0, 0 );
}

inline
const DenseLinAlgPack::DVectorSlice
AbstractLinAlgPack::dense_view( const SpVectorSlice& sv_rhs )
{
  return sv_rhs.nz()
    ? DVectorSlice( &const_cast<SpVectorSlice&>(sv_rhs).begin()->value(), sv_rhs.nz(), 2 )
    : DVectorSlice( NULL, 0, 0 );
}

#endif // SP_VECTOR_CLASS_H
