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

#ifndef GEN_PERM_MATRIX_SLICE_OP_H
#define GEN_PERM_MATRIX_SLICE_OP_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"

namespace AbstractLinAlgPack {

/** @name Operations for GenPermMatrixSlice.
 *
 * ToDo: Finish documentation!
 */
//@{

/** \brief <tt>sv_lhs = alpha * op(P_rhs1) * vs_rhs2</tt>.
 * 
 * This function will resize the sparse vector lhs and only the
 * resultant nonzero elements will be added.
 * 
 * If <tt>op(P_rhs1) is sorted by row (i.e. <tt>op(P_rhs1) = P_rhs1</tt> sorted by row
 * or <tt>op(P_rhs1) = P_rhs1'</tt> sorted by column) then <tt>sv_lhs->assume_sorted(true)</tt>
 * is called and <tt>sv_lhs->is_sorted()==true</tt> on output.
 * 
 * This function will execute in <tt>O(P_rhs1.nz())</tt> time.
 */ 
void V_StMtV(
  SpVector* sv_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
  ,BLAS_Cpp::Transp P_rhs1_trans, const DVectorSlice& vs_rhs2
  );

inline
/// <tt>sv_lhs = op(P_rhs1) * vs_rhs2</tt>
void V_MtV(
  SpVector* sv_lhs, const GenPermMatrixSlice& P_rhs1
  ,BLAS_Cpp::Transp P_rhs1_trans, const DVectorSlice& vs_rhs2
  )
{
  V_StMtV(sv_lhs,1.0,P_rhs1,P_rhs1_trans,vs_rhs2);
}

/** \brief <tt>sv_lhs = alpha * op(P_rhs1) * sv_rhs2</tt>.
 * 
 * This function will resize the sparse vector lhs and add only the
 * nonzero elements in the rhs.
 * 
 * If <tt>op(P_rhs1)</tt> is sorted by row (i.e. <tt>op(P_rhs1) = P_rhs1</tt> sorted by row
 * or <tt>op(P_rhs1) = P_rhs1'</tt> sorted by column) then <tt>sv_lhs->assume_sorted(true)</tt>
 * is called.
 * 
 * Let's assume that <tt>op(P_rhs1)</tt> is not sorted by column and <tt>sv_rhs2</tt> is
 * sorted.  In this case a linear search will have to be performed to match up
 * elements.  <tt>P_rhs1</tt> will be iterated through sequentially and
 * the corresponding nonzero element in <tt>sv_rhs2</tt> searched for (binary search).
 * Therefore, the runtime in this case will be:
 \verbatim

    O( P_rhs1.nz() * log(sv_rhs2.nz()) )
 \endverbatim
 * If <tt>P_rhs1</tt> and <tt>sv_rhs2</tt> are unsorted, then the runtime will be:
 \verbatim
 
    O( P_rhs1.nz() * sv_rhs2.nz() )
 \endverbatim
 * If <tt>op(P_rhs1)</tt> is sorted by column and <tt>sv_rhs2</tt> is also 
 * sorted then the runtime will be:
 \verbatim

    O( max( P_rhs1.nz(), sv_rhs2.nz() ) )
 \endverbatim 
 * Of course if <tt>op(P_rhs1)</tt> is not sorted by row then the output vector will
 * not be assumed sorted.
 */ 
void V_StMtV(
  SpVector* sv_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
  ,BLAS_Cpp::Transp P_rhs1_trans, const SpVectorSlice& sv_rhs2
  );

inline
/// sv_lhs = op(P_rhs1) * sv_rhs2
void V_MtV(
  SpVector* sv_lhs, const GenPermMatrixSlice& P_rhs1
  ,BLAS_Cpp::Transp P_rhs1_trans, const SpVectorSlice& sv_rhs2
  )
{
  V_StMtV(sv_lhs,1.0,P_rhs1,P_rhs1_trans,sv_rhs2);
}


/** \brief <tt>sv_lhs += alpha * op(P_rhs1) * vs_rhs2</tt>.
 * 
 * This function will not resize the sparse vector lhs and will add
 * new elements for the nonzero elements in the rhs.  Therefore it is
 * up to the client to ensure that there is sufficient storage for
 * these elements.  This function will not check to see if elements
 * with duplicate indexes are added. It is up to the client to determine
 * that.  If <tt>sv_lhs</tt> is sorted on input and <tt>op(P_rhs1)</tt>
 * is sorted by row, then <tt>sv_lhs->is_sorted() == true</tt> on output.
 * 
 * This function will execute in <tt>O(P_rhs1.nz())</tt> time.
 */ 
void Vp_StMtV(
  SpVector* sv_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
  ,BLAS_Cpp::Transp P_rhs1_trans, const DVectorSlice& vs_rhs2
  );

inline
/** \brief <tt>sv_lhs += op(P_rhs1) * vs_rhs2</tt>.
 */ 
void Vp_MtV(
  SpVector* sv_lhs, const GenPermMatrixSlice& P_rhs1
  ,BLAS_Cpp::Transp P_rhs1_trans, const DVectorSlice& vs_rhs2
  )
{
  Vp_StMtV(sv_lhs,1.0,P_rhs1,P_rhs1_trans,vs_rhs2);
}

/// <tt>vs_lhs = alpha * op(P_rhs1) * vs_rhs2 + beta * vs_lhs</tt>
void Vp_StMtV(
  DVectorSlice* vs_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
  ,BLAS_Cpp::Transp P_rhs1_trans, const DVectorSlice& vs_rhs2, value_type beta = 1.0
  );

/// <tt>vs_lhs = alpha * op(P_rhs1) * sv_rhs2 + beta * vs_lhs</tt>
void Vp_StMtV(
  DVectorSlice* vs_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
  ,BLAS_Cpp::Transp P_rhs1_trans, const SpVectorSlice& sv_rhs2, value_type beta = 1.0)
  ;

/** \brief Find the intersection between two GenPermMatrixSlice objects.
 *
 * This subroutine has two modes.  In the first mode (<tt>Q_max_nz == 0</tt>) it just
 * computes the number of nonzero entries in the matrix:
 *
 * <tt>Q = op(P1)*op(P1)</tt>
 *
 * In the second mode (<tt>Q_max_nz > 0 && Q_row_i != NULL && Q_col_j != NULL</tt>) it also
 * computes the row and column arrays for the resultant matrix <tt>Q</tt>.  In addition
 * if <tt>Q != NULL</tt> then a GenPermMatrixSlice object will be setup with as:
 *
 * <tt>Q->initialize( op(P1).rows(), op(P2).cols(), Q_nz, 0, 0, Q_ordered_by, Q_row_i, Q_col_j, false )</tt>
 *
 * Above <tt>Q_ordered_by</tt> will be determined of the fly.
 *
 * This operation might require  O(min(op(P1).cols(),op(P2).rows()) temporary storage
 * but will always be executed in:
 *
 * O(P1.nz()) + O(P2.nz()) + O(min(op(P1).cols(),op(P2).rows())
 *
 * @param  P1         [in] Right hand side permutation matrix.
 * @param  P1_trans   [in] If no_trans then op(P1) = P1 otherwise op(P1) = P1'.
 * @param  P1         [in] Left hand side permutation matrix.
 * @param  P1_trans   [in] If no_trans then op(P2) = P2 otherwise op(P2) = P2'.
 * @param  Q_nz       [out] On return will contain the number of nonzeros in the
 *                    resultant matrix <tt>Q</tt>
 * @param  Q_max_nz   [in] If <tt>Q_max_nz > 0</tt> then the resultant row <tt>Q_row_i</tt> and column <tt>Q_col_j</tt>
 *                    indices will be set.    If it turns out that <tt>Q_nz</tt> will be larger than
 *                    <tt>Q_max_nz</tt> then the exception <tt>std::length_error</tt> will be thrown and <tt>Q_row_i</tt>
 *                    and <tt>Q_col_j</tt> may be left in an inconsistent state. If <tt>Q_max_nz == 0</tt> then the
 *                    rest of the return arguments are ignored and the resultant matrix will not be returned.
 * @param  Q_row_i    [out] Array (length <tt>Q_max_nz</tt>): If <tt>Q_max_nz > 0</tt> then out return <tt>Q_row_i[k], k=0,,,Q_nz-1</tt>
 *                    will contain the row indices for the resultant matrix <tt>Q</tt>.  If <tt>Q_max_nz == 0</tt> then
 *                    <tt>Q_row_i</tt> can be <tt>NULL</tt>.
 * @param  Q_row_i    [out] Array (length <tt>Q_max_nz</tt>): If <tt>Q_max_nz > 0</tt> then out return <tt>Q_col_j[k], k=0,,,Q_nz-1</tt>
 *                    will contain the column indices for the resultant matrix <tt>Q</tt>.  If <tt>Q_max_nz == 0</tt> then
 *                    <tt>Qd_col_j</tt> can be <tt>NULL</tt>.
 * @param  Q          [out] If <tt>Q_max_nz > 0 && Q != NULL</tt> then <tt>Q</tt> will be initialized as described above.
 *                    It is allowed for <tt>Q == NULL</tt>.
 */
void intersection(
  const GenPermMatrixSlice     &P1
  ,BLAS_Cpp::Transp            P1_trans
  ,const GenPermMatrixSlice    &P2
  ,BLAS_Cpp::Transp            P2_trans
  ,size_type                   *Q_nz
  ,const size_type             Q_max_nz     = 0
  ,size_type                   Q_row_i[]    = NULL
  ,size_type                   Q_col_j[]    = NULL
  ,GenPermMatrixSlice          *Q           = NULL
  );

//@}

} // end namespace AbstractLinAlgPack

#endif   // GEN_PERM_MATRIX_SLICE_OP_H
