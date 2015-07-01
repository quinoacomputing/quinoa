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

#ifndef ALAP_MATRIX_LOAD_SPARSE_ELEMENTS_H
#define ALAP_MATRIX_LOAD_SPARSE_ELEMENTS_H

#include "AbstractLinAlgPack_MatrixBase.hpp"

namespace AbstractLinAlgPack {

/** \brief Mix-in interface for loading nonzero elements into a sparse matrix data structure.
 *
 * The formats supported are:
 *
 * Coordiante:
 \verbatim

     Aval[k], Arow[k], Acol[k], k = 0..num_nonzeros(...)-1
 \endverbatim
 * Compressed Row (Column):
 \verbatim

    Aval[k], Acol[k], k = 0..num_nonzeros(...)-1
    Arow_start[j], j = 0..rows()-1
 \endverbatim
 * This is meant to be the do-all interface for clients to use to add nonzero elements
 * for sparse matrices.
 *
 * ToDo: Discuss element uniqueness!
 *
 * ToDo: Finish documentation!
 */
class MatrixLoadSparseElements
  : public virtual AbstractLinAlgPack::MatrixBase
{
public:

  /** \brief . */
  enum EAssumeElementUniqueness {
    ELEMENTS_ASSUME_UNIQUE             ///< Entries assumed have unique row and column indexes (client must enforce this!)
    ,ELEMENTS_ASSUME_DUPLICATES_SUM    ///< Entries allowed with duplicate row and column indexes with the understanding that the values are summed
  };

  /** \brief Resize the matrix and reserve space for nonzero elements to be added.
   *
   * All of the nonzeros in the current matrix are discarded and we start fresh.
   *
   * ToDo: Finish documentation!
   */
  virtual void reinitialize(
    size_type                  rows
    ,size_type                 cols
    ,size_type                 max_nz
    ,EAssumeElementUniqueness  element_uniqueness = ELEMENTS_ASSUME_DUPLICATES_SUM
    ) = 0;

  /** \brief Reinitialize internal counter to load new nonzero values.
   *
   * The row and column index arrays are preserved from the last setup and here
   * the client only wants to set the nonzero values for the same matrix
   * structure.
   *
   * ToDo: Finish documentation!
   */
  virtual void reset_to_load_values() = 0;

  /** \brief Get pointers to buffers to add nonzero elements.
   *
   * @param  max_nz_load
   *                  [in] Maximum number of nonzero elements that will be set
   *                  in the returned buffers.  If \c reset_to_load_values() 
   * @param  val      [out] On output <tt>*val</tt> is set to a pointer to an contiguous array
   *                  of memory of at least \c max_nz_load entries for which the values of the
   *                  nonzero elements to add are to be set.
   * @param  row_i    [out] On output <tt>*row_i</tt> is set to a pointer to an contiguous array
   *                  of memory of at least \c max_nz_load entries for which the row indexes of the
   *                  nonzero elements to add are to be set.  If <tt>row_i == NULL</tt> then
   *                  no buffer is allocated for the row indexes.
   * @param  col_j    [out] On output <tt>*col_J</tt> is set to a pointer to an contiguous array
   *                  of memory of at least \c max_nz_load entries for which the column indexes of the
   *                  nonzero elements to add are to be set.  If <tt>col_j == NULL</tt> then
   *                  no buffer is allocated for the column indexes.
   *
   * Preconditions:<ul>
   * <li> If <tt>reset_to_load_values()</tt> was called to setup nonzero elements, then
   *      \c row_i and \c col_j must be \c NULL or a <tt>std::logic_error</tt> exception
   *      is thrown.
   * </ul>
   *
   * After entries in the arrays \c (*val)[], \c (*row_i)[] and \c (*col_j)[] are set, the client must
   * call <tt>this->commit_load_nonzeros_buffers()</tt> to commit the nonzero entries that are set
   * in these buffers.
   */
  virtual void get_load_nonzeros_buffers(
    size_type      max_nz_load
    ,value_type    **val
    ,index_type    **row_i
    ,index_type    **col_j
    ) = 0;

  /** \brief Commit nonzeros in buffers obtained from \c get_load_nonzeros_buffers().
   *
   * @param  nz_commit
   *                  [in] Number of nonzero elements to be loaded fron the buffers.  Note that
   *                  <tt>nz_commit <= max_nz_load</tt> on the previous call to <tt>get_load_nonzeros_buffers()</tt>.
   * @param  val      [in/out] On input <tt>(*val)[]</tt> contains an array of \c nz_commit value entries for
   *                  nonzero elements to load.  This must point to the same buffer returned from the last call to
   *                  \c get_load_nonzero_buffers(). On output <tt>*val</tt> is set to \c NULL>
   * @param  row_i    [in/out] On input <tt>(*row_i)[]</tt> contains an array of \c nz_commit row index entries for
   *                  nonzero elements to load.  This must point to the same buffer returned from the last call to
   *                  \c get_load_nonzero_buffers(). If <tt>row_i == NULL</tt> then no row indexes are set.  Here
   *                  it is assumed that the row indexes from a previous load have already been set.
   *                  On output <tt>*row_i</tt> is set to \c NULL>
   * @param  col_j    [in/out] On input <tt>(*col_j)[]</tt> contains an array of \c nz_commit column index entries for
   *                  nonzero elements to load.  This must point to the same buffer returned from the last call to
   *                  \c get_load_nonzero_buffers().  If <tt>col_J == NULL</tt> then no column indexes are set.  Here
   *                  it is assumed that the column indexes from a previous load have already been set.
   *                  On output <tt>*col_j</tt> is set to \c NULL>
   *
   * Preconditions:<ul>
   * <li> If <tt>reset_to_load_values()</tt> was called to setup nonzero elements, then
   *      \c row_i and \c col_j must be \c NULL or a <tt>std::logic_error</tt> exception
   *      is thrown.
   * </ul>
   *
   * ToDo: Finish documentation!
   */
  virtual void commit_load_nonzeros_buffers(
    size_type      nz_commit
    ,value_type    **val
    ,index_type    **row_i
    ,index_type    **col_j
    ) = 0;

  /** \brief To be called when the matrix construction is finally finished after all
   * of the nonzero entries have been added.
   *
   * If \c reset_to_load_values() was called to initialize this set of loads then
   * the number of nonzeros added must be exactly the same as the original load
   * or an <tt>std::logic_error</tt> will be thrown with an appropriate error
   * message.
   *
   * @param  test_setup  [in] If true, then the setup will be checked (ToDo: elaborate)
   *
   * Postconditions:<ul>
   * <li> <tt>this->nz()</tt> returns the sum of all of <tt>nz_commit</tt> in all previous calls to
   *      <tt>commit_load_nonzeros_buffers()</tt> since the last call to <tt>reinitialize()</tt>.
   * </ul>
   */
  virtual void finish_construction( bool test_setup ) = 0;

};	// end class MatrixLoadSparseElements

}	// end namespace AbstractLinAlgPack 

#endif	// ALAP_MATRIX_LOAD_SPARSE_ELEMENTS_H
