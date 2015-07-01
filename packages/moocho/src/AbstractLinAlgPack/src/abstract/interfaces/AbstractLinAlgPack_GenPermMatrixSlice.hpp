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

#ifndef GEN_PERM_MATRIX_SLICE_H
#define GEN_PERM_MATRIX_SLICE_H

#include "AbstractLinAlgPack_GenPermMatrixSliceIterator.hpp"

namespace AbstractLinAlgPack {

/** \brief Concrete matrix type to represent general permutation (mapping) matrices.
 *
 * These are matrices who's rows or columns represent eta vectors
 * (i.e. only one nonzero element with the value 1).  These matrices
 * can be rectangular and have one or more zero rows & columns.  Therefore, these
 * matrices can be used to represent gathering and scattering operations
 * on other vectors and matrices.
 *
 * This is only a view type.  The client specifies the mapping arrays and then
 * this class provides a clean encapsulation for the mapping.  Objects of this
 * type can also represent the identity matrix which is constructed with
 * the initialize_identity(...) function.
 *
 * The default copy constructor is allowd but the default 
 * assignment operator function is not.
 */
class GenPermMatrixSlice {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  enum EIdentityOrZero { IDENTITY_MATRIX, ZERO_MATRIX };

  /** \brief . */
  typedef GenPermMatrixSliceIteratorPack::EOrderedBy		EOrderedBy;

  /** \brief . */
  typedef GenPermMatrixSliceIteratorPack::row_col_iterator<const index_type>
                      const_iterator;
  /** \brief . */
  typedef	ptrdiff_t						difference_type;

  //@}

  /// Construct to an uninitialzied, unsized matrix.
  GenPermMatrixSlice();

  /// Construct to a matrix intialized to identity or zero (see initialize(,,,)).
  GenPermMatrixSlice( index_type rows, index_type cols, EIdentityOrZero type );

  /** \brief Initialize an identity or zero permutation.
   *
   * If <tt>type == IDENTITY_MATRIX</tt> then after this function is called <tt>this</tt> will
   * represent <tt>Q = [ I; 0 ]</tt> if <tt>rows > cols</tt> or <tt>Q = [ I, 0 ]</tt>
   * if <tt>rows < cols</tt> or <tt>Q = I</tt> if <tt>rows == cols</tt>.  If <tt>type == ZERO_MATRIX]</tt>
   * then <tt>this</tt> will represent a <tt>rows</tt> x <tt>cols</tt> zero matrix.
   *
   * Postconditions:<ul>
   * <li> <tt>this->rows() == rows</tt>
   * <li> <tt>this->cols() == cols</tt>
   * <li> [<tt>type == IDENTITY_MATRIX</tt>] <tt>this->nz() == min(rows,cols)</tt>
   * <li> [<tt>type == ZERO_MATRIX</tt>]     <tt>this->nz() == 0</tt>
   * <li> [<tt>type == IDENTITY_MATRIX</tt>] <tt>this->is_identity() == true</tt>
   * <li> [<tt>type == ZERO_MATRIX</tt>]     <tt>this->is_identity() == false</tt>
   * <li> <tt>this->ordered_by() == BY_ROW_AND_COL</tt>
   * </ul>
   */
  void initialize( index_type rows, index_type cols, EIdentityOrZero type );

  /** \brief Initialize.
   *
   *	@param	rows	[in] Number of rows in matrix
   *	@param	cols	[in] Number of columns in matrix
   *	@param	nz		[in] Number of nonzero elements in the matrix
   *	@param	row_off	[in] Row offsets for row_i[]
   *	@param	col_off	[in] Column offsets for col_i[]
   *	@param	ordered_by
   *					[in] The ordering of the nonzero elements
   *	@param	row_i	[in] Array (size nz):  If nz == 0 then
   *					row_i can be NULL
   *	@param	col_j	[in] Array (size nz): If nz == 0 then
   *					col_j can be NULL
   *	@param	test_setup
   *					[in] If true then all of the preconditions for
   *					the input arguments will be checked.
   *
   * This function sets up general permutation view.
   * The client must setup the arrays row_i[] and col_j[] to define
   * the mapping.  If nz == 0 then row_i and col_j can be NULL.
   * Otherwise, row_i != NULL and col_j != NULL.  If nz > 0 and
   * ordered_by == BY_ROW then row_i[k] must be sorted in assending order
   * else if ordered_by == BY_COL or col_j[k] must be sorted in assnding order
   * else if ordered_by == UNORDERED then no ordering for row_i[k] or
   * col_j[k] is required.  If nz == 1 then the value of ordered_by
   * is not really significant at it is ordered by row and by column.
   * 
   * It is required that if nz > 0 then:
   * 1 <= row_i[k] + row_off <= rows, for k = 1...nz
   * 1 <= col_j[k] + col_off <= cols, for k = 1...nz
   * 
   * All of these preconditions will be checked if test_setup == true.
   *
   * After setup, the memory pointed to by row_i[] and col_j[] must not
   * be altered since this object does not make an independent copy
   * of this data.
   * 
   * After construction the nonzero elements of this matrix are:
   * M(row_i[k]+row_off,col_j[k]+col_off) = 1.0, for k = 1...nz.
   *
   */
  void initialize(
    index_type			rows
    ,index_type			cols
    ,index_type			nz
    ,difference_type	row_off
    ,difference_type	col_off
    ,EOrderedBy			ordered_by
    ,const index_type	row_i[]
    ,const index_type	col_j[]
    ,bool				test_setup = false
    );

  /** \brief Initialize and sort.
    *
    * This is the same as the initialize(...) function except that
    * this function will actually sort the entries by row or
    * by column or not at all..
    *
    * ToDo: Finish documentation.
    *
    *	@param	rows	[in] Number of rows in matrix
    *	@param	cols	[in] Number of columns in matrix
    *	@param	nz		[in] Number of nonzero elements in the matrix
    *	@param	row_off	[in] Row offsets for row_i[]
    *	@param	col_off	[in] Column offsets for col_i[]
    *	@param	ordered_by
    *					[in] The ordering of the nonzero elements
    *	@param	row_i	[in/out] Array (size nz):  If nz == 0 then
    *					row_i can be NULL.  On output it will be
    *					sorted according to ordered_by
    *	@param	col_j	[in/out] Array (size nz): If nz == 0 then
    *					col_j can be NULL.  On output it will be
    *					sorted according to ordered_by
    *	@param	test_setup
    *					[in] If true then all of the preconditions for
    *					the input arguments will be checked.
    */
  void initialize_and_sort(
    index_type			rows
    ,index_type			cols
    ,index_type			nz
    ,difference_type	row_off
    ,difference_type	col_off
    ,EOrderedBy			ordered_by
    ,index_type			row_i[]
    ,index_type			col_j[]
    ,bool				test_setup = false
    );
    
  /** \brief Bind the view of another GenPermMatrixSlice object.
    *
    * After construction these objects will share points to
    * the same row_i[] and col_j[] arrays.
    */
  void bind( const GenPermMatrixSlice& gpms );

  /** \brief . */
  index_type rows() const;
  /** \brief . */
  index_type cols() const;
  /** \brief . */
  index_type nz() const;
  /** \brief . */
  EOrderedBy ordered_by() const;
  /** \brief . */
  bool is_identity() const;

  /** \brief Lookup the ith row index for the nonzero entry in the jth column
   * if it exists.
   *
   * This function will return 0 if the index is not found.  If 
   * <tt>this->ordered_by() == BY_COL || this->ordered_by() == BY_ROW_AND_COL</tt>
   * then this function will be executed in O(log(<tt>this->nz()</tt>)) time.
   * Otherwise it will execute in O(<tt>this->nz()</tt>) time.
   *
   * Preconditions:<ul>
   * <li> <tt>(1 <= col_j && col_j <= this->cols())</tt> (throw <tt>std::out_of_range</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>(1 <= return && return <= this->rows()) || return == 0</tt>
   * </ul>
   */
  index_type lookup_row_i(index_type col_j) const;

  /** \brief Lookup the jth column index for the nonzero entry in the ith row
   * if it exists.
   *
   * This function will return 0 if the index is not found.  If 
   * <tt>this->ordered_by() == BY_ROW || this->ordered_by() == BY_ROW_AND_COL</tt>
   * then this function will be executed in O(log(<tt>this->nz()</tt>)) time.
   * Otherwise it will execute in O(<tt>this->nz()</tt>) time.
   *
   * Preconditions:<ul>
   * <li> <tt>(1 <= row_i && row_i <= this->rows())</tt> (throw <tt>std::out_of_range</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>(1 <= return && return <= this->cols()) || return == 0</tt>
   * </ul>
   */
  index_type lookup_col_j(index_type row_i) const;

  /** @name Iterator Access. */
  //@{

  /** \brief Return a random access iterator for accessing which row and column that each
   * nonzero 1.0 entry is located at.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_identity() == false</tt> (throw <tt>???</tt>)
   * </ul>
   *
   * If <tt>this->is_identity() == true</tt> then these iterators are obvoisly unneccesary
   * and will throw exceptions.
   \verbatim
   for( GenPermMatrixSlice::const_iterator itr = gpms.begin(); itr != gpms.end(); ++itr )
   {
       std::cout << "row_i = " << itr->row_i();
       std::cout << "col_j = " << itr->col_j();
   }
   \endverbatim
   *
   * You can also take a difference between iterators.
   */
  const_iterator begin() const;

  /// Return the end of <tt>this->const_iterator_begin()</tt>.
  const_iterator end() const;

  //@}

  /** \brief Create a submatrix by row, by column.
    *
    * If nz > 1 and this->ordered_by() == BY_ROW then ordered_by
    * must also equal BY_ROW or if this->ordered_by() == BY_COL
    * then ordered_by must also equal BY_COL
    * or an std::logic_error exception will be thrown.  If nz == 1, then
    * obviously the nozeros are ordered by row and by column.
    * This function should not be called if this->is_identity() == true.
    * 
    * The argument rng must be explicitly sized (rng.full_range() != true) or the
    * exception std::logic_error will be thrown.  Also, the range argument
    * must obey rng.lbound() >= 1, and rng.ubound() <= this->rows()
    * if ordered_by == BY_ROW and rng.ubound() <= this->cols()
    * if ordered_by == BY_COL.  The argument ordered_by == UNORDERED
    * is not allowed and this function can not be called
    * if this->ordered_by() == UNORDERED.  This operation just does not
    * make any sense.
    * 
    * The returned submatrix will contain all the entries for the designated
    * rows if ordered_by == BY_ROW or columns if ordered_by == BY_CO.
    * 
    * ToDo: Spell out the behavior of this operation more carefully.
    */
  const GenPermMatrixSlice create_submatrix( const Range1D& rng
    , EOrderedBy ordered_by ) const;

private:

  // //////////////////////////////
  // Private data members

  index_type			rows_;
  index_type			cols_;
  index_type			nz_;
  difference_type		row_off_;
  difference_type		col_off_;
  EOrderedBy			ordered_by_;
  const index_type		*row_i_;
  const index_type		*col_j_;

  // //////////////////////////////
  // Private static data members

  // ToDo: We could allocate a class-wide array, initialize
  // it to [1,2,3 ...] and then use it for the iterators
  // when is_idenity() == true!  This would make implementing
  // a lot of code a lot easier if we don't care about a little
  // inefficiency!  We could just allocate a large chunk
  // of memory by default (or client could do this for us)
  // and then construct it when needed.  If a client ever
  // requested an iterator when not enough storage was avalible
  // then we would throw an exception.

  // //////////////////////////////
  // Private member functions

  // Validate the input data (not the ordering!)
  static void validate_input_data(
    index_type			rows
    ,index_type			cols
    ,index_type			nz
    ,difference_type	row_off
    ,difference_type	col_off
    ,EOrderedBy			ordered_by
    ,const index_type	row_i[]
    ,const index_type	col_j[]
    ,std::ostringstream &omsg
    );

  /** \brief . */
  void validate_not_identity() const;
  
  /// not defined and not to be called
  GenPermMatrixSlice& operator=( const GenPermMatrixSlice& );

};	// end class GenPermMatrixSlice

// //////////////////////////////////////////////////////////
// Inline members for GenPermMatrixSlice

inline
GenPermMatrixSlice::GenPermMatrixSlice( index_type rows, index_type cols, EIdentityOrZero type )
{
  initialize(rows,cols,type);
}

inline
index_type GenPermMatrixSlice::rows() const
{
  return rows_;
}

inline
index_type GenPermMatrixSlice::cols() const
{
  return cols_;
}

inline
index_type GenPermMatrixSlice::nz() const
{
  return nz_;
}

inline
bool GenPermMatrixSlice::is_identity() const
{
  return nz_ > 0 && row_i_ == NULL;
}

inline
GenPermMatrixSlice::EOrderedBy GenPermMatrixSlice::ordered_by() const
{
  return ordered_by_;
}

}	// end namespace AbstractLinAlgPack

#endif   // GEN_PERM_MATRIX_SLICE_H
