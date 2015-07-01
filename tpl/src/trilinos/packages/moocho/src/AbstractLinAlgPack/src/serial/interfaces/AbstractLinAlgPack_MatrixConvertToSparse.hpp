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

#ifndef MATRIX_CONVERT_TO_SPARSE_H
#define MATRIX_CONVERT_TO_SPARSE_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixBase.hpp"

namespace AbstractLinAlgPack {

/** \brief Mix-in interface for extracing explicit elements from a sparse matrix
 * in one of several Fortran compatible formats.
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
 * This is meant to be the simplest interface for clients to use
 * to extract nonzero elements from matrices.
 *
 * The following matrix is used in the documentation for the following extraction
 * functions.
 *
 \verbatim

    [	1.1				1.3		1.4	]	1
    [			2.2				2.4	]	2
  A =	[			3.2					]	3
    [	4.1				4.3		4.4	]	4
    [	5.1				5.3		5.4	]	5
        1		  2		  3		  4
 \endverbatim
 *
 * Note that above, A has:<ul>
 * <li> <tt>A_nz == this->num_nonzeros(EXTRACT_FULL_MATRIX,ELEMENTS_FORCE_UNIQUE) == 12</tt>
 * <li> <tt>A_up_nz = this->num_nonzeros(EXTRACT_UPPER_TRIANGULAR,ELEMENTS_FORCE_UNIQUE) == 6</tt>
 * <li> <tt>A_lo_nz = this->num_nonzeros(EXTRACT_LOWER_TRIANGULAR,ELEMENTS_FORCE_UNIQUE) == 9</tt>
 * </ul>
 */
class MatrixConvertToSparse
  : virtual public AbstractLinAlgPack::MatrixBase // doxygen needs full path
{
public:

  /** \brief . */
  enum EExtractRegion {
    EXTRACT_FULL_MATRIX       ///< Extract the nonzero elements from the full matrix
    ,EXTRACT_UPPER_TRIANGULAR ///< Extract the nonzero elements from the upper triangular region including the diagonal
    ,EXTRACT_LOWER_TRIANGULAR ///< Extract the nonzero elements from the lower triangular region including the diagonal
  };
  /** \brief . */
  enum EElementUniqueness {
    ELEMENTS_FORCE_UNIQUE             ///< Entries must have unique row and column indexes
    ,ELEMENTS_ALLOW_DUPLICATES_SUM    ///< Entries with duplicate row and column indexes are allowed with the understanding that the values are summed
  };

  /** \brief Returns the number of nonzeros in the matrix.
   *
   * @param  extract_region
   *               [in] Determines what matix region to extract:<ul>
   *               <li> \c EXTRACT_FULL_MATRIX: Extract the nonzero elements
   *                    for the full matrix.
   *               <li> \c EXTRACT_UPPER_TRIANGULAR: Extract the nonzero
   *                    elements for the upper triangular region
   *                    (including the diagonal).
   *               <li> \c EXTRACT_LOWER_TRIANGULAR: Extract the nonzero
   *                    elements for the lower triangular region
   *                    (including the diagonal).
   *               </ul>
   * @param  element_uniqueness
   *                [in] Determines if element row and column indexes must be unique.<ul>
   *                <li> \c ELEMENTS_FORCE_UNIQUE: The row and column indexes must be unique.
   *                <li> \c ELEMENTS_ALLOW_DUPLICATES_SUM: Entries with duplicate row and
   *                     column indexes are allowed with the understanding that the values
   *                     will be summed.
   *                </ul>
   */
  virtual index_type num_nonzeros(
    EExtractRegion        extract_region
    ,EElementUniqueness   element_uniqueness
    ) const = 0;

  /** \brief Extract sparse elements in a coordinate data structure.
   *
   * The client can extract the structure in \c Arow[] and \c Acol[] and/or
   * the nonzero elements in \c Aval[].  If the client wants to extract
   * the structure it sets\
   * <tt>len_Aij = this->num_nonzeros(extract_region,element_uniqueness)</tt>
   * and then \c Arow[] and \c Acol[] are filled with the row and column indexes.
   * If the client wants the nonzero values theit sets
   * <tt>len_Aval = this->num_nonzeros(extract_region,element_uniqueness)</tt>
   * and then \c Aval[] will be set on output.
   *
   * The client can choose to extract the nonzeros for the full matrix
   * (<tt>extract_region == EXTRACT_FULL_MATRIX</tt>) or the upper
   * (<tt>extract_region == EXTRACT_UPPER_TRIANGULAR</tt>) or lower
   * (<tt>extract_region == EXTRACT_LOWER_TRIANGULAR</tt>) triangular regions.
   *
   * The client can force the extracted nonzero elements to have unique row
   * and column indexes (<tt>element_uniqueness == ELEMENTS_FORCE_UNIQUE</tt>)
   * or can allow elements with duplicate row and column indexes with the
   * understanding that the values will be summed
   * (<tt>element_uniqueness == ELEMENTS_ALLOW_DUPLICATES_SUM</tt>).
   *
   * @param  extract_region
   *              [in] Determines what matix region to extract:<ul>
   *              <li> \c EXTRACT_FULL_MATRIX: Extract the nonzero elements
   *                   for the full matrix.
   *              <li> \c EXTRACT_UPPER_TRIANGULAR: Extract the nonzero
   *                   elements for the upper triangular region
   *                   (including the diagonal).
   *              <li> \c EXTRACT_LOWER_TRIANGULAR: Extract the nonzero
   *                   elements for the lower triangular region
   *                   (including the diagonal).
   *              </ul>
   * @param  element_uniqueness
   *              [in] Determines if element row and column indexes must be unique.<ul>
   *              <li> \c ELEMENTS_FORCE_UNIQUE: The row and column indexes must be unique.
   *              <li> \c ELEMENTS_ALLOW_DUPLICATES_SUM: Entries with duplicate row and
   *                   column indexes are allowed with the understanding that the values
   *                   will be summed.
   *              </ul>
   * @param  len_Aval
   *              [in] If the client wants to extract the nonzero values
   *              of the matrix in \c Avar[] then \c len_Aval must be set
   *              to the return value from <tt>this->num_nonzeros(extract_region)</tt>
   *              Otherwise \c len_Avar should be set to zero.
   * @param  Aval
   *              [out] If <tt>len_Aval > 0</tt> then \c Aval must point to the begining
   *              of an array of length \c len_Aval.  If <tt>len_Aval == 0</tt> then \c Aval
   *              must be \c NULL.  On output, \c Aval[k] is the nonzero value of the kth
   *              nonzero element, for <tt>k = 0...len_Aval-1</tt>.
   * @param  len_Aij
   *              [in] If the client wants to extract the structure of the matrix in
   *              \c Arow[] and \c Acol[] then this must be set to
   *              <tt>len_Aij = this->num_nonzeros(extract_region)</tt>. Otherwise
   *              \c len_Aij should be set to zero.
   * @param  Arow
   *              [out] If <tt>len_Aij > 0</tt> then \c Arow must point to the begining
   *              of an array of length \c len_Aij.  If <tt>len_Aij == 0</tt> then \c Arow
   *              must be \c NULL.  On output, \c Arow[k] is the row index of the kth
   *              nonzero element, for <tt>k = 0...len_Aij-1</tt>.
   * @param  Acol
   *              [out] If <tt>len_Aij > 0</tt> then \c Acol must point to the begining
   *              of an array of length \c len_Aij.  If <tt>len_Aij == 0</tt> then \c Acol
   *              must be \c NULL.  On output, \c Acol[k] is the column index of the kth
   *              nonzero element, for <tt>k = 0...len_Aij-1</tt>.
   * @param  row_offset
   *              [in] The row index in \c Arow[] are offset by <tt>+row_offset</tt>
   *              on output.  To leave the first row index as 1 set <tt>row_offset = 0</tt>.
   *              This is to allow the client to place this matrix as a submatrix into a larger
   *              matrix stored in the coordinate format.  Default value = 0.	  
   * @param  col_offset
   *              [in] The column index in \c Acol[] are offset by <tt>+col_offset</tt>
   *              on output.  To leave the first column index as 1 set <tt>col_offset = 0</tt>.
   *              This is to allow the client to place this matrix as a submatrix into a larger
   *              matrix stored in the coordinate format.  Default value = 0.
   *
   * Preconditions<ul>
   * <li> [<tt>len_Aval != 0</tt>] <tt>len_Aval == this->num_nonzeros(etract_region,element_uniqueness)</tt>
   *      (throw <tt>std::invalid_argument</tt>).
   * <li> [<tt>len_Aij != 0</tt>] <tt>len_Aij == this->num_nonzeros(etract_region,element_uniqueness)</tt>
   *      (throw <tt>std::invalid_argument</tt>).
   * <li> [<tt>len_Aval == 0</tt>] <tt>Aval == NULL</tt> (throw <tt>std::invalid_argument</tt>).
   * <li> [<tt>len_Aij == 0</tt>] <tt>Arow == NULL</tt> (throw <tt>std::invalid_argument</tt>).
   * <li> [<tt>len_Aij == 0</tt>] <tt>Acol == NULL</tt> (throw <tt>std::invalid_argument</tt>).
   * </ul>
   *
   * To illustrate the behavior of this function consider the example matix \c A
   * (shown above) in the following examples:
   *
   * (1) To extract all of the nonzero elements and structure for the full matrix with the
   * row and column indexes offset by 2 and 4 respectively the following
   * input would yield the following output.
   *
   * Input:
   \verbatim

      extract_region     = EXTRACT_FULL_MATRIX
            element_uniqueness = ELEMENTS_FORCE_UNIQUE
      len_Aval           = A_nz = 12
      len_Aij            = A_nz = 12
      row_offset         = 2
      col_offset         = 4
   \endverbatim
   * Output:
   \verbatim

      k   Aval[k-1]  Arow[k-1]  Acol[k-1]
      -	---------  ---------  ---------
      1         1.1          3          5
      2         4.1          5          5
      3         5.1          7          5
      4         2.2          4          6
      5         3.2          5          6
      6         1.3          3          7
      7         4.3          6          7
      8         5.3          7          7
      9         1.4          3          8
      10        2.4          4          8
      11        4.4          6          8
      12        5.4          7          8
   \endverbatim
   * Note that in the example above that the elements are sorted
   * by column but they do not have to be on output.
   *
   * (2) To extract all of the nonzero elements and structure for the
   * upper triangular region of the matrix with the
   * row and column indexes offset by 0 and 0 respectively the following
   * input would yield the following output.
   *
   * Input:
   \verbatim

      extract_region     = EXTRACT_UPPER_TRIANGULAR
            element_uniqueness = ELEMENTS_FORCE_UNIQUE
      len_Aval           = A_up_nz = 6
      len_Aij            = A_up_nz = 6
      row_offset         = 0
      col_offset         = 0
   \endverbatim
   * Output:
   \verbatim

      k  Aval[k-1]  Arow[k-1]  Acol[k-1]
      -  ---------  ---------  ---------
      1        1.1          1          1
      2        2.2          2          2
      3        1.3          1          3
      4        1.4          1          4
      5        2.4          2          4
      6        4.4          4          4
   \endverbatim
   * Note that in the example above that the elements are sorted
   * by column but they do not have to be on output.
   *
   * (3) To extract all of the nonzero elements and structure for the
   * lower triangular region of the matrix with the row and column indexes
   * offset by 0 and 0 respectively the following input would yield the
   * following output.
   *
   * Input:
   \verbatim

      extract_region     = EXTRACT_LOWER_TRIANGULAR
            element_uniqueness = ELEMENTS_FORCE_UNIQUE
      len_Aval           = A_lo_nz = 9
      len_Aij            = A_lo_nz = 9
      row_offset         = 0
      col_offset         = 0
   \endverbatim
   * Output:
   \verbatim

      k  Aval[k-1]  Arow[k-1]  Acol[k-1]
      -  ---------  ---------  ---------
      1        1.1          1          1
      2        4.1          4          1
      3        5.1          5          1
      4        2.2          2          2
      5        3.2          3          2
      6        4.3          4          3
      7        5.3          5          3
      8        4.4          4          4
      9        5.4          5          4
   \endverbatim
   * Note that in the example above that the elements are sorted
   * by column but they do not have to be on output.
   */
  virtual void coor_extract_nonzeros(
    EExtractRegion                extract_region
    ,EElementUniqueness           element_uniqueness
    ,const index_type             len_Aval
    ,value_type                   Aval[]
    ,const index_type             len_Aij
    ,index_type                   Arow[]
    ,index_type                   Acol[]
    ,const index_type             row_offset = 0
    ,const index_type             col_offset = 0
    ) const = 0;
  
  // Todo: Include extraction routines for Compressed Row (Column) when needed.

};	// end class MatrixConvertToSparse

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_CONVERT_TO_SPARSE_H
