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

#ifndef ALAP_MATRIX_EXTRACT_SPARSE_ELEMENTS_H
#define ALAP_MATRIX_EXTRACT_SPARSE_ELEMENTS_H

#include "AbstractLinAlgPack_MatrixConvertToSparse.hpp"

namespace AbstractLinAlgPack {

/** \brief Interface for extracting nonzero elements from a banded subregion
 * of a permuted sparse matrix in one of several Fortran compatible formats.
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
 * This is meant to be the do-all interface for clients to use to extract nonzero elements
 * for sparse matrices.
 *
 * The idea is that given a matrix \a A, row and column permutations \c P and \c Q,
 * a range of rows and columns <tt>(rl,ru)</tt> and <tt>(cl,cu)</tt> defining a square submatrix
 * \a M and a range of lower and upper bands \c dl and \c du within the submatrix \a M, this
 * interface is used to extract those nonzero elements.
 *
 * In Matlab like terms we have:
 *
 * let <tt>M = (P'*A*Q)(rl:ru,cl:cu)</tt>
 *
 * This interface extracts nonzero elements from \a M within the banded region
 * <tt>[dl,du]</tt> where <tt>d = 0</tt> is the diagonal of \a M, <tt>d < 0</tt>
 * is below the diagonal and <tt>d > 0</tt> is above the diagonal.
 *
 * The following matrix is used in the documentation for the following
 * extraction functions.
 \verbatim

    [	1				6		9	]	1
    [			4				10	]	2
  A =	[			5					]	3
    [	2				7		11	]	4
    [	3				8		12	]	5

      1		2		3		 4

 \endverbatim
 *
 * Note that above, A has:<ul>
 * <li> <tt>A_nz == this->num_nonzeros(EXTRACT_FULL_MATRIX,ELEMENTS_FORCE_UNIQUE) == 12</tt>
 * <li> <tt>A_up_nz = this->num_nonzeros(EXTRACT_UPPER_TRIANGULAR,ELEMENTS_FORCE_UNIQUE) == 6</tt>
 * <li> <tt>A_lo_nz = this->num_nonzeros(EXTRACT_LOWER_TRIANGULAR,ELEMENTS_FORCE_UNIQUE) == 9</tt>
 * </ul>
 */
class MatrixExtractSparseElements
  : public virtual MatrixConvertToSparse
{
public:

  /** \brief Returns the number of nonzeros in the banded submatrix of the permuted matrix.
   *
   * Let \a B be the banded submatrix that is being specifed and let \a A be this
   * full matrix object.  Then the (dense) elemements of \a B are defined by:
   \verbatim

              / 0                                    : for dl <= (j-i) <= du
    B(i,j) =  |
              \ A(row_perm(i+rl-1),col_perm(j+cl-1)) : for (j-i) < dl || du < (j-i)

    for i = 1..ru-rl+1, j = 1..cu-cl+1
   \endverbatim
   * Above <tt>rl = row_rng.lbound()</tt>, ru = row_rng.ubound()</tt>,
   * <tt>cl = col_rng.lbound()</tt> and <tt>cu = col_rng.ubound()</tt>.
   *
   * Preconditions:<ul>
   *	<li> <tt>1 <= row_perm(i) <= this->rows(), i = 1..rows()</tt> (throw <tt>std::out_of_bounds</tt>)
   *	<li> <tt>1 <= col_perm(i) <= this->cols(), i = 1..cols()</tt> (throw <tt>std::out_of_bounds</tt>)
   *	<li> <tt>row_rng.ubound() <= this->rows()</tt> (throw <tt>std::out_of_bounds</tt>)
   *	<li> <tt>col_rng.ubound() <= this->cols()</tt> (throw <tt>std::out_of_bounds</tt>)
   *	<li> <tt>du >= dl</tt> (throw <tt>std::range_error</tt>)
   *	<li> <tt>-(ru-rl) <= dl && du <= (cu-cl)</tt> (throw <tt>std::out_of_bounds</tt>)
   *	</ul>
   *
   * To illustrate the behavior of this function consider the example matix \a A (see intro)
   * in the following example:
   *
   * Define the permutations as:
   \verbatim

   row_perm = { 2, 4, 1, 3, 5 }
   col_perm = { 3, 2, 1, 4 }
   \endverbatim
   * The permuted matrix <tt>(P'*A*Q)</tt> would then be:
   \verbatim

          2	|	1	[			4				10	]
          4	|	2	[	7				2		11	]
          1	|	3	[	6				1		9	]
    (P'*A*Q) =	3	|	4	[			5					]
          5	|	5	[	8				3		12	]

                  1		2		3		4
                  -		-		-		-
                  3		2		1		4
   \endverbatim
   * Now define the square submatrix in the range:
   \verbatim

   row_rng = [ rl, ru ] = [ 2, 5 ]
   col_rng = [ cl, cu ] = [ 2, 4 ]
   \endverbatim
   * The square submatrix is then:
   \verbatim

                4	|	1	[			2		11	]
    (P'*A*Q)(r1:r2,cl:cu) =	1	|	2	[			1		9	]
                3	|	3	[	5					]

                        1		2		3
                        -		-		-
                        2		1		4
   \endverbatim
   * Now define the range of diagonals as:
   \verbatim

   dl = -1, du = 1
   \endverbatim
   * This finally gives us our matrix that we wish to extract nonzeros
   * from as (the x's show elemements that where excluded out of the 
   * diagonal range:
   \verbatim

      4	|	1	[			4		x	]
    B =	1	|	2	[			1		9	]
      3	|	3	[	x					]

              1		2		3
              -		-		-
              2		1		4
   \endverbatim
   * Now we can see that for this example that this->count_nonzeros() would
   * return 3.
   *
   * In summary, for the example A shown above we have:
   *
   * Input:
   \verbatim
        element_uniqueness = ELEMENTS_FORCE_UNIQUE
    row_perm[]         = { 2, 4, 1, 3, 5 }
    col_perm[]         = { 3, 2, 1, 4 }
    row_rng            = [ 2, 5 ]
    col_rng            = [ 2, 4 ]
    dl                 = -1
    du                 = +1
   \endverbatim
   * Output:
   \verbatim

    3 <- count_nonzeros(element_uniqueness,row_perm,col_perm,row_rng,col_rng,dl,du)
   \endverbatim
   *
   * @param  element_uniqueness
   *              [in] Determines if element row and column indexes must be unique.<ul>
   *              <li> \c ELEMENTS_FORCE_UNIQUE: The row and column indexes must be unique.
   *              <li> \c ELEMENTS_ALLOW_DUPLICATES_SUM: Entries with duplicate row and
   *                   column indexes are allowed with the understanding that the values
   *                   will be summed.
   *              </ul>
   * @param  inv_row_perm
   *              [in] Arary (length \c this->rows()) Defines the row permutation \c P.
   *              Row \c i of \c A is row \c inv_row_perm[i-1] of \c <tt>(P'*A)</tt> .
   *              If <tt>inv_row_perm==NULL</tt> on input then the identity permutation
   *              <tt>P = I</tt> is used.
   * @param  inv_col_perm
   *              [in] Arary (length \c this->cols()) Defines the column permutation \c Q.
   *              Column \c j of \c A is column \c inv_col_perm[j-1] of <tt>(A*Q)</tt>.
   *              If <tt>inv_col_perm==NULL</tt> on input then the identity permutation
   *              <tt>Q = I</tt> is used.
   * @param  row_rng
   *              [in] Defines the range of rows <tt>[rl,ru]</tt> that the submatrix
   *              of <tt>(P'*A*Q)</tt> is taken from (see preconditions).
   * @param  col_rng
   *              [in] Defines the range of columns <tt>[cl,cu]</tt> that the submatrix
   *              of </tt>(P'*A*Q)</tt> is taken from (see preconditions).
   * @param  dl   [in] Lower diagonal to extract elements above (see preconditions).
   * @param  du   [in] Upper diagonal to extract elements below (see preconditions).
   */
  virtual index_type count_nonzeros(
    EElementUniqueness    element_uniqueness
    ,const index_type     inv_row_perm[]
    ,const index_type     inv_col_perm[]
    ,const Range1D        &row_rng
    ,const Range1D        &col_rng
    ,index_type           dl
    ,index_type           du
    ) const = 0;

  /** \brief Extract elements in a coordinate data structure.
   *
   * Let B be the scaled banded submatrix that is being specifed
   * and let A be this matrix object.  Then the elemements
   * of B are defined by:
   \verbatim

              / 0                                          : for dl <= (j-i) <= du
    B(i,j) =  |
              \ alpha*A(row_perm(i+rl-1),col_perm(j+cl-1)) : for (j-i) < dl || du < (j-i)

    for i = 1..ru-rl+1, j = 1..cu-cl+1
   \endverbatim
   * were <tt>rl = row_rng.lbound()</tt>, <tt>ru = row_rng.ubound()</tt>,
   * <tt>cl = col_rng.lbound()</tt> and <tt>cu = col_rng.ubound()</tt>.
   *
   * The client can extract the structure in \c Arow[] and \c Acol[] and/or
   * the nonzero elements in \c Aval[].  If the client wants to extract
   * the structure it sets
   * <tt>len_Aij = this->count_nonzeros()</tt>
   * and then \c Arow[] and \c Acol[] are filled with the row and column indexes.
   * If the client wants the nonzero values it sets
   * <tt>len_Aval = this->count_nonzeros()</tt>
   * and then \c Aval[] will be set on output.
   *
   * The input arguments passed to <tt>this->count_nonzeros()</tt> must correspond to the
   * same arguments passed to this function.
   *
   * To illustrate the behavior of this function consider the same example as
   * outlined in the documentation for \c count_nonzeros().  Let's assume
   * that we want to scale the example banded matrix by <tt>alpha = 2.0</tt>,
   * make the row and column indexes start at (4,7) and extract the
   * structure (<tt>len_Aij > 0</tt>) and the nonzero values (<tt>len_Aval</tt>).
   * The input and output from this function for this example would then be:
   *
   * Input:
   \verbatim
        elements_uniqueness = ELEMENTS_FORCE_UNIQUE
    row_perm            = { 2, 4, 1, 3, 5 }
    col_perm            = { 3, 2, 1, 4 }
    row_rng             = [ 2, 5 ]
    col_rng             = [ 2, 4 ]
    dl                  = -1
    du                  = +1
    alpha               = 2.0
    len_Aval            = count_nonzeros(...) = 3
    len_Aij             = count_nonzeros(...) = 3
    row_offset          = 3
    col_offset          = 6 
   \endverbatim
   * Output:
   \verbatim

    k  A(i,j)  B(i,j)  Avar  Arow  Acol
    -  ------  ------  ----  ----  ----
    1  4(4,1)  4(1,2)     8     4     8
    2  1(1,1)  1(2,2)     2     5     8
    3  9(1,5)  9(2,3)    18     5     9
   \endverbatim
   * Some of the most common uses of this interface are:<ul>
   * <li> Extract rectuangular submatrices: <tt>(bl = -(ru-rl+1) && bu = (cu-cl+1))</tt>
   * <li> Extract upper triangular region:  <tt>(bl = 0 && bu = (cu-cl+1))</tt>
   * <li> Extract lower triangular region:  <tt>(bl = -(ru-rl+1) && bu = 0)</tt>
   * </ul>
   *
   * @param  element_uniqueness
   *                  [in] Same as passed to \c count_nonzeros(...).
   * @param  row_perm [in] Same as passed to \c count_nonzeros(...).
   * @param  col_perm [in] Same as passed to \c count_nonzeros(...).
   * @param  row_rng  [in] Same as passed to \c count_nonzeros(...).
   * @param  col_rng  [in] Same as passed to \c count_nonzeros(...).
   * @param  dl       [in] Same as passed to \c count_nonzeros(...).
   * @param  du       [in] Same as passed to \c count_nonzeros(...).
   * @param  alpha    [in] Scaling parameter (see above)
   * @param  len_Aval [in] If the client wants to extract the nonzero values
   *                  of the matrix in \c Aval[] then this should be set to
   *                  <tt>len_Aval = this->count_nonzeros(...)</tt>.
   *                  Otherwise \c len_Avar should be set to zero.
   * @param  Aval     [out] If <tt>len_Aval > 0</tt> then \c Aval must point to
   *                  the begining of an array of length \c len_Aval.
   *                  If <tt>len_Aval == 0</tt> then \c Aval may be NULL.
   *                  On output, \c Aval[k] is the nonzero value of the kth
   *                  nonzero element, for <tt>k = 0...len_Aval-1</tt>.
   * @param  len_Aij  [in] If the client wants to extract the structure
   *                  of the matrix in \c Arow[] and \c Acol[] then this
   *                  should be set to <tt>len_Aij = this->num_nonzeros(...)</tt>.
   *                  Otherwise \c len_Aij should be set to zero.
   * @param  Arow     [out] If <tt>len_Aij > 0</tt> then \c Arow must point to
   *                  the begining of an array of length \c len_Aij.
   *                  If <tt>len_Aij == 0</tt> then Arow must be \c NULL.
   *                  On output, \c Arow[k] is the row index of the kth
   *                  nonzero element, for <tt>k = 0...len_Aij-1</tt>.
   * @param  Acol     [out] If <tt>len_Aij > 0</tt> then \c Acol must point
   *                  to the begining of an array of length \c len_Aij.
   *                  If <tt>len_Aij == 0</tt> then \c Acol must be \c NULL.
   *                  On output, \c Acol[k] is the column index of the kth
   *                  nonzero element, for <tt>k = 0...len_Aij-1</tt>.
   * @param  row_offset
   *                  [in] The row indexes in \c Arow[] are offset by
   *                  <tt>+row_offset</tt> on output.  To leave the first
   *                  row index as 1 set <tt>row_offset = 0</tt>.  This is to
   *                  allow the client to place this matrix as a submatrix
   *                  into a larger matrix stored in the coordinate format.
   *                  Default value = 0.	  
   * @param  col_offset
   *                  [in] The column indexes in \c Acol[] are offset by
   *                  <tt>+col_offset</tt> on output.  To leave the first
   *                  column index as 1 set <tt>col_offset = 0</tt>.  This is to
   *                  allow the client to place this matrix as a submatrix
   *                  into a larger matrix stored in the coordinate format.
   *                  Default value = 0.	  
   */
  virtual void coor_extract_nonzeros(
    EElementUniqueness    element_uniqueness
    ,const index_type     inv_row_perm[]
    ,const index_type     inv_col_perm[]
    ,const Range1D        &row_rng
    ,const Range1D        &col_rng
    ,index_type           dl
    ,index_type           du
    ,value_type           alpha
    ,const index_type     len_Aval
    ,value_type           Aval[]
    ,const index_type     len_Aij
    ,index_type           Arow[]
    ,index_type           Acol[]
    ,const index_type     row_offset = 0
    ,const index_type     col_offset = 0
    ) const = 0;

  // ToDo: Add methods for extracting compressed row (column) elements!

  /** @name Overridden from MatrixConvertToSparse */
  //@{

  /** \brief . */
  index_type num_nonzeros(
    EExtractRegion        extract_region
    ,EElementUniqueness   element_uniqueness
    ) const;
  /** \brief . */
  void coor_extract_nonzeros(
    EExtractRegion                extract_region
    ,EElementUniqueness           element_uniqueness
    ,const index_type             len_Aval
    ,value_type                   Aval[]
    ,const index_type             len_Aij
    ,index_type                   Arow[]
    ,index_type                   Acol[]
    ,const index_type             row_offset
    ,const index_type             col_offset
    ) const;

  //@}

private:

  /** \brief . */
  void get_dl_du( EExtractRegion extract_region, index_type* dl, index_type* du ) const;

};	// end class MatrixExtractSparseElements

}	// end namespace AbstractLinAlgPack 

#endif	// ALAP_MATRIX_EXTRACT_SPARSE_ELEMENTS_H
