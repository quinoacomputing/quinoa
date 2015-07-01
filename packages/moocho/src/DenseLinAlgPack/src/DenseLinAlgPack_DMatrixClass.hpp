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
// General matrix and matrix region (slice) classes

#ifndef GEN_MATRIX_CLASS_H
#define GEN_MATRIX_CLASS_H

#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixAssign.hpp"

/* * @name {\bf Dense 2-D Rectangular Matrix Absractions}.
  *
  * The class DMatrix is a storage class for 2-D matrices while the class DMatrixSlice
  * is used to represent rectangular regions of a DMatrix object.
  */
// @{
//		begin General Rectangular 2-D Matrices scope

namespace DenseLinAlgPack {

class DMatrix;

/* * @name {\bf General Matrix Classes}. */
// @{

// ////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenMatrixClass
//

/** \brief . */
/* * 2-D General Rectangular Matrix Subregion Class (Slice) (column major).
  *
  * This class is used to represent a rectangular matrix.  It uses a BLAS-like
  * slice of a raw C++ array.  Objects of this class can represent
  * an entire matrix or any rectangular subregion of a matrix.
  */
class DMatrixSlice {
public:
  /* * @name {\bf Nested Member Types (STL)}.
    *
    * These nested types give the types used in the interface to the class.
    *
    * \begin{description}
    *	<li>[#value_type#]				- type being stored in the underlying valarray<>			
    *	<li>[#size_type#]				- type for the number rows and coluns
    *	<li>[#difference_type#]		- type for the stride between elements (in a row or column)
    *	<li>[#reference#]				- value_type&
    *	<li>[#const_reference#]		- const value_type&
    * \end{description}
    */
  // @{
  // @}
  typedef DenseLinAlgPack::value_type					value_type;
  typedef DenseLinAlgPack::size_type					size_type;
  typedef ptrdiff_t								difference_type;
  typedef value_type&								reference;
  typedef const value_type&						const_reference;

  /* * @name {\bf Constructors}.
    *
    * These constructors are used by the other entities in the library
    * to create DMatrixSlices.  In general user need not use these
    * constructors directly.  Instead, the general user should use the 
    * subscript operators to create subregions of DMatrix and DMatrixSlice
    * objects.
    * 
    * The default copy constructor is used and is therefore not shown here.
    */

  // @{

  /// 
  /* * Construct an empty view.
    *
    * The client can then call bind(...) to bind the view.
    */
  DMatrixSlice();

  /// 
  /* * Construct a veiw of a matrix from a raw C++ array.
    *
    * The DMatrixSlice constructed represents a 2-D matrix whos elements are stored
    * in the array starting at #ptr#.  This is how the BLAS represent general rectangular
    * matrices.
    * The class can be used to provide a non-constant view the elements (#DMatrix#)
    * or a constant view (#const DMatrixSlice#).  Here is an example of how to
    * create a constant view:
    *
    \verbatim
    const DMatrixSlice::size_type m = 4, n = 4;
    const DMatrixSlice::value_type ptr[m*n] = { ... };
    const GenMatrixslice mat(cosnt_cast<DMatrixSlice::value_type*>(ptr),m*n,m,m,n);
    \endverbatim
    *
    * The #const_cast<...># such as in the example above is perfectly safe to use
    * when constructing  #const# veiw of #const# data.  On the other hand casting
    * away #const# and then non-#const# use is not safe in general since the original
    * #const# data may reside in ROM for example.  By using a non-#const# pointer in the
    * constructor you avoid accidentally making a non-#const# view of #const# data.
    *
    * Preconditions: <ul>
    *		<li> #rows <= max_rows# (throw out_of_range)
    *		<li> #start + rows + max_rows * (cols - 1) <= v.size()# (throw out_of_range)
    *		</ul>
    *
    * @param	ptr			pointer to the storage of the elements of the matrix (column oriented).
    * @param	size		total size of the storage pointed to by #ptr# (for size checking)
    * @param	max_rows	number of rows in the full matrix this sub-matrix was taken from.
    *						This is equivalent to the leading dimension argument (LDA) in
    *						the BLAS.
    * @param	rows		number of rows this matrix represents
    *	@param	cols		number of columns this matrix represents
    */
  DMatrixSlice( value_type* ptr, size_type size
    , size_type max_rows, size_type rows, size_type cols );

  /// 
  /* * Construct a submatrix region of another DMatrixSlice.
    *
    * This constructor simplifies the creation of subregions using the subscript
    * operators.
    *
    * I and J must be bounded ranges (full_range() == false).
    *
    * Preconditions: <ul>
    *		<li> #I.full_range() == false# (throw out_of_range)
    *		<li> #J.full_range() == false# (throw out_of_range)
    *		<li> #I.ubound() <= gms.rows()# (throw out_of_range)
    *		<li> #J.ubound() <= gms.cols()# (throw out_of_range)
    *		</ul>
    */
  DMatrixSlice( DMatrixSlice& gms, const Range1D& I
    , const Range1D& J );

  // @}

  /** \brief . */
  /* * Set to the view of the input DMatrixSlice.
    *
    *
    */
  void bind( DMatrixSlice gms );

  /* * @name {\bf Dimensionality, Misc}. */
  // @{

  /// Return the number of rows
  size_type		rows() const;
  /// Return the number of columns
  size_type		cols() const;

  /// 
  /* * Returns the degree of memory overlap of #this# and the DMatrixSlice object #gms#.
    *
    * @return 
    *		\begin{description}
    *		<li>[NO_OVERLAP]	There is no memory overlap between #this# and #gms#.
    *		<li>[SOME_OVERLAP]	There is some memory locations that #this# and #gms# share.
    *		<li>[SAME_MEM]		The DMatrixSlice objects #this# and #gms# share the exact same memory locations.
    *		\end{description}
    */
  EOverLap overlap(const DMatrixSlice& gms) const;

  // @}

  /* * @name {\bf Individual Element Access Subscripting (lvalue)}. */
  // @{

  /// Return element at row i, col j (i,j) (1-based) (throws std::out_of_range if i, j are out of bounds)
  reference		operator()(size_type i, size_type j);
  /// Return element at row i, col j (i,j) (1-based) (throws std::out_of_range if i, j are out of bounds)
  const_reference	operator()(size_type i, size_type j) const;

  // @}

  /* * @name {\bf Subregion Access (1-based)}.
    *
    * These member functions allow access to subregions of the DMatrixSlice object.
    * The functions, row(i), col(j), and diag(k) return DVectorSlice objects while
    * the subscripting operators opeator()(I,J) return DMatrixSlice objects for
    * rectangular subregions.
    */
  // @{

  /// Return DVectorSlice object representing the ith row (1-based; 1,2,..,#this->rows()#, or throw std::out_of_range)
  DVectorSlice			row(size_type i);
  /// Same as above
  const DVectorSlice	row(size_type i) const;
  /// Return DVectorSlice object representing the jth column (1-based; 1,2,..,#this->cols()#, or throw std::out_of_range)
  DVectorSlice			col(size_type j);
  /// Same as above
  const DVectorSlice	col(size_type j) const;
  /// 
  /* * Return DVectorSlice object representing a diagonal.
    *
    * Passing k == 0 returns the center diagonal.  Values of k < 0 are the lower diagonals
    * (k = -1, -2, ..., -#this->rows()# + 1).  Values of k > 0 are the upper diagonals
    * (k = 1, 2, ..., #this->cols()# - 1).
    *
    * Preconditions: <ul>
    *		<li> #[k < 0] k <= this->rows() + 1# (throw out_of_range)
    *		<li> #[k > 0] k <= this->cols() + 1# (throw out_of_range)
    *		</ul>
    */
  DVectorSlice			diag(difference_type k = 0);
  /// Same as above.
  const DVectorSlice	diag(difference_type k = 0) const;
  /// 
  /* * Extract a rectangular subregion containing rows I, and columns J.
    *
    * This operator function returns a DMatrixSlice that represents a
    * rectangular region of this DMatrixSlice.  This submatrix region
    * represents the rows I.lbound() to I.ubound() and columns J.lbound()
    * to J.lbound().  If I or J is unbounded (full_range() == true, constructed
    * with Range1D()), the all of the rows or columns respectively will be
    * selected.  For example. To select all the rows and the first 5 columns of
    * a matrix #A# you would use #A(Range1D(),Range1D(1,5))#.
    *
    * Preconditions: <ul>
    *		<li> #[I.full_range() == false] I.ubound() <= this->rows()# (throw out_of_range)
    *		<li> #[J.full_range() == false] J.ubound() <= this->cols()# (throw out_of_range)
    *		</ul>
    */
  DMatrixSlice operator()(const Range1D& I, const Range1D& J);
  /// Same as above.
  const DMatrixSlice operator()(const Range1D& I, const Range1D& J) const;
  /// 
  /* * Extract a rectangular subregion containing rows i1 to i2, and columns j1 to j2.
    *
    * This operator function returns a DMatrixSlice that represents a
    * rectangular region of this DMatrixSlice.  This submatrix region
    * represents the rows i1 to 12 and colunms j1 to j2.
    * 
    * Preconditions: <ul>
    *		<li> #i1 <= i2# (throw out_of_range)
    *		<li> #i2 <= this->rows()# (throw out_of_range)
    *		<li> #j1 <= j2# (throw out_of_range)
    *		<li> #j2 <= this->cols()# (throw out_of_range)
    *		</ul>
    */
  DMatrixSlice operator()(size_type i1, size_type i2, size_type j1
    , size_type j2);
  /// Same as above.
  const DMatrixSlice operator()(size_type i1, size_type i2, size_type j1
    , size_type j2) const;
  /// Allow the address to be taken of an rvalue of this object.
  DMatrixSlice* operator&() {
    return this;
  }
  /** \brief . */
  const DMatrixSlice* operator&() const {
    return this;
  }
  /// Return reference of this.  Included for iniformity with DMatrix
  DMatrixSlice& operator()();
  /// Same as above
  const DMatrixSlice& operator()() const;

  // @}

  /* * @name {\bf Assignment operators}. */
  // @{

  /** \brief . */
  /* * Sets all elements = alpha
    *
    * If the underlying valarray is unsized (#this->v().size() == 0#) the matrix is sized to 1 x 1
    * and the single element is set to alpha.
    *
    * Postcondtions: <ul>
    *		<li> #this->operator()(i,j) == alpha#, i = 1,2,...,#this->rows()#, j = 1,2,...,#this->cols()#
    *		</ul>
    */
  DMatrixSlice& operator=(value_type alpha);
  /** \brief . */
  /* *  Copies all of the elements of the DMatrixSlice, #rhs#, into the elements of #this#.
    *
    * If the underlying valarray is unsized (#this->v().size() == 0#) the matrix is sized to
    * the size of the rhs matrix.
    *
    * Precondtions: <ul>
    *		<li> #this->rows() == gms_rhs.rows()# (throw length_error)
    *		<li> #this->cols() == gms_rhs.cols()# (throw length_error)
    *		</ul>
    *
    * Postcondtions: <ul>
    *		<li> #this->operator()(i,j) == gms_rhs(i,j)#, i = 1,2,...,#this->rows()#, j = 1,2,...,#this->cols()#
    *		</ul>
    */
  DMatrixSlice& operator=(const DMatrixSlice& gms_rhs);

  // @}

  /* * @name {\bf Raw data access}.
    */
  // @{

  /// Return the number of rows in the full matrix. Equivalent to BLAS LDA argument.
  size_type			max_rows() const;
  /** \brief . */
  /* * Return pointer to the first element in the underlying array the jth
     * col (j is 1-based here [1,cols]).  If unsized col_ptr(1) returns zero if unsized.
     */
  value_type*			col_ptr(size_type j);
  /// Same as above.
  const value_type*	col_ptr(size_type j) const;

  // @}

private:
  value_type		*ptr_;		// contains the data for the matrix region
  size_type		max_rows_,	// the number of rows in the full matrix
          rows_,		// the number of rows in this matrix region
          cols_;		// the number of cols in this matrix region

  // Assert the row subscript is in bounds (1-based), (throw std::out_of_range)
  void validate_row_subscript(size_type i) const;
  // Assert the column subscript is in bounds (1-based), (throw std::out_of_range)
  void validate_col_subscript(size_type j) const;
  // Assert that a constructed DMatrixSlice has a valid range, (throw std::out_of_range)
  void validate_setup(size_type size) const;
  
  // Get a diagonal
  DVectorSlice p_diag(difference_type k) const;

};	// end class DMatrixSlice

/** \brief . */
/* * 2-D General Rectangular Matrix (column major) Storage Class.
  *
  * This class provides the storage for 2-D rectangular matrices.
  */
class DMatrix {
public:
  /* * @name {\bf Nested Member Types (STL)}.
    *
    * These nested types give the types used in the interface to the class.
    *
    * \begin{description}
    *	<li>[#value_type#]				type being stored in the underlying valarray<>			
    *	<li>[#size_type#]				type for the number rows and coluns
    *	<li>[#difference_type#]		type for the stride between elements (in a row or column)
    *	<li>[#reference#]				value_type&
    *	<li>[#const_reference#]		const value_type&
    * \end{description}
    */
  // @{
  // @}
  
  typedef DenseLinAlgPack::value_type					value_type;
  typedef DenseLinAlgPack::size_type					size_type;
  typedef ptrdiff_t								difference_type;
  typedef value_type&								reference;
  typedef const value_type&						const_reference;
  typedef std::valarray<value_type>				valarray;

  /* * @name {\bf Constructors}.
    *
    * The general user uses these constructors to create a matrix.  
    *
    * The default constructor is used and is therefore not shown here.
    */
  // @{

  /// Construct a matrix with rows = cols = 0
  DMatrix();
  /// Construct an uninitialied rectangular matrix (rows x cols) 
  explicit DMatrix(size_type rows, size_type cols);
  /** \brief . */
  /* * Construct rectangular matrix (rows x cols) with elements initialized to val.
    *
    * Postconditions: <ul>
    *		<li> #this->operator()(i,j) == val#, i = 1,2,...,#rows#, j = 1,2,...,#cols#  
    *		</ul>
    */
  explicit DMatrix(value_type val, size_type rows, size_type cols);
  /// 
  /* * Construct rectangular matrix (rows x cols) initialized to elements of p (by column).
    *
    * Postconditions: <ul>
    *		<li> #this->operator()(i,j) == p[i-1 + rows * (j - 1)]#, i = 1,2,...,#rows#, j = 1,2,...,#cols#  
    *		</ul>
    */
  explicit DMatrix(const value_type* p, size_type rows, size_type cols);
  /** \brief . */
  /* * Construct a matrix from the elements in another DMatrixSlice, #gms#.
    *
    * Postconditions: <ul>
    *		<li> #this->operator()(i,j) == gms(i,j)#, i = 1,2,...,#rows#, j = 1,2,...,#cols#  
    *		</ul>
    */
  DMatrix(const DMatrixSlice& gms);

  // @}

  /* * @name {\bf Memory Management, Dimensionality, Misc}. */
  // @{
  
  /// Resize matrix to a (rows x cols) matrix and initializes any added elements by val
  void resize(size_type rows, size_type cols, value_type val = value_type());

  /// frees memory and leaves a (0 x 0) matrix
  void free();
  
  /// Return the number of rows
  size_type		rows() const;

  /// Return the number of columns
  size_type		cols() const;

  /// 
  /* * Returns the degree of memory overlap of #this# and the DMatrixSlice object #gms#.
    *
    * @return 
    *		\begin{description}
    *		<li>[NO_OVERLAP]	There is no memory overlap between #this# and #gms#.
    *		<li>[SOME_OVERLAP]	There is some memory locations that #this# and #gms# share.
    *		<li>[SAME_MEM]		The DMatrixSlice objects #this# and #gms# share the exact same memory locations.
    *		\end{description}
    */
  EOverLap overlap(const DMatrixSlice& gms) const;

  // @}

  /* * @name {\bf Individual Element Access Subscripting (lvalue)}. */
  // @{

  /// Return element at row i, col j (i,j) (1-based)
  reference		operator()(size_type i, size_type j);

  /// Return element at row i, col j (i,j) (1-based)
  const_reference	operator()(size_type i, size_type j) const;

  // @}

  /* * @name {\bf Subregion Access (1-based)}.
    *
    * These member functions allow access to subregions of the DMatrix object.
    * The functions, row(i), col(j), and diag(k) return DVectorSlice objects while
    * the subscripting operators opeator()(I,J) return DMatrixSlice objects for
    * rectangular subregions.
    */
  // @{

  /// Return DVectorSlice object representing the ith row (1-based; 1,2,..,#this->rows()#)
  DVectorSlice			row(size_type i);

  /** \brief . */
  const DVectorSlice	row(size_type i) const;

  /// Return DVectorSlice object representing the jth column (1-based; 1,2,..,#this->cols()#)
  DVectorSlice			col(size_type j);

  /** \brief . */
  const DVectorSlice	col(size_type j) const;

  /// 
  /* * Return DVectorSlice object representing a diagonal.
    *
    * Passing k == 0 returns the center diagonal.  Values of k < 0 are the lower diagonals
    * (k = -1, -2, ..., #this->rows()# - 1).  Values of k > 0 are the upper diagonals
    * (k = 1, 2, ..., #this->cols()# - 1).
    *
    * Preconditions: <ul>
    *		<li> #[k < 0] k <= this->rows() + 1# (throw out_of_range)
    *		<li> #[k > 0] k <= this->cols() + 1# (throw out_of_range)
    *		</ul>
    */
  DVectorSlice			diag(difference_type k = 0);

  /** \brief . */
  const DVectorSlice	diag(difference_type k = 0) const;

  /// 
  /* * Extract a rectangular subregion containing rows I, and columns J.
    *
    * This operator function returns a DMatrixSlice that represents a
    * rectangular region of this DMatrixSlice.  This submatrix region
    * represents the rows I.lbound() to I.ubound() and columns J.lbound()
    * to J.lbound().  If I or J is unbounded (full_range() == true, constructed
    * with Range1D()), the all of the rows or columns respectively will be
    * selected.  For example. To select all the rows and the first 5 columns of
    * a matrix #A# you would use #A(Range1D(),Range1D(1,5))#.
    *
    * Preconditions: <ul>
    *		<li> #[I.full_range() == false] I.ubound() <= this->rows()# (throw out_of_range)
    *		<li> #[J.full_range() == false] J.ubound() <= this->cols()# (throw out_of_range)
    *		</ul>
    */
  DMatrixSlice operator()(const Range1D& I, const Range1D& J);

  /** \brief . */
  const DMatrixSlice operator()(const Range1D& I, const Range1D& J) const;

  /// 
  /* * Extract a rectangular subregion containing rows i1 to i2, and columns j1 to j2.
    *
    * This operator function returns a DMatrixSlice that represents a
    * rectangular region of this DMatrixSlice.  This submatrix region
    * represents the rows i1 to 12 and colunms j1 to j2.
    * 
    * Preconditions: <ul>
    *		<li> #i1 <= i2# (throw out_of_range)
    *		<li> #i2 <= this->rows()# (throw out_of_range)
    *		<li> #j1 <= j2# (throw out_of_range)
    *		<li> #j2 <= this->cols()# (throw out_of_range)
    *		</ul>
    */
  DMatrixSlice operator()(size_type i1, size_type i2, size_type j1
    , size_type j2);

  /** \brief . */
  const DMatrixSlice operator()(size_type i1, size_type i2, size_type j1
    , size_type j2) const;

  /// Return a DMatrixSlice that represents this entire matrix.
  DMatrixSlice operator()();

  /** \brief . */
  const DMatrixSlice operator()() const;

  // @}

  /* * @name {\bf Implicit conversion operators}.
    *
    * These functions allow for the implicit converstion from a DMatrix to a DMatrixSlice.
    * This implicit converstion is important for the proper usage of much of the
    * libraries functionality.
    */
  // @{

  /** \brief . */
  operator DMatrixSlice();
  /** \brief . */
  operator const DMatrixSlice() const;

  // @}

  /* * @name {\bf Assignment Operators}. */
  // @{
  
  /** \brief . */
  /* * Sets all elements = alpha
    *
    * If the underlying valarray is unsized (#this->v().size() == 0#) the matrix is sized to 1 x 1
    * and the single element is set to alpha.
    *
    * Postcondtions: <ul>
    *		<li> #this->operator()(i,j) == alpha#, i = 1,2,...,#this->rows()#, j = 1,2,...,#this->cols()#
    */
  DMatrix& operator=(value_type rhs);
  /** \brief . */
  /* * Copies all of the elements of the DMatrixSlice, #rhs#, into the elements of #this#.
    *
    * If #this# is not the same size as gms_rhs the #this# is resized.
    *
    * Postcondtions: <ul>
    *		<li> #this->operator()(i,j) == gms_rhs(i,j)#, i = 1,2,...,#this->rows()#, j = 1,2,...,#this->cols()#
    */
  DMatrix& operator=(const DMatrixSlice& gms_rhs);
  /// Same as above.  Needed to override the default assignment operator.
  DMatrix& operator=(const DMatrix& rhs);

  // @}

  /* * @name {\bf Raw data access}.
    */
  // @{

  /// Return the number of rows in the full matrix. Equivalent to BLAS LDA argument.
  size_type			max_rows() const;
  /** \brief . */
  /* * Return pointer to the first element in the underlying array the jth
     * col (j is 1-based here [1,cols]).  If unsized col_ptr(1) returns zero if unsized.
     */
  value_type*			col_ptr(size_type j);
  /// Same as above.
  const value_type*	col_ptr(size_type j) const;

  // @}

private:
  std::valarray<value_type>	v_;
  size_type					rows_;

  // Assert the row subscript is in bounds (1-based), (throw std::out_of_range)
  void validate_row_subscript(size_type i) const;
  // Assert the column subscript is in bounds (1-based), (throw std::out_of_range)
  void validate_col_subscript(size_type j) const;

  // Get a diagonal, (throw std::out_of_range)
  DVectorSlice p_diag(difference_type k) const;

};	// end class DMatrix

//		end General Matix Classes scope
// @}

// ///////////////////////////////////////////////////////////////////////////////
// Non-member function declarations												//
// ///////////////////////////////////////////////////////////////////////////////

/* * @name {\bf DMatrix / DMatrixSlice Associated Non-Member Functions}. */
// @{
//		begin non-member functions scope

inline 
///  
/* * Explicit conversion function from DMatrix to DMatrixSlice.
  *
  * This is needed to allow a defered evaluation class (TCOL) to be evaluated using its
  * implicit conversion operator temp_type() (which returns DMatrix for DMatrixSlice
  * resulting expressions).
  */
//DMatrixSlice EvaluateToDMatrixSlice(const DMatrix& gm)
//{	return DMatrixSlice(gm); }

/// Assert two matrices are the same size and throws length_error if they are not (LINALGPACK_CHECK_RHS_SIZES).
void assert_gms_sizes(const DMatrixSlice& gms1, BLAS_Cpp::Transp trans1, const DMatrixSlice& gms2
  , BLAS_Cpp::Transp trans2);

inline 
/// Assert a matrix is square and throws length_error if it is not (LINALGPACK_CHECK_SLICE_SETUP).
void assert_gms_square(const DMatrixSlice& gms) {
#ifdef LINALGPACK_CHECK_SLICE_SETUP
  if(gms.rows() != gms.cols())
    throw std::length_error("Matrix must be square");
#endif
} 

inline 
/** \brief . */
/* * Utility to check if a lhs matrix slice is the same size as a rhs matrix slice.
  *
  * A DMatrixSlice can not be resized since the rows_ property of the
  * DMatrix it came from will not be updated.  Allowing a DMatrixSlice
  * to resize from unsized would require that the DMatrixSlice carry
  * a reference to the DMatrix it was created from.  If this is needed
  * then it will be added.
  */
void assert_gms_lhs(const DMatrixSlice& gms_lhs, size_type rows, size_type cols
  , BLAS_Cpp::Transp trans_rhs = BLAS_Cpp::no_trans)
{
  if(trans_rhs == BLAS_Cpp::trans) std::swap(rows,cols);
  if(gms_lhs.rows() == rows && gms_lhs.cols() == cols) return; // same size
  // not the same size so is an error
  throw std::length_error("assert_gms_lhs(...):  lhs DMatrixSlice dim does not match rhs dim");
}

/* * @name Return rows or columns from a possiblly transposed DMatrix or DMatrixSlice. */
// @{

inline 
/** \brief . */
DVectorSlice row(DMatrixSlice& gms, BLAS_Cpp::Transp trans, size_type i) {
  return (trans ==  BLAS_Cpp::no_trans) ? gms.row(i) : gms.col(i);
} 

inline 
/** \brief . */
DVectorSlice col(DMatrixSlice& gms, BLAS_Cpp::Transp trans, size_type j) {
  return (trans ==  BLAS_Cpp::no_trans) ? gms.col(j) : gms.row(j);
} 

inline 
/** \brief . */
const DVectorSlice row(const DMatrixSlice& gms, BLAS_Cpp::Transp trans, size_type i) {
  return (trans ==  BLAS_Cpp::no_trans) ? gms.row(i) : gms.col(i);
} 

inline 
/** \brief . */
const DVectorSlice col(const DMatrixSlice& gms, BLAS_Cpp::Transp trans, size_type j) {
  return (trans ==  BLAS_Cpp::no_trans) ? gms.col(j) : gms.row(j);
} 

inline 
/** \brief . */
DVectorSlice row(DMatrix& gm, BLAS_Cpp::Transp trans, size_type i) {
  return (trans ==  BLAS_Cpp::no_trans) ? gm.row(i) : gm.col(i);
} 

inline 
/** \brief . */
DVectorSlice col(DMatrix& gm, BLAS_Cpp::Transp trans, size_type j) {
  return (trans ==  BLAS_Cpp::no_trans) ? gm.col(j) : gm.row(j);
} 

inline 
/** \brief . */
const DVectorSlice row(const DMatrix& gm, BLAS_Cpp::Transp trans, size_type i) {
  return (trans ==  BLAS_Cpp::no_trans) ? gm.row(i) : gm.col(i);
} 

inline 
/** \brief . */
const DVectorSlice col(const DMatrix& gm, BLAS_Cpp::Transp trans, size_type j) {
  return (trans ==  BLAS_Cpp::no_trans) ? gm.col(j) : gm.row(j);
} 

// @}

inline 
/// Utility to resize a DMatrix to the size of a rhs matrix.
void resize_gm_lhs(DMatrix* gm_rhs, size_type rows, size_type cols
  , BLAS_Cpp::Transp trans_rhs)
{
  if(trans_rhs == BLAS_Cpp::trans) std::swap(rows,cols);
  gm_rhs->resize(rows,cols);
}

//		end non-member functions scope
// @}

//		end General Rectangular 2-D Matrices scope
// @}

// ////////////////////////////////////////////////////////////////////////////////
// Inline definitions of computationally independent member function definitions //
// ////////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// DMatrixSlice inline member function definitions

// Private utilities

#ifndef LINALGPACK_CHECK_RANGE
inline
void DMatrixSlice::validate_row_subscript(size_type i) const
{}
#endif

#ifndef LINALGPACK_CHECK_RANGE
inline
void DMatrixSlice::validate_col_subscript(size_type j) const
{}
#endif

#ifndef LINALGPACK_CHECK_SLICE_SETUP
inline
void DMatrixSlice::validate_setup(size_type size) const
{}
#endif

// Constructors

inline
DMatrixSlice::DMatrixSlice()
  : ptr_(0), max_rows_(0), rows_(0), cols_(0)
{}

inline
DMatrixSlice::DMatrixSlice( value_type* ptr, size_type size
    , size_type max_rows, size_type rows, size_type cols )
  : ptr_(ptr), max_rows_(max_rows), rows_(rows), cols_(cols)
{	
  validate_setup(size);
}

inline
DMatrixSlice::DMatrixSlice( DMatrixSlice& gms, const Range1D& I
    , const Range1D& J)
  : ptr_( gms.col_ptr(1) + (I.lbound() - 1) + (J.lbound() - 1) * gms.max_rows() )
  , max_rows_(gms.max_rows())
  , rows_(I.size())
  , cols_(J.size())
{	
  gms.validate_row_subscript(I.ubound());
  gms.validate_col_subscript(J.ubound());
}

inline
void DMatrixSlice::bind(DMatrixSlice gms) {
  ptr_		= gms.ptr_;
  max_rows_	= gms.max_rows_;
  rows_		= gms.rows_;
  cols_		= gms.cols_;
}

// Size / Dimensionality

inline
DMatrixSlice::size_type DMatrixSlice::rows() const {
  return rows_;
}

inline
DMatrixSlice::size_type DMatrixSlice::cols() const {
  return cols_;
}

// Misc

// Element access

inline
DMatrixSlice::reference DMatrixSlice::operator()(size_type i, size_type j)
{	
  validate_row_subscript(i);
  validate_col_subscript(j);
  return ptr_[(i-1) + (j-1) * max_rows_];
}

inline
DMatrixSlice::const_reference	DMatrixSlice::operator()(size_type i, size_type j) const
{
  validate_row_subscript(i);
  validate_col_subscript(j);
  return ptr_[(i-1) + (j-1) * max_rows_];
}

// Subregion access (validated by constructor for DMatrixSlice)

inline
DVectorSlice  DMatrixSlice::row(size_type i) {
  validate_row_subscript(i);
  return DVectorSlice( ptr_ + (i-1), cols(), max_rows() );
} 

inline
const DVectorSlice DMatrixSlice::row(size_type i) const {
  validate_row_subscript(i);
  return DVectorSlice( const_cast<value_type*>(ptr_) + (i-1), cols(), max_rows() );
} 

inline
DVectorSlice	DMatrixSlice::col(size_type j) {
  validate_col_subscript(j);
  return DVectorSlice( ptr_ + (j-1)*max_rows(), rows(), 1 );
} 

inline
const DVectorSlice DMatrixSlice::col(size_type j) const {
  validate_col_subscript(j);
  return DVectorSlice( const_cast<value_type*>(ptr_) + (j-1)*max_rows(), rows(), 1 );
} 

inline
DVectorSlice DMatrixSlice::diag(difference_type k) {
  return p_diag(k);
}

inline
const DVectorSlice DMatrixSlice::diag(difference_type k) const {
  return p_diag(k);
}

inline
DMatrixSlice DMatrixSlice::operator()(const Range1D& I, const Range1D& J) {
  return DMatrixSlice(*this, RangePack::full_range(I, 1, rows()), RangePack::full_range(J,1,cols()));
}

inline
const DMatrixSlice DMatrixSlice::operator()(const Range1D& I, const Range1D& J) const {
  return DMatrixSlice( const_cast<DMatrixSlice&>(*this)
    , RangePack::full_range(I, 1, rows()), RangePack::full_range(J,1,cols()) );
}

inline
DMatrixSlice DMatrixSlice::operator()(size_type i1, size_type i2, size_type j1
  , size_type j2)
{
  return DMatrixSlice(*this, Range1D(i1,i2), Range1D(j1,j2));
}

inline
const DMatrixSlice DMatrixSlice::operator()(size_type i1, size_type i2, size_type j1
  , size_type j2) const
{
  return DMatrixSlice( const_cast<DMatrixSlice&>(*this), Range1D(i1,i2)
    , Range1D(j1,j2) );
}

inline
DMatrixSlice& DMatrixSlice::operator()() {
  return *this;
}

inline
const DMatrixSlice& DMatrixSlice::operator()() const {
  return *this;
}

// Assignment operators

inline
DMatrixSlice& DMatrixSlice::operator=(value_type alpha) {
  assign(this, alpha);
  return *this;
}

inline
DMatrixSlice& DMatrixSlice::operator=(const DMatrixSlice& rhs) {
  assign(this, rhs, BLAS_Cpp::no_trans);
  return *this;
}

// Raw data access

inline
DMatrixSlice::size_type DMatrixSlice::max_rows() const
{	return max_rows_; }

inline
DMatrixSlice::value_type* DMatrixSlice::col_ptr(size_type j) {
  if( ptr_ )
    validate_col_subscript(j);
  return ptr_ + (j-1) * max_rows();	// will be 0 if not bound to a view.
}

inline
const DMatrixSlice::value_type* DMatrixSlice::col_ptr(size_type j) const {
  if( ptr_ )
    validate_col_subscript(j);
  return ptr_ + (j-1) * max_rows();	// will be 0 if not bound to a view.
}

// ////////////////////////////////////////////////////////////////////////////////////////
// DMatrix inline member function definitions

// Private utilities

#ifndef LINALGPACK_CHECK_RANGE
inline
void DMatrix::validate_row_subscript(size_type i) const
{}
#endif

#ifndef LINALGPACK_CHECK_RANGE
inline
void DMatrix::validate_col_subscript(size_type j) const
{}
#endif

// constructors

inline
DMatrix::DMatrix() : v_(), rows_(0)
{}

inline
DMatrix::DMatrix(size_type rows, size_type cols)
  : v_(rows*cols), rows_(rows)
{}

inline
DMatrix::DMatrix(value_type val, size_type rows, size_type cols)
  : v_(val,rows*cols), rows_(rows)
{}

inline
DMatrix::DMatrix(const value_type* p, size_type rows, size_type cols)
  : v_(rows*cols), rows_(rows)
{
// 6/7/00: valarray<> in libstdc++-2.90.7 has a bug in v_(p,size) so we do not
// use it.  This is a hack until I can find the time to remove valarray all
// together.
  std::copy( p, p + rows*cols, &v_[0] );
}

inline
DMatrix::DMatrix(const DMatrixSlice& gms)
  : v_(gms.rows() * gms.cols()), rows_(gms.rows())
{	
  assign(this, gms, BLAS_Cpp::no_trans);
}

// Memory management

inline
void DMatrix::resize(size_type rows, size_type cols, value_type val)
{
  v_.resize(rows*cols,val);
  v_ = val;
  rows_ = rows;
}

inline
void DMatrix::free() {
  v_.resize(0);
  rows_ = 0;
}

// Size / Dimensionality

inline
DMatrix::size_type DMatrix::rows() const {
  return rows_;
}

inline
DMatrix::size_type DMatrix::cols() const {
  return rows_ > 0 ? v_.size() / rows_ : 0;
}

// Element access

inline
DMatrix::reference DMatrix::operator()(size_type i, size_type j)
{	 
  validate_row_subscript(i); validate_col_subscript(j);
  return v_[(i-1) + (j-1) * rows_];
}

inline
DMatrix::const_reference DMatrix::operator()(size_type i, size_type j) const
{
  validate_row_subscript(i); validate_col_subscript(j);
  return (const_cast<std::valarray<value_type>&>(v_))[(i-1) + (j-1) * rows_];
}

// subregion access (range checked by constructors)

inline
DVectorSlice DMatrix::row(size_type i)
{
  validate_row_subscript(i);
  return DVectorSlice( col_ptr(1) + (i-1), cols(), rows() );
} 

inline
const DVectorSlice DMatrix::row(size_type i) const
{
  validate_row_subscript(i);
  return DVectorSlice( const_cast<value_type*>(col_ptr(1)) + (i-1), cols(), rows() );
} 

inline
DVectorSlice	DMatrix::col(size_type j)
{
  validate_col_subscript(j);
  return DVectorSlice( col_ptr(1) + (j-1) * rows(), rows(), 1 );
} 

inline
const DVectorSlice DMatrix::col(size_type j) const
{
  validate_col_subscript(j);
  return DVectorSlice( const_cast<value_type*>(col_ptr(1)) + (j-1) * rows(), rows(), 1 ) ;
} 

inline
DVectorSlice DMatrix::diag(difference_type k)
{
  return p_diag(k);
}	

inline
const DVectorSlice DMatrix::diag(difference_type k) const
{
  return p_diag(k);
}	

inline
DMatrixSlice DMatrix::operator()(const Range1D& I, const Range1D& J)
{
  Range1D Ix = RangePack::full_range(I,1,rows()), Jx = RangePack::full_range(J,1,cols());
  return DMatrixSlice( col_ptr(1) + (Ix.lbound() - 1) + (Jx.lbound() - 1) * rows()
    , max_rows() * cols(), max_rows(), Ix.size(), Jx.size() );
}

inline
const DMatrixSlice DMatrix::operator()(const Range1D& I, const Range1D& J) const
{
  Range1D Ix = RangePack::full_range(I,1,rows()), Jx = RangePack::full_range(J,1,cols());
  return DMatrixSlice( const_cast<value_type*>(col_ptr(1)) + (Ix.lbound() - 1) + (Jx.lbound() - 1) * rows()
    , max_rows() * cols(), max_rows(), Ix.size(), Jx.size() );
}

inline
DMatrixSlice DMatrix::operator()(size_type i1, size_type i2, size_type j1
  , size_type j2)
{
  return DMatrixSlice( col_ptr(1) + (i1 - 1) + (j1 - 1) * rows()
    , max_rows() * cols(), max_rows(), i2 - i1 + 1, j2 - j1 + 1 );
}

inline
const DMatrixSlice DMatrix::operator()(size_type i1, size_type i2, size_type j1
  , size_type j2) const
{
  return DMatrixSlice( const_cast<value_type*>(col_ptr(1)) + (i1 - 1) + (j1 - 1) * rows()
    , max_rows() * cols(), max_rows(), i2 - i1 + 1, j2 - j1 + 1 );
}

inline
DMatrixSlice DMatrix::operator()()
{
  return DMatrixSlice( col_ptr(1), max_rows() * cols(), max_rows(), rows(), cols() );
}

inline
const DMatrixSlice DMatrix::operator()() const
{
  return DMatrixSlice( const_cast<value_type*>(col_ptr(1)), max_rows() * cols(), max_rows()
    , rows(), cols() );
}

// Implicit conversion operators

inline
DMatrix::operator DMatrixSlice() {
  return (*this)();
}

inline
DMatrix::operator const DMatrixSlice() const
{
  return (*this)();
}

// Assignment operators

inline
DMatrix& DMatrix::operator=(value_type alpha)
{
  assign(this, alpha);
  return *this;
}

inline
DMatrix& DMatrix::operator=(const DMatrix& rhs)
{
  assign(this, rhs, BLAS_Cpp::no_trans);
  return *this;
}

inline
DMatrix& DMatrix::operator=(const DMatrixSlice& rhs)
{
  assign(this, rhs, BLAS_Cpp::no_trans);
  return *this;
}

// Raw data access

inline
DMatrix::size_type DMatrix::max_rows() const
{	return rows_; }

inline
DMatrix::value_type* DMatrix::col_ptr(size_type j)
{
  if( v_.size() ) {
    validate_col_subscript(j);
    return &v_[ (j-1) * max_rows() ];
  }
  else {
    return 0;
  }
}

inline
const DMatrix::value_type* DMatrix::col_ptr(size_type j) const 
{
  if( v_.size() ) {
    validate_col_subscript(j);
    return &const_cast<valarray&>(v_)[ (j-1) * max_rows() ];
  }
  else {
    return 0;
  }
}

}	// end namespace DenseLinAlgPack

#endif	// GEN_MATRIX_CLASS_H
