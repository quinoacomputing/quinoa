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

#ifndef VECTOR_CLASS_TMPL_H
#define VECTOR_CLASS_TMPL_H

#include <vector>

#include "DenseLinAlgPack_Types.hpp"
#include "StrideIterPack_StrideIter.hpp"

namespace DenseLinAlgPack{

// ////////////////////////////////////////////////////////////////////////////////
/* * @name {\bf Dense 1-D DVector Abstractions}.
  *
  * These are classes that abstract 1-D vectors. The class \Ref{DVector} is a storage class
  * for vectors while the class \Ref{DVectorSlice} is used to represent regions of vectors
  * , for rows, columns or diagonals of matrices (see \Ref{DMatrix}
  * and \Ref{DMatrixSlice}).
  */

// @{

/* * @name {\bf DVector Classes}. */
// @{
  
/** \brief . */
/* * Slice of a 1-D sequential C++ array treated as a vector.
  *
  * Objects of this class represent regions of vectors (continuous), rows of matrices
  * , columns of matrices or diagonals of matrices.
  * The underlying representation is of a continuous C++ array with non unit stride.
  * It uses the same convention that the BLAS use where a vector is represented as
  * the first element of
  * in an array, the stride between elements in that array and the number of elements.
  * Changes to elements through a DVectorSlice object result in changes to the elements
  * in the underlying value_type* data.
  *
  * DVectorSlice provides many STL compliant features such as typedef type members
  *, iterator returning functions
  * and the dim() function.  It also provides access to individual elements (lvalue)
  * through 0-based
  * and 1-based subscripting with operator[](i) and operator()(i) respectively.
  *  In addition and subregions can be created
  * by subscripting (operator()()) with an \Ref{Range1D} object or using lower (>0)
  * and upper bounds of the region.
  */
template<class T>
class VectorSliceTmpl {
public:
  
  /* * @name {\bf Nested Member Types (STL)}.
    *
    * These nested types give the types used in the interface to the class.
    *
    * \begin{description}
    *	<li>[#value_type#]				- type being stored in the underlying C++ array			
    *	<li>[#size_type#]				- type used as an index and for the number of elements
    *										in the vector
    *	<li>[#difference_type#]		- type for the distance between elements and the stride
    *	<li>[#iterator#]				- type for the forward non-constant iterator
    *	<li>[#const_iterator#]			- type for the forward constant iterator (can't change elements)
    *	<li>[#reverse_iterator#]		- type for the reverse non-constant iterator
    *	<li>[#const_reverse_iterator#]	- type for the reverse constant iterator (can't change elements)
    *	<li>[#reference#]				- type returned from subscripting, iterator deferencing etc.
    *	<li>[#const_reference#]		- "" "" for const vector slice objects
    * \end{description}
    */

  // @{
  // @}

  typedef T												value_type;
  typedef DenseLinAlgPack::size_type							size_type;
  typedef ptrdiff_t										difference_type;
  typedef StrideIterPack::stride_iter<value_type*
    , value_type, value_type&, value_type*
    , difference_type>									iterator;
  typedef StrideIterPack::stride_iter<const value_type*
    , value_type, const value_type&, const value_type*
    , difference_type>									const_iterator;
#if defined(_INTEL_CXX) || defined (_INTEL_CXX)
  typedef std::reverse_iterator<iterator, value_type
    , value_type&, value_type*, difference_type>		reverse_iterator;
  typedef std::reverse_iterator<const_iterator
    , value_type, const value_type&, const value_type*
    , difference_type>									const_reverse_iterator;
#else
  typedef std::reverse_iterator<iterator>					reverse_iterator;
  typedef std::reverse_iterator<const_iterator>			const_reverse_iterator;
#endif
  typedef value_type&										reference;
  typedef const value_type&								const_reference;

  /* * @name {\bf Constructors}.
    *
    * The user usually does not need not call any of these constructors
    * explicitly to create a vector slice.
    * These
    * constructors are used by the classes in the library to construct VectorSliceTmpl objects.
    * Instead, users create VectorSliceTmpl objects by indexing (\Ref{Range1D}) a \Ref{DVector}
    * , or \Ref{VectorSliceTmpl}
    * object or calling row(...), col(...) or diag(...) on a \Ref{DMatrix} or
    * \Ref{DMatrixSlice} object.
    * The default C++ copy constructor is used, and is therefore not show here.
    *
    * Constructors are also included for creating views of raw C++ arrays.
    * These constructors take non-#const# pointers.  However you can savely
    * create a #const# view of a #const# C++ array by using a constant cast.
    * For example:
    *
    \verbatim
    const VectorSliceTmpl<T>::size_type n = 5;
    const VectorSliceTmpl<T>::value_type ptr[n] = { 1.0, 2.0, 3.0, 4.0, 5,0 };
    const VectorSliceTmpl vec(const_cast<VectorSliceTmpl<T>::value_type*>(ptr),n);
    \endverbatim
    */

  // @{

  /** \brief . */
  /* * Creates an empty view.
    *
    * You must use bind(...) to bind to a view to initialize after construction.
    */
  VectorSliceTmpl();
  /** \brief . */
  /* * Creates a VectorSice object that represents a non-continous region of a raw C++ array.
    *
    * Of course the sequence of elements #ptr[stride * i]# for #i# = 0, 1, ..., #size#-1
    * must yield valid properly allocated regions of memory.  It is up the the user to insure
    * that they do.
    *
    *	@param	ptr		Pointer to the first element in the raw C++ array
    * @param	size	number of elements in the vector slice
    * @param	stride	the distance (may be negative) between each successive element (default = 1)
    */
  VectorSliceTmpl( value_type* ptr, size_type size, difference_type stride = 1 );
  /** \brief . */
  /* * Creates a VectorSliceTmpl object that represents a continous region of a raw C++ array.
    *
    * The VectorSliceTmpl Object represents the following elements of raw array:
    * 
    *      #ptr[rng.lbound()-1+i]#, for #i# = 0, 1, ..., #rng.ubound()-1#
    * 
    * Preconditions: <ul>
    *		<li> #rng.ubound() + 1 <= n# (throw std::out_of_range)
    *		</ul>
    *
    *	@param	ptr		Pointer to the first element in the raw C++ array
    * @param	size	number of elements in the vector slice
    * @param	rng		index range (1-based) of the region being represented.  
    *					Here rng.full_range() can be true.
    */
  VectorSliceTmpl( value_type* ptr, size_type size, const Range1D& rng );
  /// 
  /* * Create a VectorSliceTmpl that represents a continous region of the existing VectorSliceTmpl object, vs.
    *
    * The index, rng, is relative to the VectorSliceTmpl object, vs.
    * For example rng = [1,3] would create a VectorSliceTmpl object
    * representing the elements 2, 4 and 6.  The following
    * shows the elements represented by each of the participating objects.
    \verbatim

    vs =   [2, 4, 6, 8, 10]
    this = [2, 4, 6]

    \endverbatim

    * Preconditions: <ul>
    *		<li> rng.full_range() == false (throw #std::out_of_range#)
    *		<li> rng.dim() <= vs.dim() (throw #std::out_of_range#) 
    *		</ul>
    * 
    * @param	vs		VectorSliceTmpl object that this VectorSliceTmpl object is being created from
    * @param	rng		Range1D range of the vector slice being created.
    */
  VectorSliceTmpl( VectorSliceTmpl<value_type>& vs, const Range1D& rng );

  // @}

  /// Bind to the view of another VectorSliceTmpl
  void bind(VectorSliceTmpl<value_type> vs);

  /* * @name {\bf STL Iterator Access Functions}.
    *
    * These member functions return valid STL random access iterators to the elements in the
    * VectorSliceTmpl object.
    *
    * The forward iterators returned by begin() and end() iterator sequentialy from the first
    * element (same element as returned by operator()(1)) to the last
    * element (same element as returned by operator()(dim()).  This goes for reverse 
    * (stride() < 0) VectorSliceTmpl objects as well.  The reverse iterators returned by
    * rbegin() and rend() iterate in the reverse sequence.
    *
    * Warning! Beware of using iterations in a reverse vector slice (stride() < 0). 
    * In a reverse vector slice end() returns a slice iterator which is the current
    * element is one before the first allocated element.  Strictly speaking this is
    * not allowed so use iterators with reversed VectorSliceTmpl objects at own risk.
    */

  // @{

  /// 
  iterator begin();
  /** \brief . */
  iterator end();
  /** \brief . */
  const_iterator begin() const;
  /** \brief . */
  const_iterator end() const;
  /** \brief . */
  reverse_iterator rbegin();
  /** \brief . */
  reverse_iterator rend();
  /** \brief . */
  const_reverse_iterator rbegin() const;
  /** \brief . */
  const_reverse_iterator rend() const;

  // @}

  /* * @name {\bf Individual Element Access Subscripting (lvalue)}.
    *
    * These operator functions allow access (lvalue) to the individual elements
    * of the VectorSliceTmpl object.
    * 
    * The subscript i must be, 1 <= i <= this->dim(), for the 1-based element access
    * operators and, 0 <= i <= this->dim() - 1, for the 0-based element access operators.
    * If they are not then an #std::out_of_range# exception will be thrown.
    */

  // @{

  /// 1-based element access (lvalue)
  reference operator()(size_type i);
  /// 1-based element access (rvalue)
  const_reference operator()(size_type i) const;
  /// 1-based element access (lvalue)
  reference operator[](size_type i);
  /// 0-based element access (rvalue)
  const_reference operator[](size_type i) const;

  // @}

  /* * @name {\bf Subvector Access Operators}.
    *
    * These operator functions are used to create views of continous regions of the VectorSliceTmpl.
    * Each of them returns a VectorSliceTmpl object for the region.  Constant (const) VectorSliceTmpl objects
    * are returned for a const VectorSliceTmpl.  This means that the elements can not be changed
    * as should be the case.
    * 
    * Beware!  VC++ is returning non-const VectorSliceTmpl objects for the 
    * #VectorSliceTmpl operator()(...) const;# member functions and therefore a const \Ref{DVector} or
    * \Ref{VectorSliceTmpl} can be modifed my subsetting it.  Hopefully this problem will
    * be fixed in future versions of the compiler or I when will get another compiler.
    */

  // @{

  /// 
  /* * Returns a VectorSliceTmpl object representing the entire vector slice.
    *
    * Included for uniformity with vector.
    */
  /// Allow the address to be taken of an rvalue of this object.
  VectorSliceTmpl<value_type>* operator&() {
    return this;
  }
  /** \brief . */
  const VectorSliceTmpl<value_type>* operator&() const {
    return this;
  }
  VectorSliceTmpl<value_type>& operator()();
  /// Same as above.
  const VectorSliceTmpl<value_type>& operator()() const;
  /// 
  VectorSliceTmpl<value_type> operator()(const Range1D& rng);
  /// 
  /* * Returns a continous subregion of the VectorSliceTmpl object.
    *
    * The returned VectorSliceTmpl object represents the range of the rng argument.
    *
    * Preconditions: <ul>
    *		<li> #rng.ubound() - 1 <= this->dim()# (throw #out_of_range#)
    *		</ul>
    *
    * @param	rng		Indece range [lbound,ubound] of the region being returned.
    */
  const VectorSliceTmpl<value_type> operator()(const Range1D& rng) const;
  /// 
  /* * Returns a VectorSliceTmpl object for the continous subregion [ubound, lbound].
    * 
    * Preconditions: <ul>
    *		<li> #lbound > 1# (throw out_of_range)
    *		<li> #lbound < ubound# (throw out_of_range)
    *		<li> #ubound <= this->dim()# (throw out_of_range)
    *		</ul>
    *
    * @param	rng		Range [lbound,ubound] of the region being returned.
    */
  VectorSliceTmpl<value_type> operator()(size_type lbound, size_type ubound);
  /// Same as above.
  const VectorSliceTmpl<value_type> operator()(size_type lbound, size_type ubound) const;
  /// 
  /* * Return a const VectorSliceTmpl object the reverse of this VectorSliceTmpl.
    *
    * In the reverse VectorSliceTmpl,
    * the first element becomes the last element and visa-versa.  For example, for 
    * #VectorSliceTmpl r = x.rev()#, #&x(1) == &z(z.dim())# and #&x(x.dim()) == &z(1)# are both true.
    * The iterators returned by \Ref{begin()} iterate from the first conceptual element to the last.
    */
  VectorSliceTmpl<value_type> rev();
  /// Same as above.
  const VectorSliceTmpl<value_type> rev() const;

  // @}
  
  /* * @name {\bf Assignment operators}. */

  // @{

  /** \brief . */
  /* * vs = alpha (Sets all the elements to the constant alpha).
    *
    * Preconditions: <ul>
    *		<li> #this->dim() > 0# (throw #std::length_error#)
    *		</ul>
    *
    * Postconditions: <ul>
    *		<li> #this->operator()(i) == alpha#, i = 1, 2, ... , #this->dim()#
    *		</ul>
    */
  VectorSliceTmpl<value_type>& operator=(value_type alpha);
  /** \brief . */
  /* * vs = rhs (Copies the elements of rhs into the elements of this).
    *
    * Preconditions: <ul>
    *		<li> #this->dim() == rhs.dim()# (throw #out_of_range#)
    *		<li> #rhs.dim() > 0# (throw #out_of_range#)
    *		</ul>
    *
    * Postconditions: <ul>
    *		<li> #this->operator()(i) == rhs(i)#, i = 1, 2, ..., #this->dim()#
    *		</ul>
    */
  VectorSliceTmpl<value_type>& operator=(const VectorSliceTmpl<value_type>& rhs);

  // @}

  /* * @name {\bf Misc. Member Functions}. */

  // @{

  /// Returns the number of elements of the VectorSliceTmpl.
  size_type dim() const;
  /// 
  /* * Returns the degree of memory overlap of the two VectorSliceTmpl objects this and vs.
    *
    * @return 
    *		\begin{description}
    *		<li>[NO_OVERLAP]	There is no memory overlap between this and vs
    *		<li>[SOME_OVERLAP]	There is some memory locations that this and vs share
    *		<li>[SAME_MEM]		The VectorSliceTmpl objects this and vs share the exact same memory locations.
    *		\end{description}
    */
  EOverLap overlap(const VectorSliceTmpl<value_type>& vs) const;

  // @}

  /* * @name {\bf Raw data access}.
    *
    * Provides access to underlying raw data.
    */

  // @{

  /// Return a pointer to the address of the first memory location of underlying array.
  value_type*			raw_ptr();
  /** \brief . */
  const value_type*	raw_ptr() const;
  /// Return a pointer to the conceptual first element in the underlying array.
  value_type*			start_ptr();
  /** \brief . */
  const value_type*	start_ptr() const;
  /// Return the distance (+,-) (in units of elements) between adjacent elements in the underlying array.
  difference_type		stride() const;

  // @}

private:
  value_type					*ptr_;	// Pointer to first element in array.
  size_type					size_;  // # elements represented in v_
  difference_type				stride_;// # positions to skip between elements. Must be positive
  // Above the sequence represented is:
  //  ptr_[ i * stride_ ], for i = 0, ..., size_ - 1

}; // end class VectorSliceTmpl<T>

// /////////////////////////////////////////////////////////////////////////////////////////
// DVector
//

/** \brief . */
/* * 1-D DVector Abstraction Storage Class.
  *
  * Holds the storage space for a 1-D vector of element type value_type.  The storage space class
  * used in a standard vector<> private member.  DVector provides much of the
  * same functionaliy of a VectorSliceTmpl object accept that DVector object can be resized at any time by
  * either explicitly calling #resize(...)# or to match an assignment to a rhs linear algebra expression.
  */
template<class T>
class VectorTmpl {
public:
  /* * @name {\bf Nested Member Types (STL)}.
    *
    * These nested types give the types used in the interface to the class.
    *
    * \begin{description}
    *	<li>[#value_type#]				- type being stored in the underlying vector<>			
    *	<li>[#size_type#]				- type for the number of elements in the vector<>
    *	<li>[#difference_type#]		- type for the distance between elements
    *	<li>[#iterator#]				- type for the forward non-constant iterator
    *	<li>[#const_iterator#]			- type for the forward constant iterator (can't change elements)
    *	<li>[#reverse_iterator#]		- type for the reverse non-constant iterator
    *	<li>[#const_reverse_iterator#]	- type for the reverse constant iterator (can't change elements)
    *	<li>[#reference#]				- type returned from subscripting, iterator deferencing etc.
    *	<li>[#const_reference#]		- "" "" for const vector slice objects
    * \end{description}
    */

  // @{
  // @}

  typedef T										value_type;
  typedef DenseLinAlgPack::size_type					size_type;
  typedef ptrdiff_t								difference_type;
  typedef value_type*								iterator;
  typedef const value_type*						const_iterator;
#if 0 /* defined(_INTEL_CXX) || defined(_WINDOWS) */
  typedef std::reverse_iterator<iterator, value_type
    , value_type&, value_type*, difference_type>		reverse_iterator;
  typedef std::reverse_iterator<const_iterator
    , value_type, const value_type&, const value_type*
    , difference_type>									const_reverse_iterator;
#else
  typedef std::reverse_iterator<iterator>					reverse_iterator;
  typedef std::reverse_iterator<const_iterator>			const_reverse_iterator;
#endif
  typedef value_type&								reference;
  typedef const value_type&						const_reference;
  typedef std::vector<value_type>					valarray;

  /* * @name {\bf Constructors}.
    * 
    * These constructors allocate and may initialize the elements of a 1-D vector.
    * The default C++ copy constructor is used and is therefore not show here.
    */

  // @{

  /// Constructs a vector with 0 elements (this->dim()==0).
  VectorTmpl();
  /// Constructs a vector with n elements of initialized memory.
  VectorTmpl(size_type n);
  /// Constructs a vector with n elements initialized to val.
  VectorTmpl(value_type val, size_type n);
  /** \brief . */
  /* * Constructs a vector with n elements and intializes elements to those of an array.
    *
    * Postconditions: <ul>
    *		<li> #this->operator[](i) == p[i]#, i = 0, 1, ... n
    *		</ul>
    */  
  VectorTmpl(const value_type* p, size_type n);
  /** \brief . */
  /* * Constructs a DVector object fron a VectorSliceTmpl object.
    *
    * Postconditions: <ul>
    *		<li> #this->dim() == vs.dim()#
    *		<li> #this->operator[](i) == vs[i]#, i = 0, 1, ... n
    *		</ul>
    */  
  VectorTmpl(const VectorSliceTmpl<value_type>& vs);

  // @}

  /* * @name {\bf Memory Management / Misc}. */

  // @{

  /** \brief . */
  /* * Resize the vector to hold n elements.
    *
    * Any new elements added are initialized to val.
    *
    * Postconditions: <ul>
    *		<li> #this->dim() == n#
    *		</ul>
    */  
  void resize(size_type n, value_type val = value_type());
  /** \brief . */
  /* * Free memory and resize DVector to this->dim() == 0.
    *
    * Postconditions: <ul>
    *		<li> #this->dim() == 0#
    *		</ul>
    */  
  void free();
  /// Returns the number of elements of the DVector.
  size_type dim() const;
  /// 
  /* * Returns the degree of memory overlap of this and the VectorSliceTmpl object vs.
    *
    * @return 
    *		\begin{description}
    *		<li>[NO_OVERLAP]	There is no memory overlap between this and vs
    *		<li>[SOME_OVERLAP]	There is some memory locations that this and vs share
    *		<li>[SAME_MEM]		The VectorSliceTmpl objects this and vs share the exact same memory locations.
    *		\end{description}
    */
  EOverLap overlap(const VectorSliceTmpl<value_type>& vs) const;
  /// Conversion operator for implicit conversions from DVector to VectorSliceTmpl.
  operator VectorSliceTmpl<value_type>();
  /// Conversion operator for implicit conversions from const DVector to const VectorSliceTmpl.
  operator const VectorSliceTmpl<value_type>() const;

  // @}
  

  /* * @name {\bf STL Iterator Access functions}.
    *
    * The iterators returned are valid STL random access iterators.
    * The forward iterators returned iterate from the first element to the last element.
    * The reverse iterators returned iterate from the last element to the first element.
    */

  // @{

  /** \brief . */
  iterator begin();
  /** \brief . */
  iterator end();
  /** \brief . */
  const_iterator begin() const;
  /** \brief . */
  const_iterator end() const;
  /** \brief . */
  reverse_iterator rbegin();
  /** \brief . */
  reverse_iterator rend();
  /** \brief . */
  const_reverse_iterator rbegin() const;
  /** \brief . */
  const_reverse_iterator rend() const;

  // @}

  /* * @name {\bf Individual Element Access Subscripting (lvalue)}.
    *
    * These operator functions allow access (lvalue) to the individual elements
    * of the DVector object.
    * 
    * The subscript i must be, 1 <= i <= this->dim(), for the 1-based element access
    * operators and, 0 <= i <= this->dim() - 1, for the 0-based element access operators.
    * If they are not then an #std::out_of_range# exception will be thrown.
    */

  // @{

  /// 1-based element access (lvalue)
  reference operator()(size_type i);
  /// 1-based element access (rvalue)
  const_reference operator()(size_type i) const;
  /// 1-based element access (lvalue)
  reference operator[](size_type i);
  /// 0-based element access (rvalue)
  const_reference operator[](size_type i) const;

  // @}

  /* * @name {\bf Subvector Access Operators}.
    *
    * These operator functions are used to create views of continous regions of the DVector.
    * Each of them returns a VectorSliceTmpl object for the region.  Constant (const) VectorSliceTmpl objects
    * are returned for a const DVector.  This means that the elements can not be changed
    * as should be the case.
    * 
    * Beware!  VC++ is returning non-const VectorSliceTmpl objects for the 
    * #VectorSliceTmpl operator()(...) const;# member functions and therefore a const \Ref{DVector} or
    * \Ref{VectorSliceTmpl} can be modifed my subsetting it.  Hopefully this problem will
    * be fixed in future versions of the compiler or I when will get another compiler.
    */

  // @{

  /// 
  /* * Returns a VectorSliceTmpl object representing the entire DVector. 
    * 
    * Call this member function to force a type conversion to VectorSliceTmpl.  Using the
    * VectorSliceTmpl of a DVector for algebraic expressions used with the TCOL allows a for simplier
    * implementaion of those operations by cutting down on the number combinations.  This is
    * especialy true for longer optimized expression.
    */
  VectorSliceTmpl<value_type> operator()();
  /// Same as above
  const VectorSliceTmpl<value_type> operator()() const;
  /// 
  /* * Returns a continous subregion of the DVector object.
    *
    * The returned VectorSliceTmpl object represents the range of the rng argument.
    *
    * Preconditions: <ul>
    *		<li> #rng.ubound() - 1 <= this->dim()# (throw #out_of_range#)
    *		</ul>
    *
    * @param	rng		Indece range [lbound,ubound] of the region being returned.
    */
  VectorSliceTmpl<value_type> operator()(const Range1D& rng);
  /// Same as above
  const VectorSliceTmpl<value_type> operator()(const Range1D& rng) const;
  /// 
  /* * Returns a VectorSliceTmpl object for the continous subregion [ubound, lbound].
    * 
    * Preconditions: <ul>
    *		<li> #lbound > 1# (throw #out_of_range#)
    *		<li> #lbound < ubound# (throw #out_of_range#)
    *		<li> #ubound <= this->dim()# (throw #out_of_range#)
    *		</ul>
    *
    * @param	rng		Range [lbound,ubound] of the region being taken.
    */
  VectorSliceTmpl<value_type> operator()(size_type lbound, size_type ubound);
  /// Same as above.
  const VectorSliceTmpl<value_type> operator()(size_type lbound, size_type ubound) const;
  /// 
  /* * Return a VectorSliceTmpl object the reverse of this DVector.
    *
    * In the reverse VectorSliceTmpl,
    * the first element becomes the last element and visa-versa.  For example, for 
    * #VectorSliceTmpl r = x.rev()#, #&x(1) == &z(z.dim())# and #&x(x.dim()) == &z(1)# are both true.
    * The iterators returned by \Ref{begin()} iterate from the first conceptual element to the last.
    */
  VectorSliceTmpl<value_type> rev();
  /// Same as above.
  const VectorSliceTmpl<value_type> rev() const;

  // @}

  /* * @name {\bf Assignment Operators}. */

  // @{

  /** \brief . */
  /* * vs = alpha (Sets all the elements to the constant alpha).
    *
    * Preconditions: <ul>
    *		<li> #this->dim() > 0# (throw #std::length_error#)
    *		</ul>
    *
    * Postconditions: <ul>
    *		<li> #this->operator()(i) == alpha#, i = 1, 2, ... , #this->dim()#
    *		</ul>
    */
  VectorTmpl<value_type>& operator=(value_type alpha);
  /** \brief . */
  /* * vs = rhs (Copies the elements of rhs into the elements of this).
    *
    * Preconditions: <ul>
    *		<li> #this->dim() == rhs.dim()# (throw #out_of_range#)
    *		<li> #rhs.dim() > 0# (throw #out_of_range#)
    *		</ul>
    *
    * Postconditions: <ul>
    *		<li> #this->operator()(i) == rhs(i)#, i = 1, 2, ..., #this->dim()#
    *		</ul>
    */
  VectorTmpl<value_type>& operator=(const VectorSliceTmpl<value_type>& rhs);
  /** \brief . */
  /* * Needed to override the default assignment operator.
    */
  VectorTmpl<value_type>& operator=(const VectorTmpl<value_type>& rhs);

  // @}

  /* * @name {\bf Implementation Access}.
    *
    * Provides access to underlying raw data.
    */

  // @{

  /// Return a pointer to the address of the first memory location of underlying array.
  value_type*			raw_ptr();
  /** \brief . */
  const value_type*	raw_ptr() const;
  /// Return a pointer to the conceptual first element in the underlying array.
  value_type*			start_ptr();
  /** \brief . */
  const value_type*	start_ptr() const;
  /// Return the distance (+,-) (in units of elements) between adjacent elements in the underlying array.
  difference_type		stride() const;

  // @}
  
private:
  valarray v_;

}; // end class VectorTmpl<T>

//		end DVector Classes scope
// @}

// ///////////////////////////////////////////////////////////////////////////////
// Non-member function declarations												//
// ///////////////////////////////////////////////////////////////////////////////


/* * @name {\bf Non-Member Functions}. */

// @{
//		begin non-member functions scope

//
size_type vector_validate_sized(size_type size);
//
void vector_validate_range(size_type ubound, size_type max_ubound);
//
void vector_validate_subscript(size_type size, size_type i);
/** \brief . */
/* * Utility for checking the sizes of two VectorSliceTmpl objects and throwing an exception
  * if the sizes are not the same.
  */
void assert_vs_sizes(size_type size1, size_type size2);

/** \brief . */
/* * Create a general vector slice.
  */
template<class T>
inline
VectorSliceTmpl<T> gen_vs( VectorSliceTmpl<T>& vs, size_type start, size_type size, ptrdiff_t stride )
{
  return VectorSliceTmpl<T>( vs.start_ptr() + vs.stride() * (start-1), size, vs.stride() * stride );
}

/** \brief . */
template<class T>
inline
const VectorSliceTmpl<T> gen_vs( const VectorSliceTmpl<T>& vs, size_type start, size_type size
  , ptrdiff_t stride )
{
  return VectorSliceTmpl<T>( const_cast<typename VectorSliceTmpl<T>::value_type*>(vs.start_ptr()) + vs.stride() * (start-1)
    , size, vs.stride() * stride );
}

//		end non-member functions scope
// @}

//		end Vectors scope
// @}

} // end namespace DenseLinAlgPack

// ////////////////////////////////////////////////////////////////////////////////
// Inline definitions of member function definitions							 //
// ////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////
// Non-member functions / utilities

#ifndef LINALGPACK_CHECK_SLICE_SETUP
inline
DenseLinAlgPack::size_type DenseLinAlgPack::vector_validate_sized(size_type size)
{
  return size;
}
#endif

#ifndef LINALGPACK_CHECK_RANGE
inline
void DenseLinAlgPack::vector_validate_range(size_type ubound, size_type max_ubound)
{}
#endif

#ifndef LINALGPACK_CHECK_RANGE
inline
void DenseLinAlgPack::vector_validate_subscript(size_type size, size_type i)
{}
#endif

#ifndef LINALGPACK_CHECK_RHS_SIZES
inline
void DenseLinAlgPack::assert_vs_sizes(size_type size1, size_type size2)
{}
#endif

namespace DenseLinAlgPack {

// /////////////////////////////////////////////////////////////////////////////
// VectorSliceTmpl inline member function definitions

// Constructors.  Use default copy constructor

template<class T>
inline
VectorSliceTmpl<T>::VectorSliceTmpl()
  : ptr_(0)
  , size_(0)
  , stride_(0)
{}

template<class T>
inline
VectorSliceTmpl<T>::VectorSliceTmpl( value_type* ptr, size_type size, difference_type stride)
  : ptr_(ptr)
  , size_(size)
  , stride_(stride)
{}

template<class T>
inline
VectorSliceTmpl<T>::VectorSliceTmpl( value_type* ptr, size_type size, const Range1D& rng )
  : ptr_( ptr + rng.lbound() - 1 )
  , size_( rng.full_range() ?	vector_validate_sized(size) : rng.size() )
  , stride_(1)
{
  vector_validate_range( rng.full_range() ? size : rng.ubound(), size );
}

template<class T>
inline
VectorSliceTmpl<T>::VectorSliceTmpl( VectorSliceTmpl<T>& vs, const Range1D& rng )
  : ptr_( vs.start_ptr() + (rng.lbound() - 1) * vs.stride() )
  , size_( rng.full_range() ?	vector_validate_sized(vs.dim()) : rng.size() )
  , stride_( vs.stride() )
{	vector_validate_range(  rng.full_range() ? vs.dim() : rng.ubound(), vs.dim() ); }

template<class T>
inline
void VectorSliceTmpl<T>::bind(VectorSliceTmpl vs)
{
  ptr_	= vs.ptr_;
  size_	= vs.size_;
  stride_	= vs.stride_;
}

// Iterator functions
template<class T>
inline
typename VectorSliceTmpl<T>::iterator	VectorSliceTmpl<T>::begin()
{	return iterator(start_ptr(), stride()); }

template<class T>
inline
typename VectorSliceTmpl<T>::iterator	VectorSliceTmpl<T>::end()
{	return iterator(start_ptr() + dim() * stride(), stride()); }

template<class T>
inline
typename VectorSliceTmpl<T>::const_iterator VectorSliceTmpl<T>::begin() const
{	return const_iterator(start_ptr(), stride()); }

template<class T>
inline
typename VectorSliceTmpl<T>::const_iterator VectorSliceTmpl<T>::end() const
{	return const_iterator(start_ptr() + dim() * stride(), stride()); }

template<class T>
inline
typename VectorSliceTmpl<T>::reverse_iterator	VectorSliceTmpl<T>::rbegin()
{	return reverse_iterator(end()); }

template<class T>
inline
typename VectorSliceTmpl<T>::reverse_iterator	VectorSliceTmpl<T>::rend()
{	return reverse_iterator(begin()); }

template<class T>
inline
typename VectorSliceTmpl<T>::const_reverse_iterator VectorSliceTmpl<T>::rbegin() const
{	return const_reverse_iterator(end()); }

template<class T>
inline
typename VectorSliceTmpl<T>::const_reverse_iterator VectorSliceTmpl<T>::rend() const
{	return const_reverse_iterator(begin()); }

// Element access
template<class T>
inline
typename VectorSliceTmpl<T>::reference VectorSliceTmpl<T>::operator()(size_type i) // 1 based
{
  vector_validate_subscript(dim(),i);
  return ptr_[(i-1)*stride_];
}

template<class T>
inline
typename VectorSliceTmpl<T>::const_reference VectorSliceTmpl<T>::operator()(size_type i) const
{
  vector_validate_subscript(dim(),i);
  return ptr_[(i-1)*stride_];
}

template<class T>
inline
typename VectorSliceTmpl<T>::reference VectorSliceTmpl<T>::operator[](size_type i) // 0 based		
{
  vector_validate_subscript(dim(),i+1);
  return ptr_[(i)*stride_];
}

template<class T>
inline
typename VectorSliceTmpl<T>::const_reference VectorSliceTmpl<T>::operator[](size_type i) const
{
  vector_validate_subscript(dim(),i+1);
  return ptr_[(i)*stride_];
}

// Subregion Access.  Let the constructors of VectorSliceTmpl validate the ranges
template<class T>
inline
VectorSliceTmpl<T>& VectorSliceTmpl<T>::operator()() 
{	return *this; }

template<class T>
inline
const VectorSliceTmpl<T>& VectorSliceTmpl<T>::operator()() const 
{	return *this; }

template<class T>
inline
VectorSliceTmpl<T> VectorSliceTmpl<T>::operator()(const Range1D& rng) 
{	return VectorSliceTmpl(*this, RangePack::full_range(rng,1,dim())); }

template<class T>
inline
const VectorSliceTmpl<T> VectorSliceTmpl<T>::operator()(const Range1D& rng) const
{	return VectorSliceTmpl(const_cast<VectorSliceTmpl<T>&>(*this), RangePack::full_range(rng,1,dim())); }

template<class T>
inline
VectorSliceTmpl<T> VectorSliceTmpl<T>::operator()(size_type lbound, size_type ubound)
{	return VectorSliceTmpl(*this, Range1D(lbound, ubound)); }

template<class T>
inline
const VectorSliceTmpl<T> VectorSliceTmpl<T>::operator()(size_type lbound, size_type ubound) const
{	return VectorSliceTmpl(const_cast<VectorSliceTmpl<T>&>(*this), Range1D(lbound, ubound)); }

template<class T>
inline
VectorSliceTmpl<T> VectorSliceTmpl<T>::rev()
{	return VectorSliceTmpl( start_ptr() + stride() * (dim()-1), dim(), - stride() ); }

template<class T>
inline
const VectorSliceTmpl<T> VectorSliceTmpl<T>::rev() const
{	return VectorSliceTmpl( const_cast<value_type*>(start_ptr()) + stride() * (dim()-1), dim(), - stride() ); }

// Assignment Operators
template<class T>
inline
VectorSliceTmpl<T>& VectorSliceTmpl<T>::operator=(value_type alpha)
{
  std::fill(begin(),end(),alpha);
  return *this;
}

template<class T>
inline
VectorSliceTmpl<T>& VectorSliceTmpl<T>::operator=(const VectorSliceTmpl<T>& rhs) 
{
  assert_vs_sizes(this->dim(),rhs.dim());
  std::copy(rhs.begin(),rhs.end(),begin());
  return *this;
}

// Misc. member functions

template<class T>
inline
typename VectorSliceTmpl<T>::size_type VectorSliceTmpl<T>::dim() const
{	return size_; }

// Raw pointer access

template<class T>
inline
typename VectorSliceTmpl<T>::value_type*	VectorSliceTmpl<T>::raw_ptr()
{	return stride() > 0 ? start_ptr() : start_ptr() + stride() * (dim() - 1); }

template<class T>
inline
const typename VectorSliceTmpl<T>::value_type* VectorSliceTmpl<T>::raw_ptr() const
{	return stride() > 0 ? start_ptr() : start_ptr() + stride() * (dim() - 1); }

template<class T>
inline
typename VectorSliceTmpl<T>::value_type*	VectorSliceTmpl<T>::start_ptr()
{	return ptr_; }

template<class T>
inline
const typename VectorSliceTmpl<T>::value_type* VectorSliceTmpl<T>::start_ptr() const
{	return ptr_; }

template<class T>
inline
typename VectorSliceTmpl<T>::difference_type VectorSliceTmpl<T>::stride() const
{	return stride_; }

// /////////////////////////////////////////////////////////////////////////////
// DVector inline member function definitions

// Constructors
template<class T>
inline
VectorTmpl<T>::VectorTmpl()
{}	// used to shut satisfy compiler

template<class T>
inline
VectorTmpl<T>::VectorTmpl(size_type n) 
  : v_(n)
{}

template<class T>
inline
VectorTmpl<T>::VectorTmpl(value_type val, size_type n) 
  : v_(n)
{
  std::fill(begin(),end(),val);
}

template<class T>
inline
VectorTmpl<T>::VectorTmpl(const value_type* p, size_type n)
  : v_(n)
{
  std::copy(p,p+n,begin());
}

template<class T>
inline
VectorTmpl<T>::VectorTmpl(const VectorSliceTmpl<T>& vs)
  : v_(vs.dim())
{  
  std::copy(vs.begin(),vs.end(),begin());
}

// Memory management
template<class T>
inline
void VectorTmpl<T>::resize(size_type n, value_type val)
{
  v_.resize(n);
  std::fill(begin(),end(),val);
 }

template<class T>
inline
void VectorTmpl<T>::free()
{
  v_.resize(0);
}

// Size
template<class T>
inline
typename VectorTmpl<T>::size_type VectorTmpl<T>::dim() const
{	return v_.size(); }

// Iterator functions
template<class T>
inline
typename VectorTmpl<T>::iterator VectorTmpl<T>::begin()
{	return start_ptr(); }

template<class T>
inline
typename VectorTmpl<T>::iterator VectorTmpl<T>::end()
{	return start_ptr() + dim(); }

template<class T>
inline
typename VectorTmpl<T>::const_iterator VectorTmpl<T>::begin() const
{	return start_ptr(); }

template<class T>
inline
typename VectorTmpl<T>::const_iterator VectorTmpl<T>::end() const 
{	return start_ptr() + dim(); }

template<class T>
inline
typename VectorTmpl<T>::reverse_iterator	VectorTmpl<T>::rbegin()
{	return reverse_iterator(end()); }

template<class T>
inline
typename VectorTmpl<T>::reverse_iterator	VectorTmpl<T>::rend()
{	return reverse_iterator(begin()); }

template<class T>
inline
typename VectorTmpl<T>::const_reverse_iterator VectorTmpl<T>::rbegin() const
{	return const_reverse_iterator(end()); }

template<class T>
inline
typename VectorTmpl<T>::const_reverse_iterator VectorTmpl<T>::rend() const 
{	return const_reverse_iterator(begin()); }

// Element access
template<class T>
inline
typename VectorTmpl<T>::reference VectorTmpl<T>::operator()(size_type i)
{
  vector_validate_subscript(dim(),i);
  return start_ptr()[i-1];
}

template<class T>
inline
typename VectorTmpl<T>::const_reference VectorTmpl<T>::operator()(size_type i) const
{
  vector_validate_subscript(dim(),i);
  return start_ptr()[i-1];
}

template<class T>
inline
typename VectorTmpl<T>::reference VectorTmpl<T>::operator[](size_type i)
{
  vector_validate_subscript(dim(),i+1);
  return start_ptr()[i];
}

template<class T>
inline
typename VectorTmpl<T>::const_reference VectorTmpl<T>::operator[](size_type i) const
{
  vector_validate_subscript(dim(),i+1);
  return start_ptr()[i];
}

// Subregion Access.  Leave validation to the VectorSliceTmpl constructors.
template<class T>
inline
VectorSliceTmpl<T> VectorTmpl<T>::operator()() 
{	return VectorSliceTmpl<T>(start_ptr(),dim()); }

template<class T>
inline
const VectorSliceTmpl<T> VectorTmpl<T>::operator()() const
{	return VectorSliceTmpl<T>(const_cast<value_type*>(start_ptr()),dim()); }

template<class T>
inline
VectorSliceTmpl<T> VectorTmpl<T>::operator()(const Range1D& rng) 
{	return VectorSliceTmpl<T>(start_ptr(),dim(),rng); }

template<class T>
inline
const VectorSliceTmpl<T> VectorTmpl<T>::operator()(const Range1D& rng) const
{	return VectorSliceTmpl<T>(const_cast<value_type*>(start_ptr()),dim(),rng); }

template<class T>
inline
VectorSliceTmpl<T> VectorTmpl<T>::operator()(size_type lbound, size_type ubound)
{	return VectorSliceTmpl<T>(start_ptr(), dim(), Range1D(lbound, ubound)); }

template<class T>
inline
const VectorSliceTmpl<T> VectorTmpl<T>::operator()(size_type lbound, size_type ubound) const
{	return VectorSliceTmpl<T>(const_cast<value_type*>(start_ptr()), dim(), Range1D(lbound, ubound)); }

template<class T>
inline
VectorSliceTmpl<T> VectorTmpl<T>::rev()
{	return VectorSliceTmpl<T>( start_ptr() + dim() - 1, dim(), -1 ); }

template<class T>
inline
const VectorSliceTmpl<T> VectorTmpl<T>::rev() const
{	return VectorSliceTmpl<T>( const_cast<value_type*>(start_ptr()) + dim() - 1, dim(), -1 ); }

// Conversion operators
template<class T>
inline
VectorTmpl<T>::operator VectorSliceTmpl<T>()
{	return VectorSliceTmpl<T>(start_ptr(), dim()); }

template<class T>
inline
VectorTmpl<T>::operator const VectorSliceTmpl<T>() const
{	return VectorSliceTmpl<T>(const_cast<value_type*>(start_ptr()), dim()); }

// Assignment Operators
template<class T>
inline
VectorTmpl<T>& VectorTmpl<T>::operator=(value_type alpha) 
{
  if(!dim()) resize(1);
  std::fill(begin(),end(),alpha);
  return *this;
}

template<class T>
inline
VectorTmpl<T>& VectorTmpl<T>::operator=(const VectorTmpl<T>& rhs) 
{
  resize(rhs.dim());
  std::copy(rhs.begin(),rhs.end(),begin());
  return *this;
}

template<class T>
inline
VectorTmpl<T>& VectorTmpl<T>::operator=(const VectorSliceTmpl<T>& rhs) 
{
  resize(rhs.dim());
  std::copy(rhs.begin(),rhs.end(),begin());
  return *this;
}

// Raw pointer access

template<class T>
inline
typename VectorTmpl<T>::value_type*	VectorTmpl<T>::raw_ptr()
{	return start_ptr(); }

template<class T>
inline
const typename VectorTmpl<T>::value_type* VectorTmpl<T>::raw_ptr() const
{	return start_ptr(); }

template<class T>
inline
typename VectorTmpl<T>::value_type*	VectorTmpl<T>::start_ptr()
{	return dim() ? &(v_)[0] : 0; }

template<class T>
inline
const typename VectorTmpl<T>::value_type* VectorTmpl<T>::start_ptr() const
{	return &const_cast<valarray&>((v_))[0]; }

template<class T>
inline
typename VectorTmpl<T>::difference_type VectorTmpl<T>::stride() const
{	return 1; }

// //////////////////////////////////////////////////
// Non-inlined members

// Assume that the vector slices are the rows, cols or diag of a 2-D matrix.
template<class T>
EOverLap VectorSliceTmpl<T>::overlap(const VectorSliceTmpl<value_type>& vs) const {

  const typename VectorSliceTmpl<T>::value_type
    *raw_ptr1 = ( stride() > 0 ? start_ptr() : start_ptr() + (dim()-1)*stride() ),
    *raw_ptr2 = ( vs.stride() > 0 ? vs.start_ptr() : vs.start_ptr() + (vs.dim()-1)*vs.stride() );
  typename VectorSliceTmpl<T>::size_type
    size1 = dim(),
    size2 = vs.dim();
  typename VectorSliceTmpl<T>::difference_type
    stride1 = std::abs(stride()),
    stride2 = std::abs(vs.stride());

  // Establish a frame of reference where raw_ptr1 < raw_ptr2
  if(raw_ptr1 > raw_ptr2) {
    std::swap(raw_ptr1,raw_ptr2);
    std::swap(stride1,stride2);
    std::swap(size1,size2);
  }

  if( raw_ptr1 + stride1 * (size1 - 1) < raw_ptr2 ) {
    return NO_OVERLAP; // can't be any overlap
  }

  typename VectorSliceTmpl<T>::size_type
    start1 = 0,
    start2 = raw_ptr2 - raw_ptr1;

  if(start1 == start2 && stride1 == stride2 && size1 == size2)
    return SAME_MEM;
//	else if(start1 == start2)
//		return SOME_OVERLAP;	// First elements are the same
//	else if(stride1 + (size1 - 1) * stride1 == stride2 + (size2 - 1) * stride2)
//		return SOME_OVERLAP;	// Last elements are the same
  else if(stride1 == stride2) {
    if(!((start2 - start1) % stride1))
      return SOME_OVERLAP;
    else
      return NO_OVERLAP; // a different row, col or diag of a matrix
  }
  else {
    if(stride1 == 1 || stride2 == 1) {
      // One of them is a column vector.
      // Make vs1 the column vector.
      bool switch_them = (stride2 == 1);
      if(switch_them) {
        std::swap(start1,start2);
        std::swap(stride1,stride2);
        std::swap(size1,size2);
      }

      // Determine if the other vector could be row vector
      // or must be a diag vector.  If using stride2 makes
      // the first and last elements of the column vector
      // on different rows, then the other vector must be a diagonal.
      // col_first = start1/stride2, col_last = (start1+size1-1)/stride2
      // if(col_last - col_first > 0) then vs2 must be a diagonal vector
      // with max_rows = stride2 - 1.
      size_t max_rows = (start1+size1-1)/stride2 - start1/stride2 > 0 ? stride2 - 1 : stride2;

      // find the index (0-based) of vs2 that intersects this column.
      size_t vs2_col_i = start1/max_rows - start2/max_rows;
      
      // See if the first element of the column is above the element in vs2
      // and the last element of the column is below the element.  If it is
      // then we conclude that there is an itersection.
      size_t vs2_col_rng = start2 + vs2_col_i * stride2;
      if(start1 <= vs2_col_rng && vs2_col_rng <= start1+size1-1)
        return SOME_OVERLAP;
      else
        return NO_OVERLAP;
    }
    // They are not the same and nether is a column vector so one is a row vector
    // and the other is a diagonal vector.
    // Nether is a column vector so choose as vs1 the row vector (the smaller stride).
    bool switch_them = stride2 < stride1;
    if(switch_them) {
      std::swap(start1,start2);
      std::swap(stride1,stride2);
      std::swap(size1,size2);
    }

    size_t max_rows = stride1;
    // Determine the first and last columns (0-based) in the original
    // matrix where there vs1 and vs2 intersect.
    size_t	sec_first_col = (start1 > start2) ? start1/max_rows : start2/max_rows,
        last1 = start1 + (size1 - 1) * stride1,
        last2 = start2 + (size2 - 1) * stride2,
        sec_last_col = (last1 < last2) ? last1/max_rows : last2/max_rows;
    // Determine the vector indexes (0-based) of vs1 and vs2 for the start and end
    // in this region
    size_t	vs1_first_col = start1 / max_rows,
        vs2_first_col = start2 / max_rows;
    
    // Determine the indexes in the valarray of the two vectors for the two ends
    size_t	vs1_first_col_i = sec_first_col - vs1_first_col,
        vs1_last_col_i = sec_last_col - vs1_first_col,
        vs2_first_col_i = sec_first_col - vs2_first_col,
        vs2_last_col_i = sec_last_col - vs2_first_col;

    // Compare the indexes in the valarray at the two ends.  If they cross then
    // there must be an element of overlap.  Uses equivalent of the intermediate
    // value therorm.
    // Must cast to an int that can hold a negative value
    ptrdiff_t	diff1 = (start1 + vs1_first_col_i * stride1) 
            - static_cast<ptrdiff_t>((start2 + vs2_first_col_i * stride2)),
          diff2 = (start1 + vs1_last_col_i * stride1) 
            - static_cast<ptrdiff_t>((start2 + vs2_last_col_i * stride2));
    if(diff1 * diff2 > 0 )
      return NO_OVERLAP;		// they do not cross
    else
      return SOME_OVERLAP;	// they share an element
  }
}

template<class T>
EOverLap VectorTmpl<T>::overlap(const VectorSliceTmpl<value_type>& vs) const {

  const typename VectorSliceTmpl<T>::value_type
    *raw_ptr1 = ( stride() > 0 ? start_ptr() : start_ptr() + (dim()-1)*stride() ),
    *raw_ptr2 = ( vs.stride() > 0 ? vs.start_ptr() : vs.start_ptr() + (vs.dim()-1)*vs.stride() );
  typename VectorSliceTmpl<T>::size_type
    size1 = dim(),
    size2 = vs.dim();

  if( raw_ptr1 <= raw_ptr2 && raw_ptr2 + size2 <= raw_ptr1 + size1 ) {
    if( raw_ptr1 == raw_ptr2 && size1 == size2 && 1 == vs.stride() )
      return SAME_MEM;
    return SOME_OVERLAP;
  }
  return NO_OVERLAP;
}

} // end namespace DenseLinAlgPack

#endif	// end VECTOR_CLASS_TMPL_H
