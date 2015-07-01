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

#ifndef SPARSE_VECTOR_CLASS_DECL_H
#define SPARSE_VECTOR_CLASS_DECL_H

#include <assert.h>

#include <vector>
#include <sstream>

#include "AbstractLinAlgPack_SpVecIndexLookupClass.hpp"

namespace AbstractLinAlgPack {

namespace SparseVectorUtilityPack {

void assert_is_sorted(bool is_sorted);

/** \brief . */
class DoesNotExistException : public std::logic_error
{public: DoesNotExistException(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class NotSortedException : public std::logic_error
{public: NotSortedException(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class DuplicateIndexesException : public std::logic_error
{public: DuplicateIndexesException(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class OutOfRoomException : public std::logic_error
{public: OutOfRoomException(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class UnsizedException : public std::logic_error
{public: UnsizedException(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class NoNonZeroElementsException : public std::logic_error
{public: NoNonZeroElementsException(const std::string& what_arg) : std::logic_error(what_arg) {}};

}	// end namespace SparseVectorUtilityPack

// /////////////////////////////////////////////////////////////////////////////////////
/** @name Nonmember untility functions */
//@{

/** \brief Return a sparse vector slice.
  *
  * Preconditions:<ul>
  * <li> index_lookup.validate_state() is called
  * <li> [rng.full_range() != true] rng.ubound() <= size (throw #std::out_of_range#)
  * </ul>
  *
  * Postconditions:<ul>
  * <li> Fill them in latter
  * </ul>
  */
template<class T_Element>
SparseVectorSlice<T_Element> create_slice(
  const SparseVectorUtilityPack::SpVecIndexLookup<T_Element>& index_lookup
  , size_type size, Range1D rng);

//@}

template <class T_Element>
class SparseVectorSlice;

// ///////////////////////////////////////////////////////////////////////
/** \brief Sparse Vector class template.
  *
  * This is a class for abstracting a sparse vector of a templated element
  * type.  All of the operations are based on the element type.  Access
  * to the elements is provided by iterators.  It is also templated by 
  * an allocator that is that is used to allocate memory for the nonzero
  * elements.
  *
  * The templated type T_Element must support the following interface
  * (SparseElementTemplateInterface):
  * \begin{description}
  * <li>[value_type]			public typedef for the stored value of the element
  * <li>[index_type]			public typedef for the stored index of the element
  * <li>[value_type& value()]	function returning a lvalue for the value of the element
  * <li>[value_type value() const] const function returning a rvalue for the value of the element
  * <li>[index_type index() const] const function returning a rvalue for the index
  *										of the element.
  * <li>[T_Element& operator=(const T_Element&)] assignment operator
  * <li>[T_Element(const T_Element&)] copy constructor
  * \end{description}
  */
template <class T_Element, class T_Alloc = std::allocator<T_Element> >
class SparseVector {
public:
  /** @name Public Types. */
  //@{

  /** \brief . */
  typedef T_Alloc											allocator_type;
  /** \brief . */
  typedef T_Element										element_type;
  /** \brief . */
  typedef AbstractLinAlgPack::size_type					size_type;
  /** \brief . */
  typedef ptrdiff_t										difference_type;
  /** \brief . */
  typedef element_type*									iterator;
  /** \brief . */
  typedef const element_type*								const_iterator;

#if 0 /* defined(_WINDOWS) || defined(_INTEL_CXX) */

  typedef std::reverse_iterator<iterator, element_type
    , element_type&, element_type*, difference_type>	reverse_iterator;

  typedef std::reverse_iterator<const_iterator
    , element_type, const element_type&
    , const element_type*, difference_type>				const_reverse_iterator;

#else

  /** \brief . */
  typedef std::reverse_iterator<iterator>					reverse_iterator;
  /** \brief . */
  typedef std::reverse_iterator<const_iterator>			const_reverse_iterator;

#endif

  /** \brief . */
  typedef SparseVectorUtilityPack::DoesNotExistException		DoesNotExistException;
  /** \brief . */
  typedef SparseVectorUtilityPack::NotSortedException			NotSortedException;
  /** \brief . */
  typedef SparseVectorUtilityPack::DuplicateIndexesException	DuplicateIndexesException;
  /** \brief . */
  typedef SparseVectorUtilityPack::OutOfRoomException			OutOfRoomException;
  /** \brief . */
  typedef SparseVectorUtilityPack::UnsizedException			UnsizedException;
  /** \brief . */
  typedef SparseVectorUtilityPack::NoNonZeroElementsException	NoNonZeroElementsException;

  //@}

  /** @name Constuctors */
  //@{

  /** \brief Constructs a sparse vector with no elements (#nz() == dim() == 0#) and
    * assumes the elements are not sorted.
    */
  SparseVector(const allocator_type& alloc = allocator_type());

  /// Constructs a sparse vector with no elements (#nz() == dim() == 0#).
  SparseVector(bool assume_sorted, const allocator_type& alloc = allocator_type());

  /// Constructs a sparse vector of size #size# with storage for #max_nz# elements (#nz() == 0#)
  SparseVector(size_type size, size_type max_nz, difference_type offset = 0
    , bool assume_sorted = false, const allocator_type& alloc = allocator_type());

  /** \brief Constructs a sparse vector from another sparse vector.
    *
    * Copies the complete state including the same max_nz() but a fresh copy
    * of the elements are made.
    */
  SparseVector(const SparseVector<T_Element,T_Alloc>& sp_vec);

  /// Constructs a sparse vector of from a sparse vector slice.
  SparseVector( SparseVectorSlice<T_Element> sp_vec_slc
    , const allocator_type& alloc = allocator_type());

  /// Destructor (frees storage for elements).
  ~SparseVector();

  //@}

  /** \brief Assignment operator.
    *
    * If #max_nz() > sp_vec.nz()# then no new allocation takes place
    * otherwise #this# will will be resized to #sp_vec.nz()#.
    */
  SparseVector<T_Element,T_Alloc>& operator=(const SparseVector<T_Element,T_Alloc>& sp_vec);

  /** \brief Assignment operator.
    *
    * If #max_nz() > sp_vec_slc.nz()# then no new allocation takes place
    * otherwise #this# will will be resized to #sp_vec_slc.nz()#.
    */
  SparseVector<T_Element,T_Alloc>& operator=(const SparseVectorSlice<T_Element>& sp_vec_slc);

  /// 
  /** Returns the degree of memory overlap of this SparseVector and a SparseVectorSlice.
    *
    * @return 
    *		\begin{description}
    *		<li>[NO_OVERLAP]	There is no memory overlap between this and sv
    *		<li>[SOME_OVERLAP]	There is some memory locations that this and sv share
    *		<li>[SAME_MEM]		The DVectorSlice objects this and sv share the exact same memory locations.
    *		\end{description}
    */
  EOverLap overlap(const SparseVectorSlice<T_Element>& sv) const;

  /** @name SparseVectorTemplateInterface for linear algebra operations */
  //@{

  /// Return the number of elements in the full vector
  size_type dim() const;

  /// Return the number of non-zero elements
  size_type nz() const;

  /** \brief Return the offset for the indexes (ith real index = #begin()[i-1]->index() + offset()#
    * , for i = 1,,,#nz()#)
    */
  difference_type offset() const;

  /** \brief Return true if the sequence is sorted.
    *
    * If sorted() was called prior to this then it is garrented to be sorted
    * and if assume_sorted(true) was called then a client is assumed to be
    * responcible for it being sorted by it can not be garrented to be sorted.
    */
  bool is_sorted() const;

  /** @name Iterator access to elements.
    *
    * These functions return random access iterators that yield
    * SparseElementTemplateInterface objects when dereferenced.
    * This is required for the template argument.
    */
  //@{

  /** \brief Returns iterator that iterates forward through the nonzero elements.
    *
    * If #is_sorted() == true# then the elements will be forward iterated in accending
    * indexes.
    */
  iterator begin();

  /** \brief . */
  const_iterator begin() const;

  /// 
  iterator end();

  /** \brief . */
  const_iterator end() const;

  /** \brief Returns iterator that iterates backward through the nonzero elements.
    *
    * If #is_sorted() == true# then the elements will be forward iterated in deaccending
    * indexes.
    */
  reverse_iterator rbegin();

  /** \brief . */
  const_reverse_iterator rbegin() const;

  /// 
  reverse_iterator rend();

  /** \brief . */
  const_reverse_iterator rend() const;
  
  //	end Iterator access to elements
  //@}

  //	end SparseVectorTemplateInterface
  //@}

  /** @name Element setup and modification */
  //@{

  /** \brief Resize to #size# with a maximum of #max_nz# non-zero elements.
    *
    * This does not preserve the existing elements already in the sparse vector.
    * If you pass in #size == 0# or #max_nz == 0# then the storage will be deallocated
    * and no storage will be reallocated. 
    */
  void resize(size_type size, size_type max_nz, difference_type offset = 0);

  /** \brief Resize to #size# with a #max_nz# uninitialized non-zero elements.
    *
    * This function has the same basic behavior as #resize(...)# accept
    * on return #nz()# will equal #nz#.  The elements are initialized
    * to garbage so it is imparative that the client initialize
    * the elements before the sparse vector is used.
    */
  void uninitialized_resize(size_type size, size_type nz, size_type max_nz, difference_type offset = 0);

  /// Return the max number of elements that can be held without resizing
  size_type max_nz() const;

  /** \brief Add an unsorted element.
    *
    * If #nz() = max_nz()#) then the exception OutOfRoomException will be thrown.
    *
    * If you want to add more elements than you have reserved space for in the
    * construction or resize operation then you have to mannually copy the
    * elements in the sparse vector, resize the sparse vector, and then
    * readd the elements including the extra ones you want to add.
    */
  void add_element(element_type ele);

  /** \brief Add an element into a sorted sequence.
    *
    * If #nz() = max_nz()#) then the exception OutOfRoomException will be thrown.
    *
    * If you want to add more elements than you have reserved space for in the
    * construction or resize operation then you have to mannually copy the
    * elements in the sparse vector, resize the sparse vector, and then
    * readd the elements including the extra ones you want to add.
    */
  void insert_element(element_type ele);

  /** \brief Called by the client to inform this sparse vector object that the elements
    * be assumed to be in sequence and it is the clients responcibiliy to make sure
    * that it is.
    */
  void assume_sorted(bool assume_is_sorted);

  /// Sort the elements into assending order by index.
  void sort();

  /** \brief Assert that sparse vector is sorted.
    *
    * This function will throw an exception if any of the following are not true:
    * \begin{enumerate}
    * <li> The sequence is not sorted by index (#NotSortedException#)
    * <li> There are duplicate indexes (#DuplicateIndexesException#)
    * <li> The indexes are out of range (#std::out_of_range#)
    * \end{enumerate}
    *
    * This function will throw an exception for the first error it finds.
    */
  void assert_valid_and_sorted() const;

  //@}

  /** @name Lookup an element.
    *
    * If element v(i) exists, then a pointer to the element will
    * be returned.  If v(i) does not exist, then the NULL pointer
    * will be returned.
    *
    * If i is out of range then a std::out_of_range exception will be
    * thrown.
    * 
    * If the elements are sored then this operation is O(log(nz))
    * for a binary search.  Otherwise, it requries a O(nz) linear
    * search.
    */
  //@{

  /** \brief . */
  element_type* lookup_element(size_type i);
  /** \brief . */
  const element_type* lookup_element(size_type i) const;

  //@}

  /** @name Creating a slice (subregion) of the sparse vector.
    *
    * If the vector is not sorted (#is_sorted() == false#) then all of these
    * functions will throw an exception (#NotSortedException#).
    *
    * ** Say something about the cost of these operations! **
    */
  //@{

  /** \brief Allow an implicit conversion to a SparseVectorSlice<> object.
    *
    * This is a very cheap operation.
    */
  operator SparseVectorSlice<T_Element>();

  /** \brief . */
  operator const SparseVectorSlice<T_Element>() const;

  /** \brief Returns a SparseVectorSlice representing the entire sparse vector.
    *
    * It is used to provide a quick, explicit conversion so that
    * the SparseVector object can be used in functions that
    * expect a SparseVectorSlice object.
    */
  SparseVectorSlice<T_Element> operator()();

  /** \brief . */
  const SparseVectorSlice<T_Element> operator()() const;

  /// 
  /** Returns a continous subregion of the SparseVector object.
    *
    * The returned SparseVectorSlice object represents the range of the rng argument.
    *
    * Preconditions: <ul>
    *		<li> #rng.ubound() - 1 <= this->dim()# (throw #out_of_range#)
    *		<li> #dim() > 0#	(throw #UnsizedException#)
    *		</ul>
    *
    * Postconditions: <ul>
    *		<li> returned#.dim() == rng.ubound() - rng.lbound() + 1#
    *     <li> contains all of the elements in the range.
    *		</ul>
    *
    * @param	rng		Index range [lbound,ubound] of the region being returned.
    */
  SparseVectorSlice<T_Element> operator()(const Range1D& rng);

  /** \brief . */
  const SparseVectorSlice<T_Element> operator()(const Range1D& rng) const;

  /// 
  /** Returns a SparseVectorSlice object for the continous subregion [ubound, lbound].
    * 
    * Preconditions: <ul>
    *		<li> #lbound > 1# (throw #out_of_range#)
    *		<li> #lbound < ubound# (throw #out_of_range#)
    *		<li> #ubound <= this->dim()# (throw #out_of_range#)
    *		</ul>
    *
    * Postconditions: <ul>
    *		<li> returned#.dim() == ubound() - lbound() + 1#
    *     <li> contains all of the elements in the range.
    *		</ul>
    *
    * @param	lbound		Lower bound of range [lbound,ubound] of the region being returned.
    * @param	ubound		Upper bound of range [lbound,ubound] of the region being returned.
    */
  SparseVectorSlice<T_Element> operator()(size_type lbound, size_type ubound);

  /// Same as above.
  const SparseVectorSlice<T_Element> operator()(size_type lbound, size_type ubound) const;

  //@}

private:

  // /////////////////////////////////////////////////////////////////////////
  // Private types

  /** \brief . */
  typedef SparseVectorUtilityPack::SpVecIndexLookup<element_type> SpVecIndexLookup;

  // /////////////////////////////////////////////////////////////////////////
  // Private data members

  allocator_type			alloc_;			// allocator used to allocate memory
  size_type				size_;			// the number of elements in the full vector
  size_type				max_nz_;		// the amount of storage that has been allocated
//  commented out because of problems with MS Visual C++ 5.0
//	std::vector<element_type, allocator_type>	ele_;
  SpVecIndexLookup		index_lookup_;	// Acts as storage for elements and caching of searches.
  bool					assume_sorted_;	// true if the client said that you can assume sorted.
  bool					know_is_sorted_; // true if it has been varified that is sorted.

  // //////////////////////////
  // Private member functions

  // Throw a NotSortedException of is_sorted() == false
  void assert_is_sorted() const {
    SparseVectorUtilityPack::assert_is_sorted(is_sorted());
  }

  /// Assert (#OutOfRoom#) that there is room for n elements.
  void assert_space(size_type n) const {
#ifdef LINALGPACK_CHECK_SLICE_SETUP
    if(index_lookup_.nz() + n > max_nz_)
      throw OutOfRoomException("SparseVector<T_Element,T_Alloc>::assert_space():  There is not storage for this many elements");
#endif
  }

  /// Assert #dim() > 0# (#UnsizedException#) and #index_lookup_.ele() != 0# (#NoNonZeroElementsException#)
  void assert_sized_with_mem_set() const {
    if(!dim())
      throw UnsizedException("SparseVector<...>::assert_sized_with_mem_set() : "
        "Error: The sparse vector is unsized");
    if(!index_lookup_.ele()) {
      throw NoNonZeroElementsException("SparseVector<...>::assert_sized_with_mem_set() : "
        "Error: There is no memory set.");
    }
  }

  /// Return the entire vector slice
  SparseVectorSlice<T_Element> get_whole_sp_vec() {
    return SparseVectorSlice<T_Element>(index_lookup_.ele(), index_lookup_.nz()
          , index_lookup_.offset(), size_, is_sorted());
  }

  /** \brief . */
  const SparseVectorSlice<T_Element> get_whole_sp_vec() const {
    return SparseVectorSlice<T_Element>(index_lookup_.ele(), index_lookup_.nz()
          , index_lookup_.offset(), size_, is_sorted());
  }

  /// Return a SparseVectorSlice (inplementation for indexing operators)
  SparseVectorSlice<T_Element> get_slice(const Range1D& rng) const {
    assert_is_sorted();
    return create_slice(index_lookup_, size_, rng);
  }

};	// end class SparseVector

// ///////////////////////////////////////////////////////////////////////
/** \brief Sparse Vector Slice class template.
  *
  * This is a class for abstracting a region of a sparse vector stored
  * as an array of elements of a templated type.  The required inteface
  * for the type of these elements is given in the SparseVector documentation.
  *
  * Here if nz() == 0 then begin() == end() so it is safe to set up loops
  * in the form of:
  *
  *  for(SparseVectorSlice<T>::const_iterator itr = sv.begin(); itr != sv.end(); ++itr)
  *    // some access of *itr.
  *
  * Note that if nz()==0 then begin() may not point to a valid object so don't do it.
  *
  * The default copy constructor is allowed but the default constructor and
  * assignment operator functions are not.
  */
template <class T_Element>
class SparseVectorSlice {
public:
  /** @name Public types. */
  //@{

  /** \brief . */
  typedef T_Element										element_type;
  /** \brief . */
  typedef AbstractLinAlgPack::size_type					size_type;
  /** \brief . */
  typedef ptrdiff_t										difference_type;
  /** \brief . */
  typedef element_type*									iterator;
  /** \brief . */
  typedef const element_type*								const_iterator;

#if 0 /* defined(_WINDOWS) || defined(_INTEL_CXX) */

  typedef std::reverse_iterator<iterator, element_type
    , element_type&, element_type*, difference_type>	reverse_iterator;

  typedef std::reverse_iterator<const_iterator
    , element_type, const element_type&
    , const element_type*, difference_type>				const_reverse_iterator;

#else

  /** \brief . */
  typedef std::reverse_iterator<iterator>					reverse_iterator;
  /** \brief . */
  typedef std::reverse_iterator<const_iterator>			const_reverse_iterator;

#endif

  /** \brief . */
  typedef SparseVectorUtilityPack::DoesNotExistException	DoesNotExistException;
  /** \brief . */
  typedef SparseVectorUtilityPack::NotSortedException		NotSortedException;

  //@}

  /** @name Constuctors
    *
    * The default copy constructor is allowed since it has the proper semantics.
    */
  //@{

  /** \brief Constructs a sparse vector slice from an array of elements.
    *
    * Here a pointer to an array of elements is used instead of a
    * pointer to std::vector<T_Ele,T_Alloc> in order to insulated
    * this class from the type of allocator used since this information
    * is not needed.
    *
    * A sparse vector slice with no nonzero elements can be constructed by
    * setting nz == 0;
    *
    * Preconditions: <ul>
    *		<li> #ele != 0#
    *		<li> #size >= nz#
    *		</ul>
    *
    * @param	ele		pointer to array of elements (length #nz#)
    * @param	offset	offset for the indexes of the elements. index = ele[i].index() + offset
    * @param	size	number of elements in the full vector
    * @param	nz		number of non-zero elements in vector
    */
  SparseVectorSlice(element_type ele[], size_type nz, difference_type offset, size_type size
    , bool assume_sorted = false);

  //@}

  /** \brief Constructs a sparse vector slice view from another sparse vector slice.
    */
  void bind(SparseVectorSlice svs);

  /// 
  /** Returns the degree of memory overlap of this SparseVector and a SparseVectorSlice.
    *
    * @return 
    *		\begin{description}
    *		<li>[NO_OVERLAP]	There is no memory overlap between this and sv
    *		<li>[SOME_OVERLAP]	There is some memory locations that this and sv share
    *		<li>[SAME_MEM]		The DVectorSlice objects this and sv share the exact same memory locations.
    *		\end{description}
    */
  EOverLap overlap(const SparseVectorSlice<T_Element>& sv) const;

  /** @name Sparse Vector Templated interface for linear algebra operations */
  //@{

  /// Return the number of elements in the full vector
  size_type dim() const;

  /// Return the number of non-zero elements
  size_type nz() const;

  /** \brief Return the offset for the indexes (ith real index = #begin()[i-1]->index() + offset()#
    * , for i = 1,,,#nz()#)
    */
  difference_type offset() const;

  /** \brief Return true if the sequence is assumed sorted.
    */
  bool is_sorted() const;

  /** \brief . */
  iterator begin();

  /** \brief . */
  const_iterator begin() const;

  /// 
  iterator end();

  /** \brief . */
  const_iterator end() const;

  /** \brief . */
  reverse_iterator rbegin();

  /** \brief . */
  const_reverse_iterator rbegin() const;

  /// 
  reverse_iterator rend();

  /** \brief . */
  const_reverse_iterator rend() const;

  //@}

  /** @name Lookup an element.
    *
    * If element v(i) exists, then a pointer to the element will
    * be returned.  If v(i) does not exist, then the NULL pointer
    * will be returned.
    *
    * If i is out of range then a std::out_of_range exception will be
    * thrown.
    * 
    * If the elements are sored then this operation is O(log(nz))
    * for a binary search.  Otherwise, it requries a O(nz) linear
    * search.
    */
  //@{

  /** \brief . */
  element_type* lookup_element(size_type i);
  /** \brief . */
  const element_type* lookup_element(size_type i) const;

  //@}

  /** @name Creating a slice (subregion) of the sparse vector */
  //@{

  /// 
  /** Returns a SparseVectorSlice<> reference to this object.
    *
    * It is included for uniformity with SparseVector.
    */
  SparseVectorSlice<T_Element>& operator()();

  /** \brief . */
  const SparseVectorSlice<T_Element>& operator()() const;

  /// Allow address to be taken of an rvalue of this object
  SparseVectorSlice* operator&()
  {    return this; }

  const SparseVectorSlice* operator&() const
  {    return this; }

  /// 
  /** Returns a continous subregion of the SparseVector object.
    *
    * The returned SparseVectorSlice object represents the range of the rng argument.
    *
    * Preconditions: <ul>
    *		<li> #rng.ubound() - 1 <= this->dim()# (throw #out_of_range#)
    *		<li> #dim() > 0#	(throw #UnsizedException#)
    *		</ul>
    *
    * Postconditions: <ul>
    *		<li> returned#.dim() == rng.ubound() - rng.lbound() + 1#
    *     <li> contains all of the elements in the range.
    *		</ul>
    *
    * @param	rng		Index range [lbound,ubound] of the region being returned.
    */
  SparseVectorSlice<T_Element> operator()(const Range1D& rng);

  /** \brief . */
  const SparseVectorSlice<T_Element> operator()(const Range1D& rng) const;

  /// 
  /** Returns a SparseVectorSlice object for the continous subregion [ubound, lbound].
    * 
    * Preconditions: <ul>
    *		<li> #lbound > 1# (throw #out_of_range#)
    *		<li> #lbound < ubound# (throw #out_of_range#)
    *		<li> #ubound <= this->dim()# (throw #out_of_range#)
    *		</ul>
    *
    * Postconditions: <ul>
    *		<li> returned#.dim() == ubound() - lbound() + 1#
    *     <li> contains all of the elements in the range.
    *		</ul>
    *
    * @param	lbound		Lower bound of range [lbound,ubound] of the region being returned.
    * @param	ubound		Upper bound of range [lbound,ubound] of the region being returned.
    */
  SparseVectorSlice<T_Element> operator()(size_type lbound, size_type ubound);

  /** \brief . */
  const SparseVectorSlice<T_Element> operator()(size_type lbound, size_type ubound) const;

  //@}

private:
  // /////////////////////////////////////////////////////////////////////////
  // Private types

  /** \brief . */
  typedef SparseVectorUtilityPack::SpVecIndexLookup<element_type> index_lookup_type;

  // /////////////////////////////////////////////////////////////////////////
  // Private data members

  index_lookup_type		index_lookup_;	// Acts as storage and cacheing
  size_type				size_;			// size of the full vector
  bool					assume_sorted_;	// true if the client said that you can assume sorted.


  // /////////////////////////////////////////////////////////////////////////
  // Private member functions

  // Throw a NotSortedException of is_sorted() == false
  void assert_is_sorted() const {
    SparseVectorUtilityPack::assert_is_sorted(is_sorted());
  }

  /// Return a SparseVectorSlice (inplementation for indexing operators)
  SparseVectorSlice<T_Element> get_slice(const Range1D& rng) const {
    assert_is_sorted();
    return create_slice(index_lookup_, size_, rng);
  }

  /// Not defined and not to be called
  SparseVectorSlice();
  /// Not defined and not to be called
  SparseVectorSlice<element_type>& operator=(const SparseVectorSlice<element_type>&);

};	// end class SparseVectorSlice

// ////////////////////////////////////////////////
// Non-member non-public utility functions

namespace SparseVectorUtilityPack {

/** \brief Lookup an element.
  *
  * If the element does not exist, then NULL will be returned.
  */
template< class T_Element >
inline const T_Element* lookup_element( const SpVecIndexLookup<T_Element>& index_lookup
  , typename SpVecIndexLookup<T_Element>::index_type index, bool is_sorted )
{
  size_type poss;	
  return ( ( poss = index_lookup.find_element(index,is_sorted) ) < index_lookup.nz() )
    ? index_lookup.ele() + poss
    : NULL;
}

}	// end namespace SparseVectorUtilityPack

// /////////////////////////////////////////////////////////////////////////////////////
// Inline members for SparseVector<>

// constructors

template <class T_Element, class T_Alloc>
inline SparseVector<T_Element,T_Alloc>::SparseVector(const allocator_type& alloc)
  : alloc_(alloc), size_(0), max_nz_(0), assume_sorted_(false), know_is_sorted_(false)
{}

template <class T_Element, class T_Alloc>
inline SparseVector<T_Element,T_Alloc>::SparseVector(bool assume_sorted,const allocator_type& alloc)
  : alloc_(alloc), size_(0), max_nz_(0), assume_sorted_(assume_sorted), know_is_sorted_(false)
{}

template <class T_Element, class T_Alloc>
inline SparseVector<T_Element,T_Alloc>::SparseVector(size_type size, size_type max_nz
    , difference_type offset, bool assume_sorted, const allocator_type& alloc)
  : alloc_(alloc), size_(0), max_nz_(0), assume_sorted_(assume_sorted), know_is_sorted_(false)
{
  resize(size,max_nz,offset);
}

template <class T_Element, class T_Alloc>
inline SparseVector<T_Element,T_Alloc>::~SparseVector() {
  resize(0,0);
}

// SparseVectorTemplateInterface for linear algebra operations

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::size_type SparseVector<T_Element,T_Alloc>::dim() const {
  return size_;
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::size_type SparseVector<T_Element,T_Alloc>::nz() const {
  return index_lookup_.nz();
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::difference_type SparseVector<T_Element,T_Alloc>::offset() const {
  return index_lookup_.offset();
}

template <class T_Element, class T_Alloc>
inline bool SparseVector<T_Element,T_Alloc>::is_sorted() const {
  return nz() <= 1 || assume_sorted_ || know_is_sorted_;
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::iterator SparseVector<T_Element,T_Alloc>::begin() {
  return index_lookup_.nz() ? index_lookup_.ele() : NULL;
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::const_iterator SparseVector<T_Element,T_Alloc>::begin() const {
  return index_lookup_.nz() ? index_lookup_.ele() : NULL;
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::iterator SparseVector<T_Element,T_Alloc>::end() {
  return index_lookup_.nz() ? index_lookup_.ele() + index_lookup_.nz() : NULL;
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::const_iterator SparseVector<T_Element,T_Alloc>::end() const {
  return index_lookup_.nz() ? index_lookup_.ele() + index_lookup_.nz() : NULL;
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::reverse_iterator SparseVector<T_Element,T_Alloc>::rbegin() {
  return reverse_iterator(end());
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::const_reverse_iterator SparseVector<T_Element,T_Alloc>::rbegin() const {
  return const_reverse_iterator(end());
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::reverse_iterator SparseVector<T_Element,T_Alloc>::rend() {
  return reverse_iterator(begin());
}

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::const_reverse_iterator SparseVector<T_Element,T_Alloc>::rend() const {
  return const_reverse_iterator(begin());
}

// Element setup and modification

template <class T_Element, class T_Alloc>
inline typename SparseVector<T_Element,T_Alloc>::size_type SparseVector<T_Element,T_Alloc>::max_nz() const {
  return max_nz_;
}

template <class T_Element, class T_Alloc>
inline void SparseVector<T_Element,T_Alloc>::add_element(element_type ele) {
  assert_space(1);
  assume_sorted_ = know_is_sorted_ = false;
#ifdef _PG_CXX
  new (index_lookup_.ele() + index_lookup_.nz()) element_type;
#else
  alloc_.construct(index_lookup_.ele() + index_lookup_.nz(), ele);
#endif
  index_lookup_.incr_nz();
}

template <class T_Element, class T_Alloc>
inline void SparseVector<T_Element,T_Alloc>::assume_sorted(bool assume_is_sorted) {
  assume_sorted_ = assume_is_sorted;
}

// Lookup an element

template <class T_Element, class T_Alloc>
inline
typename SparseVector<T_Element,T_Alloc>::element_type*
SparseVector<T_Element,T_Alloc>::lookup_element(size_type i)
{
  return const_cast<element_type*>(SparseVectorUtilityPack::lookup_element(index_lookup_,i,assume_sorted_));
}

template <class T_Element, class T_Alloc>
inline
const typename SparseVector<T_Element,T_Alloc>::element_type*
SparseVector<T_Element,T_Alloc>::lookup_element(size_type i) const
{
  return SparseVectorUtilityPack::lookup_element(index_lookup_,i,assume_sorted_);
}

// Creating a slice (subregion) of the sparse vector

template <class T_Element, class T_Alloc>
inline SparseVector<T_Element,T_Alloc>::operator SparseVectorSlice<T_Element>() {
  return get_whole_sp_vec();
}

template <class T_Element, class T_Alloc>
inline SparseVector<T_Element,T_Alloc>::operator const SparseVectorSlice<T_Element>() const {
  return get_whole_sp_vec();
}

template <class T_Element, class T_Alloc>
inline SparseVectorSlice<T_Element> SparseVector<T_Element,T_Alloc>::operator()() {
  return get_whole_sp_vec();
}

template <class T_Element, class T_Alloc>
inline const SparseVectorSlice<T_Element> SparseVector<T_Element,T_Alloc>::operator()() const {
  return get_whole_sp_vec();
}

template <class T_Element, class T_Alloc>
inline SparseVectorSlice<T_Element> SparseVector<T_Element,T_Alloc>::operator()(const Range1D& rng) {
  return get_slice(rng);
}

template <class T_Element, class T_Alloc>
inline const SparseVectorSlice<T_Element> SparseVector<T_Element,T_Alloc>::operator()(const Range1D& rng) const {
  return get_slice(rng);
}

template <class T_Element, class T_Alloc>
inline SparseVectorSlice<T_Element> SparseVector<T_Element,T_Alloc>::operator()(size_type lbound, size_type ubound) {
  return get_slice(Range1D(lbound,ubound));
}

template <class T_Element, class T_Alloc>
inline const SparseVectorSlice<T_Element> SparseVector<T_Element,T_Alloc>::operator()(size_type lbound, size_type ubound) const {
  return get_slice(Range1D(lbound,ubound));
}

// /////////////////////////////////////////////////////////////////////////////////////
// Inline members for SparseVectorSlice<>

// Constuctors

template <class T_Element>
inline SparseVectorSlice<T_Element>::SparseVectorSlice(element_type ele[], size_type nz
    , difference_type offset, size_type size, bool assume_sorted)
  : index_lookup_(ele,nz,offset), size_(size), assume_sorted_(assume_sorted)
{}

template <class T_Element>
inline void SparseVectorSlice<T_Element>::bind(SparseVectorSlice svs)
{
  index_lookup_	= svs.index_lookup_;
  size_			= svs.size_;
  assume_sorted_	= svs.assume_sorted_;
}

// Sparse Vector Templated interface for linear algebra operations

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::size_type SparseVectorSlice<T_Element>::dim() const {
  return size_;
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::size_type SparseVectorSlice<T_Element>::nz() const {
  return index_lookup_.nz();
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::difference_type SparseVectorSlice<T_Element>::offset() const {
  return index_lookup_.offset();
}

template <class T_Element>
inline bool SparseVectorSlice<T_Element>::is_sorted() const {
  return nz() <= 1 || assume_sorted_;
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::iterator SparseVectorSlice<T_Element>::begin() {
  return index_lookup_.ele();
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::const_iterator SparseVectorSlice<T_Element>::begin() const {
  return index_lookup_.ele();
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::iterator SparseVectorSlice<T_Element>::end() {
  return index_lookup_.ele() + index_lookup_.nz();
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::const_iterator SparseVectorSlice<T_Element>::end() const {
  return index_lookup_.ele() +  index_lookup_.nz();
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::reverse_iterator SparseVectorSlice<T_Element>::rbegin() {
  return reverse_iterator(end());
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::const_reverse_iterator SparseVectorSlice<T_Element>::rbegin() const {
  return const_reverse_iterator(end());
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::reverse_iterator SparseVectorSlice<T_Element>::rend() {
  return reverse_iterator(begin());
}

template <class T_Element>
inline typename SparseVectorSlice<T_Element>::const_reverse_iterator SparseVectorSlice<T_Element>::rend() const {
  return const_reverse_iterator(begin());
}

// Lookup an element

template <class T_Element>
inline
typename SparseVectorSlice<T_Element>::element_type*
SparseVectorSlice<T_Element>::lookup_element(size_type i)
{
  return const_cast<element_type*>(SparseVectorUtilityPack::lookup_element(index_lookup_,i,assume_sorted_));
}

template <class T_Element>
inline
const typename SparseVectorSlice<T_Element>::element_type*
SparseVectorSlice<T_Element>::lookup_element(size_type i) const
{
  return SparseVectorUtilityPack::lookup_element(index_lookup_,i,assume_sorted_);
}

// Creating a slice (subregion) of the sparse vector

template <class T_Element>
inline SparseVectorSlice<T_Element>& SparseVectorSlice<T_Element>::operator()() {
  return *this;
}

template <class T_Element>
inline const SparseVectorSlice<T_Element>& SparseVectorSlice<T_Element>::operator()() const {
  return *this;
}

template <class T_Element>
inline SparseVectorSlice<T_Element> SparseVectorSlice<T_Element>::operator()(const Range1D& rng) {
  return get_slice(rng);
}

template <class T_Element>
inline const SparseVectorSlice<T_Element> SparseVectorSlice<T_Element>::operator()(const Range1D& rng) const {
  return get_slice(rng);
}

template <class T_Element>
inline SparseVectorSlice<T_Element> SparseVectorSlice<T_Element>::operator()(size_type lbound, size_type ubound) {
  return get_slice(Range1D(lbound,ubound));
}

template <class T_Element>
inline const SparseVectorSlice<T_Element> SparseVectorSlice<T_Element>::operator()(size_type lbound, size_type ubound) const {
  return get_slice(Range1D(lbound,ubound));
}

} // end namespace AbstractLinAlgPack 

#endif // SPARSE_VECTOR_CLASS_DECL_H
