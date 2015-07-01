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

#ifndef SPVEC_INDEX_LOOKUP_CLASS_DECL_H
#define SPVEC_INDEX_LOOKUP_CLASS_DECL_H

#include <stdexcept>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

namespace SparseVectorUtilityPack {

// ///////////////////////////////////////////////////////////////////////
/** \brief Sparse Vector Index Lookup and Caching class.
  *
  * This class is used to perform a lookup for elements in a sparse vector
  * stored as an array of nonzero elements of a templated type T_Element.
  * The type T_Element must conform to the SparseElementTemplateInterface
  * specification.  These elements must be sorted in accending order.
  *
  * The default C++ assignment operator and copy constructor are allowed.
  */
template <class T_Element>
class SpVecIndexLookup {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef T_Element							element_type;
  /** \brief . */
  typedef typename element_type::index_type	index_type;
  /** \brief . */
  typedef ptrdiff_t							difference_type;
  /** \brief . */
  class NoSpVecSetException : public std::logic_error
  {public: explicit NoSpVecSetException(const std::string& what_arg) : std::logic_error(what_arg) {}};
  /** \brief . */
  class InvalidInternalStateException : public std::logic_error
  {public: explicit InvalidInternalStateException(const std::string& what_arg) : std::logic_error(what_arg) {}};
  /** \brief . */
  enum UpperLower { UPPER_ELE, LOWER_ELE };
  /** \brief . */
  enum ElementRelation { BEFORE_ELE, AFTER_ELE, EQUAL_TO_ELE };
  /// Struct with members: size_type poss; ElementRelation rel;
  struct poss_type {
    poss_type() : poss(0), rel(EQUAL_TO_ELE) {} 
    poss_type(size_type _poss, ElementRelation _rel) : poss(_poss), rel(_rel) {} 
    size_type			poss;
    ElementRelation		rel;
  };

  //@}

  /** @name Constructors */
  //@{

  /** \brief Construct uninitialized with not sparse vector set (#ele() == 0#) */
  SpVecIndexLookup()
    : ele_(0), nz_(0), offset_(0), index_cached_(0)
  {}

  /** \brief Construct initialize with a sparse vector */
  SpVecIndexLookup(element_type* ele, size_type nz, difference_type offset)
    : ele_(ele), nz_(nz), offset_(offset), index_cached_(0)
  {}

  //@}

  /** @name Sparse vector representation setup */
  //@{

  /** \brief Set a new sparse vector.
    *
    * This will wipe out any cache that may be stored.
    */
  void set_sp_vec(element_type* ele, size_type nz, difference_type offset) {
    ele_ = ele;		nz_ = nz;		offset_ = offset;
    sp_vec_was_modified();
  }

  /// Increment nz only
  void incr_nz() {
    nz_++;
  }

  //@}

  /** @name Sparse vector representation access */
  //@{

  /** \brief . */
  element_type*		ele() const {
    return ele_;
  }
  /** \brief . */
  size_type			nz() const {
    return nz_;
  }
  /** \brief . */
  difference_type		offset() const {
    return offset_;
  }

  //@}

  /** @name Element lookup and caching */
  //@{

  /** \brief Lookup an element and cache the result if a binary search was performed.
    *
    * This function should only be used if it can be assumed that the elements
    * are sorted in assending order.
    *
    * If #index# is the same as the previously cached lookup then this function
    * will execute in O(1) time, otherwise a O(log(nz)) binary search will be
    * performed to find the element and the result of the lookup will be cached.
    *
    * To be able to utilize a previously cached search this function must know
    * if an upper element or lower element is to be found.\\
    *
    * Preconditions:<ul>
    *	<li> #ele() > 0# (throw #NoSpVecSetException#)
    * </ul>
    *
    * Postconditions:<ul>
    *	<li> [uplow == lower_ele] #index <= ele()[return.poss].index() + offset()# 
    *	<li> [uplow == upper_ele] #ele()[return.poss].index() + offset() <= index# 
    * </ul>
    *
    * @return		#poss_type# object where #return.poss# gives a position in tye
    *				underlying array and #return.rel# gives where the element with
    *				#index# is in relation to this possition.
    *									[ #BEFORE_ELE#		if #ele()[return.poss].index() + offset() < index#	]
    *					#return.rel# =	[ #EQUAL_TO_ELE#	if #ele()[return.poss].index() + offset() == index#	]
    *									[ #AFTER_ELE#		if #ele()[return.poss].index() + offset() > index#	]
    */
  poss_type find_poss(index_type index, UpperLower uplow) const;

  /** \brief Lookup an element.
    *
    * Lookup an exact match for an element.  If the element is not found, the
    * end of the sequence will be returned.
    * 
    * If is_sorted == true then a binary search will be used (O(log(nz)).  If is_sorted==false
    * a sequential search will be used (O(nz)).  No result is cached here.
    * 
    * Preconditions:<ul>
    *	<li> #ele() > 0# (throw #NoSpVecSetException#)
    * </ul>
    *
    * Postconditions:<ul>
    *	<li> [return == nz()] No element exits with this index 
    *	<li> [return < nz()] #index == ele()[return].index() + offset()# 
    * </ul>
    */
  size_type find_element( index_type index, bool is_sorted ) const;
  
  //@}

  /** @name State management */
  //@{

  /** \brief Called by client to inform the object that the sparse vector was modified
    * so that the cache may be invalidated.
    */
  void sp_vec_was_modified() {
    index_cached_ = 0;
  }

  /** \brief Called by client to ensure that the internal state is valid.
    *
    * If #ele() != 0# but #nz == 0# or #ele()[0].index() + offset() < 1# then this
    * is an invalid sparse vector and the function will throw a #NoSpVecSetException#
    * exception.  It is up to the client to ensure that a valid sparse vector is set.
    *
    * If there in some other problem with internal state then an exception
    * #InvalidInternalStateException# will be thrown.  The error message will be
    * meant for a developer to inspect.
    */
  void validate_state() const;

  //@}

private:
  // ///////////////////////////////////////////////////////////////////////////////
  // Private types

  // //////////////////////////////////////////////////////////////////////////////
  // Private data members

  element_type*		ele_;		// pointer to array of elements
  size_type			nz_;		// number of nonzero elements in ele_
  difference_type		offset_;	// offset for each element in ele_.  The actuall index for
                  // each element i is ele_[i].index() + offset().
  mutable index_type		index_cached_;	// index of last binary search
  mutable size_type		poss_cached_;	// possition last looked up
  mutable ElementRelation	ele_rel_cached_;// Specifies where the element looked up is in relation
                      // to the element at poss_cached:
                      //		before_ele:	zero element before poss_cached_
                      //		after_ele:	zero element after poss_cached_
                      //		equal_to_ele: nonzero element at poss_cached_

  // ////////////////////////
  // Private member functions

  /// Assert that a sparse vector has been set up
  void assert_sp_vec_setup() const {
    if(!ele_)
      throw NoSpVecSetException("The sparse vector (ele, nz, offset) has not been set yet");
  }

  /// Adjust the cached possition
  size_type adjust_cached_poss(UpperLower uplow) const;

  /** \brief Perform a binary search for an element in the sparse vector.
    *
    * @param	index	[I]	index begin looked for
    * @param	uplow	[I]	whether it is an upper (#UPPER_ELE#) or lower (#LOWER_ELE#) element needed
    * @param	poss	[O]	possition where the element was found.  If #uplow == LOWER_ELE# then
    *						#poss# will be the largest possible integer that satisfies:
    *						#index <= ele_[poss].index()#.  If #uplow == UPPER_ELE# then #poss#
    *						will be the smallest possible integer that satisfies:
    *						#ele_[poss].index() <= index#
    *	@param	ele_rel	[O]	Where the element with #index# is in relation to the element at #poss#
    *						returned.  There are three possible values.
    *						#BEFORE_ELE# :	The element saught is a zero element that is before
    *										the element at #poss#
    *						#AFTER_ELE# :	The element saught is a zero element that is after the 
    *										element at #poss#
    *						#EQUAL_TO_POSS#: This is a nonzero elment that exists at possition #poss#
    */
  poss_type binary_ele_search(index_type index, UpperLower uplow) const;


};	// end class SpVecIndexLookup

}	// namespace SparseVectorUtilityPack

} // end namespace AbstractLinAlgPack

#endif // SPVEC_INDEX_LOOKUP_CLASS_DECL_H
