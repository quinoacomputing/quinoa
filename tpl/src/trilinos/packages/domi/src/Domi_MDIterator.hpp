// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef DOMI_MDITERATOR_HPP
#define DOMI_MDITERATOR_HPP

// Standard include
#include <iterator>

// Domi includes
#include "Domi_Exceptions.hpp"
#include "Domi_Utils.hpp"

namespace Domi
{

/** \brief Iterator class suitable for multi-dimensional arrays
 *
 * This iterator is specialized for multi-dimensional arrays.  It can
 * be used like a standard iterator, without knowledge of the number
 * of dimensions of the array.  It is designed to iterate only over
 * valid elements of the array.  It is possible, especially with
 * sliced views into a parent array, for the data buffer to have
 * stride gaps that do not belong to the array.  This iterator avoids
 * those stride gaps by keeping an internal record of the multi-
 * dimensional index, and only iterating over valid indexes.
 *
 * To apply MDIterator to all three multi-dimensional array types
 * (MDArray, MDArrayView and MDArrayRCP), the class is templated on
 * parameter class MDARRAY, which is intended to be any one of the
 * three MDArray types.  The MDARRAY class is expected to support the
 * the typedefs size_type, dim_type, and value_type, the numDims()
 * method and the _ptr, _strides, _dimensions, and _layout attributes.
 *
 * It is intended that the array class that will use the MDIterator
 * will declare the MDIterator to be a friend and to typedef the fully
 * qualified class.  For example, within the MDArray< T > class:
 *
 *   \code
 *   friend class MDIterator< MDArray< T > >;
 *   friend class MDIterator< MDArray< const T > >;
 *
 *   typedef MDIterator< MDArray< T > > iterator;
 *   typedef MDIterator< MDArray< const T > > const_iterator;
 *   \endcode
 *
 * and declare begin(), end() and cbegin(), cend() methods:
 *
 *   \code
 *   iterator begin() { return iterator(*this      ); }
 *   iterator end()   { return iterator(*this, true); }
 *
 *   const_iterator begin() const { return const_iterator(*this      ); }
 *   const_iterator end()   const { return const_iterator(*this, true); }
 *
 *   const_iterator cbegin() const { return const_iterator(*this      ); }
 *   const_iterator cend()   const { return const_iterator(*this, true); }
 *   \endcode
 */
template< class MDARRAY >
class MDIterator : public std::iterator< std::bidirectional_iterator_tag,
                                         typename MDARRAY::value_type >
{
public:

  /** \name MDARRAY typedefs */
  //@{

  /** \brief Value type */
  typedef typename MDARRAY::value_type value_type;

  /** \brief Pointer type */
  typedef typename MDARRAY::pointer pointer;

  //@}

  /** \name Constructors and Destructor */
  //@{

  /** \brief MDIterator constructor
   *
   *  \param mdarray [in] The multi-dimensional array object on which
   *         the iterator will act upon
   *
   *  \param end_index [in] If true, set the internal index to the
   *         MDARRAY end() index.  If false, set the internal index to
   *         the MDARRAY begin() index.  Default false.
   *
   *  Produces an iterator with index corresponding to either the
   *  MDARRAY <tt>begin()</tt> or <tt>end()</tt> methods, depending on
   *  the value of the <tt>end_index</tt> argument.
   */
  MDIterator(const MDARRAY & mdarray,
             bool end_index = false);

  /** \brief Index constructor
   *
   *  \param mdarray [in] The multi-dimensional array object on which
   *         the iterator will act upon
   *
   *  \param index [in] A Teuchos::ArrayView that specifies where the
   *         internal index of the iterator should start.
   *
   *  Produces an iterator with index corresponding to the given
   *  index.
   */
  MDIterator(const MDARRAY & mdarray,
             const Teuchos::ArrayView< dim_type > & index);

  /** \brief Copy constructor
   *
   * \param source [in] The source <tt>MDIterator</tt> to be copied
   */
  MDIterator(const MDIterator< MDARRAY > & source);

  /** \brief Destructor
   */
  ~MDIterator();

  //@}

  /** \name Standard Operators */
  //@{

  /** \brief Assignment operator
   *
   *  \param source [in] Source iterator for assignment
   */
  MDIterator< MDARRAY > & operator=(const MDIterator< MDARRAY > & source);

  /** \brief Equality operator
   *
   *  \param other [in] Iterator to be compared to
   */
  bool operator==(const MDIterator< MDARRAY > & other) const;

  /** \brief Inequality operator
   *
   *  \param other [in] Iterator to be compared to
   */
  bool operator!=(const MDIterator< MDARRAY > & other) const;

  /** \brief Dereferencing operator */
  inline value_type & operator*();

  /** \brief Dereferencing arrow operator */
  inline pointer operator->() const;

  /** \brief Prefix increment operator */
  MDIterator< MDARRAY > & operator++();

  /** \brief Postfix increment operator */
  MDIterator< MDARRAY > operator++(int);

  /** \brief Prefix decrement operator */
  MDIterator< MDARRAY > & operator--();

  /** \brief Postfix decrement operator */
  MDIterator< MDARRAY > operator--(int);

  //@}

  /** \brief Return the current index value along the given axis
   *
   *  \param axis [in] Requested axis for index value
   */
  inline dim_type index(int axis) const;

  /** \brief Stream output operator
   */
  template< typename T2 >
  friend std::ostream & operator<<(std::ostream & os, const MDIterator< T2 > & a);

private:

  // A copy of the dimensions of the multi-dimensional array being
  // iterated
  const Teuchos::Array< dim_type > _dimensions;

  // A copy of the strides of the multi-dimensional array being
  // iterated
  const Teuchos::Array< size_type > _strides;

  // A pointer to the data buffer of the multi-dimensional array
  // being iterated
  value_type * _ptr;

  // A copy of the storage order of the multi-dimensional array being
  // iterated
  Layout _layout;

  // The multi-dimensional index of the current iterate
  Teuchos::Array< dim_type > _index;

  // A temporary value used to indicate the axis of the index
  // currently being incremented or decremented
  mutable int _axis;

  // A temporary value used to indicate whether an increment or
  // decrement operation is complete
  mutable bool _done;

  // We need an index that is recognized as the end index.  It must
  // not be a valid index for the MDARRAY.  Since there are a nearly
  // infinite number of indexes that could serve as the end index,
  // this method should always be used to assign the index of an end
  // iterator.
  void assign_end_index();

  // Assert that the given index is valid for the given axis
  void assert_index(dim_type i, int axis) const;

};

/////////////////////
// Implementations //
/////////////////////

template< class MDARRAY >
MDIterator< MDARRAY >::MDIterator(const MDARRAY & mdarray,
                                  bool end_index) :
  _dimensions(mdarray._dimensions),
  _strides(mdarray._strides),
  _ptr(mdarray._ptr),
  _layout(mdarray._layout),
  _index(mdarray.numDims())
{
  if (end_index)
    assign_end_index();
  else
  {
    if (computeSize(_dimensions) == 0)
      assign_end_index();
    else
      _index.assign(_dimensions.size(), 0);
  }
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
MDIterator< MDARRAY >::
MDIterator(const MDARRAY & mdarray,
           const Teuchos::ArrayView< dim_type > & index) :
  _dimensions(mdarray._dimensions),
  _strides(mdarray._strides),
  _ptr(mdarray._ptr),
  _layout(mdarray._layout),
  _index(index)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != _index.size()), RangeError,
    "Input array has " << _dimensions.size() << " dimensions, while index "
    "has " << _index.size());
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  for (_axis = 0; _axis < _index.size(); ++_axis)
    assert_index(_index[_axis], _axis);
#endif
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
MDIterator< MDARRAY >::MDIterator(const MDIterator< MDARRAY > & source) :
  _dimensions(source._dimensions),
  _strides(source._strides),
  _ptr(source._ptr),
  _layout(source._layout),
  _index(source._index)
{
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
MDIterator< MDARRAY >::~MDIterator()
{
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
MDIterator< MDARRAY > &
MDIterator< MDARRAY >::operator=(const MDIterator< MDARRAY > & source)
{
  _dimensions = source._dimensions;
  _strides    = source._strides;
  _ptr        = source._ptr;
  _layout     = source._layout;
  _index      = source._index;
  return *this;
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
bool
MDIterator< MDARRAY >::operator==(const MDIterator< MDARRAY > & other) const
{
  // If underlying MDARRAYs are different, then return not equal
  if (_ptr != other._ptr) return false;
  // If any of the current index values differ, then return not equal 
  for (_axis = 0; _axis < _index.size(); _axis++)
    if (_index[_axis] != other._index[_axis]) return false;
  // Return equal
  return true;
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
bool
MDIterator< MDARRAY >::operator!=(const MDIterator< MDARRAY > & other) const
{
  return !(*this == other);
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
typename MDIterator< MDARRAY >::value_type &
MDIterator< MDARRAY >::operator*()
{
  size_type offset = 0;
  for (_axis=0; _axis < _index.size(); ++_axis)
    offset += _index[_axis] * _strides[_axis];
  return _ptr[offset];
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
typename MDIterator< MDARRAY >::pointer
MDIterator< MDARRAY >::operator->() const
{
  return &operator*();
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
MDIterator< MDARRAY > &
MDIterator< MDARRAY >::operator++()
{
  if (_layout == FIRST_INDEX_FASTEST)
  {
    _axis = 0;
    _done = false;
    while (not _done)
    {
      _index[_axis]++;
      _done = (_index[_axis] < _dimensions[_axis]);
      if (not _done)
      {
        _index[_axis] = 0;
        _axis++;
        if (_axis >= _index.size())
        {
          _done = true;
          assign_end_index();
        }
      }
    }
  }
  else
  {
    _axis = _dimensions.size() - 1;
    _done = false;
    while (not _done)
    {
      _index[_axis]++;
      _done = (_index[_axis] < _dimensions[_axis]);
      if (not _done)
      {
        _index[_axis] = 0;
        _axis--;
        if (_axis < 0)
        {
          _done = true;
          assign_end_index();
        }
      }
    }
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
MDIterator< MDARRAY >
MDIterator< MDARRAY >::operator++(int)
{
  MDIterator< MDARRAY > result(*this);
  ++(*this);
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
MDIterator< MDARRAY > &
MDIterator< MDARRAY >::operator--()
{
  if (_layout == FIRST_INDEX_FASTEST)
  {
    _axis = 0;
    _done = false;
    while (not _done)
    {
      _index[_axis]--;
      _done = (_index[_axis] >= 0);
      if (not _done)
      {
        _index[_axis] = _dimensions[_axis] - 1;
        _axis++;
        if (_axis >= _index.size())
        {
          _done = true;
          assign_end_index();
        }
      }
    }
  }
  else
  {
    _axis = _dimensions.size() - 1;
    _done = false;
    while (not _done)
    {
      _index[_axis]--;
      _done = (_index[_axis] >= 0);
      if (not _done)
      {
        _index[_axis] = _dimensions[_axis] - 1;
        _axis--;
        if (_axis < 0)
        {
          _done = true;
          assign_end_index();
        }
      }
    }
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
MDIterator< MDARRAY >
MDIterator< MDARRAY >::operator--(int)
{
  MDIterator< MDARRAY > result(*this);
  --(*this);
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
dim_type
MDIterator< MDARRAY >::
index(int axis) const
{
  return _index[axis];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
std::ostream & operator<<(std::ostream & os, const MDIterator< T > & a)
{
  os << &(*a);
  return os;
}

/////////////////////////////
// Private implementations //
/////////////////////////////

template< class MDARRAY >
void
MDIterator< MDARRAY >::assign_end_index()
{
  // We choose the end index to be equal to the MDARRAY dimensions,
  // where each index value is one greater than the largest valid
  // index for that axis.
  for (int axis = 0; axis < _index.size(); ++axis)
    _index[axis] = _dimensions[axis];
}

////////////////////////////////////////////////////////////////////////

template< class MDARRAY >
void
MDIterator< MDARRAY >::
assert_index(dim_type i,
             int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= i && i < _dimensions[axis]), RangeError,
    "MDIterator<MDARRAY>::assert_index(i=" << i << ",axis=" << axis << "): out"
    << " of range i in [0, " << _dimensions[axis] << ")"
  );
}

}  // namespace Domi

#endif
