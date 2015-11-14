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

#ifndef GEN_PERM_MATRIX_SLICE_ITERATOR_H
#define GEN_PERM_MATRIX_SLICE_ITERATOR_H

#include <assert.h>

#include <iterator>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

namespace GenPermMatrixSliceIteratorPack {

/** \brief . */
enum EOrderedBy { BY_ROW, BY_COL, BY_ROW_AND_COL, UNORDERED };

/** \brief External storage of a row and column indice.
  * This is required for creating a temporary in an assignment operation
  * in a sorting algorithm (like std::sort(...)).
  */
template< class T >
class external_row_col_value_type {
public:
  /** \brief . */
  typedef T			index_type;
  /** \brief . */
  typedef ptrdiff_t	difference_type;
  /** \brief . */
  external_row_col_value_type(
      difference_type	row_off
    , difference_type	col_off
    , index_type		row_i
    , index_type		col_j
    )
  :
    row_off_(row_off), col_off_(col_off), row_i_(row_i), col_j_(col_j)
  {}
  difference_type	row_off_;
  difference_type	col_off_;
  index_type		row_i_;
  index_type		col_j_;
};

/** \brief Internal storage for the iterator of the
  * row and column indices.
  */
template< class T >
class row_col_value_type {
public:
  /** \brief . */
  typedef T			index_type;
  /** \brief . */
  typedef ptrdiff_t	difference_type;
  /** \brief . */
  row_col_value_type( 
      difference_type	row_off
    , difference_type	col_off
    , index_type		row_i[]
    , index_type		col_j[]
    , size_type			nz
    );
  /** \brief . */
  void bind_view( const row_col_value_type<T>& val );
  /** \brief . */
  void increment(difference_type);
  /** \brief . */
  index_type	row_i() const;
  /** \brief . */
  index_type	col_j() const;
  /// May be NULL
  index_type* row_i_ptr() const;
  /** \brief . */
  row_col_value_type<T>& operator=( const row_col_value_type<T>& val );
  /** \brief . */
  operator const external_row_col_value_type<T>() const
  {
    return external_row_col_value_type<T>(row_off_,col_off_,*row_i_,*col_j_);
  }
  /** \brief . */
  row_col_value_type<T>& operator=( const external_row_col_value_type<T>& val )
  {
    TEUCHOS_TEST_FOR_EXCEPT( !(  row_off_ == val.row_off_  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  col_off_ == val.col_off_  ) );
    *row_i_ = val.row_i_;
    *col_j_ = val.col_j_;
    return *this;
  }
  
private:
  difference_type	row_off_;
  difference_type	col_off_;
  index_type		*row_i_;
  index_type		*col_j_;
  size_type		nz_;
  int				k_;	// zero based
  /** \brief . */
  void assert_in_range() const;
  /// Not defined and not to be called
  row_col_value_type();

};	// end class row_col_value_type

/// Swap row_col_value_type<T> objects
template< class T >
inline
void swap( row_col_value_type<T>& v1, row_col_value_type<T>& v2 )
{
  row_col_value_type<T> tmp = v1;
  v1 = v2;
  v2 = tmp;
}

/** \brief This is a full random access iterator for accessing row and colunmn
  * indices.
  */
template< class T >
class row_col_iterator
#if defined(_WINDOWS) || defined(_INTEL_CXX) || defined(_PG_CXX) 
  : public std::iterator< std::random_access_iterator_tag, row_col_value_type<T>, ptrdiff_t >
#endif
{
public:
  /** \brief . */
  typedef T								index_type;
  /** \brief . */
  typedef	std::random_access_iterator_tag	iterator_category;
  /** \brief . */
  typedef	row_col_value_type<T>	value_type;
  /** \brief . */
  typedef row_col_value_type<T>&			reference;
  /** \brief . */
  typedef row_col_value_type<T>*			pointer;
  /** \brief . */
  typedef	ptrdiff_t						difference_type;
  /// Null pointer!
  row_col_iterator();
  /** \brief . */
  row_col_iterator(
     difference_type	row_off
    ,difference_type	col_off
    ,index_type		row_i[]
    ,index_type		col_j[]
    ,size_type			nz			// Number of elements in row_i[] and col_j[]
    );
  /** \brief . */
  row_col_iterator<T>& operator=( const row_col_iterator<T>& itr );
  /** \brief . */
  reference operator*();
  /** \brief . */
  reference operator*() const;
  /** \brief . */
  pointer operator->() const;
  /// itr + a
  row_col_iterator<T>	operator+(difference_type) const;
  /// itr - a
  row_col_iterator<T>	operator-(difference_type);
  /// itr += a
  row_col_iterator<T>&	operator+=(difference_type);
  /// itr -= a
  row_col_iterator<T>&	operator-=(difference_type);
  /// ++itr
  row_col_iterator<T>&	operator++();
  /// itr++
  const row_col_iterator<T> operator++(int);
  /// --itr
  row_col_iterator<T>&	operator--();
  /// itr--
  const row_col_iterator<T> operator--(int);
  /// Difference
  difference_type operator-(const row_col_iterator<T>& itr) const;
  /// itr1 < itr2
  bool operator<( const row_col_iterator<T>& itr) const;
  /// itr1 <= itr2
  bool operator<=( const row_col_iterator<T>& itr) const;
  /// itr1 > itr 2
  bool operator>( const row_col_iterator<T>& itr) const;
  /// itr1 >= itr2
  bool operator>=( const row_col_iterator<T>& itr) const;
  /// itr1 == itr2
  bool operator==( const row_col_iterator<T>& itr) const;
  /// itr1 != itr2
  bool operator!=( const row_col_iterator<T>& itr) const;
  /// !itr (check for null)
  bool operator!() const;
  
private:
  mutable row_col_value_type<T>	value_;
  
};	// end class row_col_iterator<T>

// //////////////////////////////////////////////////////////
// Inline members for row_col_value_type<T>

template<class T>
inline
row_col_value_type<T>::row_col_value_type( 
      difference_type	row_off
    , difference_type	col_off
    , index_type		row_i[]
    , index_type		col_j[]
    , size_type			nz
    )
  :
    row_off_(row_off)
    ,col_off_(col_off)
    ,row_i_(row_i)
    ,col_j_(col_j)
    ,nz_(nz)
    ,k_(0)
{}

template<class T>
inline
void row_col_value_type<T>::bind_view( const row_col_value_type<T>& val )
{
    row_off_	= val.row_off_;
    col_off_	= val.col_off_;
    row_i_		= val.row_i_;
    col_j_		= val.col_j_;
    nz_			= val.nz_;
    k_			= val.k_;
}

template< class T >
inline
void row_col_value_type<T>::increment(difference_type d)
{
  row_i_	+= d;
  col_j_	+= d;
  k_ 		+= d;
}

template< class T >
inline
typename row_col_value_type<T>::index_type row_col_value_type<T>::row_i() const
{
  assert_in_range();
  return *row_i_ + row_off_;
}

template< class T >
inline
typename row_col_value_type<T>::index_type row_col_value_type<T>::col_j() const
{
  assert_in_range();
  return *col_j_ + col_off_;
}

template< class T >
inline
typename row_col_value_type<T>::index_type* row_col_value_type<T>::row_i_ptr() const
{
  return row_i_;
}

template< class T >
inline
row_col_value_type<T>& row_col_value_type<T>::operator=(
  const row_col_value_type<T>& val )
{
  *row_i_ = *val.row_i_;
  *col_j_ = *val.col_j_;
  return *this;
}

template< class T >
inline
void row_col_value_type<T>::assert_in_range() const
{
  // ToDo: Finish this!
  TEUCHOS_TEST_FOR_EXCEPT( !(  0 <= k_ && k_ < nz_  ) );
}

/// Assert not null
void GPMS_row_col_iterator_assert_not_null(const void* p);

// //////////////////////////////////////////////////////////
// Inline members for row_col_iterator<T>

template< class T >
inline
row_col_iterator<T>::row_col_iterator()
  :
    value_(0,0,NULL,NULL,0)
{}

template< class T >
inline
row_col_iterator<T>::row_col_iterator(
     difference_type	row_off
    ,difference_type	col_off
    ,index_type         row_i[]
    ,index_type         col_j[]
    ,size_type			nz
    )
  :
    value_(row_off,col_off,row_i,col_j,nz)
{}

template< class T >
inline
row_col_iterator<T>& row_col_iterator<T>::operator=( const row_col_iterator<T>& itr )
{
  value_.bind_view( itr.value_ );
  return *this;
}

template< class T >
inline
typename row_col_iterator<T>::reference
row_col_iterator<T>::operator*()
{
  GPMS_row_col_iterator_assert_not_null(value_.row_i_ptr());
  return value_;
}


template< class T >
inline
typename row_col_iterator<T>::reference
row_col_iterator<T>::operator*() const
{
  GPMS_row_col_iterator_assert_not_null(value_.row_i_ptr());
  return value_;
}

template< class T >
inline
typename row_col_iterator<T>::pointer
row_col_iterator<T>::operator->() const
{
  GPMS_row_col_iterator_assert_not_null(value_.row_i_ptr());
  return &value_;
}

template< class T >
inline
row_col_iterator<T>
row_col_iterator<T>::operator+(difference_type d) const
{
  row_col_iterator<T> itr = *this;
  itr.value_.increment(d);
  return itr;
}

template< class T >
inline
row_col_iterator<T>
row_col_iterator<T>::operator-(difference_type d)
{
  row_col_iterator<T> itr = *this;
  itr.value_.increment(-d);
  return itr;
}

template< class T >
inline
row_col_iterator<T>&
row_col_iterator<T>::operator+=(difference_type d)
{
  value_.increment(d);
  return *this;
}

template< class T >
inline
row_col_iterator<T>&
row_col_iterator<T>::operator-=(difference_type d)
{
  value_.increment(-d);
  return *this;
}

template< class T >
inline
row_col_iterator<T>&
row_col_iterator<T>::operator++()
{
  value_.increment(1);
  return *this;
}

template< class T >
inline
const row_col_iterator<T>
row_col_iterator<T>::operator++(int)
{
  row_col_iterator<T> itr = *this;
  value_.increment(1);
  return itr;
}

template< class T >
inline
row_col_iterator<T>&
row_col_iterator<T>::operator--()
{
  value_.increment(-1);
  return *this;
}

template< class T >
inline
const row_col_iterator<T>
row_col_iterator<T>::operator--(int)
{
  row_col_iterator<T> itr = *this;
  value_.increment(-1);
  return itr;
}

template< class T >
inline
typename row_col_iterator<T>::difference_type
row_col_iterator<T>::operator-(const row_col_iterator<T>& itr) const
{
  return value_.row_i_ptr() - itr.value_.row_i_ptr();
}

template< class T >
inline
bool row_col_iterator<T>::operator<( const row_col_iterator<T>& itr) const
{
  return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
      &&  ( value_.row_i_ptr() < itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator<=( const row_col_iterator<T>& itr) const
{
  return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
      &&  ( value_.row_i_ptr() <= itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator>( const row_col_iterator<T>& itr) const
{
  return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
      &&  ( value_.row_i_ptr() > itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator>=( const row_col_iterator<T>& itr) const
{
  return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
      &&  ( value_.row_i_ptr() >= itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator==( const row_col_iterator<T>& itr) const
{
  return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
      &&  ( value_.row_i_ptr() == itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator!=( const row_col_iterator<T>& itr) const
{
  return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
      &&  ( value_.row_i_ptr() != itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator!() const
{
  return  value_.row_i_ptr() == NULL;
}

}	// end namespace GenPermMatrixSliceIteratorPack

}	// end namespace AbstractLinAlgPack

#endif   // GEN_PERM_MATRIX_SLICE_ITERATOR_H
