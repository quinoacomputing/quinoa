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

#ifndef IVECTOR_H
#define IVECTOR_H

#include <assert.h>

#include <valarray>

#include "DenseLinAlgPack_Types.hpp"
#include "Teuchos_Assert.hpp"

namespace DenseLinAlgPack {
/** \brief . */
/* * Fortran compatable integer vector for holding the  pivot information for
 * the elements of a vector, or the rows or columns of a matrix.
 */
class IVector : public std::valarray<DenseLinAlgPack::size_type> {
public:

  // STL typedefs
  typedef DenseLinAlgPack::index_type		value_type;
  typedef DenseLinAlgPack::size_type		size_type;
  typedef value_type&					reference;
  typedef const value_type&			const_reference;
  typedef value_type*					iterator;
  typedef const value_type*			const_iterator;
  typedef std::valarray<size_type>	valarray;

  // constructors

  /** \brief . */
  IVector();
  /** \brief . */
  IVector(size_type n);
  /** \brief . */
  IVector(const value_type& val, size_type n);
  /** \brief . */
  IVector(const value_type* p, size_type n);

  /// Resize on assignment
  IVector& operator=(const IVector&);

  /// 1-based element access (range checked if TEUCHOS_DEBUG is defined)
  reference operator()(size_type i);
  /// 1-based element access (range checked if TEUCHOS_DEBUG is defined)
  const_reference operator()(size_type i) const;

  /// STL iterator
  iterator begin();
  /// STL iterator
  const_iterator begin() const;
  /// STL iterator
  iterator end();
  /// STL iterator
  const_iterator end() const;

}; // end class IVector

// Inline definitions

inline IVector::IVector() : std::valarray<size_type>()
{}

inline IVector::IVector(size_type n) : std::valarray<size_type>(n)
{}

inline IVector::IVector(const value_type& val, size_type n) : std::valarray<size_type>(val,n)
{}

inline IVector::IVector(const value_type* p, size_type n) : std::valarray<size_type>(p,n)
{}

inline IVector& IVector::operator=(const IVector& iv)
{
  this->resize(iv.size());
  std::valarray<DenseLinAlgPack::size_type>::operator=(iv);
  return *this;
}

inline IVector::reference IVector::operator()(size_type i)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= i && i <= static_cast<size_type>(size())  ) );
#endif
  return operator[](i-1);
}

inline IVector::const_reference IVector::operator()(size_type i) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= i && i <= static_cast<size_type>(size())  ) );
#endif
  return const_cast<IVector*>(this)->operator[](i-1);
}

inline IVector::iterator IVector::begin()
{	return &operator[](0); }

inline IVector::const_iterator IVector::begin() const
{	return &(const_cast<IVector*>(this)->operator[](0)); }

inline IVector::iterator IVector::end()
{	return begin() + size(); }

inline IVector::const_iterator IVector::end() const
{	return begin() + size(); }

}	// end namespace DenseLinAlgPack

#endif // IVECTOR_H
