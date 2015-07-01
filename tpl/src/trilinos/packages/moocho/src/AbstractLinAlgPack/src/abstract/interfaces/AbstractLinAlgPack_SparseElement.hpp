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

#ifndef SPARSE_ELEMENT_H
#define SPARSE_ELEMENT_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Sparse storage element type.
  *
  * This class abstracts a sparse element of a templated
  * type.  It is ment to be used in a sparse vector.  Objects of
  * this type are designed so that the size of the object is
  * the same at least two value_type objects.
  *
  * The default assignment operator and copy constructor
  * are allowed.
  */
template <class T_Index, class T_Value>
class SparseElement {
public:
  /** @name Public Typedefs. */
  //@{

  /** \brief . */
  typedef T_Value						value_type;
  /** \brief . */
  typedef T_Index						index_type;

  //@}

  /** @name Constructors */
  //@{

  /// Construct uninitialized (#value() == 0.0#, #index() == 0#).
  SparseElement()
  {
    idx_pad_.index_  = 0;
    value_           = 0.0;
  }

  /// Construct with a value and index set
  SparseElement(index_type index, value_type value)
  {
    idx_pad_.index_ = index;
    value_          = value;
  }
  
  //@}

  /** @name Value and index access */
  //@{ 

  /** \brief . */
  value_type& value()
  {
    return value_;
  }
  /** \brief . */
  const value_type& value() const
  {
    return value_;
  }
  /** \brief . */
  const index_type& index() const
  {
    return idx_pad_.index_;
  }
  /// Initialize
  void initialize(index_type index, value_type value) {
    idx_pad_.index_ = index;
    value_          = value;
  }	
  /// Change the index
  void change_index(index_type index)
  {
    idx_pad_.index_ = index;
  }

  //@}

private:
  union index_and_padding {
    value_type  dummy;     // This is just included for alignment
    index_type index_;   // so that sizeof(this) == 2*sizeof(value_type)
  };
  index_and_padding		idx_pad_;
  value_type				value_;

};	// end class SparseElement

} // end namespace AbstractLinAlgPack 

#endif // SPARSE_ELEMENT_H
