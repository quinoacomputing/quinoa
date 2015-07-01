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

#ifndef SPARSE_COO_PTR_ELEMENT_H
#define SPARSE_COO_PTR_ELEMENT_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Sparse pointer element type for a COO matrix (val, ivect, jvect).
 *
 * This class abstracts a sparse element of a templated
 * type from a coordinate matrix. It
 * has a pointer to the value of the element.
 *
 * The default assignment operator and copy constructor
 * are allowed.
 */
template <class T_Index, class T_Value>
class SparseCOOPtrElement {
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

  /// Construct uninitialized (poiner to value set to zero) (#index() == 0#).
  SparseCOOPtrElement() : pvalue_(0), row_i_(0), col_j_(0)
  {}

  /// Construct with a pointer to the value and index set
  SparseCOOPtrElement(value_type* pvalue, index_type row_i, index_type col_j)
    : pvalue_(pvalue), row_i_(row_i), col_j_(col_j)
  {}

  /// Initialize
  void initialize(value_type* pvalue, index_type row_i, index_type col_j) {
    pvalue_	= pvalue;
    row_i_	= row_i;
    col_j_	= col_j;
  }
  
  //@}

  /** @name Value and index access */
  //@{ 

  /** \brief . */
  value_type& value()
  {
    return *pvalue_;
  }
  /** \brief . */
  value_type value() const
  {
    return *pvalue_;
  }
  /** \brief . */
  index_type row_i() const
  {
    return row_i_;
  }
  /** \brief . */
  index_type col_j() const
  {
    return col_j_;
  }
  /// Change the indexs
  void change_indexes(index_type row_i, index_type col_j)
  {
    row_i_ = row_i;
    col_j_ = col_j;
  }
  /// Change the element pointer
  void change_value_ptr(value_type* pvalue)
  {
    pvalue_ = pvalue;
  }

  //@}
private:
  value_type*				pvalue_;
  index_type				row_i_, col_j_;

};	// end class SparseCOOPtrElement

} // end namespace AbstractLinAlgPack 

#endif // SPARSE_COO_PTR_ELEMENT_H
