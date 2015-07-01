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

#ifndef SLAP_MATRIX_COOR_TMPL_ITFC_H
#define SLAP_MATRIX_COOR_TMPL_ITFC_H

#include <stdexcept>

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_Assert.hpp"

namespace AbstractLinAlgPack {

template<class T_Scalar, class T_Index>
class MatrixCOORTmplItfcItrEleView;

template<class T_Scalar, class T_Index>
class MatrixCOORTmplItfcItr;

/** \brief Templated class that supports the \c COOMatrixTemplateInterface
 * template interface.
 */
template<class T_Scalar, class T_Index>
class MatrixCOORTmplItfc {
public:
  typedef T_Index                                          size_type;
  typedef ptrdiff_t                                        difference_type;
  typedef MatrixCOORTmplItfcItrEleView<T_Scalar,T_Index>   element_type;
  typedef T_Scalar                                         value_type;
  typedef T_Index                                          index_type;
  typedef MatrixCOORTmplItfcItr<T_Scalar,T_Index>          const_iterator;
  MatrixCOORTmplItfc(
    size_type rows, size_type cols, size_type nz
    ,difference_type row_offset, difference_type col_offset
    ,const T_Scalar *values, const T_Index* row_i, const T_Index* col_j
    )
    :rows_(rows), cols_(cols), nz_(nz), row_offset_(row_offset), col_offset_(col_offset)
    ,values_(values), row_i_(row_i), col_j_(col_j)
  {}
  size_type       rows()       const { return rows_;       }
  size_type       cols()       const { return cols_;       }
  size_type       nz()         const { return nz_;         }
  difference_type row_offset() const { return row_offset_; }	
  difference_type col_offset() const { return col_offset_; }
  const_iterator  begin() const;
  const_iterator  end() const;
private:
  size_type        rows_;
    size_type        cols_;
  size_type        nz_;
  difference_type  row_offset_;
  difference_type  col_offset_;
  const T_Scalar   *values_;
  const T_Index    *row_i_;
  const T_Index    *col_j_;
  // Not defined and not to be called
  MatrixCOORTmplItfc();
}; // end class MatrixCOORTmplItfc


// ///////////////////////////////////////////////
// Implementatioins, not of the user to look at!

/** \brief Templated class for objects that support the
 * \c SparseCOOElementTemplatInterface specification.
 */
template<class T_Scalar, class T_Index>
class MatrixCOORTmplItfcItrEleView {
public:
  typedef T_Scalar   value_type;
  typedef T_Index    index_type;
  MatrixCOORTmplItfcItrEleView(
    const T_Scalar* value, const T_Index* row_i, const T_Index* col_j
    )
    : value_(value), row_i_(row_i), col_j_(col_j)
  {}
  void increment()       { ++value_; ++row_i_; ++col_j_; }
  bool operator!=(const MatrixCOORTmplItfcItrEleView<T_Scalar,T_Index>& ele) const
  {   return value_ != ele.value_ || row_i_ != ele.row_i_ || col_j_ != ele.col_j_; }
  T_Scalar value() const { return *value_;  }
  T_Index  row_i() const { return *row_i_;  }
  T_Index  col_j() const { return *col_j_;  }
private:
  const T_Scalar  *value_;
  const T_Index   *row_i_;
  const T_Index   *col_j_;
  // Not defined and not to be called
  MatrixCOORTmplItfcItrEleView();
}; // end class MatrixCOORTmplItfcItrEleView

/** \brief Templated class for iterator returning objects that support the
 * \c SparseCOOElementTemplatInterface specification.
 */
template<class T_Scalar, class T_Index>
class MatrixCOORTmplItfcItr {
public:
  MatrixCOORTmplItfcItr( const T_Scalar* value, const T_Index* row_i, const T_Index* col_j, T_Index nz )
    : ele_(value,row_i,col_j)
#ifdef TEUCHOS_DEBUG
    , nz_left_(nz)
#endif
  {}
  void operator++() {
#ifdef TEUCHOS_DEBUG
    --nz_left_;
#endif
    ele_.increment();
  }
  bool operator!=(const MatrixCOORTmplItfcItr<T_Scalar,T_Index>& itr) const
  {   return ele_ != itr.ele_; }
  const MatrixCOORTmplItfcItrEleView<T_Scalar,T_Index>* operator->() const
  {	assert_nz(); return &ele_; }
private:
  MatrixCOORTmplItfcItrEleView<T_Scalar,T_Index>  ele_;
#ifdef TEUCHOS_DEBUG
  T_Index   nz_left_;
  void assert_nz() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      nz_left_ <= 0, std::logic_error
      ,"MatrixCOORTmplItfcItr<>::assert_nz: Error, trying to access past storage!" );
  }
#else
  void assert_nz() const {}
#endif
  // Not defined and not to be called
  MatrixCOORTmplItfcItr();
}; // end class MatrixCOORTmplItfcItr

// ///////////////////////////
// Inline members

// MatrixCOORTmplItfc

template<class T_Scalar, class T_Index>
inline
typename MatrixCOORTmplItfc<T_Scalar,T_Index>::const_iterator
MatrixCOORTmplItfc<T_Scalar,T_Index>::begin() const
{
  return MatrixCOORTmplItfcItr<T_Scalar,T_Index>(values_,row_i_,col_j_,nz_);
}

template<class T_Scalar, class T_Index>
inline
typename MatrixCOORTmplItfc<T_Scalar,T_Index>::const_iterator
MatrixCOORTmplItfc<T_Scalar,T_Index>::end() const
{
  return MatrixCOORTmplItfcItr<T_Scalar,T_Index>(values_+nz_,row_i_+nz_,col_j_+nz_,0);
}

} // end namespace AbstractLinAlgPack

#endif // SLAP_MATRIX_COOR_TMPL_ITFC_H
