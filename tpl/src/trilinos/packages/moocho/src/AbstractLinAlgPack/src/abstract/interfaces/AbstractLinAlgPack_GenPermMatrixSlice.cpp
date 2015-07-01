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

#include <assert.h>

#include <sstream>
#include <functional>
#include <algorithm>

#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "Teuchos_Assert.hpp"

#ifdef _WINDOWS

namespace std {

// Help some compilers lookup the function.
inline
void swap(
    AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::row_col_value_type<
      DenseLinAlgPack::size_type>& v1
  , AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::row_col_value_type<
    DenseLinAlgPack::size_type>& v2
  )
{
  AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::swap(v1,v2);
}

}	// end namespace std

#endif // _WINDOWS

namespace {

//
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }

// Return a string with the same name as the enumeration
const char* ordered_by_str(
  AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::EOrderedBy ordered_by )
{
  switch( ordered_by ) {
    case AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::BY_ROW:
      return "BY_ROW";
    case AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::BY_COL:
      return "BY_COL";
    case AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::BY_ROW_AND_COL:
      return "BY_ROW_AND_COL";
    case AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::UNORDERED:
      return "UNORDERED";
  }
  TEUCHOS_TEST_FOR_EXCEPT(true);	// should never be executed
  return 0;
}

// Define function objects for comparing by row and by column

template<class T>
struct imp_row_less
   : public std::binary_function<
     AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::row_col_value_type<T>
     ,AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::row_col_value_type<T>
     ,bool
     >
{
  bool operator()(
      const AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::external_row_col_value_type<T>& v1
    , const AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::external_row_col_value_type<T>& v2
    )
  {
    return v1.row_i_ < v2.row_i_;
  }
};

template<class T>
struct imp_col_less
   : public std::binary_function<
     AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::row_col_value_type<T>
     ,AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::row_col_value_type<T>
     ,bool
     >
{
  bool operator()(
      const AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::external_row_col_value_type<T>& v1
    , const AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::external_row_col_value_type<T>& v2
    )
  {
    return v1.col_j_ < v2.col_j_;
  }
};

}	// end namespace

//#ifdef _WINDOWS
//namespace std {
//	using imp_row_less;
//}
//#endif // _WINDOWS

namespace AbstractLinAlgPack {

GenPermMatrixSlice::GenPermMatrixSlice()
  : rows_(0), cols_(0), nz_(0)
{}

void GenPermMatrixSlice::initialize( size_type rows, size_type cols, EIdentityOrZero type )
{
  rows_       = rows;
  cols_       = cols;
  nz_         = type == IDENTITY_MATRIX ? my_min(rows,cols) : 0;
  ordered_by_ = GenPermMatrixSliceIteratorPack::BY_ROW_AND_COL;
  row_i_      = NULL;
  col_j_      = NULL;   // Don't need to be set but might as well for safely
  row_off_    = 0;      // ...
  col_off_    = 0;      // ...
}

void GenPermMatrixSlice::initialize(
  size_type			rows
  ,size_type			cols
  ,size_type			nz
  ,difference_type	row_off
  ,difference_type	col_off
  ,EOrderedBy			ordered_by
  ,const size_type	row_i[]
  ,const size_type	col_j[]
  ,bool				test_setup
  )
{
  namespace GPMSIP = GenPermMatrixSliceIteratorPack;

  if( test_setup ) {
    std::ostringstream omsg;
    omsg << "\nGenPermMatrixSlice::initialize(...) : Error: ";
    // Validate the input data (at least partially)
    validate_input_data(rows,cols,nz,row_off,col_off,ordered_by,row_i,col_j,omsg);
    // Validate the ordering and uniquness
    const size_type *ordered_sequence = NULL;
    if( ordered_by == GPMSIP::BY_ROW || ordered_by == GPMSIP::BY_ROW_AND_COL ) {
      for( size_type k = 1; k < nz; ++k ) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          row_i[k-1] >= row_i[k], std::invalid_argument
          ,"GenPermMatrixSlice::initialize(...) : Error: "
          "row_i[" << k-1 << "] = " << row_i[k-1]
          << " >= row_i[" << k << "] = " << row_i[k]
          << "\nThis is not sorted by row!" );
      }
    }
    if( ordered_by == GPMSIP::BY_COL || ordered_by == GPMSIP::BY_ROW_AND_COL ) {
      for( size_type k = 1; k < nz; ++k ) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          col_j[k-1] >= col_j[k], std::invalid_argument
          ,"GenPermMatrixSlice::initialize(...) : Error: "
          "col_j[" << k-1 << "] = " << col_j[k-1]
          << " >= col_j[" << k << "] = " << col_j[k]
          << "\nThis is not sorted by column!" );
      }
    }
  }
  // Set the members
  rows_		= rows;
  cols_		= cols;
  nz_			= nz;
  row_off_	= row_off;
  col_off_	= col_off;
  ordered_by_	= ordered_by;
  row_i_		= nz ? row_i : NULL;
  col_j_		= nz ? col_j : NULL;
}

void GenPermMatrixSlice::initialize_and_sort(
  size_type			rows
  ,size_type			cols
  ,size_type			nz
  ,difference_type	row_off
  ,difference_type	col_off
  ,EOrderedBy			ordered_by
  ,size_type			row_i[]
  ,size_type			col_j[]
  ,bool				test_setup
  )
{
  namespace GPMSIP = GenPermMatrixSliceIteratorPack;
  TEUCHOS_TEST_FOR_EXCEPTION(
    ordered_by == GPMSIP::BY_ROW_AND_COL, std::invalid_argument
    ,"GenPermMatrixSlice::initialize_and_sort(...) : Error, "
    "ordered_by == GPMSIP::BY_ROW_AND_COL, we can not sort by row and column!" );
  if( test_setup ) {
    std::ostringstream omsg;
    omsg << "\nGenPermMatrixSlice::initialize_and_sort(...) : Error:\n";
    // Validate the input data (at least partially)
    validate_input_data(rows,cols,nz,row_off,col_off,ordered_by,row_i,col_j,omsg);
  }

  // Sort by row or column
  typedef GPMSIP::row_col_iterator<size_type> row_col_itr_t;
  row_col_itr_t
    row_col_itr = row_col_itr_t( row_off, col_off, row_i, col_j, nz );
  if( ordered_by == GPMSIP::BY_ROW ) {
    std::stable_sort( row_col_itr, row_col_itr + nz
      , imp_row_less<size_type>() );
  }
  else if( ordered_by == GPMSIP::BY_COL ) {
    std::stable_sort( row_col_itr, row_col_itr + nz
      , imp_col_less<size_type>() );
  }

  initialize(rows,cols,nz,row_off,col_off,ordered_by,row_i,col_j,test_setup);
}
  
void GenPermMatrixSlice::bind( const GenPermMatrixSlice& gpms )
{
  this->initialize( gpms.rows_, gpms.cols_, gpms.nz_, gpms.row_off_
    , gpms.col_off_, gpms.ordered_by_, gpms.row_i_, gpms.col_j_
    , false );
}

size_type GenPermMatrixSlice::lookup_row_i(size_type col_j) const
{
  namespace QPMSIP = GenPermMatrixSliceIteratorPack;
  if( col_j < 1 || cols_ < col_j )
    std::out_of_range(
      "GenPermMatrixSlice::lookup_row_i(col_j) : Error, "
      "col_j is out of bounds" );
  if(!nz_)
    return 0;
  if(is_identity())
    return col_j <= nz_ ? col_j : 0;
  const size_type
    *itr = NULL;
  if( ordered_by() == QPMSIP::BY_COL || ordered_by() == QPMSIP::BY_ROW_AND_COL )
    itr = std::lower_bound( col_j_, col_j_ + nz_, col_j );
  else
    itr = std::find( col_j_, col_j_ + nz_, col_j );
  return (itr != col_j_ + nz_ && *itr == col_j) ? row_i_[itr - col_j_] : 0;
}

size_type GenPermMatrixSlice::lookup_col_j(size_type row_i) const
{
  namespace QPMSIP = GenPermMatrixSliceIteratorPack;
  if( row_i < 1 || rows_ < row_i )
    std::out_of_range(
      "GenPermMatrixSlice::lookup_col_j(row_i) : Error, "
      "row_i is out of bounds" );
  if(!nz_)
    return 0;
  if(is_identity())
    return row_i <= nz_ ? row_i : 0;
  const size_type
    *itr = NULL;
  if( ordered_by() == QPMSIP::BY_ROW || ordered_by() == QPMSIP::BY_ROW_AND_COL )
    itr = std::lower_bound( row_i_, row_i_ + nz_, row_i );
  else
    itr = std::find( row_i_, row_i_ + nz_, row_i );
  return (itr != row_i_ + nz_ && *itr == row_i) ? col_j_[itr - row_i_] : 0;
}

GenPermMatrixSlice::const_iterator GenPermMatrixSlice::begin() const
{
  validate_not_identity();
  return const_iterator(row_off_,col_off_,row_i_,col_j_,nz_);
}

GenPermMatrixSlice::const_iterator GenPermMatrixSlice::end() const
{
  validate_not_identity();
  return const_iterator(row_off_,col_off_,row_i_+nz_,col_j_+nz_,0);
}

const GenPermMatrixSlice GenPermMatrixSlice::create_submatrix(
  const Range1D& rng, EOrderedBy ordered_by ) const
{
  namespace GPMSIP = GenPermMatrixSliceIteratorPack;

  validate_not_identity();

  // Validate the input
  TEUCHOS_TEST_FOR_EXCEPTION(
    ordered_by == GPMSIP::BY_ROW_AND_COL, std::invalid_argument
    ,"GenPermMatrixSlice::initialize_and_sort(...) : Error, "
    "ordered_by == GPMSIP::BY_ROW_AND_COL, we can not sort by row and column!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    rng.full_range(), std::logic_error,
    "GenPermMatrixSlice::create_submatrix(...) : Error, "
    "The range argument can not be rng.full_range() == true" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    ordered_by == GPMSIP::BY_ROW && rng.ubound() > rows(), std::logic_error
    ,"GenPermMatrixSlice::create_submatrix(...) : Error, "
    "rng.ubound() can not be larger than this->rows()" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    ordered_by == GPMSIP::BY_COL && rng.ubound() > cols(), std::logic_error
    ,"GenPermMatrixSlice::create_submatrix(...) : Error, "
    "rng.ubound() can not be larger than this->cols()" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    ordered_by == GPMSIP::UNORDERED, std::logic_error
    ,"GenPermMatrixSlice::create_submatrix(...) : Error, "
    "You can have ordered_by == GPMSIP::UNORDERED" );

  // Find the upper and lower k for row[k], col[k] indices
  size_type
    k_l = 0,		// zero based 
    k_u = nz() + 1;	// zero based (== nz + 1 flag that there are no nonzeros to even search)
  size_type
    rows = 0,
    cols = 0;
  difference_type
    row_off = 0,
    col_off = 0;
  switch( ordered_by ) {
    case GPMSIP::BY_ROW:
    case GPMSIP::BY_COL:
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        this->ordered_by() != GPMSIP::BY_ROW_AND_COL
        && ( nz() > 1 && ordered_by != this->ordered_by() )
        ,std::logic_error
        ,"GenPermMatrixSlice::create_submatrix(...) : Error, "
        << "nz = " << nz() << " > 1 and "
        << "ordered_by = " << ordered_by_str(ordered_by)
        << " != this->ordered_by() = "
        << ordered_by_str(this->ordered_by()) );
      // Search the rows or the columns.
      const size_type
        *search_k = NULL;
      difference_type
        search_k_off;
      if( ordered_by == GPMSIP::BY_ROW ) {
        search_k = row_i_;	// may be null
        search_k_off = row_off_;
        rows = rng.size();
        cols = this->cols();
        row_off = row_off_ - (difference_type)(rng.lbound() - 1);
        col_off = col_off_;
      }
      else {	// BY_COL
        search_k = col_j_;	// may be null
        search_k_off = col_off_;
        rows = this->rows();
        cols = rng.size();
        row_off = row_off_;
        col_off = col_off_ - (difference_type)(rng.lbound() - 1);;
      }
      if( search_k ) {
        const size_type
          *l = std::lower_bound( search_k, search_k + nz()
              , rng.lbound() - search_k_off );
        k_l = l - search_k;
        // If we did not find the lower bound in the range, don't even bother
        // looking for the upper bound.
        if( k_l != nz() ) {
          const size_type
            *u = std::upper_bound( search_k, search_k + nz()
                , rng.ubound() - search_k_off );
          k_u = u - search_k;
          // Here, if there are any nonzero elements in this range then
          // k_u - k_l > 0 will always be true!
        }
      }
      break;
    }
    case GPMSIP::UNORDERED:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  GenPermMatrixSlice gpms;
  if( k_u - k_l > 0 && k_u != nz() + 1 ) {
    gpms.initialize(
        rows, cols
      , k_u - k_l
      , row_off, col_off
      , ordered_by
      , row_i_ + k_l
      , col_j_ + k_l
      );
  }
  else  {
    gpms.initialize(
        rows, cols
      , 0
      , row_off, col_off
      , ordered_by
      , NULL
      , NULL
      );
  }
  return gpms;
}

// static members

void GenPermMatrixSlice::validate_input_data(
  size_type			rows
  ,size_type			cols
  ,size_type			nz
  ,difference_type	row_off
  ,difference_type	col_off
  ,EOrderedBy			ordered_by
  ,const size_type	row_i[]
  ,const size_type	col_j[]
  ,std::ostringstream &omsg
  )
{
  namespace GPMSIP = GenPermMatrixSliceIteratorPack;

  TEUCHOS_TEST_FOR_EXCEPTION(
    nz > rows * cols, std::invalid_argument
    ,omsg.str() << "nz = " << nz << " can not be greater than rows * cols = "
    << rows << " * " << cols << " = " << rows * cols );
  
  // First see if everything is in range.
  for( size_type k = 0; k < nz; ++k ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      row_i[k] + row_off < 1 || rows < row_i[k] + row_off, std::invalid_argument
      ,omsg.str() << "row_i[" << k << "] + row_off = " << row_i[k] << " + " << row_off
      << " = " << (row_i[k] + row_off)
      << " is out of range [1,rows] = [1," << rows << "]" );
    TEUCHOS_TEST_FOR_EXCEPTION(
      col_j[k] + col_off < 1 || cols < col_j[k] + col_off, std::invalid_argument
      ,omsg.str() << "col_j[" << k << "] + col_off = " << col_j[k] << " + " << col_off
      << " = " << (col_j[k] + col_off)
      << " is out of range [1,cols] = [1," << cols << "]" );
  }
  
  // ToDo: Technically, we need to validate that the nonzero elements row_i[k], col_j[k] are
  // unique but this is much harder to do!
  
}

// private

void GenPermMatrixSlice::validate_not_identity() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    is_identity(), std::logic_error
    ,"GenPermMatrixSlice::validate_not_identity() : "
    "Error, this->is_identity() is true" );
}

}	// end namespace AbstractLinAlgPack
