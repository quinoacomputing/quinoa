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

#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

void AbstractLinAlgPack::V_StMtV(
    SpVector* y, value_type a, const GenPermMatrixSlice& P
  , BLAS_Cpp::Transp P_trans, const DVectorSlice& x )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  using DenseLinAlgPack::MtV_assert_sizes;

  MtV_assert_sizes( P.rows(), P.cols(), P_trans, x.dim() );

  y->resize( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans ), P.nz() );

  typedef SpVector::element_type ele_t;

  if( P.is_identity() ) {
    for( size_type i = 1; i <= P.nz(); ++i ) {
      const value_type x_j = x(i);
      if( x_j != 0.0 )
        y->add_element( ele_t( i, a * x_j ) );
    }
  }		
  else if( P_trans == no_trans ) {
    for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
      const size_type
        i = itr->row_i(),
        j = itr->col_j();
      const value_type x_j = x(j);
      if( x_j != 0.0 )
        y->add_element( ele_t( i, a * x_j ) );
    }
  }
  else {
    for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
      const size_type
        j = itr->row_i(),
        i = itr->col_j();
      const value_type x_j = x(j);
      if( x_j != 0.0 )
        y->add_element( ele_t( i, a * x_j ) );
    }
  }
  if( P.ordered_by() == GPMSIP::BY_ROW_AND_COL
    || ( P_trans == no_trans	&& P.ordered_by() == GPMSIP::BY_ROW )
    ||	( P_trans == trans		&& P.ordered_by() == GPMSIP::BY_COL )	)
  {
    y->assume_sorted(true);
  }
}

void AbstractLinAlgPack::V_StMtV(
    SpVector* y, value_type a, const GenPermMatrixSlice& P
  , BLAS_Cpp::Transp P_trans, const SpVectorSlice& x )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  using DenseLinAlgPack::MtV_assert_sizes;
  MtV_assert_sizes( P.rows(), P.cols(), P_trans, x.dim() );

  y->resize( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans ), P.nz() );

  typedef SpVector::element_type ele_t;
  const SpVectorSlice::element_type *ele_ptr;

  if( P.is_identity() ) {
    AbstractLinAlgPack::add_elements( y, 1.0, x(1,P.nz()) );
    AbstractLinAlgPack::Vt_S( &(*y)(), a );
  }		
  else if( x.is_sorted() ) {
    if( P_trans == no_trans ) {
      for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
        const size_type
          i = itr->row_i(),
          j = itr->col_j();
        if( ele_ptr = x.lookup_element(j) )
          y->add_element( ele_t( i, a * ele_ptr->value() ) );
      }
    }
    else {
      for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
        const size_type
          j = itr->row_i(),
          i = itr->col_j();
        if( ele_ptr = x.lookup_element(j) )
          y->add_element( ele_t( i, a * ele_ptr->value() ) );
      }
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: Implement the other cases!
  }
  if(	 P.ordered_by() == GPMSIP::BY_ROW_AND_COL
    || 	( P_trans == no_trans	&& P.ordered_by() == GPMSIP::BY_ROW )
    ||	( P_trans == trans		&& P.ordered_by() == GPMSIP::BY_COL )	)
  {
    y->assume_sorted(true);
  }
}

void AbstractLinAlgPack::Vp_StMtV(
    SpVector* y, value_type a, const GenPermMatrixSlice& P
  , BLAS_Cpp::Transp P_trans, const DVectorSlice& x )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  using DenseLinAlgPack::Vp_MtV_assert_sizes;

  Vp_MtV_assert_sizes( y->dim(), P.rows(), P.cols(), P_trans, x.dim() );

  typedef SpVector::element_type ele_t;

  const bool was_sorted = y->is_sorted();

  if( P.is_identity() ) {
    for( size_type i = 1; i <= P.nz(); ++i ) {
      const value_type x_j = x(i);
      if( x_j != 0.0 )
        y->add_element( ele_t( i, a * x_j ) );
    }
  }		
  else if( P_trans == no_trans ) {
    for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
      const size_type
        i = itr->row_i(),
        j = itr->col_j();
      const value_type x_j = x(j);
      if( x_j != 0.0 )
        y->add_element( ele_t( i, a * x_j ) );
    }
  }
  else {
    for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
      const size_type
        j = itr->row_i(),
        i = itr->col_j();
      const value_type x_j = x(j);
      if( x_j != 0.0 )
        y->add_element( ele_t( i, a * x_j ) );
    }
  }
  if( was_sorted &&
    ( P.ordered_by() == GPMSIP::BY_ROW_AND_COL
      || ( P_trans == no_trans	&& P.ordered_by() == GPMSIP::BY_ROW )
      ||	( P_trans == trans		&& P.ordered_by() == GPMSIP::BY_COL ) )	)
  {
    y->assume_sorted(true);
  }
}

void AbstractLinAlgPack::Vp_StMtV(
    DVectorSlice* y, value_type a, const GenPermMatrixSlice& P
  , BLAS_Cpp::Transp P_trans, const DVectorSlice& x, value_type b )
{
  using DenseLinAlgPack::Vt_S;
  using DenseLinAlgPack::Vp_MtV_assert_sizes;
  Vp_MtV_assert_sizes( y->dim(), P.rows(), P.cols(), P_trans, x.dim() );
  // y = b*y
  if( b == 0.0 )
    *y = 0.0;
  else
    Vt_S(y,b);	
  // y += a*op(P)*x
  if( P.is_identity() ) {
    if( b == 0.0 )
      *y = 0.0;
    else
      DenseLinAlgPack::Vt_S( y, b );
    DenseLinAlgPack::Vp_StV( &(*y)(1,P.nz()), a, x(1,P.nz()) );
  }		
  else if( P_trans == BLAS_Cpp::no_trans ) {
    for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
      const size_type
        i = itr->row_i(),
        j = itr->col_j();
      (*y)(i) += a * x(j);
    }
  }
  else {
    for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
      const size_type
        j = itr->row_i(),
        i = itr->col_j();
      (*y)(i) += a * x(j);
    }
  }
}

void AbstractLinAlgPack::Vp_StMtV(
    DVectorSlice* y, value_type a, const GenPermMatrixSlice& P
  , BLAS_Cpp::Transp P_trans, const SpVectorSlice& x, value_type b )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  using DenseLinAlgPack::Vt_S;
  using DenseLinAlgPack::Vp_MtV_assert_sizes;
  
  Vp_MtV_assert_sizes( y->dim(), P.rows(), P.cols(), P_trans, x.dim() );
  // y = b*y
  if( b == 0.0 )
    *y = 0.0;
  else
    Vt_S(y,b);
  // y += a*op(P)*x
  if( P.is_identity() ) {
    DenseLinAlgPack::Vt_S( y, b ); // takes care of b == 0.0 and y == NaN
    AbstractLinAlgPack::Vp_StV( &(*y)(1,P.nz()), a, x(1,P.nz()) );
  }		
  else if( x.is_sorted() ) {
    const SpVectorSlice::difference_type x_off = x.offset();
    if( P_trans == no_trans && P.ordered_by() == GPMSIP::BY_COL ) {
      TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: implement this!
    }
    else if( ( P_trans == trans && P.ordered_by() == GPMSIP::BY_ROW )
      || P.ordered_by() == GPMSIP::BY_ROW_AND_COL )
    {
      GenPermMatrixSlice::const_iterator
        P_itr = P.begin(),
        P_end = P.end();
      SpVectorSlice::const_iterator
        x_itr = x.begin(),
        x_end = x.end();
      while( P_itr != P_end && x_itr != x_end ) {
        const size_type
          i = rows(P_itr->row_i(),P_itr->col_j(),P_trans),
          j = cols(P_itr->row_i(),P_itr->col_j(),P_trans);
        if( j < x_itr->index() + x_off ) {
          ++P_itr;
          continue;
        }
        else if( j > x_itr->index() + x_off ) {
          ++x_itr;
          continue;
        }
        else {	// they are equal
          (*y)(i) += a * x_itr->value();
          ++P_itr;
          ++x_itr;
        }
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement what ever this case is?	
    }
  }
  else {
    // Since things do not match up we will have to create a
    // temporary dense copy of x to operate on.
    TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: Implement this!
  }
}

namespace {

AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::EOrderedBy
ordered_by(
  AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::EOrderedBy P_ordered_by
  , BLAS_Cpp::Transp P_trans
  )
{
  using BLAS_Cpp::no_trans;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  GPMSIP::EOrderedBy
    opP_ordered_by;
  switch( P_ordered_by ) {
      case GPMSIP::BY_ROW_AND_COL:
      opP_ordered_by = GPMSIP::BY_ROW_AND_COL;
      break;
      case GPMSIP::BY_ROW:
      opP_ordered_by = P_trans == no_trans ? GPMSIP::BY_ROW : GPMSIP::BY_COL;
      break;
      case GPMSIP::BY_COL:
      opP_ordered_by = P_trans == no_trans ? GPMSIP::BY_COL : GPMSIP::BY_COL;
      break;
      case GPMSIP::UNORDERED:
      opP_ordered_by = GPMSIP::UNORDERED;
      break;
        default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never happen
  }
  return opP_ordered_by;
}

} // end namespace

void AbstractLinAlgPack::intersection(
  const GenPermMatrixSlice     &P1
  ,BLAS_Cpp::Transp            P1_trans
  ,const GenPermMatrixSlice    &P2
  ,BLAS_Cpp::Transp            P2_trans
  ,size_type                   *Q_nz
  ,const size_type             Q_max_nz
  ,size_type                   Q_row_i[]
  ,size_type                   Q_col_j[]
  ,GenPermMatrixSlice          *Q
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  //
  // Q = op(P1)*op(P2)
  //
  // There are several different possibilities for how to compute this
  // intersection.
  //
  DenseLinAlgPack::MtM_assert_sizes(
    P1.rows(), P1.cols() , P1_trans, P2.rows(), P2.cols() , P2_trans );
  //
  const size_type
    opP1_rows = rows(P1.rows(),P1.cols(),P1_trans),
    opP1_cols = cols(P1.rows(),P1.cols(),P1_trans),
    opP2_rows = rows(P2.rows(),P2.cols(),P2_trans),
    opP2_cols = cols(P2.rows(),P2.cols(),P2_trans);
  GPMSIP::EOrderedBy
    opP1_ordered_by = ordered_by(P1.ordered_by(),P1_trans),
    opP2_ordered_by = ordered_by(P2.ordered_by(),P2_trans);
  //
  *Q_nz = 0;
  // Either is zero?
  if( !P1.nz() || !P2.nz() ) {
    *Q_nz = 0;
    if(Q)
      Q->initialize(opP1_rows,opP2_cols,GenPermMatrixSlice::ZERO_MATRIX);
    return;
  }
  // Both are identity?
  if( P1.is_identity() && P2.is_identity() ) {
    *Q_nz = P1.nz(); // Should be the same as P2.nz();
    TEUCHOS_TEST_FOR_EXCEPT( !(  P1.nz() == P2.nz()  ) );
    if(Q)
      Q->initialize(opP1_rows,opP2_cols,GenPermMatrixSlice::IDENTITY_MATRIX);
    return;
  }
  // One is identity?
  if( P1.is_identity() || P2.is_identity() ) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement this but it is a little tricking?
    return;
  }
  //
  // Both are not identity or zero
  //
  if( ( opP1_ordered_by == GPMSIP::BY_COL || opP1_ordered_by == GPMSIP::BY_ROW_AND_COL )
    && ( opP2_ordered_by == GPMSIP::BY_ROW || opP2_ordered_by == GPMSIP::BY_ROW_AND_COL ) )
  {
    // This is great!  Both of the matrices are sorted and we don't need any temparary storage
    if( Q_max_nz ) {
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement initializing Q_row_i, Q_col_j
    }
    else {
      GenPermMatrixSlice::const_iterator
        P1_itr = P1.begin(), // Should not throw exception!
        P1_end = P1.end(),
        P2_itr = P2.begin(), // Should not throw exception!
        P2_end = P2.end();
      while( P1_itr != P1_end && P2_itr != P2_end ) {
        const size_type
          opP1_col_j = cols(P1_itr->row_i(),P1_itr->col_j(),P1_trans),
          opP2_row_i = rows(P2_itr->row_i(),P2_itr->col_j(),P2_trans);
        if( opP1_col_j < opP2_row_i ) {
          ++P1_itr;
          continue;
        }
        if( opP1_col_j > opP2_row_i ) {
          ++P2_itr;
          continue;
        }
        ++(*Q_nz);
        ++P1_itr;
        ++P2_itr;
      }
    }
  }
  else {
    //
    // We will load the row indices of op(P1) or the column indices op(P1)
    // indexed by column or row indices (whichever is smaller)
    // into a temporary sorted buffer and then loop through the nonzeros of the other.
    //
    // First let's get reorder P1 and P2 so that we put the rows of P2 into a buffer
    //
    const GenPermMatrixSlice
      &oP1      = opP1_cols > opP2_rows ? P1 : P2,
      &oP2      = opP1_cols > opP2_rows ? P2 : P1;
    const BLAS_Cpp::Transp
      oP1_trans = opP1_cols > opP2_rows ? P1_trans : trans_not(P1_trans),
      oP2_trans = opP1_cols > opP2_rows ? P2_trans : trans_not(P2_trans);
    // Load the column indices of op(op(P2)) into a lookup array
    typedef std::vector<size_type> oP2_col_j_lookup_t;      // Todo: use tempoarary workspace
    oP2_col_j_lookup_t oP2_col_j_lookup(rows(oP2.rows(),oP2.rows(),oP2_trans));
    std::fill( oP2_col_j_lookup.begin(), oP2_col_j_lookup.end(), 0 );
    {
      GenPermMatrixSlice::const_iterator
        itr     = oP2.begin(), // Should not throw exception!
        itr_end = oP2.end();
      while( itr != itr_end ) {
        oP2_col_j_lookup[rows(itr->row_i(),itr->col_j(),oP2_trans)]
          = cols(itr->row_i(),itr->col_j(),oP2_trans);
      }
    }
    // Loop through the columns of op(op(P1)) and look for matches
    if( Q_max_nz ) {
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement initializing Q_row_i, Q_col_j
    }
    else {
      GenPermMatrixSlice::const_iterator
        itr     = oP1.begin(), // Should not throw exception!
        itr_end = oP1.end();
      while( itr != itr_end ) {
        const size_type
          oP2_col_j = oP2_col_j_lookup[cols(itr->row_i(),itr->col_j(),oP1_trans)];
        if(oP2_col_j)
          ++(*Q_nz);
      }
    }

  }
  // Setup Q
  TEUCHOS_TEST_FOR_EXCEPT( !( Q == NULL ) ); // ToDo: Implement initializing Q
}
