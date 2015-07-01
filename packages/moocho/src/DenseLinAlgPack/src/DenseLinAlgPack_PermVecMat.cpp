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

#include "DenseLinAlgPack_PermVecMat.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_IVector.hpp"

#ifdef TEUCHOS_DEBUG   // Debug only!
bool DenseLinAlgPack::PermVecMat_print = false;
#include <iostream>
#include "DenseLinAlgPack_PermOut.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#endif

// Local assert function
namespace {
inline void i_assert_perm_size(size_t size1, size_t size2)
{
#ifdef LINALGPACK_CHECK_RHS_SIZES
  if(size1 != size2)
    throw std::length_error("The size of the permutation vector is not correct");
#endif
}

} // end namespace

void DenseLinAlgPack::identity_perm(IVector* perm) {
  if(!perm->size())
    throw std::length_error("DenseLinAlgPack::identity_perm(): perm must be sized");
  IVector::iterator	itr_perm		= perm->begin();
  for(size_type i = 1; i <= perm->size(); ++i)
    *itr_perm++ = i;
}

void DenseLinAlgPack::inv_perm(const IVector& perm, IVector* inv_perm) {
  inv_perm->resize(perm.size());
  for(size_type i = 1; i <= perm.size(); ++i)
    (*inv_perm)(perm(i)) = i;
}

void DenseLinAlgPack::perm_ele(const IVector& perm, DVectorSlice* vs)
{
  i_assert_perm_size(vs->dim(),perm.size());
  DVector tmp_v(vs->dim());
  DVector::iterator		v_itr		= tmp_v.begin(),
              v_itr_end	= tmp_v.end();
  IVector::const_iterator perm_itr = perm.begin();
  // Copy the elements in the permed order into the temp vector
  for(; v_itr != v_itr_end; ++v_itr, ++perm_itr)
    *v_itr = (*vs)(*perm_itr);
    
  // Copy the permed vector back
  (*vs) = tmp_v();
}

void DenseLinAlgPack::perm_ele(const DVectorSlice& x, const IVector& perm, DVectorSlice* y)
{
  i_assert_perm_size(x.dim(),perm.size());
  i_assert_perm_size(y->dim(),perm.size());

  IVector::const_iterator
    perm_itr	= perm.begin();
  DVectorSlice::iterator
    y_itr		= y->begin(),
    y_end		= y->end();
  while(y_itr != y_end)
    *y_itr++ = x(*perm_itr++);
}

void DenseLinAlgPack::inv_perm_ele(const DVectorSlice& y, const IVector& perm, DVectorSlice* x)
{
  i_assert_perm_size(y.dim(),perm.size());
  i_assert_perm_size(x->dim(),perm.size());
#ifdef TEUCHOS_DEBUG
  if( PermVecMat_print ) {
    std::cerr
      << "enter inv_perm_ele(y,perm,x):\n"
      << "before:\n"
      << "y =\n" << y
      << "perm =\n" << perm
      << "x =\n" << *x;
  }
#endif	
  DVectorSlice::const_iterator
    y_itr		= y.begin(),
    y_end		= y.end();
  IVector::const_iterator
    perm_itr	= perm.begin();
  while(y_itr != y_end)
    (*x)(*perm_itr++) = *y_itr++;
#ifdef TEUCHOS_DEBUG
  if( PermVecMat_print ) {
    std::cerr
      << "inv_perm_ele(y,perm,x):\n"
      << "after:\n"
      << "x =\n" << *x
      << "exit inv_perm_ele(...) ...\n";
  }
#endif	
}

void DenseLinAlgPack::perm_rows(const IVector& row_perm, DMatrixSlice* gms)
{
  i_assert_perm_size(gms->rows(),row_perm.size());
  DMatrix tmp_gm(gms->rows(),gms->cols());
  DMatrixSlice::size_type rows = gms->rows(), i;
  // Copy the rows in the correct order into the temp matrix.
  for(i = 1; i <= rows; ++i)
    tmp_gm.row(i) = gms->row(row_perm(i));
  // Copy the permed matrix back
  (*gms) = tmp_gm();
}

void DenseLinAlgPack::perm_cols(const IVector& col_perm, DMatrixSlice* gms)
{
  i_assert_perm_size(gms->cols(),col_perm.size());
  DMatrix tmp_gm(gms->rows(),gms->cols());
  DMatrixSlice::size_type cols = gms->cols(), i;
  // Copy the columns in the correct order into the temp matrix.
  for(i = 1; i <= cols; ++i)
    tmp_gm.col(i) = gms->col(col_perm(i));
  // Copy the permed matrix back
  (*gms) = tmp_gm();
}

void DenseLinAlgPack::perm_rows_cols(const IVector& row_perm, const IVector& col_perm
  , DMatrixSlice* gms)
{
  i_assert_perm_size(gms->rows(),row_perm.size());
  i_assert_perm_size(gms->cols(),col_perm.size());
  DMatrix tmp_gm(gms->rows(),gms->cols());
  DMatrixSlice::size_type rows = gms->rows(), cols = gms->cols(), i;
  // Copy the rows in the correct order into the temp matrix.
  for(i = 1; i <= rows; ++i)
    tmp_gm.row(i) = gms->row(row_perm(i));
  // Copy the columns in the correct order back into matrix.
  for(i = 1; i <= cols; ++i)
    gms->col(i) = tmp_gm.col(col_perm(i));
}
