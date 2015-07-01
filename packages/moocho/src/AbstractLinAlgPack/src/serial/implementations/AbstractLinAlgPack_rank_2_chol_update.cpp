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

#include "AbstractLinAlgPack_rank_2_chol_update.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "DenseLinAlgPack_BLAS_Cpp.hpp"

void AbstractLinAlgPack::rank_2_chol_update(
  const value_type     a
  ,DVectorSlice         *u
  ,const DVectorSlice   &v
  ,DVectorSlice         *w
  ,DMatrixSliceTriEle         *R
  ,BLAS_Cpp::Transp    R_trans
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::upper;
  using BLAS_Cpp::lower;
  using BLAS_Cpp::rotg;
  using DenseLinAlgPack::row;
  using DenseLinAlgPack::col;
  using DenseLinAlgPack::rot;
  using DenseLinAlgPack::Vp_StV;
  //
  // The idea for this routine is to perform a set of givens rotations to return
  // op(R) +a *u*v' back to triangular form.  The first set of (n-1) givens
  // rotations Q1 eleminates the elements in u to form ||u||2 * e(last_i) and transform:
  //
  //     Q1 * (op(R) + a*u*v') -> op(R1) + a*||u||2 * e(last_i) * v'
  //
  // where op(R1) is a upper or lower Hessenberg matrix.  The ordering of the rotations
  // and and wether last_i = 1 or n depends on whether op(R) is logically upper or lower
  // triangular.
  //
  // If op(R) logically lower triangular then we have:
  //
  // [ x          ]       [ x ]
  // [ x  x       ]       [ x ]
  // [ x  x  x    ]       [ x ]
  // [ x  x  x  x ] + a * [ x ] * [ x x x x ]
  // 
  //       op(R)    + a *   u   *  v'
  // 
  //  Rotations are applied to zero out u(i), for i = 1..n-1 to form
  //  (first_i = 1, di = +1, last_i = n):
  // 
  //  [ x  x       ]       [   ]
  //  [ x  x  x    ]       [   ]
  //  [ x  x  x  x ]       [   ]
  //  [ x  x  x  x ] + a * [ x ] * [ x x x x ]
  // 
  //      op(R1)     + a * ||u||2*e(n) * v'
  // 
  //  If op(R) logically upper triangular then we have:
  // 
  //  [ x  x  x  x ] + a * [ x ] * [ x x x x ]
  //  [    x  x  x ]       [ x ]
  //  [       x  x ]       [ x ]
  //  [          x ]       [ x ]
  // 
  //       op(R)     + a *   u   *  v'
  // 
  //  Rotations are applied to zero out u(i), for i = n..2 to form
  //  (first_i = n, di = -1, last_i = 1):
  // 
  //  [ x  x  x  x ] + a * [ x ] * [ x x x x ]
  //  [ x  x  x  x ]       [   ]
  //  [    x  x  x ]       [   ]
  //  [       x  x ]       [   ]
  // 
  //      op(R1)     + a * ||u||2*e(1) * v'
  // 
  //  The off diagonal elements in R created by this set of rotations is not stored
  //  in R(*,*) but instead in the workspace vector w(*).  More precisely,
  //
  //       w(i+wo) = R(i,i+di), for i = first_i...last_i-di
  //
  //       where wo = 0 if lower, wo = -1 if upper
  //
  //  This gives the client more flexibility in how workspace is handeled.
  //  Other implementations overwrite the off diagonal of R (lazyness I guess).
  // 
  //  The update a*||u||2*e(last_i)*v' is then added to the row op(R1)(last_i,:) to
  //  give op(R2) which is another upper or lower Hessenberg matrix:
  //  
  //       op(R1) + a*||u||2 * e(last_i) * v' -> op(R2)
  // 
  //  Then, we apply (n-1) more givens rotations to return op(R2) back to triangular
  //  form and we are finished:
  // 
  //       Q2 * op(R2) -> op(R_new)
  // 
  const size_type n = R->rows();
  DenseLinAlgPack::Mp_M_assert_sizes( n, n, no_trans, u->dim(), v.dim(), no_trans );
  // Is op(R) logically upper or lower triangular
  const BLAS_Cpp::Uplo
    opR_uplo = ( (R->uplo()==upper && R_trans==no_trans) || (R->uplo()==lower && R_trans==trans)
           ? upper : lower );
  // Get frame of reference
  size_type first_i, last_i;
  int di, wo;
  if( opR_uplo == lower ) {
    first_i =  1;
    di      = +1;
    last_i  =  n;
    wo      =  0;
  }
  else {
    first_i =  n;
    di      = -1;
    last_i  =  1;
    wo      = -1;
  }
  // Zero out off diagonal workspace
  if(n > 1)
    *w = 0.0;
  // Find the first u(k) != 0 in the direction of the forward sweeps
  size_type k = first_i;
  for( ; k != last_i+di; k+=di )
    if( (*u)(k) != 0.0 ) break;
  if( k == last_i+di ) return; // All of u is zero, no update to perform!
  //
  // Transform Q(.)*...*Q(.) * u -> ||u||2 * e(last_i) while applying
  // Q1*R forming R1 (upper or lower Hessenberg)
  //
  {for( size_type i = k; i != last_i; i+=di ) {
    // [  c  s ] * [ u(i+di) ] -> [ u(i+id) ]
    // [ -s  c ]   [ u(i)    ]    [ 0       ]
    value_type c, s;  // Here c and s are computed to zero out u(i)
    rotg( &(*u)(i+di), &(*u)(i), &c, &s );
    // [  c  s ] * [ op(R)(i+di,:) ] -> [ op(R)(i+di,:) ]
    // [ -s  c ]   [ op(R)(i   ,:) ]    [ op(R)(i   ,:) ]
    DVectorSlice
      opR_row_idi = row(R->gms(),R_trans,i+di),
      opR_row_i   = row(R->gms(),R_trans,i);
    if( opR_uplo == lower ) {
      // op(R)(i   ,:) = [ x  x  x  w(i+wo)     ] i
      // op(R)(i+di,:) = [ x  x  x   x          ]
      //                         i
      rot( c, s, &opR_row_idi(1  ,i  ), &opR_row_i(1,i)  ); // First nonzero columns
      rot( c, s, &opR_row_idi(i+1,i+1), &(*w)(i+wo,i+wo) ); // Last nonzero column
    }
    else {
      // op(R)(i+di,:) = [         x  x  x  x ]
      // op(R)(i   ,:) = [    w(i+wo) x  x  x ] i
      //                              i
      rot( c, s, &opR_row_idi(i-1,i-1), &(*w)(i+wo,i+wo) ); // First nonzero column
      rot( c, s, &opR_row_idi(i  ,n  ), &opR_row_i(i,n)  ); // Last nonzero columns
    }
  }}
  //
  // op(R1) + a * ||u||2 * e(last_i) * v' -> op(R2)
  //
  Vp_StV( &row(R->gms(),R_trans,last_i), a * (*u)(last_i), v );
  //
  // Apply (k-1) more givens rotations to return op(R2) back to triangular form:
  //
  // Q2 * R2 -> R_new
  //
  {for( size_type i = last_i; i != k; i-=di ) {
    DVectorSlice
      opR_row_i   = row(R->gms(),R_trans,i),
      opR_row_idi = row(R->gms(),R_trans,i-di);
    value_type
      &w_idiwo = (*w)(i-di+wo);
    // [  c  s ]   [ op(R)(i,i)    ]    [ op(R)(i,i) ]
    // [ -s  c ] * [ op(R)(i-id,i) ] -> [ 0          ]  w(i-di+wo)
    value_type &R_i_i = opR_row_i(i);
    value_type c, s;  // Here c and s are computed to zero out w(i-di+wo)
    rotg( &R_i_i, &w_idiwo, &c, &s );
    // Make sure the diagonal is positive
    value_type scale = +1.0;
    if( R_i_i < 0.0 ) {
      R_i_i    = -R_i_i;
      scale    = -1.0;
      w_idiwo  = -w_idiwo; // Must also change this stored value
    }
    // [  c  s ] * [ op(R)(i,:)   ] -> [ op(R)(i,:)    ]
    // [ -s  c ]   [ op(R)(i-id,:)]    [ op(R)(i-id,:) ]
    if( opR_uplo == lower ) {
      // op(R)(i-di,:) = [ x  x  x   0     ]      
      // op(R)(i,:)    = [ x  x  x   x     ] i
      //                             i
      rot( scale*c, scale*s, &opR_row_i(1,i-1), &opR_row_idi(1,i-1) );
    }
    else {
      // op(R)(i,:)     = [      x  x  x  x ] i
      // op(R)(i-di,:)  = [      0  x  x  x ]
      //                         i
      rot( scale*c, scale*s, &opR_row_i(i+1,n), &opR_row_idi(i+1,n) );
    }
  }}
  // When we get here note that u contains the first set of rotations Q1 and
  // w contains the second set of rotations Q2.
}
