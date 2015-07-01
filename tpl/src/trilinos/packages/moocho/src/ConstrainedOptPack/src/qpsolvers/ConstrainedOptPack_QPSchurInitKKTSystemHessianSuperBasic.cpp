#if 0

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

#include "ConstrainedOptPack_QPSchurInitKKTSystemHessianSuperBasic.hpp"
#include "ConstrainedOptPack_MatrixHessianSuperBasic.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Midynamic_cast_verbose.h"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptPack {

void QPSchurInitKKTSystemHessianSuperBasic::initialize_kkt_system(
  const DVectorSlice&    g
  ,const MatrixOp&  G
  ,value_type           etaL
  ,const SpVectorSlice& dL
  ,const SpVectorSlice& dU
  ,const MatrixOp*  F
  ,BLAS_Cpp::Transp     trans_F
  ,const DVectorSlice*   f
  ,const DVectorSlice&   d
  ,const SpVectorSlice& nu
  ,size_type*           n_R
  ,i_x_free_t*          i_x_free
  ,i_x_fixed_t*         i_x_fixed
  ,bnd_fixed_t*         bnd_fixed
  ,j_f_decomp_t*        j_f_decomp
  ,DVector*              b_X
  ,Ko_ptr_t*            Ko
  ,DVector*              fo
  ) const
{
  using BLAS_Cpp::trans;
  using Teuchos::dyn_cast;
  using LinAlgOpPack::V_mV;
  using LinAlgOpPack::V_StMtV;
  using AbstractLinAlgPack::Vp_StMtV;

  // Validate type of and convert G
  const MatrixHessianSuperBasic
    *G_super_ptr = dynamic_cast<const MatrixHessianSuperBasic*>(&G);

  if( G_super_ptr == NULL ) {
    init_kkt_full_.initialize_kkt_system(
      g,G,etaL,dL,dU,F,trans_F,f,d,nu,n_R,i_x_free,i_x_fixed,bnd_fixed
      ,j_f_decomp,b_X,Ko,fo);
    return;
  }

  const MatrixHessianSuperBasic
    &G_super = *G_super_ptr;

  // get some stuff
  const GenPermMatrixSlice
    &Q_R = G_super.Q_R(),
    &Q_X = G_super.Q_X();
  const size_type
    nd   = G_super.rows(),
    nd_R = Q_R.cols(),
    nd_X = Q_X.cols();
  TEUCHOS_TEST_FOR_EXCEPT( !(  nd_R + nd_X == nd  ) );

  // Setup output arguments

  // n_R = nd_R
  *n_R = nd_R;
  // i_x_free[(G.Q_R.begin()+l-1)->col_j()-1] = (G.Q_R.begin()+l-1)->row_i(), l = 1...nd_R
  i_x_free->resize( Q_R.is_identity() ? 0: nd_R );
  if( nd_R && !Q_R.is_identity() ) {
    GenPermMatrixSlice::const_iterator
      Q_itr = Q_R.begin();
    for( ; Q_itr != Q_R.end(); ++Q_itr ) {
      const size_type i = Q_itr->row_i();
      const size_type k = Q_itr->col_j();
      TEUCHOS_TEST_FOR_EXCEPT( !(  0 < i && i <= nd  ) );
      TEUCHOS_TEST_FOR_EXCEPT( !(  0 < k && k <= nd_R  ) );
      (*i_x_free)[k-1] = i;
    }
  }
  // i_x_fixed[]
  i_x_fixed->resize(nd_X+1);
  if(nd_X) {
    // i_x_fixed[(G.Q_X.begin()+l-1)->col_j()-1] = (G.Q_X.begin()+l-1)->row_i(), l = 1...nd_X
    GenPermMatrixSlice::const_iterator
      Q_itr = Q_X.begin();
    for( ; Q_itr != Q_X.end(); ++Q_itr ) {
      const size_type i = Q_itr->row_i();
      const size_type k = Q_itr->col_j();
      TEUCHOS_TEST_FOR_EXCEPT( !(  0 < i && i <= nd  ) );
      TEUCHOS_TEST_FOR_EXCEPT( !(  0 < k && k <= nd_X  ) );
      (*i_x_fixed)[k-1] = i;
    }
  }
  (*i_x_fixed)[nd_X] = nd+1; // relaxation is always initially active
  // bnd_fixed[]
  bnd_fixed->resize(nd_X+1);
  if(nd_X) {
    // bnd_fixed[l-1] = G.bnd_fixed[l-1], l = 1...nd_X
    typedef MatrixHessianSuperBasic MHSB;
    const MHSB::bnd_fixed_t
      &bnd_fixed_from = G_super.bnd_fixed();
    TEUCHOS_TEST_FOR_EXCEPT( !( bnd_fixed_from.size() == nd_X ) );
    MHSB::bnd_fixed_t::const_iterator
      bnd_from_itr = bnd_fixed_from.begin();
    bnd_fixed_t::iterator
      bnd_to_itr = bnd_fixed->begin();
    std::copy( bnd_from_itr, bnd_fixed_from.end(), bnd_to_itr );
  }
  (*bnd_fixed)[nd_X] = LOWER; // relaxation is always initially active
  // j_f_decomp[]
  j_f_decomp->resize(0);
  // b_X
  b_X->resize(nd_X+1);
  if(nd_X) {
    // b_X[l-1] = { dL(i) if bnd_fixed[l-1] == LOWER or EQUALITY
    //              dU(i) if bnd_fixed[l-1] == UPPER }
    //             , l = 1...nd_X
    //             (where i = i_x_fixed[l-1])
    bnd_fixed_t::const_iterator
      bnd_itr     = const_cast<const bnd_fixed_t&>(*bnd_fixed).begin(),
      bnd_itr_end = const_cast<const bnd_fixed_t&>(*bnd_fixed).begin() + nd_X;
    i_x_fixed_t::const_iterator
      i_x_itr     = const_cast<const i_x_fixed_t&>(*i_x_fixed).begin();
    DVector::iterator
      b_X_itr     = b_X->begin();
    const SpVectorSlice::element_type
      *ele = NULL;
    for( ; bnd_itr != bnd_itr_end; ++bnd_itr, ++i_x_itr, ++b_X_itr ) {
      const size_type i = *i_x_itr;
      switch(*bnd_itr) {
          case LOWER:
          case EQUALITY:
          *b_X_itr = (ele = dL.lookup_element(i))->value(); // Should not be null!
          break;
          case UPPER:
          *b_X_itr = (ele = dU.lookup_element(i))->value(); // Should not be null!
          break;
          default:
          TEUCHOS_TEST_FOR_EXCEPT(true);
      }
    }
  }
  (*b_X)[nd_X] = etaL; // relaxation is always initially active
  // Ko = G.B_RR
  *Ko = G_super.B_RR_ptr(); // now B_RR is a shared object
  // fo = - G.Q_R'*g - op(G.B_RX)*b_X(1:nd_X)
  V_StMtV( fo, -1.0, Q_R, trans, g );
  if( nd_X && G_super.B_RX_ptr().get() )
    Vp_StMtV( &(*fo)(), -1.0, *G_super.B_RX_ptr(), G_super.B_RX_trans(), (*b_X)(1,nd_X) );

}

} // end namesapce ConstrainedOptPack

#endif // 0
