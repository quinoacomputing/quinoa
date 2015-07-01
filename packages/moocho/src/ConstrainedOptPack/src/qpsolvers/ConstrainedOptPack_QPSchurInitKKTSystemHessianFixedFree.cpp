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

#include "ConstrainedOptPack_QPSchurInitKKTSystemHessianFixedFree.hpp"
#include "ConstrainedOptPack_initialize_Q_R_Q_X.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_sparse_bounds.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Midynamic_cast_verbose.h"
#include "MiWorkspacePack.h"
#include "Miprofile_hack.h"

namespace LinAlgOpPack {
    using AbstractLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptPack {

void QPSchurInitKKTSystemHessianFixedFree::initialize_kkt_system(
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
  using Teuchos::dyn_cast;
  using LinAlgOpPack::V_mV;
  namespace rcp = MemMngPack;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "QPSchurInitKKTSystemHessianFixedFree::initialize_kkt_system(...)" );
#endif

  // Validate type of and convert G
#ifdef _WINDOWS
  const MatrixSymOp&
    G_sym = dynamic_cast<const MatrixSymOp&>(G);
#else
  const MatrixSymOp&
    G_sym = dyn_cast<const MatrixSymOp>(G);
#endif

  const size_type nd = g.size();

  // Determine the number of initially fixed variables
  Workspace<EBounds> x_frfx(wss,nd);
  std::fill_n( &x_frfx[0], nd, FREE ); // make all free initially
  size_type
    num_init_fixed = 0;
  {
    const value_type inf_bnd = std::numeric_limits<value_type>::max();
    AbstractLinAlgPack::sparse_bounds_itr
      dLU_itr(
        dL.begin(), dL.end(), dL.offset(),
        dU.begin(), dU.end(), dU.offset(), inf_bnd );
    SpVectorSlice::const_iterator
      nu_itr = nu.begin(),
      nu_end = nu.end();
    const SpVector::difference_type o = nu.offset();
    while( !dLU_itr.at_end() || nu_itr != nu_end ) {
      if( dLU_itr.at_end() ) { // Add the rest of the elements in nu
        for( ; nu_itr != nu_end; ++num_init_fixed, ++nu_itr )
          x_frfx[nu_itr->indice() + o - 1] = ( nu_itr->value() > 0.0 ? UPPER : LOWER );
      }
      else { // Be carefull to add fixed dL(i) == dU(i)
        // Add elements in nu up to the current dLU_itr.indice()
        for( ; nu_itr != nu_end && nu_itr->indice() + o < dLU_itr.indice(); ++num_init_fixed, ++nu_itr )
          x_frfx[nu_itr->indice() + o - 1] = ( nu_itr->value() > 0.0 ? UPPER : LOWER );
        if( dLU_itr.lbound() == dLU_itr.ubound() ) {
          // This is a permanently fixed variable!
          x_frfx[dLU_itr.indice() - 1] = EQUALITY;
          ++num_init_fixed;
          // Don't add a duplicate entry in nu
          if( nu_itr != nu_end && nu_itr->indice() + o == dLU_itr.indice() )
            ++nu_itr;
        }
        ++dLU_itr;
      }
    }
  }
  TEUCHOS_TEST_FOR_EXCEPT( !(  nd >= num_init_fixed  ) );

  // n_R
  *n_R = nd - num_init_fixed;
  
  // Set up i_x_free[], i_x_fixed[], bnd_fixed[], and b_X
  i_x_free->resize(*n_R);
  i_x_fixed->resize(num_init_fixed+1);
  bnd_fixed->resize(num_init_fixed+1);
  b_X->resize(num_init_fixed+1);
  {
    const value_type inf_bnd = std::numeric_limits<value_type>::max();
    AbstractLinAlgPack::sparse_bounds_itr
      dLU_itr(
        dL.begin(), dL.end(), dL.offset(),
        dU.begin(), dU.end(), dU.offset(), inf_bnd );
    size_type i_R = 0, i_X = 0;
    for( size_type i = 1; i <= nd; ++i ) {
      const EBounds
        bnd_i = x_frfx[i-1];
      if( bnd_i == FREE ) {
        (*i_x_free)[i_R] = i;
        ++i_R;
      }
      else {
        (*i_x_fixed)[i_X] = i;
        (*bnd_fixed)[i_X] = bnd_i;
        TEUCHOS_TEST_FOR_EXCEPT( !(  !dLU_itr.at_end()  ) );    // find entry in b_X
        while( dLU_itr.indice() < i )
          ++dLU_itr;
        TEUCHOS_TEST_FOR_EXCEPT( !(  dLU_itr.indice() == i  ) );
        value_type b_X_val = 0.0;
        switch( bnd_i ) {
          case EQUALITY:
          case LOWER:
            b_X_val = dLU_itr.lbound();
            break;
          case UPPER:
            b_X_val = dLU_itr.ubound();
            break;
          default:
            TEUCHOS_TEST_FOR_EXCEPT(true); // Local error only?
        }
        (*b_X)[i_X] = b_X_val;
        ++i_X;
      }
    }
    (*i_x_fixed)[i_X] = nd+1;   // built-in relaxation variable
    (*bnd_fixed)[i_X] = LOWER;
    (*b_X)[i_X]       = etaL;
    ++i_X;
  }
  
  // j_f_decomp[] = empty
  j_f_decomp->resize(0);

  // Initialize temporary Q_R and Q_X (not including extra relaxation variable)
  Workspace<size_type>
    Q_R_row_i(wss,*n_R),
    Q_R_col_j(wss,*n_R),
    Q_X_row_i(wss,num_init_fixed),
    Q_X_col_j(wss,num_init_fixed);
  GenPermMatrixSlice
    Q_R, Q_X;
  initialize_Q_R_Q_X(
    *n_R,num_init_fixed,&(*i_x_free)[0],&(*i_x_fixed)[0],false
    ,&Q_R_row_i[0],&Q_R_col_j[0],&Q_R
    ,&Q_X_row_i[0],&Q_X_col_j[0],&Q_X
    );

  //
  // Create and initialize object for Ko = G_RR = Q_R'*G*Q_R
  //

  // Compute the dense matrix G_RR
  DMatrix G_RR_dense(*n_R,*n_R);
  DMatrixSliceSym sym_G_RR_dense(G_RR_dense(),BLAS_Cpp::lower);
  AbstractLinAlgPack::Mp_StPtMtP(
    &sym_G_RR_dense, 1.0, MatrixSymOp::DUMMY_ARG
    ,G_sym, Q_R, BLAS_Cpp::no_trans, 0.0 );
  // Initialize a factorization object for this matrix
  typedef Teuchos::RCP<MatrixSymPosDefCholFactor> G_RR_ptr_t;
  G_RR_ptr_t
    G_RR_ptr = new MatrixSymPosDefCholFactor();
  G_RR_ptr->initialize(sym_G_RR_dense);
  
  *Ko = Teuchos::rcp_implicit_cast<Ko_ptr_t::element_type>(G_RR_ptr); // Ko is initialized!

  // ToDo: (2001/07/05) We could be more carefull about how memory is initialized and reused
  // in the future but this implementation is just easier.

  // fo = - Q_R'*g - Q_R'*G*(Q_X*b_X)
  LinAlgOpPack::V_StMtV( fo, -1.0, Q_R, BLAS_Cpp::trans, g );
  if( num_init_fixed ) {
    SpVector b_XX;
    AbstractLinAlgPack::V_MtV( &b_XX, Q_X, BLAS_Cpp::no_trans, (*b_X)(1,num_init_fixed) );
    AbstractLinAlgPack::Vp_StPtMtV( &(*fo)(), -1.0, Q_R, BLAS_Cpp::trans, G, BLAS_Cpp::no_trans, b_XX() );
  }

}

} // end namesapce ConstrainedOptPack

#endif // 0
