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

#include "MoochoPack_ReducedHessianSecantUpdateBFGSProjected_Strategy.hpp"
#include "MoochoPack_PBFGS_helpers.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "ConstrainedOptPack/src/AbstractLinAlgPack_MatrixSymAddDelUpdateable.hpp"
#include "ConstrainedOptPack/src/AbstractLinAlgPack_BFGS_helpers.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSliceOut.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymInitDiag.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Midynamic_cast_verbose.h"
#include "MiWorkspacePack.h"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

namespace MoochoPack {

ReducedHessianSecantUpdateBFGSProjected_Strategy::ReducedHessianSecantUpdateBFGSProjected_Strategy(
  const bfgs_update_ptr_t&      bfgs_update
  ,value_type                   act_set_frac_proj_start
  ,value_type                   project_error_tol
  ,value_type                   super_basic_mult_drop_tol
  )
  : bfgs_update_(bfgs_update)
  , act_set_frac_proj_start_(act_set_frac_proj_start)
  , project_error_tol_(project_error_tol)
  , super_basic_mult_drop_tol_(super_basic_mult_drop_tol)
{}

bool ReducedHessianSecantUpdateBFGSProjected_Strategy::perform_update(
  DVectorSlice* s_bfgs, DVectorSlice* y_bfgs, bool first_update
  ,std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
  ,MatrixOp *rHL_k
  )
{
  using std::setw;
  using std::endl;
  using std::right;
  using Teuchos::dyn_cast;
  using DenseLinAlgPack::dot;
  using LinAlgOpPack::V_MtV;
  using AbstractLinAlgPack::norm_inf;
  typedef ConstrainedOptPack::MatrixHessianSuperBasic MHSB_t;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\n*** (PBFGS) Projected BFGS ...\n";
  }

#ifdef _WINDOWS
  MHSB_t &rHL_super = dynamic_cast<MHSB_t&>(*rHL_k);
#else
  MHSB_t &rHL_super = dyn_cast<MHSB_t>(*rHL_k);
#endif

  const size_type
    n    = algo->nlp().n(),
    r    = algo->nlp().r(),
    n_pz = n-r;
  const GenPermMatrixSlice
    &Q_R = rHL_super.Q_R(),
    &Q_X = rHL_super.Q_X();

  bool do_projected_rHL_RR = false;

  // Check that the current update is sufficiently p.d. before we do anything
  const value_type
    sTy = dot(*s_bfgs,*y_bfgs),
    yTy = dot(*y_bfgs,*y_bfgs);
  if( !ConstrainedOptPack::BFGS_sTy_suff_p_d(
    *s_bfgs,*y_bfgs,&sTy
    ,  int(olevel) >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL )
    && !bfgs_update().use_dampening()
    )
  {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out	<< "\nWarning!  use_damening == false so there is no way we can perform any kind BFGS update (projected or not) so we will skip it!\n";
    }
    quasi_newton_stats_(*s).set_k(0).set_updated_stats(
      QuasiNewtonStats::INDEF_SKIPED );
    return true;
  }

  // Get matrix scaling
  value_type
    rHL_XX_scale = sTy > 0.0 ? yTy/sTy : 1.0;

  // 
  // Initialize or adjust the active set before the BFGS update
  //
  if( !s->nu().updated_k(-1) ) {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out	<< "\nWarning!  nu_k(-1) has not been updated.  No adjustment to the set of superbasic variables is possible!\n";
    }
  }
  else if( Q_R.is_identity() ) {
    // Determine when to start adding and removing rows/cols form rHL_RR
    if( act_set_stats_(*s).updated_k(-1) ) {
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
        out	<< "\nDetermining if projected BFGS updating of superbasics should begin ...\n";
      }
      // Determine if the active set has calmed down enough
      const SpVector
        &nu_km1 = s->nu().get_k(-1);
      const SpVectorSlice
        nu_indep = nu_km1(s->var_indep());
      do_projected_rHL_RR = PBFGSPack::act_set_calmed_down(
        act_set_stats_(*s).get_k(-1)
        ,act_set_frac_proj_start()
        ,olevel,out
        );
      if( do_projected_rHL_RR ) {
        //
        // Determine the set of initially fixed and free independent variables.
        //
        typedef Workspace<size_type>                              i_x_t;
        typedef Workspace<ConstrainedOptPack::EBounds>   bnd_t;
        i_x_t   i_x_free(wss,n_pz);
        i_x_t   i_x_fixed(wss,n_pz);
        bnd_t   bnd_fixed(wss,n_pz);
        i_x_t   l_x_fixed_sorted(wss,n_pz);
        size_type n_pz_X = 0, n_pz_R = 0;
        value_type sRTBRRsR = 0.0, sRTyR = 0.0, sXTBXXsX = 0.0, sXTyX = 0.0;
        // Get the elements in i_x_free[] for variables that are definitely free
        // and initialize s_R'*y_R
        PBFGSPack::init_i_x_free_sRTsR_sRTyR(
          nu_indep, *s_bfgs, *y_bfgs
          , &n_pz_R, &i_x_free[0], &sRTBRRsR, &sRTyR );  // We don't really want sRTRBBsR here
        // rHL_XX = some_scaling * I
        if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
          out	<< "\nScaling for diagonal rHL_XX = rHL_XX_scale*I, rHL_XX_scale = " << rHL_XX_scale << std::endl;
        }
        Workspace<value_type> rHL_XX_diag_ws(wss,nu_indep.nz());
        DVectorSlice rHL_XX_diag(&rHL_XX_diag_ws[0],rHL_XX_diag_ws.size());
        rHL_XX_diag = rHL_XX_scale;
        // s_R'*B_RR_*s_R
        Workspace<value_type> Q_R_Q_RT_s_ws(wss,n_pz);
        DVectorSlice Q_R_Q_RT_s(&Q_R_Q_RT_s_ws[0],Q_R_Q_RT_s_ws.size());
        Q_R_Q_RT_s = 0.0;
        {for( size_type k = 0; k < n_pz_R; ++k ) {
          const size_type i = i_x_free[k];
          Q_R_Q_RT_s(i) = (*s_bfgs)(i);
        }}
        sRTBRRsR = AbstractLinAlgPack::transVtMtV( Q_R_Q_RT_s, *rHL_k, BLAS_Cpp::no_trans,  Q_R_Q_RT_s );
        // Sort fixed variables according to |s_X(i)^2*B_XX(i,i)|/|sRTBRRsR| + |s_X(i)*y_X(i)|/|sRTyR|
        // and initialize s_X'*B_XX*s_X and s_X*y_X
        PBFGSPack::sort_fixed_max_cond_viol(
          nu_indep,*s_bfgs,*y_bfgs,rHL_XX_diag,sRTBRRsR,sRTyR
          ,&sXTBXXsX,&sXTyX,&l_x_fixed_sorted[0]);
        // Pick initial set of i_x_free[] and i_x_fixed[] (sorted!)
        PBFGSPack::choose_fixed_free(
          project_error_tol(),super_basic_mult_drop_tol(),nu_indep
          ,*s_bfgs,*y_bfgs,rHL_XX_diag,&l_x_fixed_sorted[0]
          ,olevel,out,&sRTBRRsR,&sRTyR,&sXTBXXsX,&sXTyX
          ,&n_pz_X,&n_pz_R,&i_x_free[0],&i_x_fixed[0],&bnd_fixed[0] );
        //
        // Delete rows/cols from rHL_RR for fixed variables
        //
#ifdef _WINDOWS
        MatrixSymAddDelUpdateable
          &rHL_RR = dynamic_cast<MatrixSymAddDelUpdateable&>(
            const_cast<MatrixSymWithOpFactorized&>(*rHL_super.B_RR_ptr())
            );
#else
        MatrixSymAddDelUpdateable
          &rHL_RR = dyn_cast<MatrixSymAddDelUpdateable>(
            const_cast<MatrixSymWithOpFactorized&>(*rHL_super.B_RR_ptr())
            );
#endif
        if( n_pz_R < n_pz ) {
          if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
            out << "\nDeleting n_pz_X = " << n_pz_X << " rows/columns from rHL_RR for fixed independent variables...\n";
          }
          {for( size_type k = n_pz_X; k > 0; --k ) { // Delete from the largest to the smallest index (cheaper!)
            rHL_RR.delete_update( i_x_fixed[k-1], false );
          }}
          TEUCHOS_TEST_FOR_EXCEPT( !(  rHL_RR.rows() == n_pz_R  ) );
          if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
            out << "\nrHL_RR after rows/columns where removed =\n" << *rHL_super.B_RR_ptr();
          }
          // Initialize rHL_XX = rHL_XX_scale*I
#ifdef _WINDOWS
          MatrixSymInitDiag
            &rHL_XX = dynamic_cast<MatrixSymInitDiag&>(
              const_cast<MatrixSymOp&>(*rHL_super.B_XX_ptr())
              );
#else
          MatrixSymInitDiag
            &rHL_XX = dyn_cast<MatrixSymInitDiag>(
              const_cast<MatrixSymOp&>(*rHL_super.B_XX_ptr())
              );
#endif
          rHL_XX.init_identity(n_pz_X,rHL_XX_scale);
          // Reinitialize rHL for new active set
          rHL_super.initialize(
            n_pz, n_pz_R, &i_x_free[0], &i_x_fixed[0], &bnd_fixed[0]
            ,rHL_super.B_RR_ptr(),NULL,BLAS_Cpp::no_trans,rHL_super.B_XX_ptr()
            );
          if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
            out << "\nFull rHL after reinitialization but before BFGS update:\n"
              << "\nrHL =\n" << *rHL_k
              << "\nQ_R =\n" << rHL_super.Q_R()
              << "\nQ_X =\n" << rHL_super.Q_X();
          }
        }
        else {
          do_projected_rHL_RR = false;
          if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
            out << "\nWith n_pz_X = " << n_pz_X << ", there where no variables to drop from superbasis!\n";
          }
        }
      }
    }
  }
  else {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out	<< "\nAdjust the set of superbasic variables and the projected reduced Hessian rHL_RR ...\n";
    }
    //
    // Modify rHL_RR by adding and dropping rows/cols for freeded and fixed variables
    //
    const SpVectorSlice
      nu_indep = s->nu().get_k(-1)(s->var_indep());
    //
    // Determine new Q_R and Q_X
    //
    typedef Workspace<size_type>                              i_x_t;
    typedef Workspace<ConstrainedOptPack::EBounds>   bnd_t;
    i_x_t   i_x_free(wss,n_pz);
    i_x_t   i_x_fixed(wss,n_pz);
    bnd_t   bnd_fixed(wss,n_pz);
    i_x_t   l_x_fixed_sorted(wss,n_pz);
    size_type n_pz_X = 0, n_pz_R = 0;
    value_type sRTBRRsR = 0.0, sRTyR = 0.0, sXTBXXsX = 0.0, sXTyX = 0.0;
    // Get the elements in i_x_free[] for variables that are definitely free
    // and initialize s_R'*y_R.  This will be the starting point for the new Q_R.
    PBFGSPack::init_i_x_free_sRTsR_sRTyR(
      nu_indep, *s_bfgs, *y_bfgs
      , &n_pz_R, &i_x_free[0], &sRTBRRsR, &sRTyR );  // We don't really want sRTBRRsR here
    // Initialize rHL_XX_diag = some_scaling * I as though all of the currently fixed variables
    // will be left out of Q_R.  Some of these variables might already be in Q_R and B_RR
    // and may still be in Q_R and B_RR after we are finished adjusting the sets Q_R and Q_X
    // and we don't want to delete these rows/cols in B_RR just yet!
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out	<< "\nScaling for diagonal rHL_XX = rHL_XX_scale*I, rHL_XX_scale = " << rHL_XX_scale << std::endl;
    }
    Workspace<value_type> rHL_XX_diag_ws(wss,nu_indep.nz());
    DVectorSlice rHL_XX_diag(&rHL_XX_diag_ws[0],rHL_XX_diag_ws.size());
    rHL_XX_diag = rHL_XX_scale;
    // Initialize rHL_XX = rHL_XX_scale * I so that those variables in the current Q_R
    // not in the estimate i_x_free[] will have their proper value when s_R'*B_RR*s_R computed
    // for the estimate i_x_free[].  This is needed to change the behavior of *rHL_k which
    // is used below to compute s_R'*B_RR*s_R
#ifdef _WINDOWS
    MatrixSymInitDiag
      &rHL_XX = dynamic_cast<MatrixSymInitDiag&>(
        const_cast<MatrixSymOp&>(*rHL_super.B_XX_ptr())
        );
#else
    MatrixSymInitDiag
      &rHL_XX = dyn_cast<MatrixSymInitDiag>(
        const_cast<MatrixSymOp&>(*rHL_super.B_XX_ptr())
        );
#endif
    rHL_XX.init_identity(rHL_XX.rows(),rHL_XX_scale); // Don't change its size yet!
    // s_R'*B_RR_*s_R
    // This will only include those terms for the variable actually free.
    Workspace<value_type> Q_R_Q_RT_s_ws(wss,n_pz);
    DVectorSlice Q_R_Q_RT_s(&Q_R_Q_RT_s_ws[0],Q_R_Q_RT_s_ws.size());
    Q_R_Q_RT_s = 0.0;
    {for( size_type k = 0; k < n_pz_R; ++k ) {
      const size_type i = i_x_free[k];
      Q_R_Q_RT_s(i) = (*s_bfgs)(i);
    }}
    sRTBRRsR = AbstractLinAlgPack::transVtMtV( Q_R_Q_RT_s, *rHL_k, BLAS_Cpp::no_trans,  Q_R_Q_RT_s );
    // Sort fixed variables according to |s_X(i)^2*B_XX(i,i)|/|sRTBRRsR| + |s_X(i)*y_X(i)|/|sRTyR|
    PBFGSPack::sort_fixed_max_cond_viol(
      nu_indep,*s_bfgs,*y_bfgs,rHL_XX_diag,sRTBRRsR,sRTyR
      ,&sXTBXXsX,&sXTyX,&l_x_fixed_sorted[0]);
    // Pick initial set of i_x_free[] and i_x_fixed[] (sorted!)
    PBFGSPack::choose_fixed_free(
      project_error_tol(),super_basic_mult_drop_tol(),nu_indep
      ,*s_bfgs,*y_bfgs,rHL_XX_diag,&l_x_fixed_sorted[0]
      ,olevel,out,&sRTBRRsR,&sRTyR,&sXTBXXsX,&sXTyX
      ,&n_pz_X,&n_pz_R,&i_x_free[0],&i_x_fixed[0],&bnd_fixed[0] );
    // Get the changes to the set of superbasic variables
    size_type num_free_to_fixed = 0, num_fixed_to_free = 0;
    i_x_t  i_x_free_to_fixed(wss,Q_R.cols());
    i_x_t  i_x_fixed_to_free(wss,Q_X.cols());
    i_x_t  i_x_free_still(wss,Q_R.cols());             // Will be set to those indices still in Q_R
    std::fill_n( &i_x_free_still[0], Q_R.cols(), 0 );  // in the same order as in Q_R and B_RR
    {
      GenPermMatrixSlice::const_iterator
        Q_R_begin     = Q_R.begin(),
        Q_R_itr       = Q_R_begin,
        Q_R_end       = Q_R.end();
      const size_type
        *i_x_free_itr = &i_x_free[0],
        *i_x_free_end = i_x_free_itr + n_pz_R;
      for( size_type i = 1; i <= n_pz; ++i ) {
        if( Q_R_itr == Q_R_end && i_x_free_itr == i_x_free_end ) {
          break; // The rest of these variables were and still are not superbasic
        }
        else if( i_x_free_itr == i_x_free_end ) {
          // A variable that was in the superbasis now is not
          i_x_free_to_fixed[num_free_to_fixed] = Q_R_itr->row_i();
          num_free_to_fixed++;
          ++Q_R_itr;
        }
        else if( Q_R_itr == Q_R_end ) {
          // A variable that was not in the superbasis now is
          i_x_fixed_to_free[num_fixed_to_free] = *i_x_free_itr;
          num_fixed_to_free++;
          ++i_x_free_itr;
        }
        else {
          if( Q_R_itr->row_i() == *i_x_free_itr ) {
            // Both still superbasic
            i_x_free_still[Q_R_itr-Q_R_begin] = Q_R_itr->row_i();
            ++Q_R_itr;
            ++i_x_free_itr;
          }
          else if( Q_R_itr->row_i() < *i_x_free_itr ) {
            // A variable that was in the superbasis now is not
            i_x_free_to_fixed[num_free_to_fixed] = Q_R_itr->row_i();
            num_free_to_fixed++;
            ++Q_R_itr;
          }
          else {
            // A variable that was not in the superbasis now is
            i_x_fixed_to_free[num_fixed_to_free] = *i_x_free_itr;
            num_fixed_to_free++;
            ++i_x_free_itr;
          }
        }
      }
    }
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "\nThere will be " << num_fixed_to_free  << " independent variables added to the superbasis and rHL_RR";
      if( num_fixed_to_free && int(olevel) >= int(PRINT_ACTIVE_SET) ) {
        out << " and their indexes are:\n";
        for(size_type k = 0; k < num_fixed_to_free; ++k)
          out << " " << i_x_fixed_to_free[k];
        out << std::endl;
      }
      else {
        out << "\n";
      }
      out << "\nThere will be " << num_free_to_fixed  << " independent variables removed from the superbasis and rHL_RR";
      if( num_free_to_fixed && int(olevel) >= int(PRINT_ACTIVE_SET) ) {
        out << " and their indexes are:\n";
        for(size_type k = 0; k < num_free_to_fixed; ++k)
          out << " " << i_x_free_to_fixed[k];
        out << std::endl;
      }
      else {
        out << "\n";
      }
    }
    // Get reference to rHL_RR = B_RR
#ifdef _WINDOWS
    MatrixSymAddDelUpdateable
      &rHL_RR = dynamic_cast<MatrixSymAddDelUpdateable&>(
        const_cast<MatrixSymWithOpFactorized&>(*rHL_super.B_RR_ptr())
        );
#else
    MatrixSymAddDelUpdateable
      &rHL_RR = dyn_cast<MatrixSymAddDelUpdateable>(
        const_cast<MatrixSymWithOpFactorized&>(*rHL_super.B_RR_ptr())
        );
#endif
    // Remove rows/cols from rHL_RR.
    if( num_free_to_fixed ) {
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
        out << "\nDeleting " << num_free_to_fixed << " rows/columns from rHL_RR ...\n";
      }
      {for( size_type k = i_x_free_still.size(); k > 0; --k ) { // Delete from the largest to the smallest index (cheaper!)
        if( !i_x_free_still[k-1] )
          rHL_RR.delete_update( k, false );
      }}
      TEUCHOS_TEST_FOR_EXCEPT( !(  rHL_RR.rows() == n_pz_R - num_fixed_to_free  ) );
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
        out << "\nrHL_RR after rows/columns where removed =\n" << *rHL_super.B_RR_ptr();
      }
    }
    // Add new rows/cols to rHL_RR.
    if( num_fixed_to_free ) {
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
        out << "\nAppending " << num_fixed_to_free << " rows/columns to rHL_RR ...\n";
      }
      {for( size_type k = 0; k < num_fixed_to_free; ++k ) {
        rHL_RR.augment_update( NULL, rHL_XX_scale, false );
      }}
      TEUCHOS_TEST_FOR_EXCEPT( !(  rHL_RR.rows() == n_pz_R  ) );
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
        out << "\nrHL_RR after rows/columns where appended =\n" << *rHL_super.B_RR_ptr();
      }
    }
    // Resort i_x_free[] to reflect the actual order of the indices in rHL_RR
    {
      size_type tmp_n_pz_R = 0;
      {for(size_type k = 0; k < i_x_free_still.size(); ++k) {
        if( i_x_free_still[k] ) {
          i_x_free[tmp_n_pz_R] = i_x_free_still[k];
          ++tmp_n_pz_R;
        }
      }}
      {for(size_type k = 0; k < num_fixed_to_free; ++k) {
        i_x_free[tmp_n_pz_R] = i_x_fixed_to_free[k];
        ++tmp_n_pz_R;
      }}
      TEUCHOS_TEST_FOR_EXCEPT( !(  tmp_n_pz_R == n_pz_R  ) );
    }
    // Initialize rHL_XX = rHL_XX_scale * I resized to the proper dimensions
    rHL_XX.init_identity(n_pz_X,rHL_XX_scale);
    // Reinitalize rHL for new active set
    rHL_super.initialize(
      n_pz, n_pz_R, &i_x_free[0], &i_x_fixed[0], &bnd_fixed[0]
      ,rHL_super.B_RR_ptr(),NULL,BLAS_Cpp::no_trans,rHL_super.B_XX_ptr()
      );
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
      out << "\nFull rHL after reinitialization but before BFGS update:\n"
        << "\nrHL =\n" << *rHL_k
        << "\nQ_R =\n" << rHL_super.Q_R()
        << "\nQ_X =\n" << rHL_super.Q_X();
    }
    // Now we will do the PBFGS updating from now on!
    do_projected_rHL_RR = true;
  }
  //
  // Perform the BFGS update
  //
  if( do_projected_rHL_RR ) {
    // Perform BFGS update on smaller rHL_RR.
    // By the time we get here rHL_RR should be resize and ready to update
    const GenPermMatrixSlice
      &Q_R = rHL_super.Q_R(),
      &Q_X = rHL_super.Q_X();
    const size_type
      n_pz_R = Q_R.cols(),
      n_pz_X = Q_X.cols();
    TEUCHOS_TEST_FOR_EXCEPT( !(  n_pz_R + n_pz_X == n_pz  ) );
    // Get projected BFGS update vectors y_bfgs_R, s_bfgs_R
    Workspace<value_type>
      y_bfgs_R_ws(wss,Q_R.cols()),
      s_bfgs_R_ws(wss,Q_R.cols());
    DVectorSlice y_bfgs_R(&y_bfgs_R_ws[0],y_bfgs_R_ws.size());
    DVectorSlice s_bfgs_R(&s_bfgs_R_ws[0],s_bfgs_R_ws.size());
    V_MtV( &y_bfgs_R, Q_R, BLAS_Cpp::trans, *y_bfgs );  // y_bfgs_R = Q_R'*y_bfgs
    V_MtV( &s_bfgs_R, Q_R, BLAS_Cpp::trans, *s_bfgs );  // s_bfgs_R = Q_R'*s_bfgs
    // Update rHL_RR
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "\nPerform BFGS update on " << n_pz_R << " x " << n_pz_R << " projected reduced Hessian for the superbasic variables where B = rHL_RR...\n";
    }
    bfgs_update().perform_update(
      &s_bfgs_R(),&y_bfgs_R(),first_update,out,olevel,algo->algo_cntr().check_results()
      ,const_cast<MatrixOp*>(static_cast<const MatrixOp*>(rHL_super.B_RR_ptr().get()))
      ,&quasi_newton_stats_(*s).set_k(0)
      );
  }
  else {
    // Update the full reduced Hessain matrix (rHL = rHL_RR)
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "\nPerform BFGS update on the full reduced Hessian where B = rHL...\n";
    }
    bfgs_update().perform_update(
      s_bfgs,y_bfgs,first_update,out,olevel,algo->algo_cntr().check_results()
      ,const_cast<MatrixOp*>(static_cast<const MatrixOp*>(rHL_super.B_RR_ptr().get()))
      ,&quasi_newton_stats_(*s).set_k(0)
      );
  }

  return true;
}

void ReducedHessianSecantUpdateBFGSProjected_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Perform BFGS update on only free independent (superbasic) variables.\n"
    << L << "if s'*y < sqrt(macheps) * ||s||2 * ||y||2 and use_dampening == false then"
    << L << "    Skip the update and exit this step!\n"
    << L << "end\n"
    << L << "rHL_XX_scale = max(y'*y/s'*y,1.0)\n"
    << L << "do_projected_rHL_RR = false\n"
    << L << "if nu_km1 is updated then\n"
    << L << "    if rHL_k.Q_R is the identity matrix then\n"
    << L << "        *** Determine if the active set has calmed down enough\n"
    << L << "        Given (num_active_indep,num_adds_indep,num_drops_indep) from act_set_stats_km1\n"
    << L << "        fact_same =\n"
    << L << "            ( num_adds_indep== NOT_KNOWN || num_drops_indep==NOT_KNOWN\n"
    << L << "            || num_active_indep==0\n"
    << L << "                ? 0.0\n"
    << L << "                : std::_MAX((num_active_indep-num_adds_indep-num_drops_indep)\n"
    << L << "                    /num_active_indep,0.0)\n"
    << L << "            )\n"
    << L << "        do_projected_rHL_RR\n"
    << L << "            = fact_same >= act_set_frac_proj_start && num_active_indep > 0\n"
    << L << "        if do_projected_rHL_RR == true then\n"
    << L << "            Determine the sets of superbasic variables given the mapping matrix\n"
    << L << "            Q = [ Q_R, Q_X ] where pz_R = Q_R'*pz <: R^n_pz_R are the superbasic variables and\n"
    << L << "            pz_X = Q_X'*pz <: R^n_pz_X are the nonbasic variables that only contain fixed\n"
    << L << "            variables in nu_km1(indep) where the following condidtions are satisfied:\n"
    << L << "                (s'*Q_X*rHL_XX_scale*Q_X'*s)/(s'*Q_R*Q_R'*rHL_k*Q_R*Q_R'*s) <= project_error_tol\n"
    << L << "                |s'*Q_X*Q_X'*y|/|s'*Q_R*Q_R'*s| <= project_error_tol\n"
    << L << "                |Q_X'*nu_km1(indep)|/||nu_km1(indep)||inf >= super_basic_mult_drop_tol\n"
    << L << "            if n_pz_R < n-r then\n"
    << L << "                Delete rows/columns of rHL_k to form rHL_RR = Q_R'*rHL_k*Q_R\n"
    << L << "                Define new rHL_k = [ Q_R, Q_X ] * [ rHL_RR, 0; 0; rHL_XX_scale*I ] [ Q_R'; Q_X ]\n"
    << L << "            else\n"
    << L << "                do_projected_rHL_RR = false\n"
    << L << "            end\n"
    << L << "        end\n"
    << L << "    else\n"
    << L << "        Determine the new Q_n = [ Q_R_n, Q_X_n ] that satisfies:\n"
    << L << "            (s'*Q_X_n*rHL_XX_scale*Q_X_n'*s)/(s'*Q_R_n*Q_R_n'*rHL_k*Q_R_n*Q_R_n'*s) <= project_error_tol\n"
    << L << "            |s'*Q_X_n*Q_X_n'*y|/|s'*Q_R_n*Q_R_n'*s| <= project_error_tol\n"
    << L << "            |Q_X_n'*nu_km1(indep)|/||nu_km1(indep)||inf >= super_basic_mult_drop_tol\n"
    << L << "        Remove rows/cols from rHL_k.rHL_RR for variables in rHL_k.Q_R that are not in Q_R_n.\n"
    << L << "        Add digonal entries equal to rHL_XX_scale to rHL_k.rHL_RR for variables in Q_R_n\n"
    << L << "        that are not in rHL_k.Q_R\n"
    << L << "        Define new rHL_k = [ Q_R_n, Q_X_n ] * [ rHL_k.rHL_RR, 0; 0; rHL_XX_scale*I ] [ Q_R_n'; Q_X_n ]\n"
    << L << "        do_projected_rHL_RR = true\n"
    << L << "    end\n"
    << L << "end\n"
    << L << "if do_projected_rHL_RR == true then\n"
    << L << "    Perform projected BFGS update (see below): (rHL_k.rHL_RR, Q_R_n'*s, Q_R_n'*y) -> rHL_k.rHL_RR\n"
    << L << "else\n"
    << L << "    Perform full BFGS update: (rHL_k, s, y) -> rHL_k\n"
    << L << "    begin BFGS update where B = rHL_k\n";
  bfgs_update().print_step( out, L + "        " );
  out
    << L << "    end BFGS update\n"
    << L << "else\n"
    ;
}

}  // end namespace MoochoPack

#endif // 0
