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

#include "MoochoPack_ReducedHessianSecantUpdateLPBFGS_Strategy.hpp"
#include "MoochoPack_PBFGS_helpers.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "ConstrainedOptPack_MatrixSymPosDefLBFGS.hpp"
#include "ConstrainedOptPack/src/AbstractLinAlgPack_MatrixSymPosDefCholFactor.hpp"
#include "ConstrainedOptPack/src/AbstractLinAlgPack_BFGS_helpers.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymInitDiag.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Midynamic_cast_verbose.h"
#include "MiWorkspacePack.h"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

namespace MoochoPack {

ReducedHessianSecantUpdateLPBFGS_Strategy::ReducedHessianSecantUpdateLPBFGS_Strategy(
  const proj_bfgs_updater_ptr_t&  proj_bfgs_updater
  ,size_type                      min_num_updates_proj_start
  ,size_type                      max_num_updates_proj_start
  ,size_type                      num_superbasics_switch_dense
  ,size_type                      num_add_recent_updates
  )
  : proj_bfgs_updater_(proj_bfgs_updater)
  , min_num_updates_proj_start_(min_num_updates_proj_start)
  , max_num_updates_proj_start_(max_num_updates_proj_start)
  , num_superbasics_switch_dense_(num_superbasics_switch_dense)
  , num_add_recent_updates_(num_add_recent_updates)
{}

bool ReducedHessianSecantUpdateLPBFGS_Strategy::perform_update(
  DVectorSlice* s_bfgs, DVectorSlice* y_bfgs, bool first_update
  ,std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
  ,MatrixOp *rHL_k
  )
{
  using std::setw;
  using std::endl;
  using std::right;
  using Teuchos::dyn_cast;
  namespace rcp = MemMngPack;
  using Teuchos::RCP;
  using LinAlgOpPack::V_MtV;
  using DenseLinAlgPack::dot;
  using AbstractLinAlgPack::norm_inf;
  using AbstractLinAlgPack::transVtMtV;
  typedef ConstrainedOptPack::MatrixHessianSuperBasic MHSB_t;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\n*** (LPBFGS) Full limited memory BFGS to projected BFGS ...\n";
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

  bool do_projected_rHL_RR = false;

  // See if we still have a limited memory BFGS update matrix
  RCP<MatrixSymPosDefLBFGS> // We don't want this to be deleted until we are done with it
    lbfgs_rHL_RR = Teuchos::rcp_const_cast<MatrixSymPosDefLBFGS>(
      Teuchos::rcp_dynamic_cast<const MatrixSymPosDefLBFGS>(rHL_super.B_RR_ptr()) );

  if( lbfgs_rHL_RR.get() && rHL_super.Q_R().is_identity()  ) {
    //
    // We have a limited memory BFGS matrix and have not started projected BFGS updating
    // yet so lets determine if it is time to consider switching
    //
    // Check that the current update is sufficiently p.d. before we do anything
    const value_type
      sTy = dot(*s_bfgs,*y_bfgs),
      yTy = dot(*y_bfgs,*y_bfgs);
    if( !ConstrainedOptPack::BFGS_sTy_suff_p_d(
      *s_bfgs,*y_bfgs,&sTy
      ,  int(olevel) >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL )
      && !proj_bfgs_updater().bfgs_update().use_dampening()
      )
    {
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
        out	<< "\nWarning!  use_damening == false so there is no way we can perform any"
            " kind of BFGS update (projected or not) so we will skip it!\n";
      }
      quasi_newton_stats_(*s).set_k(0).set_updated_stats(
        QuasiNewtonStats::INDEF_SKIPED );
      return true;
    }
    // Consider if we can even look at the active set yet.
    const bool consider_switch  = lbfgs_rHL_RR->num_secant_updates() >= min_num_updates_proj_start();
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "\nnum_previous_updates = " << lbfgs_rHL_RR->num_secant_updates()
        << ( consider_switch ? " >= " : " < " )
        << "min_num_updates_proj_start = " << min_num_updates_proj_start()
        << ( consider_switch
           ? "\nConsidering if we should switch to projected BFGS updating of superbasics ...\n"
           : "\nNot time to consider switching to projected BFGS updating of superbasics yet!" );
    }
    if( consider_switch ) {
      // 
      // Our job here is to determine if it is time to switch to projected projected
      // BFGS updating.
      //
      if( act_set_stats_(*s).updated_k(-1) ) {
        if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
          out	<< "\nDetermining if projected BFGS updating of superbasics should begin ...\n";
        }
        // Determine if the active set has calmed down enough
        const SpVector
          &nu_km1 = s->nu().get_k(-1);
        const SpVectorSlice
          nu_indep = nu_km1(s->var_indep());
        const bool 
          act_set_calmed_down
          = PBFGSPack::act_set_calmed_down(
            act_set_stats_(*s).get_k(-1)
            ,proj_bfgs_updater().act_set_frac_proj_start()
            ,olevel,out
            ),
          max_num_updates_exceeded
          = lbfgs_rHL_RR->m_bar() >= max_num_updates_proj_start();
        do_projected_rHL_RR = act_set_calmed_down || max_num_updates_exceeded;
        if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
          if( act_set_calmed_down ) {
            out << "\nThe active set has calmed down enough so lets further consider switching to\n"
              << "projected BFGS updating of superbasic variables ...\n";
          }
          else if( max_num_updates_exceeded ) {
            out << "\nThe active set has not calmed down enough but num_previous_updates = "
              << lbfgs_rHL_RR->m_bar() << " >= max_num_updates_proj_start = "	<< max_num_updates_proj_start()
              << "\nso we will further consider switching to projected BFGS updating of superbasic variables ...\n";
          }
          else {
            out << "\nIt is not time to switch to projected BFGS so just keep performing full limited memory BFGS for now ...\n";
          }
        }
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
          // rHL = rHL_scale * I
          value_type
            rHL_scale = sTy > 0.0 ? yTy/sTy : 1.0;
          if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
            out	<< "\nScaling for diagonal intitial rHL = rHL_scale*I, rHL_scale = " << rHL_scale << std::endl;
          }
          value_type sRTBRRsR = 0.0, sRTyR = 0.0, sXTBXXsX = 0.0, sXTyX = 0.0;
          // Get the elements in i_x_free[] for variables that are definitely free
          // and initialize s_R'*B_RR*s_R and s_R'*y_R
          PBFGSPack::init_i_x_free_sRTsR_sRTyR(
            nu_indep, *s_bfgs, *y_bfgs
            , &n_pz_R, &i_x_free[0], &sRTBRRsR, &sRTyR );
          sRTBRRsR *= rHL_scale;
          Workspace<value_type> rHL_XX_diag_ws(wss,nu_indep.nz());
          DVectorSlice rHL_XX_diag(&rHL_XX_diag_ws[0],rHL_XX_diag_ws.size());
          rHL_XX_diag = rHL_scale;
          // Sort fixed variables according to |s_X(i)^2*B_XX(i,i)|/|sRTBRRsR| + |s_X(i)*y_X(i)|/|sRTyR|
          // and initialize s_X'*B_XX*s_X and s_X*y_X
          PBFGSPack::sort_fixed_max_cond_viol(
            nu_indep,*s_bfgs,*y_bfgs,rHL_XX_diag,sRTBRRsR,sRTyR
            ,&sXTBXXsX,&sXTyX,&l_x_fixed_sorted[0]);
          // Pick initial set of i_x_free[] and i_x_fixed[] (sorted!)
          PBFGSPack::choose_fixed_free(
            proj_bfgs_updater().project_error_tol()
            ,proj_bfgs_updater().super_basic_mult_drop_tol(),nu_indep
            ,*s_bfgs,*y_bfgs,rHL_XX_diag,&l_x_fixed_sorted[0]
            ,olevel,out,&sRTBRRsR,&sRTyR,&sXTBXXsX,&sXTyX
            ,&n_pz_X,&n_pz_R,&i_x_free[0],&i_x_fixed[0],&bnd_fixed[0] );
          if( n_pz_R < n_pz ) {
            //
            // We are ready to transition to projected BFGS updating!
            //
            // Determine if we are be using dense or limited memory BFGS?
            const bool
              low_num_super_basics = n_pz_R <= num_superbasics_switch_dense();
            if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
              out	<< "\nSwitching to projected BFGS updating ..."
                << "\nn_pz_R = " << n_pz_R << ( low_num_super_basics ? " <= " : " > " )
                << " num_superbasics_switch_dense = " << num_superbasics_switch_dense()
                << ( low_num_super_basics
                   ? "\nThere are not too many superbasic variables so use dense projected BFGS ...\n"
                   : "\nThere are too many superbasic variables so use limited memory projected BFGS ...\n"
                  );
            }
            // Create new matrix to use for rHL_RR initialized to rHL_RR = rHL_scale*I
            RCP<MatrixSymSecant>
              rHL_RR = NULL;
            if( low_num_super_basics ) {
              rHL_RR = new MatrixSymPosDefCholFactor(
                NULL    // Let it allocate its own memory
                ,NULL   // ...
                ,n_pz   // maximum size
                ,lbfgs_rHL_RR->maintain_original()
                ,lbfgs_rHL_RR->maintain_inverse()
                );
            }
            else {
              rHL_RR = new MatrixSymPosDefLBFGS(
                n_pz, lbfgs_rHL_RR->m()
                ,lbfgs_rHL_RR->maintain_original()
                ,lbfgs_rHL_RR->maintain_inverse()
                ,lbfgs_rHL_RR->auto_rescaling()
                );
            }
            rHL_RR->init_identity( n_pz_R, rHL_scale );
            if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
              out << "\nrHL_RR after intialized to rHL_RR = rHL_scale*I but before performing current BFGS update:\nrHL_RR =\n"
                << dynamic_cast<MatrixOp&>(*rHL_RR); // I know this is okay!
            }
            // Initialize rHL_XX = rHL_scale*I
#ifdef _WINDOWS
            MatrixSymInitDiag
              &rHL_XX = dynamic_cast<MatrixSymInitDiag&>(
                const_cast<MatrixSymOp&>(*rHL_super.B_XX_ptr()));
#else
            MatrixSymInitDiag
              &rHL_XX = dyn_cast<MatrixSymInitDiag>(
                const_cast<MatrixSymOp&>(*rHL_super.B_XX_ptr()));
#endif
            rHL_XX.init_identity( n_pz_X, rHL_scale );
            // Reinitialize rHL
            rHL_super.initialize(
              n_pz, n_pz_R, &i_x_free[0], &i_x_fixed[0], &bnd_fixed[0]
              ,Teuchos::rcp_const_cast<const MatrixSymWithOpFactorized>(
                Teuchos::rcp_dynamic_cast<MatrixSymWithOpFactorized>(rHL_RR))
              ,NULL,BLAS_Cpp::no_trans,rHL_super.B_XX_ptr()
              );
            //
            // Perform the current BFGS update first
            //
            MatrixSymOp
              &rHL_RR_op = dynamic_cast<MatrixSymOp&>(*rHL_RR);
            const GenPermMatrixSlice
              &Q_R = rHL_super.Q_R(),
              &Q_X = rHL_super.Q_X();
            // Get projected BFGS update vectors y_bfgs_R, s_bfgs_R
            Workspace<value_type>
              y_bfgs_R_ws(wss,Q_R.cols()),
              s_bfgs_R_ws(wss,Q_R.cols()),
              y_bfgs_X_ws(wss,Q_X.cols()),
              s_bfgs_X_ws(wss,Q_X.cols());
            DVectorSlice y_bfgs_R(&y_bfgs_R_ws[0],y_bfgs_R_ws.size());
            DVectorSlice s_bfgs_R(&s_bfgs_R_ws[0],s_bfgs_R_ws.size());
            DVectorSlice y_bfgs_X(&y_bfgs_X_ws[0],y_bfgs_X_ws.size());
            DVectorSlice s_bfgs_X(&s_bfgs_X_ws[0],s_bfgs_X_ws.size());
            V_MtV( &y_bfgs_R, Q_R, BLAS_Cpp::trans, *y_bfgs );  // y_bfgs_R = Q_R'*y_bfgs
            V_MtV( &s_bfgs_R, Q_R, BLAS_Cpp::trans, *s_bfgs );  // s_bfgs_R = Q_R'*s_bfgs
            V_MtV( &y_bfgs_X, Q_X, BLAS_Cpp::trans, *y_bfgs );  // y_bfgs_X = Q_X'*y_bfgs
            V_MtV( &s_bfgs_X, Q_X, BLAS_Cpp::trans, *s_bfgs );  // s_bfgs_X = Q_X'*s_bfgs
            // Update rHL_RR
            if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
              out << "\nPerform current BFGS update on " << n_pz_R << " x " << n_pz_R << " projected "
                << "reduced Hessian for the superbasic variables where B = rHL_RR...\n";
            }
            QuasiNewtonStats quasi_newton_stats;
            proj_bfgs_updater().bfgs_update().perform_update(
              &s_bfgs_R(),&y_bfgs_R(),false,out,olevel,algo->algo_cntr().check_results()
              ,&rHL_RR_op, &quasi_newton_stats );
            // Perform previous updates if possible
            if( lbfgs_rHL_RR->m_bar() ) {
              const size_type num_add_updates = std::_MIN(num_add_recent_updates(),lbfgs_rHL_RR->m_bar());
              if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
                out	<< "\nAdd the min(num_previous_updates,num_add_recent_updates) = min(" << lbfgs_rHL_RR->m_bar()
                  << "," << num_add_recent_updates() << ") = " << num_add_updates << " most recent BFGS updates if possible ...\n";
              }
              // Now add previous update vectors
              const value_type
                project_error_tol = proj_bfgs_updater().project_error_tol();
              const DMatrixSlice
                S = lbfgs_rHL_RR->S(),
                Y = lbfgs_rHL_RR->Y();
              size_type k = lbfgs_rHL_RR->k_bar();  // Location in S and Y of most recent update vectors
              for( size_type l = 1; l <= num_add_updates; ++l, --k ) {
                if(k == 0) k = lbfgs_rHL_RR->m_bar();  // see MatrixSymPosDefLBFGS
                // Check to see if this update satisfies the required conditions.
                // Start with the condition on s'*y since this are cheap to check.
                V_MtV( &s_bfgs_X(), Q_X, BLAS_Cpp::trans, S.col(k) ); // s_bfgs_X = Q_X'*s_bfgs
                V_MtV( &y_bfgs_X(), Q_X, BLAS_Cpp::trans, Y.col(k) ); // y_bfgs_X = Q_X'*y_bfgs
                sRTyR    = dot( s_bfgs_R, y_bfgs_R );
                sXTyX    = dot( s_bfgs_X, y_bfgs_X );
                bool
                  sXTyX_cond    = ::fabs(sXTyX/sRTyR) <= project_error_tol,
                  do_update     = sXTyX_cond,
                  sXTBXXsX_cond = false;
                if( sXTyX_cond ) {
                  // Check the second more expensive condition
                  V_MtV( &s_bfgs_R(), Q_R, BLAS_Cpp::trans, S.col(k) ); // s_bfgs_R = Q_R'*s_bfgs
                  V_MtV( &y_bfgs_R(), Q_R, BLAS_Cpp::trans, Y.col(k) ); // y_bfgs_R = Q_R'*y_bfgs
                  sRTBRRsR = transVtMtV( s_bfgs_R, rHL_RR_op, BLAS_Cpp::no_trans, s_bfgs_R );
                  sXTBXXsX = rHL_scale * dot( s_bfgs_X, s_bfgs_X );
                  sXTBXXsX_cond = sXTBXXsX/sRTBRRsR <= project_error_tol;
                  do_update     = sXTBXXsX_cond && sXTyX_cond;
                }
                if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
                  out << "\n---------------------------------------------------------------------"
                    << "\nprevious update " << l
                    << "\n\nChecking projection error:\n"
                    << "\n|s_X'*y_X|/|s_R'*y_R| = |" << sXTyX << "|/|" << sRTyR
                    << "| = " << ::fabs(sXTyX/sRTyR)
                    << ( sXTyX_cond ? " <= " : " > " ) << " project_error_tol = "
                    << project_error_tol;
                  if( sXTyX_cond ) {
                    out	<< "\n(s_X'*rHL_XX*s_X/s_R'*rHL_RR*s_R) = (" << sXTBXXsX
                        << ") = " << (sXTBXXsX/sRTBRRsR)
                        << ( sXTBXXsX_cond ? " <= " : " > " ) << " project_error_tol = "
                        << project_error_tol;
                  }
                  out << ( do_update
                    ? "\n\nAttemping to add this previous update where B = rHL_RR ...\n"
                    : "\n\nCan not add this previous update ...\n" );
                }
                if( do_update ) {
                  // ( rHL_RR, s_bfgs_R, y_bfgs_R ) -> rHL_RR (this should not throw an exception!)
                  try {
                    proj_bfgs_updater().bfgs_update().perform_update(
                      &s_bfgs_R(),&y_bfgs_R(),false,out,olevel,algo->algo_cntr().check_results()
                      ,&rHL_RR_op, &quasi_newton_stats );
                  }
                  catch( const MatrixSymSecant::UpdateSkippedException& excpt ) {
                    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
                      out	<< "\nOops!  The " << l << "th most recent BFGS update was rejected?:\n"
                        << excpt.what() << std::endl;
                    }
                  }
                }
              }
              if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
                out << "\n---------------------------------------------------------------------\n";
                }
              if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
                out << "\nrHL_RR after adding previous BFGS updates:\nrHL_BRR =\n"
                  << dynamic_cast<MatrixOp&>(*rHL_RR);
              }
            }
            else {
              if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
                out	<< "\nThere were no previous update vectors to add!\n";
              }
            }
            if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
              out << "\nFull rHL after complete reinitialization:\nrHL =\n" << *rHL_k;
            }
            quasi_newton_stats_(*s).set_k(0).set_updated_stats(
              QuasiNewtonStats::REINITIALIZED );
            return true;
          }
          else {
            if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
              out << "\nn_pz_R == n_pz = " << n_pz_R << ", No variables would be removed from "
                << "the superbasis so just keep on performing limited memory BFGS for now ...\n";
            }
            do_projected_rHL_RR = false;
          }
        }
      }
    }
    // If we have not switched to PBFGS then just update the full limited memory BFGS matrix
    if(!do_projected_rHL_RR) {
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
        out << "\nPerform current BFGS update on " << n_pz << " x " << n_pz << " full "
          << "limited memory reduced Hessian B = rHL ...\n";
      }
      proj_bfgs_updater().bfgs_update().perform_update(
        s_bfgs,y_bfgs,first_update,out,olevel,algo->algo_cntr().check_results()
        ,lbfgs_rHL_RR.get()
        ,&quasi_newton_stats_(*s).set_k(0)
        );
      return true;
    }
  }
  else {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out	<< "\nWe have already switched to projected BFGS updating ...\n";
    }
  }
  //
  // If we get here then we must have switched to
  // projected updating so lets just pass it on!
  //
  return proj_bfgs_updater().perform_update(
    s_bfgs,y_bfgs,first_update,out,olevel,algo,s,rHL_k);
}

void ReducedHessianSecantUpdateLPBFGS_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Perform limited memory LBFGS updating initially then switch to dense\n"
    << L << "*** projected BFGS updating when appropriate.\n"
    << L << "ToDo: Finish implementation!\n";
}

}  // end namespace MoochoPack

#endif // 0
