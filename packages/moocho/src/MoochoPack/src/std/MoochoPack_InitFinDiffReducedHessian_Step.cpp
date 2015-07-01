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

#include <math.h>

#include <ostream>

#include "MoochoPack_InitFinDiffReducedHessian_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "NLPInterfacePack_NLPObjGrad.hpp"
#include "AbstractLinAlgPack_MatrixSymInitDiag.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
//#include "AbstractLinAlgPack_SpVectorClass.hpp"
//#include "AbstractLinAlgPack/src/max_near_feas_step.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
} // end namespace

namespace MoochoPack {

InitFinDiffReducedHessian_Step::InitFinDiffReducedHessian_Step(
  EInitializationMethod   initialization_method
  ,value_type             max_cond
  ,value_type             min_diag
  ,value_type             step_scale
  )
  :initialization_method_(initialization_method)
  ,max_cond_(max_cond)
  ,min_diag_(min_diag)
  ,step_scale_(step_scale)
{}

bool InitFinDiffReducedHessian_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  using Teuchos::dyn_cast;
  using LinAlgOpPack::Vt_S;
  using LinAlgOpPack::Vp_StV;
  using LinAlgOpPack::V_MtV;
  using AbstractLinAlgPack::max_near_feas_step;

  NLPAlgo          &algo  = rsqp_algo(_algo);
  NLPAlgoState         &s     = algo.rsqp_state();
  NLPObjGrad    &nlp   = dyn_cast<NLPObjGrad>(algo.nlp());
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  // Get the iteration quantity container objects
  IterQuantityAccess<index_type>
    &num_basis_iq = s.num_basis();

  const bool new_basis     = num_basis_iq.updated_k(0);
  const int  k_last_offset = s.rHL().last_updated();

  // If the basis has changed or there is no previous matrix to use
  // then reinitialize.

  if( new_basis || k_last_offset == IterQuantity::NONE_UPDATED ) {

    // Compute a finite difference along the null space of the
    // constraints

    IterQuantityAccess<VectorMutable>
      &x_iq     = s.x(),
      &rGf_iq   = s.rGf();
    IterQuantityAccess<MatrixOp>
      &Z_iq     = s.Z();
    IterQuantityAccess<MatrixSymOp>
      &rHL_iq   = s.rHL();

    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\nReinitializing the reduced Hessain using a finite difference\n";
    }

    MatrixSymInitDiag &rHL_diag = dyn_cast<MatrixSymInitDiag>(rHL_iq.set_k(0));
    const MatrixOp    &Z_k      = Z_iq.get_k(0);
    const Vector    &x_k      = x_iq.get_k(0);
    const Vector    &rGf_k    = rGf_iq.get_k(0);

    // one vector
    VectorSpace::vec_mut_ptr_t  e = Z_k.space_rows().create_member(1.0);

    // Ze
    VectorSpace::vec_mut_ptr_t Ze = x_k.space().create_member();
    V_MtV( Ze.get(), Z_k, BLAS_Cpp::no_trans, *e );
    
    // This does not have to be an accurate finite difference so lets just
    // take step_scale/||Ze|| as the step size unless it is out of bounds
    // If we assume that our variables are scaled near
    // one (step_scale == 1?) then this will give us an appreciable step.  Beside we
    // should not be near the solution so the reduced gradient should not
    // be near zero.

    const value_type nrm_Ze = Ze->norm_inf();
    value_type u = step_scale() / nrm_Ze;

    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\n||Ze||inf = " << nrm_Ze << std::endl;
    }

    if( (int)olevel >= (int)PRINT_VECTORS ) {
      out << "\nZe =\n" << *Ze;
    }

    if( algo.nlp().num_bounded_x() ) {

      // Find the maximum step u
      // we can take along x_fd = x_k + u*Ze
      // that don't violate variable bounds by too much.
      // If a positive step can't be found then this may be a negative step.
      
      std::pair<value_type,value_type>
        u_steps = max_near_feas_step(
          x_k, *Ze
          ,nlp.xl(), nlp.xu()
          ,nlp.max_var_bounds_viol()
          );
      
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nMaximum steps ( possitive, negative ) in bounds u = ("
          << u_steps.first << "," << u_steps.second << ")\n";
      }

      if( u_steps.first < u )
        u = u_steps.first;
      if( ::fabs(u_steps.second) < u )
        u = u_steps.second;
    }

    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\nFinite difference step length u = " << u << std::endl;
    }

    // Take the finite difference from x_k to x_fd = x_k + u * Ze
    //
    // rGf_fd = ( Z_k'*Gf(x_k + u*Ze) - rGf_k ) / u
    //

    VectorSpace::vec_mut_ptr_t x_fd = x_k.space().create_member();
    Vp_StV( x_fd.get(), u, *Ze );

    // Gf_fd = Gf(x_fd)
    VectorSpace::vec_mut_ptr_t Gf_fd = x_k.space().create_member();
    nlp.unset_quantities();
    nlp.set_Gf(	Gf_fd.get() );
    nlp.calc_Gf( *x_fd );

    if( (int)olevel >= (int)PRINT_VECTORS ) {
      out << "\nGf_fd =\n" << *Gf_fd;
    }

    // rGf_fd = Z'*Gf_fd
    VectorSpace::vec_mut_ptr_t rGf_fd = Z_k.space_rows().create_member();
    V_MtV( rGf_fd.get(), Z_k, BLAS_Cpp::trans, *Gf_fd );

    // rGf_fd = rGf_fd - rGf_k
    Vp_StV( rGf_fd.get(), -1.0, rGf_k );

    // rGf_fd = rGf_fd / u
    Vt_S( rGf_fd.get(), 1.0 / u );

    const value_type
      nrm_rGf_fd = rGf_fd->norm_inf();

    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\n||(rGf_fd - rGf_k)/u||inf = " << nrm_rGf_fd << std::endl;
    }
    if( (int)olevel >= (int)PRINT_VECTORS ) {
      out << "\n(rGf_fd - rGf_k)/u =\n" << *rGf_fd;
    }

    if( nrm_rGf_fd <= min_diag() ) {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\n||(rGf_fd - rGf_k)/u||inf = " << nrm_rGf_fd
            << " < min_diag = " << min_diag() << std::endl
          << "\nScale by min_diag ... \n";
      }
      rHL_diag.init_identity(Z_k.space_rows(),min_diag());
    }
    else {
      switch( initialization_method() ) {
        case SCALE_IDENTITY: {
          if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
            out << "\nScale the identity matrix by ||(rGf_fd - rGf_k)/u||inf ... \n";
          }
          rHL_diag.init_identity(Z_k.space_rows(),nrm_rGf_fd);
          break;
        }
        case SCALE_DIAGONAL:
        case SCALE_DIAGONAL_ABS:
        {
          if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
            out << "\nScale diagonal by modified finite difference ... \n";
          }
          // In order to keep the initial reduced Hessian well conditioned we
          // will not let any diagonal element drop below
          // ||rGf_fd||inf / max_cond
        
          const value_type
            min_ele = my_max( nrm_rGf_fd / max_cond(), min_diag() );

          if( initialization_method() == SCALE_DIAGONAL )
            AbstractLinAlgPack::max_vec_scalar( min_ele, rGf_fd.get() );
          else
            AbstractLinAlgPack::max_abs_vec_scalar( min_ele, rGf_fd.get() );

          if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
            out << "\n||diag||inf = " << rGf_fd->norm_inf() << std::endl;
          }
          if( (int)olevel >= (int)PRINT_VECTORS ) {
            out << "\ndiag =\n" << *rGf_fd;
          }
          rHL_diag.init_diagonal(*rGf_fd);
          break;
        }
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true);	// only local programming error?
      }
    }
    nlp.unset_quantities();

    quasi_newton_stats_(s).set_k(0).set_updated_stats(
      QuasiNewtonStats::REINITIALIZED );

    if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
      out << "\nrHL_k =\n" << rHL_iq.get_k(0);
    }

  }

  return true;
}

void InitFinDiffReducedHessian_Step::print_step(
  const Algorithm& algo
  ,poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  ,std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Initialize the reduced Hessian using a single finite difference.\n"
    << L << "*** Where the nlp must support the NLPObjGrad interface and\n"
    << L << "*** rHL_k must support the MatrixSymInitDiag interface or exceptions\n"
    << L << "*** will be thrown.\n"
    << L << "default: num_basis_remembered = NO_BASIS_UPDATED_YET\n"
    << L << "         initialization_method = SCALE_DIAGONAL\n"
    << L << "         max_cond              = " << max_cond() << std::endl
    << L << "         min_diag              = " << min_diag() << std::endl
    << L << "         step_scale            = " << step_scale() << std::endl
    << L << "if num_basis_k is updated then\n"
    << L << "    new_basis = true\n"
    << L << "else\n"
    << L << "    new_basis = false\n"
    << L << "end\n"
    << L << "if new_basis == true or no past rHL as been updated then\n"
    << L << "    *** Reinitialize the reduced Hessian using finite differences\n"
    << L << "    Ze = Z * e\n"
    << L << "    u = step_scale / norm(Ze,inf)\n"
    << L << "    if there are bounds on the problem then\n"
    << L << "        Find the largest (in magnitude) positive (u_pos) and\n"
    << L << "        negative (u_neg) step u where the slightly relaxed variable\n"
    << L << "        bounds:\n"
    << L << "            xl - delta <= x_k + u * Ze <= xu + delta\n"
    << L << "        are strictly satisfied (where delta = max_var_bounds_viol).\n"
    << L << "        if u_pos < u then\n"
    << L << "            u = u_pos\n"
    << L << "        end\n"
    << L << "        if abs(u_neg) < u then\n"
    << L << "            u = u_neg\n"
    << L << "        end\n"
    << L << "    end\n"
    << L << "    x_fd = x_k + u * Ze\n"
    << L << "    rGf_fd = ( Z_k' * Gf(x_fd) - rGf_k ) / u\n"
    << L << "    if norm(rGf_fd,inf) <= min_diag then\n"
    << L << "        rHL_k = min_diag * eye(n-r)\n"
    << L << "    else\n"
    << L << "        if initialization_method == SCALE_IDENTITY then\n"
    << L << "            rHL_k = norm(rGf_fd,inf) * eye(n-r)\n"
    << L << "        else if initialization_method == SCALE_DIAGONAL or SCALE_DIAGONAL_ABS then\n"
    << L << "            *** Make sure that all of the diagonal elements are\n"
    << L << "            *** positive and that the smallest element is\n"
    << L << "            *** no smaller than norm(rGf_fd,inf) / max_cond\n"
    << L << "            *** So that rHL_k will be positive definite an\n"
    << L << "            *** well conditioned\n"
    << L << "            min_ele = max( norm(rGf_fd,inf)/max_cond, min_diag )\n"
    << L << "            if initialization_method == SCALE_DIAGONAL then\n"
    << L << "                for i = 1 ... n-r\n"
    << L << "                   diag(i) = max( rGf_fd(i), min_ele )\n"
    << L << "                end\n"
    << L << "            else *** SCALE_DIAGONAL_ABS\n"
    << L << "                for i = 1 ... n-r\n"
    << L << "                   diag(i) = max( abs(rGf_fd(i)), min_ele )\n"
    << L << "                end\n"
    << L << "            end\n"
    << L << "            rHL_k = diag(diag)\n"
    << L << "        end\n"
    << L << "    end\n"
    << L << "end\n";
}

} // end namespace MoochoPack
