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

#include <ostream>
#include <limits>
#include <sstream>

#include "MoochoPack_CheckConvergenceStd_Strategy.hpp"
#include "MoochoPack_NLPAlgoContainer.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

namespace MoochoPack {

CheckConvergenceStd_Strategy::CheckConvergenceStd_Strategy(
  EOptErrorCheck         opt_error_check
  ,EScaleKKTErrorBy      scale_opt_error_by
  ,EScaleKKTErrorBy      scale_feas_error_by
  ,EScaleKKTErrorBy      scale_comp_error_by
  ,bool                  scale_opt_error_by_Gf
  )
  :
  CheckConvergence_Strategy(
    opt_error_check,
    scale_opt_error_by,
    scale_feas_error_by,
    scale_comp_error_by,
    scale_opt_error_by_Gf
    )
  {}

bool CheckConvergenceStd_Strategy::Converged(
  Algorithm& _algo
  )
  {
  using AbstractLinAlgPack::assert_print_nan_inf;
  using AbstractLinAlgPack::combined_nu_comp_err;
  
  NLPAlgo        &algo = rsqp_algo(_algo);
  NLPAlgoState   &s    = algo.rsqp_state();
  NLP            &nlp  = algo.nlp();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  const size_type
    n  = nlp.n(),
    m  = nlp.m(),
    nb = nlp.num_bounded_x();

  // Get the iteration quantities
  IterQuantityAccess<value_type>
    &opt_kkt_err_iq  = s.opt_kkt_err(),
    &feas_kkt_err_iq = s.feas_kkt_err(),
      &comp_kkt_err_iq = s.comp_kkt_err();
  
  IterQuantityAccess<VectorMutable>
    &x_iq       = s.x(),
    &d_iq       = s.d(),
    &Gf_iq      = s.Gf(),
    *c_iq       = m     ? &s.c()      : NULL,
    *rGL_iq     = n > m ? &s.rGL()    : NULL,
    *GL_iq      = n > m ? &s.GL()     : NULL,
    *nu_iq      = n > m ? &s.nu()     : NULL;

  // opt_err = (||rGL||inf or ||GL||) / (||Gf|| + scale_kkt_factor)
  value_type 
    norm_inf_Gf_k = 0.0,
    norm_inf_GLrGL_k = 0.0;

  if( n > m && scale_opt_error_by_Gf() && Gf_iq.updated_k(0) ) {
    assert_print_nan_inf(
      norm_inf_Gf_k = Gf_iq.get_k(0).norm_inf(),
      "||Gf_k||inf",true,&out
      );
  }

  // NOTE:
  // The strategy object CheckConvergenceIP_Strategy assumes
  // that this will always be the gradient of the lagrangian
  // of the original problem, not the gradient of the lagrangian
  // for psi. (don't use augmented nlp info here)
  if( n > m ) {
    if( opt_error_check() == OPT_ERROR_REDUCED_GRADIENT_LAGR ) {
      assert_print_nan_inf( norm_inf_GLrGL_k = rGL_iq->get_k(0).norm_inf(),
                  "||rGL_k||inf",true,&out);
    }
    else {
      assert_print_nan_inf( norm_inf_GLrGL_k = GL_iq->get_k(0).norm_inf(),
                  "||GL_k||inf",true,&out);
    }
  }

  const value_type
    opt_scale_factor = 1.0 + norm_inf_Gf_k,
    opt_err = norm_inf_GLrGL_k / opt_scale_factor;
  
  // feas_err
  const value_type feas_err = ( ( m ? c_iq->get_k(0).norm_inf() : 0.0 ) );

  // comp_err
  value_type comp_err = 0.0;
  if ( n > m ) {
    if (nb > 0) {
      comp_err = combined_nu_comp_err(nu_iq->get_k(0), x_iq.get_k(0), nlp.xl(), nlp.xu());
    }
    if(m) {
      assert_print_nan_inf( feas_err,"||c_k||inf",true,&out);
    }
  }

  // scaling factors
  const value_type 
    scale_opt_factor = CalculateScalingFactor(s, scale_opt_error_by()),
    scale_feas_factor = CalculateScalingFactor(s, scale_feas_error_by()),
    scale_comp_factor = CalculateScalingFactor(s, scale_comp_error_by());
  
  // kkt_err
  const value_type
    opt_kkt_err_k  = opt_err/scale_opt_factor,
     feas_kkt_err_k = feas_err/scale_feas_factor,
    comp_kkt_err_k = comp_err/scale_comp_factor;

  // update the iteration quantities
  if(n > m) opt_kkt_err_iq.set_k(0) = opt_kkt_err_k;
  feas_kkt_err_iq.set_k(0) = feas_kkt_err_k;
  comp_kkt_err_iq.set_k(0) = comp_kkt_err_k;

  // step_err
  value_type step_err = 0.0;
  if( d_iq.updated_k(0) ) {
    step_err = AbstractLinAlgPack::max_rel_step(x_iq.get_k(0),d_iq.get_k(0));
    assert_print_nan_inf( step_err,"max(d(i)/max(1,x(i)),i=1...n)",true,&out);
  }

  const value_type
    opt_tol		= algo.algo_cntr().opt_tol(),
    feas_tol	= algo.algo_cntr().feas_tol(),
    comp_tol	= algo.algo_cntr().comp_tol(),
    step_tol	= algo.algo_cntr().step_tol();
  
  const bool found_solution = 
    opt_kkt_err_k < opt_tol 
    && feas_kkt_err_k < feas_tol 
    && comp_kkt_err_k < comp_tol 
    && step_err < step_tol;
  
  if( int(olevel) >= int(PRINT_ALGORITHM_STEPS) || (int(olevel) >= int(PRINT_BASIC_ALGORITHM_INFO) && found_solution) )
  {
    out	
      << "\nscale_opt_factor = " << scale_opt_factor
      << " (scale_opt_error_by = " << (scale_opt_error_by()==SCALE_BY_ONE ? "SCALE_BY_ONE"
                       : (scale_opt_error_by()==SCALE_BY_NORM_INF_X ? "SCALE_BY_NORM_INF_X"
                        : "SCALE_BY_NORM_2_X" ) ) << ")"

      << "\nscale_feas_factor = " << scale_feas_factor
      << " (scale_feas_error_by = " << (scale_feas_error_by()==SCALE_BY_ONE ? "SCALE_BY_ONE"
                       : (scale_feas_error_by()==SCALE_BY_NORM_INF_X ? "SCALE_BY_NORM_INF_X"
                        : "SCALE_BY_NORM_2_X" ) ) << ")"

      << "\nscale_comp_factor = " << scale_comp_factor
      << " (scale_comp_error_by = " << (scale_comp_error_by()==SCALE_BY_ONE ? "SCALE_BY_ONE"
                       : (scale_comp_error_by()==SCALE_BY_NORM_INF_X ? "SCALE_BY_NORM_INF_X"
                        : "SCALE_BY_NORM_2_X" ) ) << ")"
      << "\nopt_scale_factor = " << opt_scale_factor
      << " (scale_opt_error_by_Gf = " << (scale_opt_error_by_Gf()?"true":"false") << ")"
      << "\nopt_kkt_err_k    = " << opt_kkt_err_k << ( opt_kkt_err_k < opt_tol ? " < " : " > " )
      << "opt_tol  = " << opt_tol
      << "\nfeas_kkt_err_k   = " << feas_kkt_err_k << ( feas_kkt_err_k < feas_tol ? " < " : " > " )
      << "feas_tol = " << feas_tol
      << "\ncomp_kkt_err_k   = " << comp_kkt_err_k << ( comp_kkt_err_k < comp_tol ? " < " : " > " )
      << "comp_tol = " << comp_tol
      << "\nstep_err         = " << step_err << ( step_err < step_tol ? " < " : " > " )
      << "step_tol = " << step_tol
      << std::endl;
    }
  
  return found_solution;

  }

void CheckConvergenceStd_Strategy::print_step( const Algorithm& _algo, std::ostream& out, const std::string& L ) const
  {
  out
    << L << "*** Check to see if the KKT error is small enough for convergence\n"
    << L << "if scale_(opt|feas|comp)_error_by == SCALE_BY_ONE then\n"
    << L << "    scale_(opt|feas|comp)_factor = 1.0\n"
    << L << "else if scale_(opt|feas|comp)_error_by == SCALE_BY_NORM_2_X then\n"
    << L << "    scale_(opt|feas|comp)_factor = 1.0 + norm_2(x_k)\n"
    << L << "else if scale_(opt|feas|comp)_error_by == SCALE_BY_NORM_INF_X then\n"
    << L << "    scale_(opt|feas|comp)_factor = 1.0 + norm_inf(x_k)\n"
    << L << "end\n"
    << L << "if scale_opt_error_by_Gf == true then\n"
    << L << "    opt_scale_factor = 1.0 + norm_inf(Gf_k)\n"
    << L << "else\n"
    << L << "    opt_scale_factor = 1.0\n"
    << L << "end\n";
  if( opt_error_check() == OPT_ERROR_REDUCED_GRADIENT_LAGR )
    {
    out
      << L << "opt_err = norm_inf(rGL_k)/opt_scale_factor\n";
    }
  else
    {
    out
      << L << "opt_err = norm_inf(GL_k)/opt_scale_factor\n";
    }

  out
    << L << "feas_err = norm_inf(c_k)\n"
    << L << "comp_err = max(i, nu(i)*(xu(i)-x(i)), -nu(i)*(x(i)-xl(i)))\n"
    << L << "opt_kkt_err_k = opt_err/scale_opt_factor\n"
    << L << "feas_kkt_err_k = feas_err/scale_feas_factor\n"
    << L << "comp_kkt_err_k = feas_err/scale_comp_factor\n"
    << L << "if d_k is updated then\n"
    << L << "    step_err = max( |d_k(i)|/(1+|x_k(i)|), i=1..n )\n"
    << L << "else\n"
    << L << "    step_err = 0\n"
    << L << "end\n"
    << L << "if opt_kkt_err_k < opt_tol\n"
    << L << "       and feas_kkt_err_k < feas_tol\n"
    << L << "       and step_err < step_tol then\n"
    << L << "   report optimal x_k, lambda_k and nu_k to the nlp\n"
    << L << "   terminate, the solution has beed found!\n"
    << L << "end\n";
  }


value_type CheckConvergenceStd_Strategy::CalculateScalingFactor( NLPAlgoState& state, EScaleKKTErrorBy scale_by ) const
  {
  // scale_kkt_factor
  value_type scale_factor = 1.0;
  switch(scale_by) 
    {
    case SCALE_BY_ONE:
      scale_factor = 1.0;
      break;
    case SCALE_BY_NORM_2_X:
      scale_factor = 1.0 + state.x().get_k(0).norm_2();
      break;
    case SCALE_BY_NORM_INF_X:
      scale_factor = 1.0 + state.x().get_k(0).norm_inf();
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Should never be called
    }

  return scale_factor;
  }

}	// end namespace MoochoPack

