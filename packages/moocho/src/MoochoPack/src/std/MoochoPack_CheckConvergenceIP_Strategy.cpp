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
//#include <limits>
//#include <sstream>

#include "MoochoPack_CheckConvergenceIP_Strategy.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

CheckConvergenceIP_Strategy::CheckConvergenceIP_Strategy(
  EOptErrorCheck         opt_error_check
  ,EScaleKKTErrorBy      scale_opt_error_by
  ,EScaleKKTErrorBy      scale_feas_error_by
  ,EScaleKKTErrorBy      scale_comp_error_by
  ,bool                  scale_opt_error_by_Gf
  )
  :
  CheckConvergenceStd_Strategy(
    opt_error_check,
    scale_opt_error_by,
    scale_feas_error_by,
    scale_comp_error_by,
    scale_opt_error_by_Gf
    )
  {}

bool CheckConvergenceIP_Strategy::Converged(
  Algorithm& _algo
  )
  {
  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::num_bounded;
  using AbstractLinAlgPack::IP_comp_err_with_mu;

  // Calculate kkt errors and check for overall convergence
  //bool found_solution = CheckConvergenceStd_Strategy::Converged(_algo);
  bool found_solution = false;

  // Recalculate the complementarity error including mu
  
  // Get the iteration quantities
  IpState &s = dyn_cast<IpState>(*_algo.get_state());
  NLPAlgo& algo = rsqp_algo(_algo);
  NLP& nlp = algo.nlp();
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // Get necessary iteration quantities
  const value_type &mu_km1 = s.barrier_parameter().get_k(-1);
  const Vector& x_k = s.x().get_k(0);
  const VectorMutable& Gf_k = s.Gf().get_k(0);
  const Vector& rGL_k = s.rGL().get_k(0);
  const Vector& c_k = s.c().get_k(0);
  const Vector& vl_k = s.Vl().get_k(0).diag();
  const Vector& vu_k = s.Vu().get_k(0).diag();
  
  // Calculate the errors with Andreas' scaling
  value_type& opt_err = s.opt_kkt_err().set_k(0);
  value_type& feas_err = s.feas_kkt_err().set_k(0);
  value_type& comp_err = s.comp_kkt_err().set_k(0);
  value_type& comp_err_mu = s.comp_err_mu().set_k(0);

  // scaling
  value_type scale_1 = 1 + x_k.norm_1()/x_k.dim();

  Teuchos::RCP<VectorMutable> temp = Gf_k.clone();
  temp->axpy(-1.0, vl_k);
  temp->axpy(1.0, vu_k);
  value_type scale_2 = temp->norm_1();
  scale_2 += vl_k.norm_1() + vu_k.norm_1();

  *temp = nlp.infinite_bound();
  const size_type nlb = num_bounded(nlp.xl(), *temp, nlp.infinite_bound());
  *temp = -nlp.infinite_bound();
  const size_type nub = num_bounded(*temp, nlp.xu(), nlp.infinite_bound());
  scale_2 = 1 + scale_2/(1+nlp.m()+nlb+nub);

  // Calculate the opt_err
  opt_err = rGL_k.norm_inf() / scale_2;

  // Calculate the feas_err
  feas_err = c_k.norm_inf() / scale_1;
  
  // Calculate the comp_err
  if( (int)olevel >= (int)PRINT_VECTORS )
    {
    out << "\nx =\n"    << s.x().get_k(0);
    out << "\nxl =\n"   << nlp.xl();
    out << "\nvl =\n"   << s.Vl().get_k(0).diag();
    out << "\nxu =\n"   << nlp.xu();
    out << "\nvu =\n"   << s.Vu().get_k(0).diag();
    }

  comp_err = IP_comp_err_with_mu(
    0.0, nlp.infinite_bound(), s.x().get_k(0), nlp.xl(), nlp.xu()
    ,s.Vl().get_k(0).diag(), s.Vu().get_k(0).diag());

  comp_err_mu = IP_comp_err_with_mu(
    mu_km1, nlp.infinite_bound(), s.x().get_k(0), nlp.xl(), nlp.xu()
    ,s.Vl().get_k(0).diag(), s.Vu().get_k(0).diag());

  comp_err = comp_err / scale_2;
  comp_err_mu = comp_err_mu / scale_2;

  // check for convergence
  
  const value_type opt_tol = algo.algo_cntr().opt_tol();
  const value_type feas_tol = algo.algo_cntr().feas_tol();
  const value_type comp_tol = algo.algo_cntr().comp_tol();

  if (opt_err < opt_tol && feas_err < feas_tol && comp_err < comp_tol)
    {
    found_solution = true;
    }

  if( int(olevel) >= int(PRINT_ALGORITHM_STEPS) || (int(olevel) >= int(PRINT_BASIC_ALGORITHM_INFO) && found_solution) )
    {
    out	
      << "\nopt_kkt_err_k   = " << opt_err << ( opt_err < opt_tol ? " < " : " > " )
      << "opt_tol = " << opt_tol
      << "\nfeas_kkt_err_k   = " << feas_err << ( feas_err < feas_tol ? " < " : " > " )
      << "feas_tol = " << feas_tol
      << "\ncomp_kkt_err_k   = " << comp_err << ( comp_err < comp_tol ? " < " : " > " )
      << "comp_tol = " << comp_tol
      << "\ncomp_err_mu      = " << comp_err_mu
      << "\nbarrier_parameter_k (mu_km1) = " << mu_km1
      << "comp_tol = " << comp_tol
      << std::endl;
    }
    
  return found_solution;
  }

void CheckConvergenceIP_Strategy::print_step( const Algorithm& _algo, std::ostream& out, const std::string& L ) const
  {
  out 
    << L << "CheckConvergenceIP_Strategy\n"
    << L << " IP_found_solution = CheckConvergedStd_Strategy::Converged(_algo, reportFinalSolution)\n";
  
  CheckConvergenceStd_Strategy::print_step(_algo, out, L+"   ");
  
  out 
    << L << "*** recalculate comp_err\n"
    << L << "comp_err_k = 0.0"
    << L << "for all i = 1 to n\n"
    << L << "   comp_err_k = max( comp_err_k, vl_k(i)*(x_k(i)-xl_k(i))-mu_km1, vu_k(i)*(xu_k(i)-x_k(i))-mu_k )\n"
    << L << "next i\n"
    << L << "if IP_found_solution then\n"
    << L << "   IP_found_solution = false\n"
    << L << "   if (comp_err_k < comp_tol && mu_k < comp_tol) then\n"
    << L << "      IP_found_solution = true\n"
    << L << "   end\n"
    << L << "end\n"
    << L << "return IP_found_solution\n";
  }

}	// end namespace MoochoPack

