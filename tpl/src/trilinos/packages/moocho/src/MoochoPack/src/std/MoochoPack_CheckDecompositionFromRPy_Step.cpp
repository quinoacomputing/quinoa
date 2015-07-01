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

#include <typeinfo>

#include "MoochoPack_CheckDecompositionFromRPy_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace MoochoPack {

CheckDecompositionFromRPy_Step::CheckDecompositionFromRPy_Step(
  const new_decomp_strategy_ptr_t   &new_decomp_strategy
  ,value_type                       max_decomposition_cond_change_frac
  )
  :new_decomp_strategy_(new_decomp_strategy)
  ,max_decomposition_cond_change_frac_(max_decomposition_cond_change_frac)
{
  reset();
}

void CheckDecompositionFromRPy_Step::reset() {
  beta_min_ = std::numeric_limits<value_type>::max();
}

// Overridden

bool CheckDecompositionFromRPy_Step::do_step( Algorithm& _algo, poss_type step_poss
  , IterationPack::EDoStepType type, poss_type assoc_step_poss )
{
  NLPAlgo                &algo       = rsqp_algo(_algo);
  NLPAlgoState               &s          = algo.rsqp_state();
  const Range1D           equ_decomp  = s.equ_decomp();
  EJournalOutputLevel     olevel      = algo.algo_cntr().journal_output_level();
  std::ostream            &out        = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  bool select_new_decomposition = false;

  // Compute: resid = (Gc(decomp)'*Y) * py + c(decomp)
  const Vector                  &py_k       = s.py().get_k(0);
  const Vector                  &c_k        = s.c().get_k(0);
  Vector::vec_ptr_t             c_decomp_k  = c_k.sub_view(equ_decomp);
  VectorMutable::vec_mut_ptr_t  resid       = c_decomp_k->space().create_member();

  // resid = R*py + c(equ_decomp)
  LinAlgOpPack::V_MtV( resid.get(), s.R().get_k(0), BLAS_Cpp::no_trans, py_k );
  LinAlgOpPack::Vp_V( resid.get(), *c_decomp_k );

  const value_type
    small_num    = std::numeric_limits<value_type>::min(),
    epsilon      = std::numeric_limits<value_type>::epsilon(),
    nrm_resid    = resid->norm_inf(),
    nrm_c_decomp = c_decomp_k->norm_inf(),
    beta         = nrm_resid / (nrm_c_decomp+small_num);

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    out	<< "\nbeta = ||R*py_k + c_k(decomp)||inf / (||c_k(decomp)||inf + small_number)"
      << "\n     = "<<nrm_resid<<" / ("<<nrm_c_decomp<<" + "<<small_num<<")"
      << "\n     = " << beta << std::endl;
  }

  // Check to see if a new basis was selected or not
  IterQuantityAccess<index_type>
    &num_basis_iq = s.num_basis();
  if( num_basis_iq.updated_k(0) ) {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
      out	<< "\nnum_basis_k was updated so the basis changed so we will skip this check\n"
        << "    reset min ||R*py+c||/||c|| to current value + epsilon(" << epsilon << ")\n";
    beta_min_ = beta + epsilon;
    return true;
  }
  
  if( beta != 0.0 ) {
    if( beta < beta_min_ ) {
      beta_min_ = beta;
    }
    else {
      if( beta / beta_min_ > max_decomposition_cond_change_frac() ) {
        if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO ) {
          out
            << "\nbeta_change = ( ||R*py+c||/||c|| = " << beta
            << " ) / ( min ||R*py+c||/||c|| = " << beta_min_ << " )\n"
            << "              = " << (beta/beta_min_) << " > max_decomposition_cond_change_frac = "
            << max_decomposition_cond_change_frac()
            << "\n\nSelecting a new decomposition"
            << " (k = " << algo.state().k() << ") ...\n";
        }
        select_new_decomposition = true;
      }
    }
    if(select_new_decomposition) {
      reset();
      return new_decomp_strategy().new_decomposition(algo,step_poss,type,assoc_step_poss);
    }
  }
  return true;
}

void CheckDecompositionFromRPy_Step::print_step( const Algorithm& algo, poss_type step_poss
  , IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Try to detect when the decomposition is becomming illconditioned\n"
    << L << "default: beta_min = inf\n"
    << L << "         max_decomposition_cond_change_frac = " << max_decomposition_cond_change_frac() << std::endl
    << L << "beta = norm_inf(R*py_k + c_k(decomp)) / (norm_inf(c_k(decomp))+small_number)\n"
    << L << "select_new_decomposition = false\n"
    << L << "if num_basis_k is updated then\n"
    << L << "  beta_min = beta\n"
    << L << "end\n"
    << L << "if beta < beta_min then\n"
    << L << "  beta_min = beta\n"
    << L << "else\n"
    << L << "  if beta/ beta_min > max_decomposition_cond_change_frac then\n"
    << L << "        select_new_decomposition = true\n"
    << L << "    end\n"
    << L << "end\n"
    << L << "if select_new_decomposition == true then\n"
    << L << "    new decomposition selection : " << typeName(new_decomp_strategy()) << std::endl
    ;
  new_decomp_strategy().print_new_decomposition(
    rsqp_algo(algo),step_poss,type,assoc_step_poss,out, L + "    " );
  out
    << L << "    end new decomposition selection\n"
    << L << "end\n"
    ;
}

}	// end namespace MoochoPack 
