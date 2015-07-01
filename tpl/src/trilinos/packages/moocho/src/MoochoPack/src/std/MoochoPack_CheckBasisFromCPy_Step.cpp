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

#include <typeinfo>

#include "MoochoPack/src/std/CheckBasisFromCPy_Step.h"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack_DecompositionSystemVarReduct.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Midynamic_cast_verbose.h"
#include "MiWorkspacePack.h"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

namespace MoochoPack {

CheckBasisFromCPy_Step::CheckBasisFromCPy_Step(
  const new_basis_strategy_ptr_t& new_basis_strategy
  ,value_type max_basis_cond_change_frac
  )
  : new_basis_strategy_(new_basis_strategy)
  , max_basis_cond_change_frac_( max_basis_cond_change_frac )
{
  reset();
}

void CheckBasisFromCPy_Step::reset() {
  beta_min_ = std::numeric_limits<value_type>::max();
}

// Overridden

bool CheckBasisFromCPy_Step::do_step( Algorithm& _algo, poss_type step_poss
  , IterationPack::EDoStepType type, poss_type assoc_step_poss )
{
  using DenseLinAlgPack::norm_inf;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_V;
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  Range1D		decomp	= s.con_decomp();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  bool select_new_basis = false;

  // Compute: resid = (Gc(decomp)'*Y) * py + c(decomp)
  DVectorSlice                  py_k = s.py().get_k(0)();
  DVectorSlice                  c_k  = s.c().get_k(0)();
  DVectorSlice                  c_decomp_k  = c_k(decomp);
  Workspace<value_type>   resid_ws(wss,py_k.size());
  DVectorSlice                  resid(&resid_ws[0],resid_ws.size());
  {
#ifdef _WINDOWS
    DecompositionSystemVarReduct
      &decomp_sys = dynamic_cast<DecompositionSystemVarReduct&>(s.decomp_sys());
#else
    DecompositionSystemVarReduct
      &decomp_sys = dyn_cast<DecompositionSystemVarReduct>(s.decomp_sys());
#endif
    V_MtV( &resid, decomp_sys.C(), BLAS_Cpp::no_trans, py_k );
    Vp_V( &resid, c_decomp_k );
  }

  const value_type
    small_num    = std::numeric_limits<value_type>::min(),
    nrm_resid    = norm_inf(resid),
    nrm_c_decomp = norm_inf(c_decomp_k),
    beta         = nrm_resid / (nrm_c_decomp+small_num);

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    out	<< "\nbeta = ||(Gc(decomp)'*Y)*py_k + c_k(decomp)||inf / (||c_k(decomp)||inf + small_number)"
      << "\n     = "<<nrm_resid<<" / ("<<nrm_c_decomp<<" + "<<small_num<<")"
      << "\n     = " << beta << std::endl;
  }
  
  if( beta != 0.0 ) {
    if( beta < beta_min_ ) {
      beta_min_ = beta;
    }
    else {
      if( beta / beta_min_ > max_basis_cond_change_frac() ) {
        if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO ) {
          out
            << "\nbeta_change = ( ||R*py+c||/||c|| = " << beta
            << " ) / ( min ||R*py+c||/||c|| = " << beta_min_ << " )\n"
            << "              = " << (beta/beta_min_) << " > max_basis_cond_change_frac = "
            << max_basis_cond_change_frac()
            << "\n\nSelecting a new basis"
            << " (k = " << algo.state().k() << ") ...\n";
        }
        select_new_basis = true;
      }
    }
    if(select_new_basis) {
      // reset memory
      // ToDo: You really need to check to see if someone else
      //		 changed the basis also.
      beta_min_ = std::numeric_limits<value_type>::max();
      return new_basis_strategy().transistion_basis(algo,step_poss,type,assoc_step_poss);
    }
  }
  return true;
}

void CheckBasisFromCPy_Step::print_step( const Algorithm& algo, poss_type step_poss
  , IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Try to detect when the basis is becomming illconditioned\n"
    << L << "default: beta_min = inf\n"
    << L << "         max_basis_cond_change_frac = " << max_basis_cond_change_frac() << std::endl
    << L << "beta = norm_inf((Gc(decomp)'*Y)*py_k + c_k(decomp)) / (norm_inf(c_k(decomp))+small_number)\n"
    << L << "select_new_basis = false\n"
    << L << "if beta < beta_min then\n"
    << L << "    beta_min = beta\n"
    << L << "else\n"
    << L << "    if beta/ beta_min > max_basis_cond_change_frac then\n"
    << L << "        select_new_basis = true\n"
    << L << "    end\n"
    << L << "end\n"
    << L << "if select_new_basis == true then\n"
    << L << "    new basis selection : " << typeName(new_basis_strategy()) << std::endl;
  new_basis_strategy().print_method( rsqp_algo(algo),step_poss,type,assoc_step_poss,out
    ,L+"        " );
  out
    << L << "    end new basis selection\n"
    << L << "end\n";
}

}	// end namespace MoochoPack 

#endif // 0
