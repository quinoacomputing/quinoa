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

#include <limits>
#include <ostream>

#include "MoochoPack_CalcDFromYPYZPZ_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
//#include "ConstrainedOptPack_print_vector_change_stats.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_implicit_cast.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

bool MoochoPack::CalcDFromYPYZPZ_Step::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{

  using Teuchos::implicit_cast;
  using AbstractLinAlgPack::dot;
  using LinAlgOpPack::V_VpV;

  NLPAlgo &algo = rsqp_algo(_algo);
  NLPAlgoState &s = algo.rsqp_state();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  EJournalOutputLevel ns_olevel = algo.algo_cntr().null_space_journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( implicit_cast<int>(olevel) >= implicit_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  // d = Ypy + Zpz
  VectorMutable &d_k = s.d().set_k(0);
  const Vector &Ypy_k = s.Ypy().get_k(0);
  const Vector &Zpz_k = s.Zpz().get_k(0);
  V_VpV( &d_k, Ypy_k, Zpz_k );

  Range1D
    var_dep = s.var_dep(),
    var_indep = s.var_indep();
  
  if( implicit_cast<int>(olevel) >= implicit_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    const value_type very_small = std::numeric_limits<value_type>::min();
    out
      << "\n(Ypy_k'*Zpz_k)/(||Ypy_k||2 * ||Zpz_k||2 + eps)\n"
      << "  = ("<<dot(Ypy_k,Zpz_k)<<")/("<<Ypy_k.norm_2()<<" * "<<Zpz_k.norm_2()<<" + "<<very_small<<")\n"
      << "  = " << dot(Ypy_k,Zpz_k) / ( Ypy_k.norm_2() * Zpz_k.norm_2() + very_small ) << "\n";
/*
  ConstrainedOptPack::print_vector_change_stats(
  s.x().get_k(0), "x", s.d().get_k(0), "d", out );
*/
    // ToDo: Replace the above with a reduction operator!
  }

  if( implicit_cast<int>(olevel) >= implicit_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\n||d_k||inf            = " << d_k.norm_inf();
    if( var_dep.size() )
      out << "\n||d(var_dep)_k||inf   = " << d_k.sub_view(var_dep)->norm_inf();
    if( var_indep.size() )
      out << "\n||d(var_indep)_k||inf = " << d_k.sub_view(var_indep)->norm_inf();
    out << std::endl;
  }
  if( implicit_cast<int>(olevel) >= implicit_cast<int>(PRINT_VECTORS) ) {
    out << "\nd_k = \n" << d_k;
    if( var_dep.size() )
      out << "\nd(var_dep)_k = \n" << *d_k.sub_view(var_dep);
  }
  if( implicit_cast<int>(ns_olevel) >= implicit_cast<int>(PRINT_VECTORS) ) {
    if( var_indep.size() )
      out << "\nd(var_indep)_k = \n" << *d_k.sub_view(var_indep);
  }
  
  return true;

}

void MoochoPack::CalcDFromYPYZPZ_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Calculates the search direction d from Ypy and Zpz\n"
    << L << "d_k = Ypy_k + Zpz_k \n";
}
