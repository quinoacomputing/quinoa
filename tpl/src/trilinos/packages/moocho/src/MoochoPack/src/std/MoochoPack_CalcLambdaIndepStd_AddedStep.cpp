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

#include <ostream>

#include "MoochoPack_CalcLambdaIndepStd_AddedStep.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack_ComputeMinMult.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::Vp_StMtV;
}

bool MoochoPack::CalcLambdaIndepStd_AddedStep::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{

  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;

  using DenseLinAlgPack::Vt_S;
  using DenseLinAlgPack::norm_inf;

  using AbstractLinAlgPack::Vp_StMtV;
  
  using ConstrainedOptPack::min_abs;

  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::Vp_MtV;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  Range1D		decomp	= s.equ_decomp();
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  // Compute: lambda(decomp) = inv(Gc(decomp)'* Y)' * ( - Y' * (Gf + nu) - U' * lambda(undecomp) )
  // where U = Gc(undecomp)' * Y

  // Must resize lambda here explicitly since we will only be updating a region of it.
  // If lambda(undecomp) has already been updated then lambda will have been resized
  // already but lambda(decomp) will not be initialized yet.
  if( !s.lambda().updated_k(0) ) s.lambda().set_k(0).v().resize( algo.nlp().m() );
  
  DVectorSlice lambda_decomp = s.lambda().get_k(0).v()(decomp);
  
  // lambda_decomp_tmp1 = - Y' * (Gf + nu)
  if( algo.nlp().has_bounds() ) {
    // _tmp = Gf + nu
    DVector _tmp = s.Gf().get_k(0)();
    DVectorSlice _vs_tmp = _tmp;	// only create this DVectorSlice once
    Vp_V( &_vs_tmp, s.nu().get_k(0)() );
    // lambda_decomp_tmp1 = - Y' * _tmp
    V_StMtV( &lambda_decomp, -1.0, s.Y().get_k(0), trans, _vs_tmp );
  }
  else {
    // lambda_decomp__tmp1 = - Y' * Gf
    V_StMtV( &lambda_decomp, -1.0, s.Y().get_k(0), trans, s.Gf().get_k(0)() );
  }

  // lambda_decomp_tmp2 = lambda_decomp_tmp1 - U' * lambda(undecomp)
  if( algo.nlp().r() < algo.nlp().m() ) {
    Range1D undecomp = s.equ_undecomp();
    Vp_StMtV( &lambda_decomp, -1.0, s.U().get_k(0), trans, s.lambda().get_k(0).v()(undecomp) );
  }
  // else lambda(decomp)_tmp2 = lambda(decomp)_tmp1

  // lambda_decomp = inv(Gc(decomp)'* Y)' * lambda_decomp_tmp2
  s.decomp_sys().solve_transAtY( lambda_decomp, trans, &lambda_decomp );

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\nmax(|lambda_k(equ_decomp)(i)|) = " << norm_inf(lambda_decomp)
      << "\nmin(|lambda_k(equ_decomp)(i)|) = " << min_abs(lambda_decomp)  << std::endl;
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out	<< "\nlambda_k(equ_decomp) = \n" << lambda_decomp;
  }

  return true;
}

void MoochoPack::CalcLambdaIndepStd_AddedStep::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Calculate the Lagrange multipliers for the decomposed constraints\n"
    << L << "lambda_k(equ_decomp) = - inv(Gc_k(:,equ_decomp)'*Y_k)\n"
    << L << "                        * (Y_k'*(Gf_k + nu_k) + U_k'*lambda_k(equ_undecomp))\n";
}

#endif // 0
