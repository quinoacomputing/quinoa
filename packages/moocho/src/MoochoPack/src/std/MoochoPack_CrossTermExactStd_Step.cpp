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

#include "MoochoPack_CrossTermExactStd_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}

bool MoochoPack::CrossTermExactStd_Step::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{
  using LinAlgOpPack::V_MtV;
  using DenseLinAlgPack::norm_inf;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
  }

  // tmp = HL * Ypy
  DVector tmp;
  V_MtV( &tmp, s.HL().get_k(0), BLAS_Cpp::no_trans, s.Ypy().get_k(0)() );
  // w = Z' * tmp
  V_MtV( &s.w().set_k(0).v(), s.Z().get_k(0), BLAS_Cpp::trans, tmp() );

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\n||w||inf = "	<< s.w().get_k(0).norm_inf() << std::endl;
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nw_k =\n" << s.w().get_k(0)();
  }

  return true;
}

void MoochoPack::CrossTermExactStd_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Evaluate the exact reduced QP cross term\n"
    << L << "w_k = Z_k' * HL_k * Ypy_k\n";
}

#endif // 0
