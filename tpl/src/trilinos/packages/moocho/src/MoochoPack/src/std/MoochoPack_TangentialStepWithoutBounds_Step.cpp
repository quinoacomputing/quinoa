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

#include "MoochoPack_TangentialStepWithoutBounds_Step.hpp"
#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"


namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
}


namespace MoochoPack {


TangentialStepWithoutBounds_Step::TangentialStepWithoutBounds_Step()
  :max_pz_norm_(-1.0),
   num_pz_damp_iters_(0)
{}


bool TangentialStepWithoutBounds_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{

  using std::endl;
  using BLAS_Cpp::no_trans;
  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using AbstractLinAlgPack::Vt_S;
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::V_InvMtV;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_MtV;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState &s = algo.rsqp_state();
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  EJournalOutputLevel ns_olevel = algo.algo_cntr().null_space_journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  if( algo.nlp().num_bounded_x() )
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error
      ,"TangentialStepWithoutBounds_Step::do_step(...): Error, "
      "can't solve for pz for NLP with undecomposed constraints or "
      "has bounds on the variables");

  // Comupte qp_grad which is an approximation to rGf + Z' * HL * Y * py

  // qp_grad_k = rGf_k
  VectorMutable &qp_grad_k = s.qp_grad().set_k(0) = s.rGf().get_k(0);

  IterQuantityAccess<value_type> &zeta_iq = s.zeta();
  IterQuantityAccess<VectorMutable> &w_iq = s.w();
  if( w_iq.updated_k(0) && zeta_iq.updated_k(0) ) {
    // qp_grad += zeta * w
    Vp_StV( &qp_grad_k, zeta_iq.get_k(0), w_iq.get_k(0) );
  }

  // Solve the system pz = - inv(rHL) * qp_grad
  VectorMutable &pz_k = s.pz().set_k(0);
  const MatrixSymOpNonsing &rHL_k = dyn_cast<MatrixSymOpNonsing>(s.rHL().get_k(0));
  V_InvMtV( &pz_k, rHL_k, no_trans, qp_grad_k );
  Vt_S( &pz_k, -1.0 );

  // nu = 0.0
  s.nu().set_k(0) = 0.0;

  value_type pz_norm_inf = -1.0;
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    pz_norm_inf = pz_k.norm_inf();
    out	<< "\n||pz_k||inf   = " << pz_norm_inf << endl;
  }
  if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_VECTORS) )
    out << "\npz_k = \n" << pz_k << std::endl;

  // Check to see if we need to dampen pz
  const bool dampen_pz = max_pz_norm() >= 0.0 && s.k() <= num_pz_damp_iters();
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out
      << "\nChecking if we need to dampen pz:"
      << "  ( (max_pz_norm="<<max_pz_norm()<<") >= 0.0 )"
      << " && ( (k="<<s.k()<<") <= (num_pz_damp_iters="<<num_pz_damp_iters()<<") ) : "
      << ( dampen_pz ? "true, dampen pz ..." : "false, no dampen pz!" )
      << endl;
  }

  if ( dampen_pz ) {
    // pz_new = ( max_pz_norm / pz_norm_inf ) * pz
    if (pz_norm_inf < 0.0 )
      pz_norm_inf = pz_k.norm_inf();
    const value_type pz_damp_factor = ( max_pz_norm() / pz_norm_inf );
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out
        << "\npz_damp_factor = max_pz_norm / ||pz||inf = "
        << max_pz_norm() << " / " << pz_norm_inf << " = "
        << pz_damp_factor;
    }
    Vt_S( &pz_k, pz_damp_factor );
  }

  // Zpz = Z * pz
  V_MtV( &s.Zpz().set_k(0), s.Z().get_k(0), no_trans, pz_k );

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\n||pz_k||inf   = " << s.pz().get_k(0).norm_inf()
        << "\n||Zpz_k||2    = " << s.Zpz().get_k(0).norm_2()  << std::endl;
  }

  if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\npz_k = \n" << s.pz().get_k(0);
    out << std::endl;
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nnu_k = \n" << s.nu().get_k(0);
    out << "\nZpz_k = \n" << s.Zpz().get_k(0);
    out << std::endl;
  }

  if(algo.algo_cntr().check_results()) {
    assert_print_nan_inf(s.pz().get_k(0),  "pz_k",true,&out);
    assert_print_nan_inf(s.Zpz().get_k(0), "Zpz_k",true,&out);
  }

  return true;
}


void TangentialStepWithoutBounds_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Calculate the null space step by solving an unconstrainted QP\n"
    << L << "qp_grad_k = rGf_k + zeta_k * w_k\n"
    << L << "solve:\n"
    << L << "    min     qp_grad_k' * pz_k + 1/2 * pz_k' * rHL_k * pz_k\n"
    << L << "    pz_k <: R^(n-r)\n"
    << L << "Zpz_k = Z_k * pz_k\n"
    << L << "nu_k = 0\n";
}


} // end namespace MoochoPack
