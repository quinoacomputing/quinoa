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
#include <iostream>

#include "MoochoPack_TangentialStepIP_Step.hpp"
#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

namespace MoochoPack {

bool TangentialStepIP_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
  {
  using BLAS_Cpp::no_trans;
  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using LinAlgOpPack::Vt_S;
  using LinAlgOpPack::Vp_StV;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::V_InvMtV;
   using LinAlgOpPack::M_StM;
  using LinAlgOpPack::Mp_StM;
  using LinAlgOpPack::assign;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  IpState	    &s      = dyn_cast<IpState>(_algo.state());

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  // Compute qp_grad which is an approximation to rGf + Z'*(mu*(invXu*e-invXl*e) + no_cross_term
  // minimize round off error by calc'ing Z'*(Gf + mu*(invXu*e-invXl*e))

  // qp_grad_k = Z'*(Gf + mu*(invXu*e-invXl*e))
  const MatrixSymDiagStd  &invXu = s.invXu().get_k(0);
  const MatrixSymDiagStd  &invXl = s.invXl().get_k(0);
  const value_type            &mu    = s.barrier_parameter().get_k(0);
  const MatrixOp          &Z_k   = s.Z().get_k(0);

  Teuchos::RCP<VectorMutable> rhs = s.Gf().get_k(0).clone();
  Vp_StV( rhs.get(), mu,      invXu.diag() );
  Vp_StV( rhs.get(), -1.0*mu, invXl.diag() );
  
  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) 
    {
    out << "\n||Gf_k + mu_k*(invXu_k-invXl_k)||inf = " << rhs->norm_inf() << std::endl;
    }
  if( (int)olevel >= (int)PRINT_VECTORS)
    {
    out << "\nGf_k + mu_k*(invXu_k-invXl_k) =\n" << *rhs;
    }

  VectorMutable &qp_grad_k = s.qp_grad().set_k(0);
  V_MtV(&qp_grad_k, Z_k, BLAS_Cpp::trans, *rhs);
  
  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) 
    {
    out << "\n||Z_k'*(Gf_k + mu_k*(invXu_k-invXl_k))||inf = " << qp_grad_k.norm_inf() << std::endl;
    }
  if( (int)olevel >= (int)PRINT_VECTORS )
    {
    out << "\nZ_k'*(Gf_k + mu_k*(invXu_k-invXl_k)) =\n" << qp_grad_k;
    }

  // error check for cross term
  value_type         &zeta    = s.zeta().set_k(0);
  const Vector &w_sigma = s.w_sigma().get_k(0);
  
  // need code to calculate damping parameter
  zeta = 1.0;

  Vp_StV(&qp_grad_k, zeta, w_sigma);

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) 
    {
    out << "\n||qp_grad_k||inf = " << qp_grad_k.norm_inf() << std::endl;
    }
  if( (int)olevel >= (int)PRINT_VECTORS ) 
    {
    out << "\nqp_grad_k =\n" << qp_grad_k;
    }

  // build the "Hessian" term B = rHL + rHB
  // should this be MatrixSymOpNonsing
  const MatrixSymOp      &rHL_k = s.rHL().get_k(0);
  const MatrixSymOp      &rHB_k = s.rHB().get_k(0);
  MatrixSymOpNonsing &B_k   = dyn_cast<MatrixSymOpNonsing>(s.B().set_k(0));
  if (B_k.cols() != Z_k.cols())
    {
    // Initialize space in rHB
    dyn_cast<MatrixSymInitDiag>(B_k).init_identity(Z_k.space_rows(), 0.0);
    }

  //	M_StM(&B_k, 1.0, rHL_k, no_trans);
  assign(&B_k, rHL_k, BLAS_Cpp::no_trans);
  if( (int)olevel >= (int)PRINT_VECTORS ) 
    {
    out << "\nB_k = rHL_k =\n" << B_k;
    }
  Mp_StM(&B_k, 1.0, rHB_k, BLAS_Cpp::no_trans);
  if( (int)olevel >= (int)PRINT_VECTORS ) 
    {
    out << "\nB_k = rHL_k + rHB_k =\n" << B_k;
    }

  // Solve the system pz = - inv(rHL) * qp_grad
  VectorMutable   &pz_k  = s.pz().set_k(0);
  V_InvMtV( &pz_k, B_k, no_trans, qp_grad_k );
  Vt_S( &pz_k, -1.0 );

  // Zpz = Z * pz
  V_MtV( &s.Zpz().set_k(0), s.Z().get_k(0), no_trans, pz_k );

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
    {
    out	<< "\n||pz||inf   = " << s.pz().get_k(0).norm_inf()
      << "\nsum(Zpz)    = " << AbstractLinAlgPack::sum(s.Zpz().get_k(0))  << std::endl;
    }

  if( (int)olevel >= (int)PRINT_VECTORS )
    {
    out << "\npz_k = \n" << s.pz().get_k(0);
    out << "\nnu_k = \n" << s.nu().get_k(0);
    out << "\nZpz_k = \n" << s.Zpz().get_k(0);
    out << std::endl;
    }

  if(algo.algo_cntr().check_results())
    {
    assert_print_nan_inf(s.pz().get_k(0),  "pz_k",true,&out);
    assert_print_nan_inf(s.Zpz().get_k(0), "Zpz_k",true,&out);
    }

  return true;
  }

void TangentialStepIP_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
  {
  out
    << L << "*** Calculate the null space step by solving an unconstrainted QP\n"
    << L << "zeta_k = 1.0\n"
    << L << "qp_grad_k = Z_k'*(Gf_k + mu_k*(invXu_k-invXl_k)) + zeta_k*w_sigma_k\n"
    << L << "B_k = rHL_k + rHB_k\n"
    << L << "pz_k = -inv(B_k)*qp_grad_k\n"
    << L << "Zpz_k = Z_k*pz_k\n"
    ;
  }

} // end namespace MoochoPack
