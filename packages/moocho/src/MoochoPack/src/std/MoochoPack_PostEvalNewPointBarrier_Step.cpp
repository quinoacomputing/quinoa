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
#include <typeinfo>
#include <iostream>

#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_PostEvalNewPointBarrier_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"

#include "OptionsFromStreamPack_StringToBool.hpp"

#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

bool PostEvalNewPointBarrier_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
  {
  using Teuchos::dyn_cast;
  using IterationPack::print_algorithm_step;
  using AbstractLinAlgPack::inv_of_difference;
  using AbstractLinAlgPack::correct_upper_bound_multipliers;
  using AbstractLinAlgPack::correct_lower_bound_multipliers;
  using LinAlgOpPack::Vp_StV;

  NLPAlgo            &algo   = dyn_cast<NLPAlgo>(_algo);
  IpState             &s      = dyn_cast<IpState>(_algo.state());
  NLP                 &nlp    = algo.nlp();
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();
  
  if(!nlp.is_initialized())
    nlp.initialize(algo.algo_cntr().check_results());

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
    }

  IterQuantityAccess<VectorMutable> &x_iq = s.x();
  IterQuantityAccess<MatrixSymDiagStd> &Vl_iq = s.Vl();
  IterQuantityAccess<MatrixSymDiagStd> &Vu_iq = s.Vu();

  ///***********************************************************
  // Calculate invXl = diag(1/(x-xl)) 
  //  and invXu = diag(1/(xu-x)) matrices
  ///***********************************************************

  // get references to x, invXl, and invXu
  VectorMutable& x = x_iq.get_k(0);
  MatrixSymDiagStd& invXu = s.invXu().set_k(0);
  MatrixSymDiagStd& invXl = s.invXl().set_k(0);
  
  //std::cout << "xu=\n";
  //nlp.xu().output(std::cout);

  inv_of_difference(1.0, nlp.xu(), x, &invXu.diag());
  inv_of_difference(1.0, x, nlp.xl(), &invXl.diag());

  //std::cout << "invXu=\v";
  //invXu.output(std::cout);

  //std::cout << "\ninvXl=\v";
  //invXl.output(std::cout);
  
  // Check for divide by zeros - 
    using AbstractLinAlgPack::assert_print_nan_inf;
    assert_print_nan_inf(invXu.diag(), "invXu", true, &out); 
    assert_print_nan_inf(invXl.diag(), "invXl", true, &out); 
  // These should never go negative either - could be a useful check

  // Initialize Vu and Vl if necessary
  if ( /*!Vu_iq.updated_k(0) */ Vu_iq.last_updated() == IterQuantity::NONE_UPDATED )
    {
    VectorMutable& vu = Vu_iq.set_k(0).diag();		
    vu = 0;
    Vp_StV(&vu, s.barrier_parameter().get_k(-1), invXu.diag());

    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
      {
      out << "\nInitialize Vu with barrier_parameter * invXu ...\n";
      }
    }

if ( /*!Vl_iq.updated_k(0) */ Vl_iq.last_updated() == IterQuantity::NONE_UPDATED  )
    {
    VectorMutable& vl = Vl_iq.set_k(0).diag();
    vl = 0;
    Vp_StV(&vl, s.barrier_parameter().get_k(-1), invXl.diag());

    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
      {
      out << "\nInitialize Vl with barrier_parameter * invXl ...\n";
      }
    }

  if (s.num_basis().updated_k(0))
    {
    // Basis changed, reorder Vl and Vu
    if (Vu_iq.updated_k(-1))
      { Vu_iq.set_k(0,-1); }
    if (Vl_iq.updated_k(-1))
      { Vl_iq.set_k(0,-1); }
      
    VectorMutable& vu = Vu_iq.set_k(0).diag();
    VectorMutable& vl = Vl_iq.set_k(0).diag();

    s.P_var_last().permute( BLAS_Cpp::trans, &vu ); // Permute back to original order
    s.P_var_last().permute( BLAS_Cpp::trans, &vl ); // Permute back to original order

    if( (int)olevel >= (int)PRINT_VECTORS ) 
      {
      out	<< "\nx resorted vl and vu to the original order\n" << x;
      }

    s.P_var_current().permute( BLAS_Cpp::no_trans, &vu ); // Permute to new (current) order
    s.P_var_current().permute( BLAS_Cpp::no_trans, &vl ); // Permute to new (current) order

    if( (int)olevel >= (int)PRINT_VECTORS ) 
      {
      out	<< "\nx resorted to new basis \n" << x;
      }
    }

  correct_upper_bound_multipliers(nlp.xu(), +NLP::infinite_bound(), &Vu_iq.get_k(0).diag());
  correct_lower_bound_multipliers(nlp.xl(), -NLP::infinite_bound(), &Vl_iq.get_k(0).diag());
  
  if( (int)olevel >= (int)PRINT_VECTORS ) 
    {
    out << "x=\n"  << s.x().get_k(0);
    out << "xl=\n" << nlp.xl();
    out << "vl=\n" << s.Vl().get_k(0).diag();
    out << "xu=\n" << nlp.xu();
    out << "vu=\n" << s.Vu().get_k(0).diag();
    }

  // Update general algorithm bound multipliers
  VectorMutable& v = s.nu().set_k(0);
  v = Vu_iq.get_k(0).diag();
  Vp_StV(&v,-1.0,Vl_iq.get_k(0).diag());

  // Print vector information
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
    {
    out	<< "invXu_k.diag()=\n" << invXu.diag();
    out	<< "invXl_k.diag()=\n" << invXl.diag();
    out	<< "Vu_k.diag()=\n"    << Vu_iq.get_k(0).diag();
    out	<< "Vl_k.diag()=\n"    << Vl_iq.get_k(0).diag();
    out << "nu_k=\n"           << s.nu().get_k(0);
    }

  return true;
  }


void PostEvalNewPointBarrier_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
  {
  //const NLPAlgo   &algo = rsqp_algo(_algo);
  //const NLPAlgoState  &s    = algo.rsqp_state();
  out << L << "# Evaluate information specific to primal / dual barrier algorithms (Post EvalNewPoint)\n"
    << L << "invXl_k = diag(i, 1/(x(i)-xl))"
    << L << "invXu_k = diag(i, 1/(xu-x(i)))\n"
    << L << "if (Vu_k not updated) then\n"
    << L << "   Vu_k = mu_k*invXu_k\n"
    << L << "end\n"
    << L << "if (Vl_k not updated) then\n"
    << L << "   Vl_k = mu_k*invXl_k\n"
    << L << "end\n"
    << L << "nu_k_k = Vu_k.diag() - Vl_k.diag()\n";
  }

} // end namespace MoochoPack
