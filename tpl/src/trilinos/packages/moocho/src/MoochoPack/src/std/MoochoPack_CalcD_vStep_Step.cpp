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
#include <iostream>

#include "MoochoPack_CalcD_vStep_Step.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
//#include "ConstrainedOptPack_print_vector_change_stats.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"


bool MoochoPack::CalcD_vStep_Step::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
  {
  using Teuchos::dyn_cast;
  using IterationPack::print_algorithm_step;
  using AbstractLinAlgPack::ele_wise_prod;
  using AbstractLinAlgPack::lowerbound_multipliers_step;
  using AbstractLinAlgPack::upperbound_multipliers_step;

  NLPAlgo &algo = rsqp_algo(_algo);
  IpState	 &s    = dyn_cast<IpState>(_algo.state());
  NLP      &nlp  = algo.nlp();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
    }

  // Get iteration quantities
  const value_type& mu = s.barrier_parameter().get_k(0);
  const Vector &d_k = s.d().get_k(0);
  const MatrixSymDiagStd& invXl = s.invXl().get_k(0);
  const MatrixSymDiagStd& invXu = s.invXu().get_k(0);
  const MatrixSymDiagStd& Vl = s.Vl().get_k(0);
  const MatrixSymDiagStd& Vu = s.Vu().get_k(0);

  VectorMutable& dvl_k = s.dvl().set_k(0);
  VectorMutable& dvu_k = s.dvu().set_k(0);

  lowerbound_multipliers_step(mu, invXl.diag(), Vl.diag(), d_k, &dvl_k);
  upperbound_multipliers_step(mu, invXu.diag(), Vu.diag(), d_k, &dvu_k);

  /*
  // dvl = mu*invXl*e - vl - invXl*Vl*d_k
  dvl_k = 0;
  ele_wise_prod(-1.0, invXl.diag(), Vl.diag(), &dvl_k);
  ele_wise_prod(1.0, dvl_k, d_k, &dvl_k);

  std::cout << "d_k =\n" << d_k;
   std::cout << "-invXl*Vl*d_k = \n" << dvl_k;
 
  dvl_k.axpy(-1.0, Vl.diag());
  
   std::cout << "-vl-invXl*Vl*d_k = \n" << dvl_k;

  dvl_k.axpy(mu, invXl.diag());

   std::cout << "dvl_k = \n" << dvl_k;

  // dvu = mu*invXu*e - vu + invXu*Vu*d_k
  dvu_k = 0;
  ele_wise_prod(1.0, invXu.diag(), Vu.diag(), &dvu_k);
  ele_wise_prod(1.0, dvu_k, d_k, &dvu_k);

  dvu_k.axpy(-1.0, Vu.diag());
  
  dvu_k.axpy(mu, invXu.diag());
  */
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
    {
    out	<< "\nx_k = \n" << s.x().get_k(0)
      << "\nxl = \n" << nlp.xl()
      << "\nxu = \n" << nlp.xu()
      << "\ndvl_k = \n" << dvl_k
      << "\ndvu_k = \n" << dvu_k;
    }

  return true;
  }

void MoochoPack::CalcD_vStep_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Calculates the search direction for the dual variables\n"
    << L << "dvl_k = mu*invXl_k*e - vl_k - invXl_k*Vl_k*d_k\n"
    << L << "dvu_k = mu*invXu_k*e - vu_k + invXu_k*Vu_k*d_k\n";
}
