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

#include "MoochoPack_ReducedHessianSecantUpdateBFGSFull_Strategy.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_NLPAlgoState.hpp"

namespace MoochoPack {

ReducedHessianSecantUpdateBFGSFull_Strategy::ReducedHessianSecantUpdateBFGSFull_Strategy(
  const bfgs_update_ptr_t&      bfgs_update
  )
  :bfgs_update_(bfgs_update)
{}

bool ReducedHessianSecantUpdateBFGSFull_Strategy::perform_update(
  VectorMutable           *s_bfgs
  ,VectorMutable          *y_bfgs
  ,bool                   first_update
  ,std::ostream           & out
  ,EJournalOutputLevel    olevel
  ,NLPAlgo                *algo
  ,NLPAlgoState           *s
  ,MatrixSymOp            *rHL_k
  )
{
  bfgs_update().perform_update(
    s_bfgs,y_bfgs,first_update,out,olevel,algo->algo_cntr().check_results()
    ,rHL_k, &quasi_newton_stats_(*s).set_k(0)
    );
  return true;
}

void ReducedHessianSecantUpdateBFGSFull_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Perform BFGS update on full matrix where: B = rHL_k\n";
  bfgs_update().print_step(out,L);
}

}  // end namespace MoochoPack
