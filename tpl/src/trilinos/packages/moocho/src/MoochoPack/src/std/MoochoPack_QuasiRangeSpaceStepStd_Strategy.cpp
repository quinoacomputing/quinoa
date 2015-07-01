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

#include "MoochoPack_QuasiRangeSpaceStepStd_Strategy.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace MoochoPack {

bool QuasiRangeSpaceStepStd_Strategy::solve_quasi_range_space_step(
  std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
  ,const Vector& xo, const Vector& c_xo, VectorMutable* v
    )
{
  using LinAlgOpPack::V_InvMtV;
  using LinAlgOpPack::V_StMtV;
  const MatrixOpNonsing
    &R_k = s->R().get_k(0);
  VectorSpace::vec_mut_ptr_t
    vy = R_k.space_rows().create_member();
  // vy = inv(R_k) * c_xo
  V_InvMtV( vy.get(), R_k, BLAS_Cpp::no_trans, c_xo );
  // v = -Y_k*vy
  V_StMtV( v, -1.0, s->Y().get_k(0), BLAS_Cpp::no_trans, *vy );

  return true;
}

void QuasiRangeSpaceStepStd_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
  out << L << "*** Compute the approximate range space step:\n"
    << L << "vy = inv(R_k) * c_xo\n"
    << L << "v = -Y_k*vy\n";
}

} // end namespace MoochoPack
