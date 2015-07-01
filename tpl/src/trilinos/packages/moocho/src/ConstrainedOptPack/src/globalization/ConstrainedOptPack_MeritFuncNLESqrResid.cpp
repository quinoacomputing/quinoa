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

#include "ConstrainedOptPack_MeritFuncNLESqrResid.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace ConstrainedOptPack {

MeritFuncNLESqrResid::MeritFuncNLESqrResid()
  : deriv_(0.0)
{}

value_type MeritFuncNLESqrResid::calc_deriv( const Vector& c_k )
{
  using LinAlgOpPack::dot;
  return deriv_ = - dot(c_k,c_k);
}

// Overridden from MeritFuncNLP

value_type MeritFuncNLESqrResid::value(const Vector& c) const
{
  using LinAlgOpPack::dot;
  return 0.5 * dot(c,c);
}

value_type MeritFuncNLESqrResid::deriv() const
{
  return deriv_;
}

void MeritFuncNLESqrResid::print_merit_func(std::ostream& out
  , const std::string& L ) const
{
  out
    << L << "*** Define a square of constraint residuals merit funciton\n"
    << L << "*** (assumes Gc_k'*d_k + c_k = 0):\n"
    << L << "phi(c) = 1/2 * dot(c,c)\n"
    << L << "Dphi(x_k,d_k) = - dot(c_k,c_k)\n";
}

}	// end namespace ConstrainedOptPack 
