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

#include "ConstrainedOptPack_MeritFuncCalc1DQuadratic.hpp"
#include "ConstrainedOptPack_MeritFuncCalc.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "Teuchos_Assert.hpp"

namespace ConstrainedOptPack {

MeritFuncCalc1DQuadratic::MeritFuncCalc1DQuadratic(
  const MeritFuncCalc&      phi
  ,size_type                p
  ,const_VectorWithOp_ptr   d[]
  ,VectorMutable*     x
  )
  : phi_(phi), p_(p), x_(x)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(1 <= p && p <= 2 ), std::invalid_argument
    ,"MeritFuncCalc1DQuadratic::MeritFuncCalc1DQuadratic(...) : Error! "
    "p = " << p << " must be in the range 1 <= p <= 2"
    );
  for( size_type i = 0; i <= p; ++i )
    d_[i] = d[i];
}

value_type MeritFuncCalc1DQuadratic::operator()(value_type alpha) const
{
  using AbstractLinAlgPack::Vp_StV;
  *x_ = *d_[0];
  value_type alpha_i = alpha;
  for( size_type i = 1; i <= p_; ++i, alpha_i *= alpha ) {
    Vp_StV( x_, alpha_i, *d_[i] );
  }
  return phi_( *x_ );
}

value_type  MeritFuncCalc1DQuadratic::deriv() const
{
  return phi_.deriv();
}

void MeritFuncCalc1DQuadratic::print_merit_func(
  std::ostream& out, const std::string& L
  ) const
{
  out	<< L << "*** MeritFuncCalc1DQuadratic\n"
    << L << "x = xo + alpha*d[1] + alpha^2*d[2]\n";
  phi_.print_merit_func( out, L );
}

}	// end namespace ConstrainedOptPack
