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

#include "ConstrainedOptPack_MeritFuncCalcNLP.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

namespace ConstrainedOptPack {

MeritFuncCalcNLP::MeritFuncCalcNLP( const MeritFuncNLP* phi, const NLP* nlp )
  : phi_(phi), nlp_(nlp)
{}

value_type MeritFuncCalcNLP::operator()(const Vector& x) const
{
  const size_type
    m  = nlp().m(),
    ns = nlp().ns();
  nlp().calc_f(x);
  if(m)  nlp().calc_c(x,false);
  return phi().value(
    nlp().f()
    ,m  ? &nlp().c()  : NULL
    ,NULL  // h
    ,NULL  // hl
    ,NULL  // hu
    );
/* RAB: 20020112: ToDo: Get this working
  if(m)  nlp().calc_c_breve(x,false);
  if(ns) nlp().calc_h_breve(x,false);
  return phi().value(
    nlp().f()
    ,m  ? &nlp().c_breve()  : NULL
    ,ns ? &nlp().h_breve()  : NULL
    ,ns ? &nlp().hl_breve() : NULL
    ,ns ? &nlp().hu_breve() : NULL
    );
*/
}

value_type MeritFuncCalcNLP::deriv() const {
  return phi().deriv();
}

void MeritFuncCalcNLP::print_merit_func(
  std::ostream& out, const std::string& L
  ) const
{
  out	<< L << "*** MeritFuncCalcNLP\n"
    << L << "f = f(x), c = c_breve(x_breve), h = h_breve(x_breve)\n";
  phi().print_merit_func(out,L);
}

}	// end namespace ConstrainedOptPack
