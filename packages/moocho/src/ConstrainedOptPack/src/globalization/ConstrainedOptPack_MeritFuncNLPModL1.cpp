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

#include "ConstrainedOptPack_MeritFuncNLPModL1.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "Teuchos_Assert.hpp"

namespace ConstrainedOptPack {

MeritFuncNLPModL1::MeritFuncNLPModL1()
  : deriv_(0.0)
{}

// Overridden from MeritFuncNLP

value_type MeritFuncNLPModL1::value(
  value_type             f
  ,const Vector    *c
  ,const Vector    *h
  ,const Vector    *hl
  ,const Vector    *hu
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    h || hl || hu, std::logic_error
    ,"MeritFuncNLPModL1::value(...) : Error! general inequalities are not supported!" );
/*
  using DenseLinAlgPack::norm_1;
  return f + local_constr_term( mu_, c, "calc_deriv" );
*/
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Write a reduction operator for the above operation
  return 0.0;
}

value_type MeritFuncNLPModL1::deriv() const
{
  return deriv_;
}

void MeritFuncNLPModL1::print_merit_func(
  std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Define a modified L1 merit funciton that uses different\n"
    << L << "*** penalty parameters for each constriant.\n"
    << L << "*** (assumes Gc_k'*d_k + c_k = 0):\n"
    << L << "phi(f,c) = f + sum( mu(j) * abs(c(j)), j = 1,...,m )\n"
    << L << "Dphi(x_k,d_k) = Gf_k' * d_k - sum( mu(j) * abs(c(j)), j = 1,...,m )\n";
}

// Overridden from MeritFuncNLPDirecDeriv

value_type MeritFuncNLPModL1::calc_deriv(
  const Vector    &Gf_k
  ,const Vector   *c_k
  ,const Vector   *h_k
  ,const Vector   *hl
  ,const Vector   *hu
  ,const Vector   &d_k
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    h_k || hl || hu, std::logic_error
    ,"MeritFuncNLPModL1::value(...) : Error! general inequalities are not supported!" );
/*
  using DenseLinAlgPack::dot; using DenseLinAlgPack::norm_1;
  return deriv_ = dot( Gf_k, d_k ) - local_constr_term( mu_, c_k, "calc_deriv" );
*/
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Write a reduction operator for the above operation
  return 0.0;
}

// Overridden from MeritFuncPenaltyParam

void MeritFuncNLPModL1::set_space_c( const VectorSpace::space_ptr_t& space_c )
{
  mu_  = space_c->create_member();
  *mu_ = 0.0;
}

VectorMutable& MeritFuncNLPModL1::set_mu()
{
  return *mu_;
}

const Vector& MeritFuncNLPModL1::get_mu() const
{
  return *mu_;
}

}	// end namespace ConstrainedOptPack

/* ToDo: Write a reduction operator for the following!

namespace {

value_type local_constr_term( const DVector& mu, const DVectorSlice& c
  , const char func_name[] )
{
  if( mu.size() != c.size() ) {
    std::ostringstream omsg;
    omsg
      << "MeritFuncNLPModL1::" << func_name << "(...) : "
      << "Error, the sizes mu.size() == " << mu.size()
      << " != c.size() == " << c.size();
    throw ConstrainedOptPack::MeritFuncNLP::InvalidInitialization(omsg.str());
  }
  value_type r = 0.0;
  DVector::const_iterator
    mu_itr = mu.begin();
  DVectorSlice::const_iterator
    c_itr = c.begin();
  while( mu_itr != mu.end() )
    r += *mu_itr++ * ::fabs( *c_itr++ );
  return r;
}

}	// end namespace

*/
