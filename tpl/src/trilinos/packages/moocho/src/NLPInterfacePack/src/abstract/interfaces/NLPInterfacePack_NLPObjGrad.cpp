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

#include "NLPInterfacePack_NLPObjGrad.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"

namespace {
  const char name_Gf[] = "Gf";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPObjGrad::NLPObjGrad()
  : Gf_(NULL)
{}

void NLPObjGrad::initialize(bool test_setup) {
  num_Gf_evals_ = 0;
  NLP::initialize(test_setup);
}

// Information

bool NLPObjGrad::supports_Gf() const
{
  return true;
}

bool NLPObjGrad::supports_Gf_prod() const
{
  return false;
}

// <<std aggr>> members for Gf

void NLPObjGrad::set_Gf(VectorMutable* Gf)
{
  Gf_ = Gf;
}

AbstractLinAlgPack::VectorMutable* NLPObjGrad::get_Gf()
{
  return StandardCompositionRelationshipsPack::get_role_name(Gf_, false, name_Gf);
}

AbstractLinAlgPack::VectorMutable& NLPObjGrad::Gf()
{
  return StandardCompositionRelationshipsPack::role_name(Gf_, false, name_Gf);
}

const AbstractLinAlgPack::Vector& NLPObjGrad::Gf() const
{
  return StandardCompositionRelationshipsPack::role_name(Gf_, false, name_Gf);
}

void NLPObjGrad::unset_quantities()
{
  NLP::unset_quantities();
  Gf_ = NULL;
}

// calculations

void NLPObjGrad::calc_Gf(const Vector& x, bool newx) const
{
  StandardCompositionRelationshipsPack::assert_role_name_set(Gf_, "NLP::calc_Gf()", name_Gf);
  imp_calc_Gf(x,newx,obj_grad_info());
  num_Gf_evals_++;
}

value_type NLPObjGrad::calc_Gf_prod(const Vector& x, const Vector& d, bool newx) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"Error, the function calc_Gf_prod(...) is not implemented for the class "
    << typeName(*this) << "!"
    );

  //execution should never reach this point, but compilers expect a non-void
  //function to return something. So we'll create a dummy value to use in a
  //return statement.
  //(a better design would not require function bodies for unimplemented
  //functions like this...)
  value_type* dummy = NULL;
  return(*dummy);
}

size_type NLPObjGrad::num_Gf_evals() const
{
  return num_Gf_evals_;
}

}	// end namespace NLPInterfacePack 
