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

#include <assert.h>

#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Teuchos_Assert.hpp"

namespace NLPInterfacePack {

// NLPDirect

void NLPDirect::set_factories(
  const mat_sym_fcty_ptr_t             &factory_transDtD
  ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
  )
{
  factory_transDtD_ = factory_transDtD;
  factory_S_        = factory_S;
}

size_type NLPDirect::r() const
{
  return this->con_decomp().size();
}

Range1D NLPDirect::var_dep() const
{
  return Range1D(1,m());
}
Range1D NLPDirect::var_indep() const
{
  return Range1D(m()+1,n());
}
Range1D NLPDirect::con_decomp() const
{
  return Range1D(1,m());
}

Range1D NLPDirect::con_undecomp() const
{
  return Range1D::Invalid;
}

const NLPDirect::mat_fcty_ptr_t
NLPDirect::factory_GcU() const
{
  return Teuchos::null;
}

const NLPDirect::mat_fcty_ptr_t
NLPDirect::factory_Uz() const
{
  return Teuchos::null;
}

const NLPDirect::mat_fcty_ptr_t
NLPDirect::factory_GcUD() const
{
  return Teuchos::null;
}

const NLPDirect::mat_sym_fcty_ptr_t
NLPDirect::factory_transDtD() const
{
  return factory_transDtD_;
}
  
const NLPDirect::mat_sym_nonsing_fcty_ptr_t
NLPDirect::factory_S() const
{
  return factory_S_;
}

void NLPDirect::initialize(bool test_setup)
{
  NLPObjGrad::initialize(test_setup);
}

} // end namespace NLPIntefacePack
