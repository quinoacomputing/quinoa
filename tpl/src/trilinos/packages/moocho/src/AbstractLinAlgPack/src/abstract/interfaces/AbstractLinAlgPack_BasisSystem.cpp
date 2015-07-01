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

#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"

namespace AbstractLinAlgPack {

BasisSystem::BasisSystem(
  const mat_sym_fcty_ptr_t             &factory_transDtD
  ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
  )
{
  this->initialize(factory_transDtD,factory_S);
}

void BasisSystem::initialize(
  const mat_sym_fcty_ptr_t             &factory_transDtD
  ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
  )
{
  factory_transDtD_ = factory_transDtD;
  factory_S_        = factory_S;
}

Range1D BasisSystem::equ_decomp() const
{
  const size_type r = this->var_dep().size();
  return r ? Range1D(1,r) : Range1D::Invalid;
}

Range1D BasisSystem::equ_undecomp() const
{
  return Range1D::Invalid;
}

const BasisSystem::mat_fcty_ptr_t BasisSystem::factory_GcUP() const
{
  return Teuchos::null;
}

const BasisSystem::mat_sym_fcty_ptr_t
BasisSystem::factory_transDtD() const
{
  return factory_transDtD_;
}
  
const BasisSystem::mat_sym_nonsing_fcty_ptr_t
BasisSystem::factory_S() const
{
  return factory_S_;
}

} // end namespace AbstractLinAlgPack
