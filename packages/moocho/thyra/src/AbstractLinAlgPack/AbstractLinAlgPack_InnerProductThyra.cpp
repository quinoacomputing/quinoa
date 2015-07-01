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

#include <stdexcept>

#include "AbstractLinAlgPack_InnerProductThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

InnerProductThyra::InnerProductThyra()
{}

InnerProductThyra::InnerProductThyra(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc
  )
{
  this->initialize(thyra_vec_spc);
}

void InnerProductThyra::initialize(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    thyra_vec_spc.get()==NULL, std::invalid_argument
    ,"InnerProductThyra::initialize(thyra_vec_spc): Error!"
    );
  thyra_vec_spc_ = thyra_vec_spc;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > 
InnerProductThyra::set_uninitialized()
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > tmp_thyra_vec_spc = thyra_vec_spc_;
  thyra_vec_spc_ = Teuchos::null;
  return tmp_thyra_vec_spc;
}

// Overridden from InnerProduct

value_type InnerProductThyra::inner_prod(const Vector& v1, const Vector& v2) const
{
  using Teuchos::dyn_cast;
  return thyra_vec_spc_->scalarProd(
    *dyn_cast<const VectorMutableThyra>(v1).thyra_vec()
    ,*dyn_cast<const VectorMutableThyra>(v1).thyra_vec()
    );
}

} // end namespace AbstractLinAlgPack
