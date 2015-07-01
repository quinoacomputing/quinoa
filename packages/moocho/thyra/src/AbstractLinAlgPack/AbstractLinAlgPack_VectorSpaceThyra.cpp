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

#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "AbstractLinAlgPack_VectorSpaceFactoryThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "AbstractLinAlgPack_MultiVectorMutableThyra.hpp"
#include "AbstractLinAlgPack_InnerProductThyra.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

VectorSpaceThyra::VectorSpaceThyra()
{}

VectorSpaceThyra::VectorSpaceThyra(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >    &thyra_vec_spc
  ,const inner_prod_ptr_t                                                  &inner_prod
  )
{
  this->initialize(thyra_vec_spc,inner_prod);
}

void VectorSpaceThyra::initialize(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > &thyra_vec_spc,
  const inner_prod_ptr_t &inner_prod
  )
{
  namespace mmp = MemMngPack;
  TEUCHOS_TEST_FOR_EXCEPTION(
    thyra_vec_spc.get()==NULL, std::invalid_argument
    ,"VectorSpaceThyra::initialize(thyra_vec_spc): Error!"
    );
  thyra_vec_spc_ = thyra_vec_spc;
  if(inner_prod.get())
    this->inner_prod(inner_prod);
  else
    this->inner_prod(Teuchos::rcp(new InnerProductThyra(thyra_vec_spc)));
}

Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >
VectorSpaceThyra::set_uninitialized()
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > tmp_thyra_vec_spc = thyra_vec_spc_;
  thyra_vec_spc_ = Teuchos::null;
  return tmp_thyra_vec_spc;
}

// Overridden from VectorSpace

VectorSpace::space_ptr_t
VectorSpaceThyra::clone() const
{
  return Teuchos::rcp(new VectorSpaceThyra(thyra_vec_spc_->clone()));
}

bool VectorSpaceThyra::is_compatible(const VectorSpace& vec_spc ) const
{
  if( this->dim()==vec_spc.dim() && this->is_in_core() && vec_spc.is_in_core() )
    return true;
  const VectorSpaceThyra
    *thyra_vec_spc = dynamic_cast<const VectorSpaceThyra*>(&vec_spc);
  if( thyra_vec_spc->thyra_vec_spc()->isCompatible(*thyra_vec_spc_) )
    return true;
  return false;
}

bool VectorSpaceThyra::is_in_core() const
{
  return thyra_vec_spc_->hasInCoreView();
}

index_type VectorSpaceThyra::dim() const
{
  return thyra_vec_spc_->dim();
}

VectorSpace::vec_mut_ptr_t
VectorSpaceThyra::create_member() const
{
  return Teuchos::rcp(new VectorMutableThyra(Thyra::createMember(thyra_vec_spc_)));
}

VectorSpace::space_fcty_ptr_t
VectorSpaceThyra::small_vec_spc_fcty() const
{
  return Teuchos::rcp(new VectorSpaceFactoryThyra(thyra_vec_spc_->smallVecSpcFcty()));
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpaceThyra::create_members(size_type num_vecs) const
{
  return Teuchos::rcp(new MultiVectorMutableThyra(Thyra::createMembers(thyra_vec_spc_,num_vecs)));
}

} // end namespace AbstractLinAlgPack
