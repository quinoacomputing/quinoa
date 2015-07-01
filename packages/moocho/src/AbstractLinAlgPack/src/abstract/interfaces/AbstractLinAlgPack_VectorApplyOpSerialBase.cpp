/*
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
*/

// /////////////////////////////////////////////////////
// VectorApplyOpSerialBase.cpp

#include "AbstractLinAlgPack_VectorApplyOpSerialBase.hpp"
#include "AbstractLinAlgPack_apply_op_helper.hpp"

namespace AbstractLinAlgPack {

VectorApplyOpSerialBase::VectorApplyOpSerialBase()
  : in_apply_op_(false)
{}

void VectorApplyOpSerialBase::apply_op_serial(
  const RTOpPack::RTOp& op
  ,const size_t num_vecs, const Vector* vecs[]
  ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    in_apply_op_, std::logic_error
    ,"VectorApplyOpSerialBase::apply_op_serial(...): Error, this function has been entered "
    "recursively which most likely means that the explicit sub-vector access methods Vector::get_sub_vector(...), "
    "Vector::free_sub_vector(...), VectorMutable::get_sub_vector(...), VectorMutable::commit_sub_vector(...) "
    "have not been overridden correctly on this concrete class \'" << typeName(*this) << "\' to not call "
    "apply_op(...) in there implemenations."
    );
  in_apply_op_ = true;
  AbstractLinAlgPack::apply_op_serial(
    op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
    ,first_ele,sub_dim,global_offset
    );
  in_apply_op_ = false;
}

} // namespace AbstractLinAlgPack
