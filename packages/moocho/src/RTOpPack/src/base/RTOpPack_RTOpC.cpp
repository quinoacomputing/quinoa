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

// ///////////////////////////////
// RTOpPack_RTOpC.cpp

#include "RTOpPack_RTOpC.hpp"
#include "Teuchos_Workspace.hpp"


namespace RTOpPack {


RTOpC::RTOpC()
  :RTOpT<RTOp_value_type>("RTOpC") // Should be unused since op_name() if overridden here!
{
  op_.vtbl     = NULL;
  op_.obj_data = NULL;
}


RTOpC::~RTOpC()
{
  if(op_.obj_data)
    RTOp_free_op( &op_ );
}


// Overridden from RTOpT


void RTOpC::get_reduct_type_num_entries_impl(
  const Teuchos::Ptr<int> &num_values,
  const Teuchos::Ptr<int> &num_indexes,
  const Teuchos::Ptr<int> &num_chars
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    0!=RTOp_get_reduct_type_num_entries(&op_,&*num_values,&*num_indexes,&*num_chars)
    ,UnknownError
    ,"RTOpC::get_reduct_type_num_entries(...): Error, "
    "RTOp_get_reduct_type_num_entries(...) returned != 0"
    );
}


Teuchos::RCP<ReductTarget>
RTOpC::reduct_obj_create_impl() const
{
  RTOp_ReductTarget reduct_obj_raw = RTOp_REDUCT_OBJ_NULL;
  TEUCHOS_TEST_FOR_EXCEPTION(
    0!=RTOp_reduct_obj_create(&op_,&reduct_obj_raw)
    ,UnknownError
    ,"RTOpC::reduct_obj_create(...): Error, "
    "RTOp_reduct_obj_create(...) returned != 0"
    );
  return Teuchos::rcp(new ReductTargetC(op_,reduct_obj_raw));
}


void RTOpC::reduce_reduct_objs_impl(
  const ReductTarget &in_reduct_obj,
  const Teuchos::Ptr<ReductTarget> &inout_reduct_obj
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    0!=RTOp_reduce_reduct_objs( &op_, (*this)(in_reduct_obj), (*this)(*inout_reduct_obj) )
    ,UnknownError
    ,"RTOpC::reduce_reduct_objs(...): Error, "
    "RTOp_reduce_reduct_objs(...) returned != 0"
    );
}


void RTOpC::reduct_obj_reinit_impl(
  const Teuchos::Ptr<ReductTarget> &reduct_obj ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    0!=RTOp_reduct_obj_reinit( &op_, (*this)(*reduct_obj) )
    ,UnknownError
    ,"RTOpC::reduct_obj_reinit(...): Error, "
    "RTOp_reduct_obj_reinit(...) returned != 0"
    );
}


void RTOpC::extract_reduct_obj_state_impl(
  const ReductTarget &reduct_obj,
  const Teuchos::ArrayView<primitive_value_type> &value_data,
  const Teuchos::ArrayView<index_type> &index_data,
  const Teuchos::ArrayView<char_type> &char_data
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    0!=RTOp_extract_reduct_obj_state(
      &op_, (*this)(reduct_obj),
      value_data.size(), value_data.getRawPtr(),
      index_data.size(), index_data.getRawPtr(),
      char_data.size(), char_data.getRawPtr()
      ),
    UnknownError,
    "RTOpC::extract_reduct_obj_state(...): Error, "
    "RTOp_extract_reduct_obj_state(...) returned != 0"
    );
}


void RTOpC::load_reduct_obj_state_impl(
  const Teuchos::ArrayView<const primitive_value_type> &value_data,
  const Teuchos::ArrayView<const index_type> &index_data,
  const Teuchos::ArrayView<const char_type> &char_data,
  const Teuchos::Ptr<ReductTarget> &reduct_obj
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    0!=RTOp_load_reduct_obj_state(
      &op_,
      value_data.size(), value_data.getRawPtr(),
      index_data.size(), index_data.getRawPtr(),
      char_data.size(), char_data.getRawPtr(),
      (*this)(*reduct_obj)
      ),
    UnknownError,
    "RTOpC::load_reduct_obj_state(...): Error, "
    "RTOp_load_reduct_obj_state(...) returned != 0"
    );
}


bool RTOpC::coord_invariant_impl() const
{
  return false; // We have to assume this to be safe!
}


std::string RTOpC::op_name_impl() const
{
  const char* op_name = NULL;
  TEUCHOS_TEST_FOR_EXCEPTION(
    0!=RTOp_get_op_name(&op_,&op_name)
    ,UnknownError
    ,"RTOpC::get_op_name(...): Error, "
    "RTOp_op_name(...) returned != 0"
    );
  return op_name;
}


void RTOpC::apply_op_impl(
    const Teuchos::ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const Teuchos::ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Teuchos::Ptr<ReductTarget> &_reduct_obj
  ) const
{

  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss =
    Teuchos::get_default_workspace_store().get();

  const int num_vecs = sub_vecs.size();
  const int num_targ_vecs = targ_sub_vecs.size();

  RTOp_ReductTarget reduct_obj = RTOp_REDUCT_OBJ_NULL;
  if(!is_null(_reduct_obj))
    reduct_obj = (*this)(*_reduct_obj);

  int k;
  Workspace<RTOp_SubVector>        c_sub_vecs(wss,num_vecs,false);
  for( k = 0; k < num_vecs; ++k ) {
    const SubVector& v = sub_vecs[k];
    RTOp_sub_vector(v.globalOffset(),v.subDim(),v.values(),v.stride(),&c_sub_vecs[k]);
  }
  Workspace<RTOp_MutableSubVector>  c_targ_sub_vecs(wss,num_targ_vecs,false);
  for( k = 0; k < num_targ_vecs; ++k ) {
    const MutableSubVector& v = targ_sub_vecs[k];
    RTOp_mutable_sub_vector(v.globalOffset(),v.subDim(),v.values(),v.stride(),&c_targ_sub_vecs[k]);
  }

  const int err = RTOp_apply_op(
    &op_
    ,num_vecs,       num_vecs       ? &c_sub_vecs[0]      : (RTOp_SubVector*)NULL
    ,num_targ_vecs,  num_targ_vecs  ? &c_targ_sub_vecs[0] : (RTOp_MutableSubVector*)NULL
    ,reduct_obj
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    err==RTOp_ERR_INVALID_NUM_VECS, InvalidNumVecs
    ,"RTOpC::apply_op(...): Error, "
    "RTOp_apply_op(...) returned RTOp_ERR_INVALID_NUM_VECS" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    err==RTOp_ERR_INVALID_NUM_TARG_VECS, InvalidNumTargVecs
    ,"RTOpC::apply_op(...): Error, "
    "RTOp_apply_op(...) returned RTOp_ERR_INVALID_NUM_TARG_VECS" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    err!=0, UnknownError
    ,"RTOpC::apply_op(...): Error, "
    "RTOp_apply_op(...) returned != 0 with unknown meaning" );

}


} // namespace RTOpPack
