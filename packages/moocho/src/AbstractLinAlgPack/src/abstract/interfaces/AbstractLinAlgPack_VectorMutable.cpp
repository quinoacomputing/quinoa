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

#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorMutableSubView.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "RTOp_TOp_assign_scalar.h"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_TOp_axpy.h"
#include "RTOp_TOp_set_sub_vector.h"
#include "RTOpPack_RTOpC.hpp"
#include "Teuchos_Assert.hpp"

namespace {

// vector scalar assignment operator
static RTOpPack::RTOpC          assign_scalar_op;
// vector assignment operator
static RTOpPack::RTOpC          assign_vec_op;
// set element operator
static RTOpPack::RTOpC          set_ele_op;
// set sub-vector operator
static RTOpPack::RTOpC          set_sub_vector_op;
// axpy operator
static RTOpPack::RTOpC          axpy_op;

// Simple class for an object that will initialize the operator objects
class init_rtop_server_t {
public:
  init_rtop_server_t() {
    // Vector scalar assignment operator
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_construct(0.0,&assign_scalar_op.op()));
    // Vector assignment operator
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_vectors_construct(&assign_vec_op.op()));
    // Set sub-vector operator
    RTOp_SparseSubVector spc_sub_vec;
    RTOp_sparse_sub_vector_null(&spc_sub_vec);
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_set_sub_vector_construct(&spc_sub_vec,&set_sub_vector_op.op()));
    // axpy operator
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_axpy_construct(0.0,&axpy_op.op()));
  }
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace AbstractLinAlgPack {

// VectorMutable

VectorMutable& VectorMutable::operator=(value_type alpha)
{
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_set_alpha(alpha,&assign_scalar_op.op()));
  VectorMutable* targ_vecs[1] = { this };
  AbstractLinAlgPack::apply_op(assign_scalar_op,0,NULL,1,targ_vecs,NULL);
  return *this;
}

VectorMutable& VectorMutable::operator=(const Vector& vec)
{
  if( dynamic_cast<const void*>(&vec) == dynamic_cast<const void*>(this) )
    return *this; // Assignment to self!
  const Vector*   vecs[1]      = { &vec };
  VectorMutable*  targ_vecs[1] = { this };
  AbstractLinAlgPack::apply_op(assign_vec_op,1,vecs,1,targ_vecs,NULL);
  return *this;
}

VectorMutable& VectorMutable::operator=(const VectorMutable& vec)
{
  return this->operator=(static_cast<const Vector&>(vec));
}

void VectorMutable::set_ele( index_type i, value_type alpha )
{
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_assign_scalar_set_alpha(alpha,&assign_scalar_op.op()));
  VectorMutable* targ_vecs[1] = { this };
  AbstractLinAlgPack::apply_op(
    assign_scalar_op,0,NULL,1,targ_vecs,NULL
    ,i,1,0 // first_ele, sub_dim, global_offset
    );
}

VectorMutable::vec_mut_ptr_t
VectorMutable::sub_view( const Range1D& rng_in )
{
  namespace rcp = MemMngPack;
  const index_type dim = this->dim();
  const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    rng.ubound() > dim, std::out_of_range
    ,"VectorMutable::sub_view(rng): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
    "is not in the range [1,this->dim()] = [1,"<<dim<<"]" );
#endif	
  if( rng.lbound() == 1 && rng.ubound() == dim )
    return Teuchos::rcp( this, false );
  // We are returning a view that could change this vector so we had better
  // wipe out the cache
  //this->has_changed();  // I don't think this line is needed!
  return Teuchos::rcp(
    new VectorMutableSubView(
      Teuchos::rcp( this, false )
      ,rng ) );
}

void VectorMutable::zero()
{
  this->operator=(0.0);
}

void VectorMutable::axpy( value_type alpha, const Vector& x )
{
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_axpy_set_alpha(alpha,&axpy_op.op()));
  const Vector*  vecs[1]      = { &x   };
  VectorMutable* targ_vecs[1] = { this };
  AbstractLinAlgPack::apply_op(axpy_op,1,vecs,1,targ_vecs,NULL);
}

void VectorMutable::get_sub_vector( const Range1D& rng, RTOpPack::MutableSubVector* sub_vec_inout )
{
  //
  // Here we get a copy of the data for the sub-vector that the
  // client will modify.  We must later commit these changes to the
  // actual vector when the client calls commitDetachedView(...).
  // Note, this implementation is very dependent on the behavior of
  // the default implementation of constant version of
  // Vector<Scalar>::get_sub_vector(...) and the implementation of
  // Vector<Scalar>::set_sub_vector(...)!
  //
  RTOpPack::SubVector sub_vec;
  Vector::get_sub_vector( rng, &sub_vec );
  sub_vec_inout->initialize(
    sub_vec.globalOffset(),sub_vec.subDim(),const_cast<value_type*>(sub_vec.values()),sub_vec.stride());
}

void VectorMutable::commit_sub_vector( RTOpPack::MutableSubVector* sub_vec_inout )
{
  RTOpPack::SparseSubVector spc_sub_vec(
    sub_vec_inout->globalOffset(), sub_vec_inout->subDim(),
    Teuchos::arcp(sub_vec_inout->values(), 0, sub_vec_inout->stride()*sub_vec_inout->subDim(), false),
    sub_vec_inout->stride()
    );
  VectorMutable::set_sub_vector( spc_sub_vec );            // Commit the changes!
  RTOpPack::SubVector sub_vec(*sub_vec_inout);
  Vector::free_sub_vector( &sub_vec );                     // Free the memory!
  sub_vec_inout->set_uninitialized();                      // Make null as promised!
}

void VectorMutable::set_sub_vector( const RTOpPack::SparseSubVector& sub_vec )
{
  RTOp_SparseSubVector spc_sub_vec;
  if (!is_null(sub_vec.indices())) {
    RTOp_sparse_sub_vector(
      sub_vec.globalOffset(), sub_vec.subDim(), sub_vec.subNz(),
      sub_vec.values().get(), sub_vec.valuesStride(),
      sub_vec.indices().get(), sub_vec.indicesStride(),
      sub_vec.localOffset(), sub_vec.isSorted(),
      &spc_sub_vec
      );
  }
  else {
    RTOp_SubVector _sub_vec;
    RTOp_sub_vector(
      sub_vec.globalOffset(), sub_vec.subDim(),
      sub_vec.values().get(), sub_vec.valuesStride(),
      &_sub_vec
      );
    RTOp_sparse_sub_vector_from_dense( &_sub_vec, &spc_sub_vec );
  }
  RTOpPack::RTOpC  set_sub_vector_op;
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_set_sub_vector_construct(&spc_sub_vec,&set_sub_vector_op.op()));
  VectorMutable* targ_vecs[1] = { this };
  AbstractLinAlgPack::apply_op(
    set_sub_vector_op,0,NULL,1,targ_vecs,NULL
    ,sub_vec.globalOffset()+1,sub_vec.subDim(),sub_vec.globalOffset() // first_ele, sub_dim, global_offset
    );
}

void VectorMutable::Vp_StMtV(
  value_type                       alpha
  ,const GenPermMatrixSlice        &P
  ,BLAS_Cpp::Transp                P_trans
  ,const Vector                    &x
  ,value_type                      beta
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"VectorMutable::Vp_StMtV(...): Error, this function has not yet been "
    "given a default implementation and has not been overridden for the "
    "subclass \'" << typeName(*this) << "\'!"
    );
  // ToDo: Implement using reduction or transformation operators that will
  // work with any type of vector.
}

// Overridden from Vector

Vector::vec_ptr_t
VectorMutable::sub_view( const Range1D& rng ) const
{
  namespace rcp = MemMngPack;
  return const_cast<VectorMutable*>(this)->sub_view(rng);
}

} // end namespace AbstractLinAlgPack
