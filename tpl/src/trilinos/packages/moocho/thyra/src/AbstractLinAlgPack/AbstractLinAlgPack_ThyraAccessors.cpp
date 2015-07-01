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

// ////////////////////////////////////////////////////
// ThyraAccessors.cpp

#include "AbstractLinAlgPack_ThyraAccessors.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "Teuchos_dyn_cast.hpp"

void AbstractLinAlgPack::get_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,const Vector                                                  &vec
  ,Teuchos::RCP<const Thyra::VectorBase<value_type> >    *thyra_vec
  )
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( thyra_vec==NULL, std::invalid_argument, "Error!" );
#endif
  const VectorMutableThyra *vmthyra_vec = dynamic_cast<const VectorMutableThyra*>(&vec);
  if(vmthyra_vec) {
    // We can just grap the const smart pointer to the underlying object 
    *thyra_vec = vmthyra_vec->thyra_vec();
  }
  else if(vec.space().is_in_core()) {
    // We need to create a temporary copy and then copy the explicit elements
    Teuchos::RCP<Thyra::VectorBase<value_type> >
      _thyra_vec = ::Thyra::createMember(thyra_vec_spc.thyra_vec_spc());
    // Get explicit views of the elements
    RTOpPack::SubVector vec_sv;
    vec.get_sub_vector( Range1D(), &vec_sv );
    RTOpPack::SubVectorView<value_type> _thyra_vec_sv;
    _thyra_vec->acquireDetachedView( convert(Range1D()), &_thyra_vec_sv );
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vec_sv.subDim() != _thyra_vec_sv.subDim(), std::logic_error, "Error!" );
#endif
    // Copy the elements
    for( int i = 0; i < vec_sv.subDim(); ++i )
      _thyra_vec_sv[i] = vec_sv[i];
    // Free/commit the explicit views
    vec.free_sub_vector( &vec_sv );
    _thyra_vec->commitDetachedView( &_thyra_vec_sv );
    // Set the output smart pointer
    *thyra_vec = _thyra_vec;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error
      ,"AbstractLinAlgPack::get_thyra_vector(...): Error, the vector of concrete type \'"
      << typeName(vec) << "\' is not an incore vector."
      );
  }
}

void AbstractLinAlgPack::free_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,const Vector                                                  &vec
  ,Teuchos::RCP<const Thyra::VectorBase<value_type> >    *thyra_vec
  )
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( thyra_vec==NULL, std::invalid_argument, "Error!" );
#endif
  *thyra_vec = Teuchos::null;  // This works in both cases above!
}

void AbstractLinAlgPack::get_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,VectorMutable                                                 *vec
  ,Teuchos::RCP<Thyra::VectorBase<value_type> >          *thyra_vec
  )
{ 
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vec==NULL || thyra_vec==NULL, std::invalid_argument, "Error!" );
#endif
  VectorMutableThyra *vmthyra_vec = dynamic_cast<VectorMutableThyra*>(vec);
  if(vmthyra_vec) {
    // We can just directly grap the Thyra vector
    *thyra_vec = vmthyra_vec->set_uninitialized();
  }
  else if(thyra_vec_spc.is_in_core()) {
    // We need to create a temporary copy and then copy the explicit elements
    Teuchos::RCP<Thyra::VectorBase<value_type> >
      _thyra_vec = ::Thyra::createMember(thyra_vec_spc.thyra_vec_spc());
    // Get explicit views of the elements
    RTOpPack::SubVector vec_sv;
    vec->get_sub_vector( Range1D(), &vec_sv );
    RTOpPack::SubVectorView<value_type> _thyra_vec_sv;
    _thyra_vec->acquireDetachedView( convert(Range1D()), &_thyra_vec_sv );
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vec_sv.subDim() != _thyra_vec_sv.subDim(), std::logic_error, "Error!" );
#endif
    // Copy the elements
    for( int i = 0; i < vec_sv.subDim(); ++i )
      _thyra_vec_sv[i] = vec_sv[i];
    // Free/commit the explicit views
    vec->free_sub_vector( &vec_sv );
    _thyra_vec->commitDetachedView( &_thyra_vec_sv );
    // Set the output smart pointer
    *thyra_vec = _thyra_vec;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error
      ,"AbstractLinAlgPack::get_thyra_vector(...): Error, the vector of concrete type \'"
      << typeName(vec) << "\' is not an incore vector."
      );
  }
}

void AbstractLinAlgPack::commit_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,VectorMutable                                                 *vec
  ,Teuchos::RCP<Thyra::VectorBase<value_type> >          *thyra_vec_in
  )
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vec==NULL || thyra_vec_in==NULL, std::invalid_argument, "Error!" );
#endif
  Teuchos::RCP<Thyra::VectorBase<value_type> >  &thyra_vec = *thyra_vec_in;
  VectorMutableThyra *vmthyra_vec = dynamic_cast<VectorMutableThyra*>(vec);
  if(vmthyra_vec) {
    // We can just directly reset the Thyra vector
    vmthyra_vec->initialize(thyra_vec);
  }
  else if(thyra_vec_spc.is_in_core()) {
    // We need to copy back the temporary
    // Get explicit views of the elements
    RTOpPack::ConstSubVectorView<value_type> thyra_vec_sv;
    thyra_vec->acquireDetachedView( convert(Range1D()), &thyra_vec_sv );
    RTOpPack::MutableSubVector vec_sv;
    vec->get_sub_vector( Range1D(), &vec_sv );
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vec_sv.subDim() != thyra_vec_sv.subDim(), std::logic_error, "Error!" );
#endif
    // Copy the elements
    for( int i = 0; i < vec_sv.subDim(); ++i )
      vec_sv[i] = thyra_vec_sv[i];
    // Free/commit the explicit views
    thyra_vec->releaseDetachedView( &thyra_vec_sv );
    vec->commit_sub_vector( &vec_sv );
    // Set the output smart pointer
    thyra_vec = Teuchos::null;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(	true, std::logic_error, "Should never get here?."	);
  }
}
