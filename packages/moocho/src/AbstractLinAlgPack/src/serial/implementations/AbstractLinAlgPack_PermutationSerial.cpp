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

#include "AbstractLinAlgPack_PermutationSerial.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "DenseLinAlgPack_PermVecMat.hpp"
#include "DenseLinAlgPack_PermOut.hpp"
#include "Teuchos_Assert.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

PermutationSerial::PermutationSerial( size_type dim )
  :space_(dim)
{}

PermutationSerial::PermutationSerial(
  const i_vector_ptr_t      &perm
  ,const i_vector_ptr_t     &inv_perm
  ,bool                     allocate_missing_perm
  ,bool                     check_inv_perm
  )
{
  this->initialize(perm,inv_perm,allocate_missing_perm,check_inv_perm);
}
  
void PermutationSerial::initialize_identity( size_type dim )
{
  namespace rcp = MemMngPack;
  space_.initialize(dim);
  perm_     = Teuchos::null;
  inv_perm_ = Teuchos::null;
}

void PermutationSerial::initialize(
  const i_vector_ptr_t      &perm
  ,const i_vector_ptr_t     &inv_perm
  ,bool                     allocate_missing_perm
  ,bool                     check_inv_perm
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    perm.get() == NULL && inv_perm.get() == NULL, std::invalid_argument
    ,"PermutationSerial::initialize(...) : Error!" );
  if( perm.get() != NULL && inv_perm.get() != NULL ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      perm->size() != inv_perm->size(), std::invalid_argument
      ,"PermutationSerial::initialize(...) : Error!" );
    if(check_inv_perm) {
      // ToDo: Permform this check!
    }
  }
  space_.initialize( perm.get() ? perm->size() : inv_perm->size() );
  perm_     = perm;
  inv_perm_ = inv_perm;
  if( allocate_missing_perm && perm_.get() == NULL ) {
    Teuchos::RCP<IVector>
      _perm = Teuchos::rcp(new IVector(inv_perm_->size()));
    DenseLinAlgPack::inv_perm( *inv_perm_, _perm.get() );
    perm_ = _perm;
  }
  if( allocate_missing_perm && inv_perm_.get() == NULL ) {
    Teuchos::RCP<IVector>
      _inv_perm = Teuchos::rcp(new IVector(perm_->size()));
    DenseLinAlgPack::inv_perm( *perm_, _inv_perm.get() );
    inv_perm_ = _inv_perm;
  }
}

// Overridden from Permutation

const VectorSpace& PermutationSerial::space() const
{
  return space_;
}

size_type PermutationSerial::dim() const
{
  return space_.dim();
}

bool PermutationSerial::is_identity() const
{
  return perm_.get() == NULL && inv_perm_.get() == NULL;
}

std::ostream& PermutationSerial::output(std::ostream& out) const
{
  const size_type dim = this->dim();
  out << "Serial " << dim << " x " << dim << " permtutation matrix:\n";
  out << "perm =";
  if( perm_.get() )
    out << "\n" << *perm_;
  else
    out << " NULL\n";
  out << "inv_perm =";
  if( inv_perm_.get() )
    out << "\n" << *inv_perm_;
  else
    out << " NULL\n";
  return out;
}

void PermutationSerial::permute( 
  BLAS_Cpp::Transp          P_trans
  ,const Vector       &x
  ,VectorMutable      *y
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    y == NULL, std::invalid_argument
    ,"PermutationSerial::permute(P_trans,x,y) : Error!" );
#endif
#ifdef ABSTRACTLINALGPACK_ASSERT_COMPATIBILITY
  bool is_compatible;
  is_compatible = space_.is_compatible(x.space());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_compatible, std::invalid_argument
    ,"PermutationSerial::permute(P_trans,x,y) : Error, "
    "this->space().is_compatible(x.space()) returned false!" );
  is_compatible = space_.is_compatible(y->space());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_compatible, std::invalid_argument
    ,"PermutationSerial::permute(P_trans,x,y) : Error, "
    "this->space().is_compatible(y->space()) returned false!" );
#endif
  VectorDenseMutableEncap       y_d(*y);
  VectorDenseEncap              x_d(x);
  const IVector                 *p            = NULL;
  bool                          call_inv_perm = false;
  if( ( p = perm_.get() ) != NULL ) {
    if( P_trans == BLAS_Cpp::no_trans )
      call_inv_perm = false;
    else
      call_inv_perm = true;
  }
  else if( ( p = inv_perm_.get() ) != NULL ) {
    if( P_trans == BLAS_Cpp::no_trans )
      call_inv_perm = true;
    else
      call_inv_perm = false;
  }
  if( p ) {
    if( call_inv_perm )
      DenseLinAlgPack::inv_perm_ele( x_d(), *p, &y_d() );
    else
      DenseLinAlgPack::perm_ele( x_d(), *p, &y_d() );
  }
  else {
    // Just the identity permutation, nothing to do!
  }
}

void PermutationSerial::permute( 
  BLAS_Cpp::Transp          P_trans
  ,VectorMutable      *y
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    y == NULL, std::invalid_argument
    ,"PermutationSerial::permute(P_trans,y) : Error!" );
#endif
  VectorSpace::vec_mut_ptr_t
    t = y->clone();
  this->permute(P_trans,*t,y);
}

} // end namespace AbstractLinAlgPack
