// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_HPP
#define THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_HPP


#include "Thyra_DefaultMultiVectorProductVector_decl.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultMultiVectorProductVector<Scalar>::DefaultMultiVectorProductVector()
{
  uninitialize();
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::initialize(
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace_in,
  const RCP<MultiVectorBase<Scalar> > &multiVec
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(productSpace_in));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVec));
  THYRA_ASSERT_VEC_SPACES(
    "DefaultMultiVectorProductVector<Scalar>::initialize(productSpace,multiVec)",
    *multiVec->range(), *productSpace_in->getBlock(0)
    );
  TEUCHOS_ASSERT_EQUALITY( multiVec->domain()->dim(), productSpace_in->numBlocks());
#endif

  numBlocks_ = productSpace_in->numBlocks();

  productSpace_ = productSpace_in;

  multiVec_ = multiVec;

}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::initialize(
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace_in,
  const RCP<const MultiVectorBase<Scalar> > &multiVec
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(productSpace_in));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVec));
  THYRA_ASSERT_VEC_SPACES(
    "DefaultMultiVectorProductVector<Scalar>::initialize(productSpace_in,multiVec)",
    *multiVec->range(), *productSpace_in->getBlock(0)
    );
  TEUCHOS_ASSERT_EQUALITY( multiVec->domain()->dim(), productSpace_in->numBlocks() );
#endif

  numBlocks_ = productSpace_in->numBlocks();

  productSpace_ = productSpace_in;

  multiVec_ = multiVec;

}


template <class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstMultiVector()
{
  return multiVec_.getNonconstObj();
}


template <class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getMultiVector() const
{
  return multiVec_.getConstObj();
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::uninitialize()
{
  numBlocks_ = 0;
  productSpace_ = Teuchos::null;
  multiVec_.uninitialize();
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultMultiVectorProductVector<Scalar>::description() const
{
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description()
    << "{"
    << "dim="<<this->space()->dim()
    << ",numColumns = "<<numBlocks_
    << "}";
  return oss.str();
}

template<class Scalar>
void DefaultMultiVectorProductVector<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::OSTab;
  using Teuchos::describe;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << Teuchos::Describable::description() << "{"
        << "dim=" << this->space()->dim()
        << "}\n";
      OSTab tab2(out);
      *out <<  "multiVec = " << Teuchos::describe(*multiVec_.getConstObj(),verbLevel);
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// Overridden from ProductVectorBase


template <class Scalar>
RCP<VectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstVectorBlock(const int k)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( k, 0, numBlocks_ );
#endif
  return multiVec_.getNonconstObj()->col(k);
}


template <class Scalar>
RCP<const VectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getVectorBlock(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( k, 0, numBlocks_ );
#endif
  return multiVec_.getConstObj()->col(k);
}


// Overridden from ProductMultiVectorBase


template <class Scalar>
RCP<const ProductVectorSpaceBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::productSpace() const
{
  return productSpace_;
}


template <class Scalar>
bool DefaultMultiVectorProductVector<Scalar>::blockIsConst(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( k, 0, numBlocks_ );
#endif
  return multiVec_.isConst();
}


template <class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstMultiVectorBlock(const int k)
{
  return getNonconstVectorBlock(k);
}


template <class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getMultiVectorBlock(const int k) const
{
  return getVectorBlock(k);
}


// Overridden public functions from VectorBase


template <class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::space() const
{
  return productSpace_;
}


// protected


// Overridden protected functions from VectorBase


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::applyOpImpl(
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset
  ) const
{
  this->getDefaultProductVector()->applyOp(
    op, vecs, targ_vecs, reduct_obj, global_offset );
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::acquireDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  this->getDefaultProductVector()->acquireDetachedView(rng_in,sub_vec);
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::releaseDetachedVectorViewImpl(
  RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  this->getDefaultProductVector()->releaseDetachedView(sub_vec);
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::acquireNonconstDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::acquireNonconstDetachedVectorViewImpl(...)!");
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::commitNonconstDetachedVectorViewImpl(...)!");
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::setSubVectorImpl(
  const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
  )
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::setSubVector(...)!");
}


// Overridden protected functions from VectorBase


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::assignImpl(Scalar alpha)
{
  multiVec_.getNonconstObj()->assign(alpha);
}


// private


template <class Scalar>
RCP<const DefaultProductVector<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getDefaultProductVector() const
{

  // This function exists since in general we can not create views of a column
  // vectors and expect the changes to be mirrored in the mulit-vector
  // automatically.  Later, we might be able to change this once we have a
  // Thyra::MultiVectorBase::hasDirectColumnVectorView() function and it
  // returns true.  Until then, this is the safe way to do this ...

  Array<RCP<const VectorBase<Scalar> > > vecArray;
  for ( int k = 0; k < numBlocks_; ++k) {
    vecArray.push_back(multiVec_.getConstObj()->col(k));
  }

  return Thyra::defaultProductVector<Scalar>(
    productSpace_->getDefaultProductVectorSpace(),
    vecArray()
    );

}


} // namespace Thyra


#endif // THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_HPP
