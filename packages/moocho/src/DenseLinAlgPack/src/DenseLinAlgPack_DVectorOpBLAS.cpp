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

#include <algorithm>
#include <functional>
#include <cmath>

#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_BLAS_Cpp.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace {

using DenseLinAlgPack::DVector;
using DenseLinAlgPack::DVectorSlice;
typedef DenseLinAlgPack::value_type value_type;
typedef DVectorSlice::size_type size_type;
typedef DVectorSlice::difference_type difference_type;

//
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
//
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }

// Utility for copying vector slices.  Takes care of aliasing etc. but not sizes.
void i_assign(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  switch(vs_lhs->overlap(vs_rhs)) {
    case DenseLinAlgPack::SAME_MEM:	return; // assignment to self so skip the copy
    case DenseLinAlgPack::SOME_OVERLAP:		// create a temp for the copy
    {
      DVector temp(vs_rhs);
      BLAS_Cpp::copy( temp.dim(), temp.raw_ptr(), 1, vs_lhs->raw_ptr(), vs_lhs->stride() );
      return;
    }
    default:				// no overlap so just perform the copy
    {
      BLAS_Cpp::copy( vs_rhs.dim(),vs_rhs.raw_ptr(), vs_rhs.stride()
        , vs_lhs->raw_ptr(), vs_lhs->stride() );
      return;
    }
  }
}

inline value_type local_prod( const value_type &v1, const value_type &v2 ) { return v1*v2; }

}	// end namespace

// ///////////////////////
// op= operations

void DenseLinAlgPack::Vp_S(DVectorSlice* vs_lhs, value_type alpha) {
  if(vs_lhs->stride() == 1) {
    DVectorSlice::value_type
      *itr		= vs_lhs->start_ptr(),
      *itr_end	= itr + vs_lhs->dim();
    while(itr != itr_end)
      *itr++ += alpha;
  }
  else {
    DVectorSlice::iterator
      itr		= vs_lhs->begin(),
      itr_end	= vs_lhs->end();
    while(itr != itr_end)
      *itr++ += alpha;
  }
}

void DenseLinAlgPack::Vt_S(DVectorSlice* vs_lhs, value_type alpha) {
  if( alpha == 1.0 )
    return;
  else if( alpha == 0.0 )
    *vs_lhs = 0.0;
  else
    BLAS_Cpp::scal( vs_lhs->dim(), alpha, vs_lhs->raw_ptr(), vs_lhs->stride() );
}

void DenseLinAlgPack::Vp_StV(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  Vp_V_assert_sizes(vs_lhs->dim(), vs_rhs.dim());
  BLAS_Cpp::axpy( vs_lhs->dim(), alpha, vs_rhs.raw_ptr(), vs_rhs.stride()
    , vs_lhs->raw_ptr(), vs_lhs->stride());
}

// DVectorSlice as lhs

void DenseLinAlgPack::assign(DVectorSlice* vs_lhs, value_type alpha) {
  if(!vs_lhs->dim())
    throw std::length_error("DenseLinAlgPack::assign(vs_lhs,alpha): DVectorSlice must be bound and sized.");
  BLAS_Cpp::copy( vs_lhs->dim(), &alpha, 0, vs_lhs->raw_ptr(), vs_lhs->stride() );
}
void DenseLinAlgPack::assign(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs.dim() );
  i_assign(vs_lhs,vs_rhs);
}
void DenseLinAlgPack::V_VpV(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2) {
  VopV_assert_sizes( vs_rhs1.dim(), vs_rhs2.dim() );
  Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs1.dim() );
  std::transform(vs_rhs1.begin(),vs_rhs1.end(),vs_rhs2.begin(),vs_lhs->begin(),std::plus<value_type>());
}
void DenseLinAlgPack::V_VmV(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2) {
  VopV_assert_sizes( vs_rhs1.dim(), vs_rhs2.dim() );
  Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs1.dim() );
  std::transform(vs_rhs1.begin(),vs_rhs1.end(),vs_rhs2.begin(),vs_lhs->begin(),std::minus<value_type>());
}
void DenseLinAlgPack::V_mV(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs.dim() );
  (*vs_lhs) = vs_rhs;
  BLAS_Cpp::scal( vs_lhs->dim(), -1.0, vs_lhs->raw_ptr(), vs_lhs->stride() );
}
void DenseLinAlgPack::V_StV(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs)
{
  Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs.dim() );
  (*vs_lhs) = vs_rhs;
  BLAS_Cpp::scal( vs_lhs->dim(), alpha, vs_lhs->raw_ptr(), vs_lhs->stride() );
}

// Elementwise math vector functions. VC++ 5.0 is not allowing use of ptr_fun() with overloaded
// functions so I have to perform the loops straight out.
// ToDo: use ptr_fun() when you get a compiler that works this out.  For now I will just use macros.
#define UNARYOP_VECSLC(LHS, RHS, FUNC)																		\
  DenseLinAlgPack::Vp_V_assert_sizes( (LHS)->dim(), (RHS).dim() );											\
  DVectorSlice::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs;										\
  for(itr_lhs = LHS->begin(), itr_rhs = RHS.begin(); itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs)			\
  {	*itr_lhs = FUNC(*itr_rhs); }


#define BINARYOP_VECSLC(LHS, RHS1, RHS2, FUNC)																\
  DenseLinAlgPack::VopV_assert_sizes( (RHS1).dim(), (RHS2).dim() );											\
  DenseLinAlgPack::Vp_V_assert_sizes( (LHS)->dim(), (RHS1).dim() );											\
  DVectorSlice::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs1, itr_rhs2;							\
  for(itr_lhs = LHS->begin(), itr_rhs1 = RHS1.begin(), itr_rhs2 = RHS2.begin();							\
    itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)											\
  {	*itr_lhs = FUNC(*itr_rhs1, *itr_rhs2); }

#define BINARYOP_BIND1ST_VECSLC(LHS, RHS1, RHS2, FUNC)														\
  DenseLinAlgPack::Vp_V_assert_sizes( (LHS)->dim(), (RHS2).dim() );											\
  DVectorSlice::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs2;									\
  for(itr_lhs = LHS->begin(), itr_rhs2 = RHS2.begin();													\
    itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs2)														\
  {	*itr_lhs = FUNC(RHS1, *itr_rhs2); }

#define BINARYOP_BIND2ND_VECSLC(LHS, RHS1, RHS2, FUNC)														\
  DenseLinAlgPack::Vp_V_assert_sizes( (LHS)->dim(), (RHS1).dim());											\
  DVectorSlice::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs1;									\
  for(itr_lhs = LHS->begin(), itr_rhs1 = RHS1.begin();													\
    itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs1)														\
  {	*itr_lhs = FUNC(*itr_rhs1, RHS2); }

void DenseLinAlgPack::abs(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::fabs)
}
void DenseLinAlgPack::asin(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::asin)
}
void DenseLinAlgPack::acos(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::acos)
}
void DenseLinAlgPack::atan(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::atan)
}
void DenseLinAlgPack::atan2(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2) {
  BINARYOP_VECSLC(vs_lhs, vs_rhs1, vs_rhs2, ::atan2)
}
void DenseLinAlgPack::atan2(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs, value_type alpha) {
  BINARYOP_BIND2ND_VECSLC(vs_lhs, vs_rhs, alpha, ::atan2)
}
void DenseLinAlgPack::atan2(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs)
{
  BINARYOP_BIND1ST_VECSLC(vs_lhs, alpha, vs_rhs, ::atan2)
}
void DenseLinAlgPack::cos(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::cos)
}
void DenseLinAlgPack::cosh(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::cosh)
}
void DenseLinAlgPack::exp(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::exp)
}
void DenseLinAlgPack::max(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2) {
  DenseLinAlgPack::VopV_assert_sizes( vs_rhs1.dim(), vs_rhs2.dim() );
  DenseLinAlgPack::Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs1.dim() );
  DVectorSlice::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs1, itr_rhs2;
  for(itr_lhs = vs_lhs->begin(), itr_rhs1 = vs_rhs1.begin(), itr_rhs2 = vs_rhs2.begin();
    itr_lhs != vs_lhs->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)
  {	*itr_lhs = my_max(*itr_rhs1, *itr_rhs2); }
}
void DenseLinAlgPack::max(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  DenseLinAlgPack::Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs.dim() );
  DVectorSlice::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs;
  for(itr_lhs = vs_lhs->begin(), itr_rhs = vs_rhs.begin(); itr_lhs != vs_lhs->end(); ++itr_lhs, ++itr_rhs)
  {	*itr_lhs = my_max(alpha,*itr_rhs); }
}
void DenseLinAlgPack::min(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2) {
  DenseLinAlgPack::VopV_assert_sizes( vs_rhs1.dim(), vs_rhs2.dim() );
  DenseLinAlgPack::Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs1.dim() );
  DVectorSlice::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs1, itr_rhs2;
  for(itr_lhs = vs_lhs->begin(), itr_rhs1 = vs_rhs1.begin(), itr_rhs2 = vs_rhs2.begin();
    itr_lhs != vs_lhs->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)
  {	*itr_lhs = my_min(*itr_rhs1, *itr_rhs2); }
}
void DenseLinAlgPack::min(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  DenseLinAlgPack::Vp_V_assert_sizes( vs_lhs->dim(), vs_rhs.dim() );
  DVectorSlice::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs;
  for(itr_lhs = vs_lhs->begin(), itr_rhs = vs_rhs.begin(); itr_lhs != vs_lhs->end(); ++itr_lhs, ++itr_rhs)
  {	*itr_lhs = my_min(alpha,*itr_rhs); }
}
void DenseLinAlgPack::pow(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2) {
  BINARYOP_VECSLC(vs_lhs, vs_rhs1, vs_rhs2, std::pow)
}
void DenseLinAlgPack::pow(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs, value_type alpha) {
  BINARYOP_BIND2ND_VECSLC(vs_lhs, vs_rhs, alpha, std::pow)
}
void DenseLinAlgPack::pow(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs, int n) {
  BINARYOP_BIND2ND_VECSLC(vs_lhs, vs_rhs, n, std::pow)
}
void DenseLinAlgPack::pow(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  BINARYOP_BIND1ST_VECSLC(vs_lhs, alpha, vs_rhs, std::pow)
}
void DenseLinAlgPack::prod( DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2 ) {
  BINARYOP_VECSLC(vs_lhs, vs_rhs1, vs_rhs2, local_prod )
}
void DenseLinAlgPack::sqrt(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::sqrt)
}
void DenseLinAlgPack::sin(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::sin)
}
void DenseLinAlgPack::sinh(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::sinh)
}
void DenseLinAlgPack::tan(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::tan)
}
void DenseLinAlgPack::tanh(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VECSLC(vs_lhs, vs_rhs, ::tanh)
}

// DVector as lhs

void DenseLinAlgPack::assign(DVector* v_lhs, value_type alpha)
{
  if(!v_lhs->dim())
    throw std::length_error("DenseLinAlgPack::assign(v_lhs,alpha): DVector must be sized.");
  else
    BLAS_Cpp::copy( v_lhs->dim(), &alpha, 0, v_lhs->raw_ptr(), v_lhs->stride() );
}
void DenseLinAlgPack::assign(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  v_lhs->resize(vs_rhs.dim());
  i_assign( &(*v_lhs)(), vs_rhs );
}
void DenseLinAlgPack::V_VpV(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2)
{
  assert_vs_sizes(vs_rhs1.dim(), vs_rhs2.dim());
  v_lhs->resize(vs_rhs1.dim());
  std::transform(vs_rhs1.begin(),vs_rhs1.end(),vs_rhs2.begin(),v_lhs->begin(),std::plus<value_type>());
//	DVectorSlice::const_iterator
//		v1 = vs_rhs1.begin(),
//		v2 = vs_rhs2.begin();
//	DVector::iterator
//		v = v_lhs->begin(),
//		v_end = v_lhs->end();
//	while( v != v_end )
//		*v++ = *v1++ + *v2++;
}
void DenseLinAlgPack::V_VmV(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2)
{
  assert_vs_sizes(vs_rhs1.dim(), vs_rhs2.dim());
  v_lhs->resize(vs_rhs1.dim());
  std::transform(vs_rhs1.begin(),vs_rhs1.end(),vs_rhs2.begin(),v_lhs->begin(),std::minus<value_type>());
}
void DenseLinAlgPack::V_mV(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  v_lhs->resize(vs_rhs.dim());
  (*v_lhs) = vs_rhs;
  BLAS_Cpp::scal(v_lhs->dim(), -1.0, v_lhs->raw_ptr(), 1);
}
void DenseLinAlgPack::V_StV(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  v_lhs->resize(vs_rhs.dim());
  (*v_lhs) = vs_rhs;
  BLAS_Cpp::scal( v_lhs->dim(), alpha, v_lhs->raw_ptr(), 1);
}

void DenseLinAlgPack::rot( const value_type c, const value_type s, DVectorSlice* x, DVectorSlice* y )
{
  assert_vs_sizes( x->dim(), y->dim() );
  BLAS_Cpp::rot( x->dim(), x->raw_ptr(), x->stride(), y->raw_ptr(), y->stride(), c, s );
}

// Elementwise math vector functions. VC++ 5.0 is not allowing use of ptr_fun() with overloaded
// functions so I have to perform the loops straight out.
// ToDo: use ptr_fun() when you get a compiler that works this out.  For now I will just use macros.
#define UNARYOP_VEC(LHS, RHS, FUNC)																			\
  LHS->resize(RHS.dim());																				\
  DVector::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs;											\
  for(itr_lhs = LHS->begin(), itr_rhs = RHS.begin(); itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs)			\
  {	*itr_lhs = FUNC(*itr_rhs); }

#define BINARYOP_VEC(LHS, RHS1, RHS2, FUNC)																	\
  DenseLinAlgPack::assert_vs_sizes(RHS1.dim(), RHS2.dim()); LHS->resize(RHS1.dim());						\
  DVector::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs1, itr_rhs2;								\
  for(itr_lhs = LHS->begin(), itr_rhs1 = RHS1.begin(), itr_rhs2 = RHS2.begin();							\
    itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)											\
  {	*itr_lhs = FUNC(*itr_rhs1, *itr_rhs2); }

#define BINARYOP_BIND1ST_VEC(LHS, RHS1, RHS2, FUNC)															\
  LHS->resize(RHS2.dim());																				\
  DVector::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs2;											\
  for(itr_lhs = LHS->begin(), itr_rhs2 = RHS2.begin();													\
    itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs2)														\
  {	*itr_lhs = FUNC(RHS1, *itr_rhs2); }

#define BINARYOP_BIND2ND_VEC(LHS, RHS1, RHS2, FUNC)															\
  LHS->resize(RHS1.dim());																				\
  DVector::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs1;											\
  for(itr_lhs = LHS->begin(), itr_rhs1 = RHS1.begin();													\
    itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs1)														\
  {	*itr_lhs = FUNC(*itr_rhs1, RHS2); }

void DenseLinAlgPack::abs(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::fabs)
}
void DenseLinAlgPack::asin(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::asin)
}
void DenseLinAlgPack::acos(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::acos)
}
void DenseLinAlgPack::atan(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::atan)
}
void DenseLinAlgPack::atan2(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2)
{
  BINARYOP_VEC(v_lhs, vs_rhs1, vs_rhs2, ::atan2)
}
void DenseLinAlgPack::atan2(DVector* v_lhs, const DVectorSlice& vs_rhs, value_type alpha) {
  BINARYOP_BIND2ND_VEC(v_lhs, vs_rhs, alpha, ::atan2)
}
void DenseLinAlgPack::atan2(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  BINARYOP_BIND1ST_VEC(v_lhs, alpha, vs_rhs, ::atan2)
}
void DenseLinAlgPack::cos(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::cos)
}
void DenseLinAlgPack::cosh(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::cosh)
}
void DenseLinAlgPack::exp(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::exp)
}
void DenseLinAlgPack::max(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2)
{
  DenseLinAlgPack::assert_vs_sizes(vs_rhs1.dim(), vs_rhs2.dim()); v_lhs->resize(vs_rhs1.dim());
  DVector::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs1, itr_rhs2;
  for(itr_lhs = v_lhs->begin(), itr_rhs1 = vs_rhs1.begin(), itr_rhs2 = vs_rhs2.begin();
    itr_lhs != v_lhs->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)
  {	*itr_lhs = my_max(*itr_rhs1, *itr_rhs2); }
}
void DenseLinAlgPack::max(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  v_lhs->resize(vs_rhs.dim());
  DVector::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs;
  for(itr_lhs = v_lhs->begin(), itr_rhs = vs_rhs.begin(); itr_lhs != v_lhs->end(); ++itr_lhs, ++itr_rhs)
  {	*itr_lhs = my_max(alpha,*itr_rhs); }
}
void DenseLinAlgPack::min(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2) 
{
  DenseLinAlgPack::assert_vs_sizes(vs_rhs1.dim(), vs_rhs2.dim()); v_lhs->resize(vs_rhs1.dim());
  DVector::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs1, itr_rhs2;
  for(itr_lhs = v_lhs->begin(), itr_rhs1 = vs_rhs1.begin(), itr_rhs2 = vs_rhs2.begin();
    itr_lhs != v_lhs->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)
  {	*itr_lhs = my_min(*itr_rhs1, *itr_rhs2); }
}
void DenseLinAlgPack::min(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  v_lhs->resize(vs_rhs.dim());
  DVector::iterator itr_lhs; DVectorSlice::const_iterator itr_rhs;
  for(itr_lhs = v_lhs->begin(), itr_rhs = vs_rhs.begin(); itr_lhs != v_lhs->end(); ++itr_lhs, ++itr_rhs)
  {	*itr_lhs = my_min(alpha,*itr_rhs); }
}
void DenseLinAlgPack::pow(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2) {
  BINARYOP_VEC(v_lhs, vs_rhs1, vs_rhs2, std::pow)
}
void DenseLinAlgPack::pow(DVector* v_lhs, const DVectorSlice& vs_rhs, value_type alpha) {
  BINARYOP_BIND2ND_VEC(v_lhs, vs_rhs, alpha, std::pow)
}
void DenseLinAlgPack::pow(DVector* v_lhs, const DVectorSlice& vs_rhs, int n) {
  BINARYOP_BIND2ND_VEC(v_lhs, vs_rhs, n, std::pow)
}
void DenseLinAlgPack::pow(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs) {
  BINARYOP_BIND1ST_VEC(v_lhs, alpha, vs_rhs, std::pow)
}
void DenseLinAlgPack::prod( DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2 ) {
  BINARYOP_VEC(v_lhs, vs_rhs1, vs_rhs2, local_prod )
}
void DenseLinAlgPack::sqrt(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::sqrt)
}
void DenseLinAlgPack::sin(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::sin)
}
void DenseLinAlgPack::sinh(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::sinh)
}
void DenseLinAlgPack::tan(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::tan)
}
void DenseLinAlgPack::tanh(DVector* v_lhs, const DVectorSlice& vs_rhs) {
  UNARYOP_VEC(v_lhs, vs_rhs, ::tanh)
}

// //////////////////////////////
// Scalar returning functions

DenseLinAlgPack::value_type DenseLinAlgPack::dot(const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2)
{
  VopV_assert_sizes( vs_rhs1.dim(), vs_rhs2.dim() );
  return BLAS_Cpp::dot( vs_rhs1.dim(), vs_rhs1.raw_ptr(), vs_rhs1.stride()
    , vs_rhs2.raw_ptr(), vs_rhs2.stride() );
}
DenseLinAlgPack::value_type DenseLinAlgPack::max(const DVectorSlice& vs_rhs)
{	return *std::max_element(vs_rhs.begin(),vs_rhs.end()); }
DenseLinAlgPack::value_type DenseLinAlgPack::min(const DVectorSlice& vs_rhs)
{	return *std::min_element(vs_rhs.begin(),vs_rhs.end()); }
DenseLinAlgPack::value_type DenseLinAlgPack::norm_1(const DVectorSlice& vs_rhs)
{	return BLAS_Cpp::asum( vs_rhs.dim(), vs_rhs.raw_ptr(), vs_rhs.stride() ); }
DenseLinAlgPack::value_type DenseLinAlgPack::norm_2(const DVectorSlice& vs_rhs) 
{	return BLAS_Cpp::nrm2( vs_rhs.dim(), vs_rhs.raw_ptr(), vs_rhs.stride() ); }
DenseLinAlgPack::value_type DenseLinAlgPack::norm_inf(const DVectorSlice& vs_rhs) {
//	return BLAS_Cpp::iamax( vs_rhs.dim(), vs_rhs.raw_ptr(), vs_rhs.stride() );
//	For some reason iamax() is not working properly?
  value_type max_ele = 0.0;
  for(DVectorSlice::const_iterator itr = vs_rhs.begin(); itr != vs_rhs.end(); ) {
    value_type ele = ::fabs(*itr++);
    if(ele > max_ele) max_ele = ele;
  }
  return max_ele;
}

// /////////////////////
// Misc operations

void DenseLinAlgPack::swap( DVectorSlice* vs1, DVectorSlice* vs2 ) {
  if( vs1->overlap(*vs2) == SAME_MEM ) return;
  VopV_assert_sizes( vs1->dim(), vs2->dim() );
  BLAS_Cpp::swap( vs1->dim(), vs1->raw_ptr(), vs1->stride(), vs2->raw_ptr(), vs2->stride() );
}
