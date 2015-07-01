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

#include "AbstractLinAlgPack_SpVectorOp.hpp"

namespace {
// Setup some template classes to check at complile time that
// the layout of SpVectorSlice::element_type is proper.
template<int N, class T>
class assert_compile_time {
public:
assert_compile_time()
{
  // This should not compile if instantiated with a type T that
  // is not an integer.  However, if the compiler checks this
  // function without instantiating it, it can not cause an error
  // because it does not know the type of T to see if the
  // conversion is legal or not.
  T d;
  static_cast<int*>(d);
}
};
// Template specialization for no error
template<>
class assert_compile_time<0,double> {
public:
assert_compile_time()
{}
};
// Validate that there is an integer stride between values
assert_compile_time<
    ((int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
   % (int)sizeof(DenseLinAlgPack::value_type))
  , double
  >
    validate_value_stride;
// Validate that there is an integer stride between indexes
assert_compile_time<
    ((int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
   % (int)sizeof(DenseLinAlgPack::index_type))
  , double
  >
    validate_index_stride;
} // end namespace

void AbstractLinAlgPack::add_elements( SpVector* sv_lhs, value_type alpha, const DVectorSlice& vs_rhs
                   , size_type offset, bool add_zeros )
{
  typedef SpVector::element_type ele_t;
  const bool assume_sorted = !sv_lhs->nz() || ( sv_lhs->nz() && sv_lhs->is_sorted() );
  DVectorSlice::const_iterator
    itr = vs_rhs.begin();
  if(add_zeros) {
    for( size_type i = 1; i <= vs_rhs.dim(); ++i )
      sv_lhs->add_element( ele_t( i + offset, alpha * (*itr++) ) );
  }
  else {
    for( size_type i = 1; i <= vs_rhs.dim(); ++i, ++itr )
      if( *itr != 0.0 )
        sv_lhs->add_element( ele_t( i + offset, alpha * (*itr) ) );
  }
  sv_lhs->assume_sorted(assume_sorted);
}

void AbstractLinAlgPack::add_elements( SpVector* sv_lhs, value_type alpha, const SpVectorSlice& sv_rhs
                   , size_type offset, bool add_zeros )
{
  typedef SpVector::element_type ele_t;
  const bool assume_sorted = ( !sv_lhs->nz() || ( sv_lhs->nz() && sv_lhs->is_sorted() ) )
    && ( !sv_rhs.nz() || ( sv_rhs.nz() || sv_rhs.is_sorted() ) );
  if(add_zeros) {
    for( SpVectorSlice::const_iterator itr = sv_rhs.begin(); itr != sv_rhs.end(); ++itr )
      sv_lhs->add_element( ele_t( itr->index() + sv_rhs.offset() + offset, alpha * (itr->value()) ) );
  }
  else {
    for( SpVectorSlice::const_iterator itr = sv_rhs.begin(); itr != sv_rhs.end(); ++itr )
      if(itr->value() != 0.0 )
        sv_lhs->add_element( ele_t( itr->index() + sv_rhs.offset() + offset, alpha * (itr->value()) ) );
  }
  sv_lhs->assume_sorted(assume_sorted);
}
