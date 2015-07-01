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

#include "AbstractLinAlgPack_SpVectorView.hpp"

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
   % (int)sizeof(AbstractLinAlgPack::value_type))
  , double
  >
validate_value_stride;

// Validate that there is an integer stride between indexes
assert_compile_time<
    ((int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
   % (int)sizeof(AbstractLinAlgPack::index_type))
  , double
  >
validate_index_stride;

// Compute the stride between values
const int values_stride = (int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
  / (int)sizeof(AbstractLinAlgPack::value_type);

// Compute the stride between indexes
const int indices_stride = (int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
  / (int)sizeof(AbstractLinAlgPack::index_type);

} // end namespace

RTOpPack::SparseSubVector
AbstractLinAlgPack::sub_vec_view(
  const SpVectorSlice&   sv_in
  ,const Range1D&        rng_in
  )
{
  using Teuchos::null;
  const Range1D        rng = RangePack::full_range(rng_in,1,sv_in.dim());
  const SpVectorSlice  sv = sv_in(rng);

  RTOpPack::SparseSubVector  sub_vec;

  if(!sv.nz()) {
    sub_vec.initialize(
      rng.lbound()-1  // global_offset
      ,rng.size()     // sub_dim
      ,0              // nz
      ,null           // vlaues
      ,1              // values_stride
      ,null           // indices
      ,1              // indices_stride
      ,0              // local_offset
      ,1              // is_sorted
      );
  }
  else {
    SpVectorSlice::const_iterator itr = sv.begin();
    TEUCHOS_TEST_FOR_EXCEPT( !( itr != sv.end() ) );
    if( sv.dim() && sv.nz() == sv.dim() && sv.is_sorted() ) {
      const Teuchos::ArrayRCP<const value_type>  values =
        Teuchos::arcp(&itr->value(), 0, values_stride*rng.size(), false) ;
      sub_vec.initialize(
        rng.lbound()-1    // global_offset
        ,rng.size()       // sub_dim
        ,values           // values
        ,values_stride    // values_stride
        );
    }
    else {
      const Teuchos::ArrayRCP<const value_type>  values =
        Teuchos::arcp(&itr->value(), 0, values_stride*sv.nz(), false) ;
      const Teuchos::ArrayRCP<const index_type> indexes =
        Teuchos::arcp(&itr->index(), 0, indices_stride*sv.nz(), false);
      sub_vec.initialize(
        rng.lbound()-1    // global_offset
        ,sv.dim()         // sub_dim
        ,sv.nz()          // sub_nz
        ,values           // values
        ,values_stride    // values_stride
        ,indexes          // indices
        ,indices_stride   // indices_stride
        ,sv.offset()      // local_offset
        ,sv.is_sorted()   // is_sorted
        );
    }
  }
  
  return sub_vec;
}
