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

#include <iomanip>

#include "DenseLinAlgPack_DVectorClassTmpl.hpp"
#include "Teuchos_Assert.hpp"

#ifdef LINALGPACK_CHECK_SLICE_SETUP
DenseLinAlgPack::size_type DenseLinAlgPack::vector_validate_sized(size_type size)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !size, std::invalid_argument
    ,"vector_validate_sized(...) : Error, A vector region can not be created from an unsized vector.\n"
    );
  return size;
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void DenseLinAlgPack::vector_validate_range(size_type ubound, size_type max_ubound)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ubound > max_ubound, std::out_of_range
    ,"vector_validate_range(...) : The upper bound is out of range.\n");
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void DenseLinAlgPack::vector_validate_subscript(size_type size, size_type i)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    i < 1 || i > size, std::out_of_range
    ,"vector_validate_subscript(size,i) : Error, Subscript i out of bounds.\n");
}
#endif

#ifdef LINALGPACK_CHECK_RHS_SIZES
void DenseLinAlgPack::assert_vs_sizes(size_type size1, size_type size2)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    size1 != size2, std::length_error
    ,"assert_vs_sizes(...) : Error, size1 = " << size1 << " != size2 = " << size2 << "\n");
}
#endif
