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

#include "ConstrainedOptPack_VariableBoundsTester.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"

namespace ConstrainedOptPack {

// public

VariableBoundsTester::VariableBoundsTester(
  value_type   warning_tol
  ,value_type  error_tol
  )
  :warning_tol_(warning_tol)
  ,error_tol_(error_tol)
{}

bool VariableBoundsTester::check_in_bounds(
  std::ostream* out, bool print_all_warnings, bool print_vectors
  ,const Vector& xL, const char xL_name[]
  ,const Vector& xU, const char xU_name[]
  ,const Vector& x,  const char x_name[]
  )
{
  using AbstractLinAlgPack::max_near_feas_step;

  if(out)
    *out
      << "\n*** Checking that variables are in bounds\n";

  VectorSpace::vec_mut_ptr_t zero = x.space().create_member(0.0);
  std::pair<value_type,value_type>
    u = max_near_feas_step( x, *zero, xL, xU, warning_tol() );
  if(u.first < 0.0) {
    if(out)
      *out << "\nWarning! the variables " << xL_name << " <= " << x_name << " <= " << xU_name
        << " are out of bounds by more than warning_tol = "	<< warning_tol() << "\n";
    u = max_near_feas_step( x, *zero, xL, xU, error_tol() );
    if(u.first < 0.0) {
      if(out)
        *out << "\nError! the variables " << xL_name << " <= " << x_name << " <= " << xU_name
          << " are out of bounds by more than error_tol = " << error_tol() << "\n";
      return false;
    }
  }
  return true;
}

}	// end namespace ConstrainedOptPack
