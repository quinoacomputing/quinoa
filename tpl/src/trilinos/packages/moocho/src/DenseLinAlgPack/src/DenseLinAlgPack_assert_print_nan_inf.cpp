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

#include <ostream>
#include <sstream>
#include <iomanip>

#include "DenseLinAlgPack_assert_print_nan_inf.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "check_nan_inf.h"

bool DenseLinAlgPack::assert_print_nan_inf( const value_type& val
  , const std::string & name, bool throw_excpt, std::ostream* out )
{

  if( RTOp_is_nan_inf(val) ) {
    std::ostringstream omsg;
    omsg
      << "The scalar \"" << name
      << "\" = " << val << " is not a valid bounded number";
    if(out)
      *out << omsg.str() << std::endl;
    if( throw_excpt ) {
      if(out)
        out->flush();
      throw NaNInfException( "assert_print_nan_inf(...) : Error, "
        + omsg.str() );
    }
    return false;
  }
  return true;
}

bool DenseLinAlgPack::assert_print_nan_inf( const DVectorSlice& v
  , const std::string & name, bool throw_excpt, std::ostream* out )
{

  bool has_nan_or_inf = false;
  bool printed_header = false;

  for( DVectorSlice::const_iterator v_itr = v.begin(); v_itr != v.end(); ++v_itr ) {
    if( RTOp_is_nan_inf(*v_itr) ) {
      if(out) {
        if(!printed_header) {
          *out
            << "The vector \"" << name
            << "\" has the following NaN or Inf entries\n";
          printed_header = true;
        }
        *out
          << name << "(" << v_itr - v.begin() + 1 << ") = "
          << *v_itr << std::endl;
      }
      has_nan_or_inf = true;
    }
  }
  if( has_nan_or_inf && throw_excpt ) {
    if(out)
      out->flush();
    std::ostringstream omsg;
    omsg
      << "assert_print_nan_inf(...) : Error, the vector named "
      << name << " has at least one element which is NaN or Inf";
    throw NaNInfException( omsg.str() );
  }

  return !has_nan_or_inf;
}

bool DenseLinAlgPack::assert_print_nan_inf( const DMatrixSlice& m
  , const std::string & name, bool throw_excpt, std::ostream* out )
{

  bool has_nan_or_inf = false;
  bool printed_header = false;

  for( size_type j = 1; j <= m.cols(); ++j ) {
    const DVectorSlice& v = m.col(j);
    for( DVectorSlice::const_iterator v_itr = v.begin(); v_itr != v.end(); ++v_itr ) {
      if( RTOp_is_nan_inf(*v_itr) ) {
        if(out) {
          if(!printed_header) {
            *out
              << "The matrix \"" << name
              << "\" has the following NaN or Inf entries\n";
            printed_header = true;
          }
          *out
            << name << "(" << v_itr - v.begin() + 1 << "," << j << ") = "
            << *v_itr << std::endl;
        }
        has_nan_or_inf = true;
      }
    }
  }

  if( has_nan_or_inf && throw_excpt ) {
    if(out)
      out->flush();
    std::ostringstream omsg;
    omsg
      << "assert_print_nan_inf(...) : Error, the matrix named "
      << name << " has at least one element which is NaN or Inf";
    throw NaNInfException( omsg.str() );
  }

  return has_nan_or_inf;
}
