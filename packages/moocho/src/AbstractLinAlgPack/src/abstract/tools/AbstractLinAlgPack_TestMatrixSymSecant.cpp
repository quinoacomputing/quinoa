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
#include <iomanip>

#include <math.h>

#include "AbstractLinAlgPack_TestMatrixSymSecant.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_MatrixNonsing.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"

bool AbstractLinAlgPack::TestMatrixSymSecant(
  const MatrixOp        &B
  ,const Vector       &s
  ,const Vector       &y
  ,value_type               warning_tol
  ,value_type               error_tol
  ,bool                     print_all_warnings
  ,std::ostream             *out
  ,bool                     trase
  )
{
  using std::setw;
  using std::endl;
  using AbstractLinAlgPack::sum;
  using AbstractLinAlgPack::V_InvMtV;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using LinAlgOpPack::V_MtV;

  bool success = true;
  const char
    sep_line[] = "\n-----------------------------------------------------------------\n";

  // Check the secant property (B * s = y)
  {
    if( out && trase )
      *out
        << sep_line
        << "\n|(sum(B*s)-sum(y))/||y||inf| = ";
    VectorSpace::vec_mut_ptr_t
      Bs = y.space().create_member();
    V_MtV( Bs.get(), B, BLAS_Cpp::no_trans, s );
    const value_type
      sum_Bs = sum(*Bs),
      sum_y  = sum(y),
      nrm_y  = y.norm_inf(),
      err    = ::fabs( (sum_Bs - sum_y) / nrm_y );
    if( out && trase )
      *out
        <<"|("<<sum_Bs<<"-"<<sum_y<<")/"<<nrm_y<<"| = " << err << std::endl;
    if( err >= error_tol ) {
      if( out && trase )
        *out
          << "Error, above error = " << err << " >= error_tol = " << error_tol
          << "\nThe test has failed!\n";
      if(out && print_all_warnings) {
        *out
          << "\ns =\n" << s 
          << "\ny =\n" << y 
          << "\nBs =\n" << *Bs
          << endl;
      }
      success = false;
    }
    else if( err >= warning_tol ) {
      if( out && trase )
        *out
          << "Warning!, above error = " << err << " >= warning_tol = " << warning_tol << std::endl;
    }
  }
  // Check the secant property (s = inv(B)*y)
  const MatrixNonsing
    *B_nonsing = dynamic_cast<const MatrixNonsing*>(&B);
  if( B_nonsing ) {
    if( out && trase )
      *out
        << sep_line
        << "\n|(sum(inv(B)*y)-sum(s))/||s||inf| = ";
    VectorSpace::vec_mut_ptr_t
      InvBy = s.space().create_member();
    V_InvMtV( InvBy.get(), *B_nonsing, BLAS_Cpp::no_trans, y );
    const value_type
      sum_InvBy = sum(*InvBy),
      sum_s     = sum(s),
      nrm_s     = s.norm_inf(),
      err       = ::fabs( (sum_InvBy - sum_s) / nrm_s );
    if( out && trase )
      *out
        <<"|("<<sum_InvBy<<"-"<<sum_s<<")/"<<nrm_s<<"| = " << err << std::endl;
    if( err >= error_tol ) {
      if( out && trase )
        *out
          << "Error, above error = " << err << " >= error_tol = " << error_tol
          << "\nThe test has failed!\n";
      if(out && print_all_warnings) {
        *out
          << "\ns =\n" << s 
          << "\ny =\n" << y 
          << "\ninv(B)*y =\n" << *InvBy
          << endl;
      }
      success = false;
    }
    else if( err >= warning_tol ) {
      if( out && trase )
        *out
          << "Warning!, above error = " << err << " >= warning_tol = " << warning_tol << std::endl;
    }
  }
  if( out && trase )
    *out << sep_line;
  return success;
}
