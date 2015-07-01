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
#include <ostream>
#include <vector>
#include <typeinfo>

#include <math.h>

#include "DenseLinAlgPack_TestDenseLinAlgPack.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "DenseLinAlgPack_MatVecCompare.hpp"

namespace {

// 2/10/00: See TestVectorClass.cpp
inline
bool update_success( bool result, bool *success ) {
  return TestingHelperPack::update_success( result, success );
}

}	// end namespace

bool DenseLinAlgPack::TestingPack::TestVectorOp(std::ostream* out)
{

  using DenseLinAlgPack::comp;
  using DenseLinAlgPack::sqrt_eps;

  bool success = true;
  bool result, result1, result2;

  value_type
    rval;

  if(out)
    *out	<< "\n***********************************"
        << "\n*** Testing VectorOp Operations ***"
        << "\n***********************************\n"
        << std::boolalpha;

  try {

  if(out) *out << "\nLet alpha1 = 2.0, alpha2 = 0.5, v1val = 1.0, v2val = 2.0;\n";
  const value_type
    alpha1 = 2.0,
    alpha2 = 0.5,
    v1val = 1.0,
    v2val = 2.0;

  const int n = 6;

  if(out) *out << "\nVector v(n), v1(v1val,n), v2(v2val,n);\n";
  DVector v(n), v1(v1val,n), v2(v2val,n);

  // /////////////////////////////////
  // Test Algebraic Functions

  if(out) *out << "\n***\n*** Testing DVectorSlice modifying functions \n***\n";

  if(out) *out << "\nv = 1.0; Vp_S(&v()," << alpha1 << ");\n";
  v = 1.0;
  Vp_S(&v(),alpha1);
  result = update_success( comp( v, 1.0 + alpha1 ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == 1.0 + " << alpha1 << " : " << result << std::endl;

  if(out) *out << "\nv = 1.0; Vt_S(&v()," << alpha1 << ");\n";
  v = 1.0;
  Vt_S(&v(),alpha1);
  result = update_success( comp( v, alpha1 ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == " << alpha1 << " : " << result << std::endl;

  if(out) *out << "\nv = 1.0; Vp_StV(&v()," << alpha1 << ",v1());\n";
  v = 1.0;
  Vp_StV(&v(),alpha1,v1());
  result = update_success( comp( v, 1.0 + alpha1 ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<v1val<<" + "<<alpha1<<" : " << result << std::endl;

  // DVectorSlice as lhs

  if(out) *out << "\n***\n*** Testing DVectorSlice as lhs algebric operations\n***\n";

  if(out) *out << "\nv = -10.0; V_VpV(&v(),v1(),v2());\n";
  v = -10;
  V_VpV(&v(),v1(),v2());
  result = update_success( comp( v, v1val+v2val ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<v1val<<" + "<<v2val<<" : " << result << std::endl;

  if(out) *out << "\nv = -10.0; V_VmV(&v(),v1(),v2());\n";
  v = -10;
  V_VmV(&v(),v1(),v2());
  result = update_success( comp( v, v1val-v2val ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<v1val<<" - "<<v2val<<" : " << result << std::endl;

  if(out) *out << "\nv = -10.0; V_mV(&v(),v1());\n";
  v = -10;
  V_mV(&v(),v1());
  result = update_success( comp( v, -v1val ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<-v1val<<" : " << result << std::endl;

  if(out) *out << "\nv = -10.0; V_StV(&v(),"<<alpha2<<",v2());\n";
  v = -10;
  V_StV(&v(),alpha2,v2());
  result = update_success( comp( v, alpha2*v2val ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<alpha2<<" * "<<v2val<<" : " << result << std::endl;

  // DVector as lhs

  if(out) *out << "\n***\n*** Testing DVector as lhs algebric operations\n***\n";

  if(out) *out << "\nv.free(); V_VpV(&v,v1(),v2());\n";
  v.free();
  V_VpV(&v,v1(),v2());
  result = update_success( comp( v, v1val+v2val ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<v1val<<" + "<<v2val<<" : " << result << std::endl;

  if(out) *out << "\nv.free(); V_VmV(&v,v1(),v2());\n";
  v.free();
  V_VmV(&v,v1(),v2());
  result = update_success( comp( v, v1val-v2val ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<v1val<<" - "<<v2val<<" : " << result << std::endl;

  if(out) *out << "\nv.free(); V_mV(&v,v1());\n";
  v.free();
  V_mV(&v,v1());
  result = update_success( comp( v, -v1val ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<-v1val<<" : " << result << std::endl;

  if(out) *out << "\nv.free(); V_StV(&v,"<<alpha2<<",v2());\n";
  v.free();
  V_StV(&v,alpha2,v2());
  result = update_success( comp( v, alpha2*v2val ), &success );
  if(out)
    *out	<< "v =\n" << v
        << "v == "<<alpha2<<" * "<<v2val<<" : " << result << std::endl;

  // ////////////////////////////////
  // Test Elementary math functions
  if(out) *out << "\n***\n*** Testing Elementary math functions\n***\n";
  
  // ToDo: implement at some point
  if(out) *out << "\nWarning! Not Tested!\n";

  // DVectorSlice as lhs

  // DVector as lhs

  // /////////////////////////////////
  // Test Scalar Returning Functions

  if(out) *out << "\n***\n*** Testing Scalar Returning Functions\n***\n";

  if(out) *out << "\n(dot(v1(),v2()) -> ";
  rval = dot(v1(),v2());
  result = update_success( ::fabs( rval - v1val*v2val*n ) < sqrt_eps, &success );
  if(out)	*out << rval <<") == " << v1val*v2val*n << " : " << result << std::endl;

  if(out) *out << "\n(norm_1(v2()) -> ";
  rval = norm_1(v2());
  result = update_success( ::fabs( rval - v2val*n ) < sqrt_eps, &success );
  if(out)	*out << rval <<") == " << v2val*n << " : " << result << std::endl;

  if(out) *out << "\n(norm_2(v2()) -> ";
  rval = norm_2(v2());
  result = update_success( ::fabs( rval - v2val*::sqrt((value_type)n) ) < sqrt_eps, &success );
  if(out)	*out << rval <<") == " << v2val*::sqrt((value_type)n) << " : " << result << std::endl;

  if(out) *out << "\n(norm_inf(v2()) -> ";
  rval = norm_inf(v2());
  result = update_success( ::fabs( rval - v2val ) < sqrt_eps, &success );
  if(out)	*out << rval <<") == " << v2val << " : " << result << std::endl;

  if(out) *out << "\nv1(n/2) = 2*v1val;\n(max(v1()) -> ";
  v1(n/2) = 2*v1val;
  rval = max(v1());
  result = update_success( ::fabs( rval - 2*v1val ) < sqrt_eps, &success );
  if(out)	*out << rval <<") == " << 2*v1val << " : " << result << std::endl;

  if(out) *out << "\nv1(n/2+1) = -2*v1val;\n(min(v1()) -> ";
  v1(n/2+1) = -2*v1val;
  rval = min(v1());
  result = update_success( ::fabs( rval - (-2*v1val) ) < sqrt_eps, &success );
  if(out)	*out << rval <<") == " << (-2*v1val) << " : " << result << std::endl;

  } // end try
  catch( const std::exception& excpt ) {
    success = false;
    if(out)
      (*out)	<< "\nError, a standard exception was thrown: "
          << typeName(excpt) << ": "
          << excpt.what() << std::endl; 
  }
  catch(...) {
    success = false;
    if(out)
      (*out)	<< "\nError, an unknown exception was thrown\n";
  }

  if(out) {
    if(success)
      (*out)
        << "\n*** Congradulations, VectorOp operations seem to check out. ***\n";
    else
      (*out)
        << "\n*** Oops, all of the tests for VectorOp operations "
          "where not successful. ***\n";
  }


  return success;
}

