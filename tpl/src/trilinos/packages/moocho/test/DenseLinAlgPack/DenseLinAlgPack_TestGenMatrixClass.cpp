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

#include "DenseLinAlgPack_TestDenseLinAlgPack.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_MatVecCompare.hpp"

namespace {

using DenseLinAlgPack::sqrt_eps;

// Check consistency of row(), col(), diag() and operator()().
template<class M_t>
void check_access( M_t& M, typename M_t::size_type row_offset, typename M_t::size_type col_offset
  , std::ostream* out, bool* success )
{
  if(out)
    *out	<< "Checking M(i,j) == M.row(i)(j) == M.col(j)(i) == "
        << "M.diag(...)(...) == "
        << "(i + "<<row_offset<<") + 0.1*(j+"<<col_offset<<") : ";

  bool result = true;

  for( typename M_t::size_type i = 1; i <= M.rows(); ++i ) {
    for( typename M_t::size_type j = 1; j <= M.rows(); ++j ) {
      const typename M_t::value_type
        Mij = M(i,j);
      typename M_t::value_type
        val = (i+row_offset)+0.1*(j+col_offset);
      if( ::fabs(Mij-val) > sqrt_eps ) {
        result = false;
        if(out) *out << "(M("<<i<<","<<j<<") -> "<<Mij<<") != "<<val<<std::endl;
      }
      if( Mij != (val = M.row(i)(j)) ) {
        result = false;
        if(out) *out << "M("<<i<<","<<j<<") != (M.row("<<i<<")("<<j<<") -> "<<val<<")\n";
      }
      if( Mij != (val = M.col(j)(i)) ) {
        result = false;
        if(out) *out << "M("<<i<<","<<j<<") != (M.col("<<j<<")("<<i<<") -> "<<val<<")\n";
      }
      const int k = ( i > j ? -i + j : j - i );
      const typename M_t::size_type k_i = ( i > j ? j : i );
      if( Mij != (val = M.diag(k)(k_i) ) ) {
        result = false;
        if(out) *out << "M("<<i<<","<<j<<") != (M.diag("<<k<<")("<<k_i<<") -> "<<val<<")\n";
      }
    }
  }
  if(out) *out << result << std::endl;
  if(!result) *success = false;
}

// Print out a string for overlap
const char* overlap_str( DenseLinAlgPack::EOverLap overlap ) {
  switch(overlap) {
    case DenseLinAlgPack::NO_OVERLAP:
      return "NO_OVERLAP";
    case DenseLinAlgPack::SOME_OVERLAP:
      return "SOME_OVERLAP";
    case DenseLinAlgPack::SAME_MEM:
      return "SAME_MEM";
  }
  return "Invalid value for EOverLap";
}

}	// end namespace

bool DenseLinAlgPack::TestingPack::TestGenMatrixClass(std::ostream* out)
{

  using DenseLinAlgPack::comp;
  using DenseLinAlgPack::sqrt_eps;

  bool success = true;
  bool result;

  if(out)
    *out	<< "\n****************************************************"
        << "\n*** Testing DMatrix and DMatrixSlice classes ***"
        << "\n****************************************************\n"
        << std::boolalpha;

  try {

  const size_type
    m = 6,
    n = 8;

  const value_type
    ptr[m*n] =
       {	1.1,	2.1,	3.1,	4.1,	5.1,	6.1,
        1.2,	2.2,	3.2,	4.2,	5.2,	6.2,
        1.3,	2.3,	3.3,	4.3,	5.3,	6.3,
        1.4,	2.4,	3.4,	4.4,	5.4,	6.4,
        1.5,	2.5,	3.5,	4.5,	5.5,	6.5,
        1.6,	2.6,	3.6,	4.6,	5.6,	6.6,
        1.7,	2.7,	3.7,	4.7,	5.7,	6.7,
        1.8,	2.8,	3.8,	4.8,	5.8,	6.8	};

  // /////////////////////////////
  // Test Constructors

  if(out)
    *out	<< "\n***\n*** Testing constructors\n***\n";

  // DMatrixSlice
  if(out) *out << "\nGenMatrixSlice gms1;\n";
  DMatrixSlice gms1;
  if(out) *out << "gms1 =\n" << gms1;
  update_success( result = (gms1.rows() == 0 && gms1.cols() == 0 ), &success );
  if(out)
    *out	<< "((gms1.rows() -> "<<gms1.rows()
        << ") == 0 && (gms1.cols() -> "<<gms1.cols()<<") == 0 ) : "
        << result << std::endl;


  if(out) *out << "\nGenMatrixSlice gms2( const_cast<value_type*>(ptr), m*n, m, m, n );\n";
  const DMatrixSlice gms2( const_cast<value_type*>(ptr), m*n, m, m, n );
  if(out) *out << "gms2 =\n" << gms2;

  if(out) *out << "\nGenMatrixSlice gms3( const_cast<DMatrixSlice&>(gms2), Range1D(1,m), Range1D(1,n) );\n";
  const DMatrixSlice gms3( const_cast<DMatrixSlice&>(gms2), Range1D(1,m), Range1D(1,n) );
  if(out) *out << "gms3 =\n" << gms3;

  // DMatrix

  if(out) *out << "\nGenMatrix gm1;\n";
  DMatrix gm1;	
  if(out) *out << "gm1 =\n" << gm1();
  update_success( result = (gm1.rows() == 0 && gm1.cols() == 0 ), &success );
  if(out)
    *out	<< "((gm1.rows() -> "<<gm1.rows()
        << ") == 0 && (gm1.cols() -> "<<gm1.cols()<<") == 0 ) : "
        << result << std::endl;
  
  if(out) *out << "\nGenMatrix gm2(m,n);\n";
  DMatrix gm2(m,n);
  if(out) *out << "gm2 =\n" << gm2();

  if(out) *out << "\nGenMatrix gm3(1.0,m,n);\n";
  DMatrix gm3(1.0,m,n);
  if(out) *out << "gm3 =\n" << gm3();
  update_success( result = comp( gm3(), 1.0 ), &success );
  if(out) *out << "gm3 == 1.0 : " << result << std::endl;

  if(out) *out << "\nGenMatrix gm4(ptr,m,n);\n";
  DMatrix gm4(ptr,m,n);
  if(out) *out << "gm4 =\n" << gm4();

  if(out) *out << "\nGenMatrix gm5(gms2);\n";
  DMatrix gm5(gms2);
  if(out) *out << "gm5 =\n" << gm5();

  // ////////////////////////////
  // Test DMatrixSlice binding

  if(out)
    *out	<< "\n***\n*** Testing DMatrixSlice binding\n***\n";

  if(out) *out << "\ngms1.bind(gm4());\n";
  gms1.bind(gm4());
  if(out) *out << "gms1 =\n" << gms1();

  // ////////////////////////////
  // Test DMatrix resizing

  if(out)
    *out	<< "\n***\n*** Testing DMatrix resizing\n***\n";

  if(out) *out << "\ngm1.resize(m,n,1.0);\n";
  gm1.resize(m,n,1.0);
  if(out) *out << "gm1 =\n" << gm1();
  update_success( result = comp( gm1(), 1.0 ), &success );
  if(out) *out << "gm1 == 1.0 : " << result << std::endl;

  // ///////////////////////////////////////////////
  // Test row, col, diag access and element access

  // DMatrixSlice

  if(out)
    *out	<< "\n***\n*** Testing row, col, diag access and element access\n***\n";

  if(out) *out << "\nLet M = gms1\n";
  check_access( gms1, 0, 0, out, &success );

  if(out) *out << "\nLet M = const_cast<const DMatrixSlice&>(gms1)\n";
  check_access( const_cast<const DMatrixSlice&>(gms1), 0, 0, out, &success );

  // DMatrix

  if(out) *out << "\nLet M = gm4\n";
  check_access( gm4, 0, 0, out, &success );

  if(out) *out << "\nLet M = const_cast<const DMatrix&>(gm4)\n";
  check_access( const_cast<const DMatrix&>(gm4), 0, 0, out, &success );

  // ////////////////////////////
  // Test submatrix access

  if(out)
    *out	<< "\n***\n*** Testing submatrix access\n***\n";

  if(out) *out << "\nRange1D r_rng(2,m-1), c_rng(2,n-1);\n";
  Range1D r_rng(2,m-1), c_rng(2,n-1);

  // DMatrixSlice

  if(out) *out << "\nLet M = const_cast<DMatrixSlice&>(gms2)(r_rng,c_rng)\n";
  gms1.bind( const_cast<DMatrixSlice&>(gms2)(r_rng,c_rng) );
  if(out) *out << "M =\n" << gms1;
  check_access( gms1, 1, 1, out, &success );

  if(out) *out << "\nLet M = const_cast<DMatrixSlice&>(gms2(r_rng,c_rng))\n";
  gms1.bind( const_cast<DMatrixSlice&>(gms2)(r_rng,c_rng) );
  if(out) *out << "M =\n" << gms1;
  check_access( gms1, 1, 1, out, &success );

  if(out) *out << "\nLet M = const_cast<DMatrixSlice&>(gms2)(2,m-1,2,n-1)\n";
  gms1.bind(const_cast<DMatrixSlice&>(gms2)(2,m-1,2,n-1) );
  if(out) *out << "M =\n" << gms1;
  check_access( gms1, 1, 1, out, &success );

  if(out) *out << "\nLet M = const_cast<DMatrixSlice&>(gms2(2,m-1,2,n-1))\n";
  gms1.bind( const_cast<DMatrixSlice&>(gms2)(2,m-1,2,n-1) );
  if(out) *out << "M =\n" << gms1;
  check_access( gms1, 1, 1, out, &success );

  // DMatrix

  if(out) *out << "\nLet M = gm4(r_rng,c_rng)\n";
  gms1.bind( gm4(r_rng,c_rng) );
  if(out) *out << "M =\n" << gms1;
  check_access( gms1, 1, 1, out, &success );

  if(out) *out << "\nLet M = const_cast<const DMatrixSlice&>(gm4)(r_rng,c_rng)\n";
  gms1.bind( const_cast<const DMatrix&>(gm4)(r_rng,c_rng) );
  if(out) *out << "M =\n" << gms1;
  check_access( gms1, 1, 1, out, &success );

  if(out) *out << "\nLet M = gm4(2,m-1,2,n-1)\n";
  gms1.bind( gm4(2,m-1,2,n-1) );
  if(out) *out << "M =\n" << gms1;
  check_access( gms1, 1, 1, out, &success );

  if(out) *out << "\nLet M = const_cast<const DMatrixSlice&>(gm4)(2,m-1,2,n-1)\n";
  gms1.bind( const_cast<const DMatrix&>(gm4)(2,m-1,2,n-1) );
  if(out) *out << "M =\n" << gms1;
  check_access( gms1, 1, 1, out, &success );

  // ////////////////////
  // Test matrix overlap

  if(out)
    *out	<< "\n***\n*** matrix overlap\n***\n";

  EOverLap ovlap;

  // DMatrixSlice

  if(out) *out << "(gms2.overlap(gms2) -> ";
  ovlap = gms2.overlap(gms2);
  result = update_success( ovlap == SAME_MEM, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SAME_MEM : " << result << std::endl;

  if(out) *out << "(gms2.overlap(gms2(r_rng,c_rng)) -> ";
  ovlap = gms2.overlap(gms2(r_rng,c_rng));
  result = update_success( ovlap == SOME_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

  if(out) *out << "(gms2(1,m/2,1,n/2).overlap(gms2(m/2,m,n/2,n)) -> ";
  ovlap = gms2(1,m/2,1,n/2).overlap(gms2(m/2,m,n/2,n));
  result = update_success( ovlap == SOME_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

  if(out) *out << "(gms2(1,m/2,1,n/2).overlap(gms2(m/2+1,m,n/2+1,n)) -> ";
  ovlap = gms2(1,m/2,1,n/2).overlap(gms2(m/2+1,m,n/2+1,n));
  result = update_success( ovlap == NO_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == NO_OVERLAP : " << result << std::endl;

  // DMatrix

  if(out) *out << "(gm4.overlap(gm4) -> ";
  ovlap = gm4.overlap(gm4());
  result = update_success( ovlap == SAME_MEM, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SAME_MEM : " << result << std::endl;

  if(out) *out << "(gm4.overlap(gm4(r_rng,c_rng)) -> ";
  ovlap = gm4.overlap(gm4(r_rng,c_rng));
  result = update_success( ovlap == SOME_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

  // //////////////////////////////////////////////////////
  // Test vector overlap (continuation of vector testing)
  //
  // ToDo: Finish this someday once you get it figured out.

  // ///////////////////////////
  // Test assignment operators

  if(out)
    *out	<< "\n***\n*** assignment operators\n***\n";

  // DMatrixSlice

  if(out) *out << "\ngms1.bind(gm1());\n";
  gms1.bind(gm1());

  if(out) *out << "\ngms1 = 2.0;\n";
  gms1 = 2.0;
  if(out) *out << "gms1 =\n" << gms1;
  update_success( result = comp(gms1,2.0), &success );
  if(out) *out << "gms1 == 2.0 : " << result << std::endl; 

  if(out) *out << "\ngms1 = gms2;\n";
  gms1 = gms2;
  if(out) *out << "gms1 =\n" << gms1;
  update_success( result = comp(gms1,gms2), &success );
  if(out) *out << "gms1 == gms2 : " << result << std::endl; 

  // DMatrix

  if(out) *out << "\ngm1 = 3.0;\n";
  gm1 = 3.0;
  if(out) *out << "gm1 =\n" << gm1();
  update_success( result = comp(gm1(),3.0), &success );
  if(out) *out << "gm1 == 3.0 : " << result << std::endl; 

  if(out) *out << "\ngm1.resize(0,0); gm1 = gms2;\n";
  gm1.resize(0,0);
  gm1 = gms2;
  if(out) *out << "gm1 =\n" << gm1();
  update_success( result = comp(gm1(),gms2()), &success );
  if(out) *out << "gm1 == gms2 : " << result << std::endl; 

  if(out) *out << "\ngm1.resize(0,0); gm1 = gm4;\n";
  gm1.resize(0,0);
  gm1 = gm4;
  if(out) *out << "gm1 =\n" << gm1();
  update_success( result = comp(gm1(),gm4()), &success );
  if(out) *out << "gm1 == gm4 : " << result << std::endl; 

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
        << "\n*** Congradulations, DMatrix and DMatrixSlice seem to check out. ***\n";
    else
      (*out)
        << "\n*** Oops, all of the tests for DMatrix and DMatrixSlice "
          "where not successful. ***\n";
  }

  return success;
}

