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
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "DenseLinAlgPack_MatVecCompare.hpp"
#include "TestingHelperPack_update_success.hpp"

namespace {

// This template function checks that iterator and subscriping access all
// give the same results.
template<class V, class I>
void test_access( V* _v, I begin, I end, std::ostream*out, bool* success ) {
  using std::setw;
  using TestingHelperPack::update_success;

  V &v = *_v;
  bool result;

  if(out)
    *out	<< "\nbegin + i == v[i] == v(i+1), for i = 0,1,...,v.dim()\n";

  I itr;
  int i;
  if(out)
    *out	<< "\n i,  *(begin + i)  ==  v[i]  ==  v(i+1)  ==    ?  "
        << "\n--,  ------------     ------     ------     ------\n";
  for( itr = begin, i = 0; itr != end; ++itr, ++i ) {
    result = update_success( *itr == v[i] && v[i] == v(i+1), success );
    if(out)
      *out	<< setw(2)	<< i << ','
          << setw(14) << *itr
          << setw(11) << v[i]
          << setw(11) << v(i+1)
          << setw(11) << std::right << result << std::endl << std::left;
  }
  if(out) *out << std::endl;
}

// This template function checks that a subregion creates the expected view.
// Here rng must be rng.full_range() == true.
template<class V, class VS>
void test_subregion_access( V* _v, VS* _vs, const DenseLinAlgPack::Range1D& rng
  , std::ostream* out, bool* success )
{
  using std::setw;
  using TestingHelperPack::update_success;

  bool result;
  V &v = *_v;
  VS &vs = *_vs;

  if(out)
    *out	<< "\nv.begin() + i1 == vs.begin() + i2 == vs[i2] == vs(i2+1)"
        << ", for i1 = lb-1,..,ub-1, for i2 = 0,..,rng.size()-1\n";

  typename V::const_iterator itr1;
  typename VS::const_iterator itr2;
  int i1, i2;
  if(out)
    *out	<< "\ni1, i2, *(v.begin() + i1)  ==  *(vs.begin() + i2)  ==  vs[i2]  ==  vs(i2+1)  ==    ?  "
        << "\n--, --, -----------------      ------------------      ------      --------      ------\n";
  for(	i1 = rng.lbound()-1, itr1 = v.begin() + i1, i2 = 0, itr2 = vs.begin();
      i2 < rng.size();
      ++i1, ++itr1, ++i2, ++itr2												)
  {
    result = update_success( *itr1 == *itr2 && *itr2 == vs[i2] && *itr2 == vs(i2+1), success );
    if(out)
      *out	<< setw(2)	<< i1 << ','
          << setw(3)	<< i2 << ','
          << setw(18) << *itr1
          << setw(24) << *itr2
          << setw(12) << vs[i2]
          << setw(14) << vs(i2+1)
          << setw(12) << std::right << result << std::endl << std::left;
  }
  if(out) *out << std::endl;
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

bool DenseLinAlgPack::TestingPack::TestVectorClass(std::ostream* out)
{
  using DenseLinAlgPack::comp;
  using DenseLinAlgPack::sqrt_eps;
  using TestingHelperPack::update_success;

  bool success = true;
  bool result, result1, result2;

  if(out)
    *out	<< "\n**********************************************"
        << "\n*** Testing DVector and DVectorSlice classes ***"
        << "\n**********************************************\n"
        << std::boolalpha;

  try {

  if(out)
    *out	<< "\nLet vvz[i-1] = i + 0.1*i, for i = 1...,n\n";

  // std::vector<> which is starndard for comparisons
  const DVector::size_type n = 6;
  std::vector<DVector::value_type> vvz(6);
  {for(int i = 1; i <= n; ++i)
    vvz[i-1] = i + 0.1 * i;
  }

  if(out) *out << "\nLet alpha1 = 22.5\n";
  const value_type alpha1 = 22.5; 

  // ///////////////////////
  // Test Constructors

  if(out)
    *out	<< "\n***\n*** Testing constructors\n***\n";

  // DVectorSlice Constructors

  if(out) *out << "\nVectorSlice vs1\n";
  DVectorSlice	vs1;
  if(out) *out << "vs1 =\n" << vs1;

  if(out) *out << "\nVectorSlice vs2(vvz.begin(),n)\n";
  DVectorSlice	vs2(&vvz[0],n);
  if(out) *out << "vs2 =\n" << vs2;

  if(out) *out << "\nVectorSlice vs3(vvz.begin(),n,Range1D())\n";
  DVectorSlice	vs3(&vvz[0],n,Range1D());
  if(out) *out << "vs3 =\n" << vs3;

  if(out) *out << "\nVectorSlice vs4(vs3,Range1D())\n";
  DVectorSlice	vs4(vs3,Range1D());
  if(out) *out << "vs4 =\n" << vs4;

  // DVector Constructors

  if(out) *out << "\nVector v1\n";
  DVector v1;
  if(out) *out << "v1 =\n" << v1;

  if(out) *out << "\nVector v2(alpha1,n)\n";
  DVector v2(alpha1,n);
  if(out) *out << "v2 =\n" << v2;

  if(out) *out << "\nVector v3(vvz.begin(),n)\n";
  DVector v3(&vvz[0],n);
  if(out) *out << "v3 =\n" << v3;

  if(out) *out << "\nVector v4(vs4)\n";
  DVector v4(vs4);
  if(out) *out << "v4 =\n" << v4;

  // //////////////////////////////////////
  // Test Binding Views (DVectorSlice)

  if(out)
    *out	<< "\n***\n*** Testing DVectorSlice binding and "
          "conversion from DVector -> DVectorSlice\n***\n"
        << "\nvs1.bind(v2());\n";
  vs1.bind(v2());
  if(out) *out << "vs1 =\n" << vs1;

  // ///////////////////////
  // Test DVector Resizing

  if(out)
    *out	<< "\n***\n*** Testing DVector resizing\n***\n"
        << "\nv2.free();\n";
  v2.free();
  if(out) *out << "v2.dim() == 0 : " << update_success( v2.dim() == 0, &success ) << std::endl;

  if(out)
    *out	<< "\nv2.resize(n,2*alpha1);\n";
  v2.resize(n,2*alpha1);
  result1 = update_success( v2.dim() == n, &success );
  result2 = update_success( comp(v2,2*alpha1), &success );
  if(out)
    *out	<< "( (v2.dim() -> " << v2.dim() << ") == n && "
        << "(comp(v2,2*alpha1) -> " << result2 << ") ) : " << (result1 && result2) << std::endl
        << "v2 =\n" << v2;

  // //////////////////////////////////////////////////////////
  // Test Iterator Access, Subscriping and Reverse VectorSlices

  if(out) *out << "\n***\n*** Testing Iterator Access, Subscriping and Reverse VectorSlices\n***\n";

  if(out)
    *out	<< "\nLet v == v3, begin = v.begin()";
  test_access(&v3,v3.begin(),v3.end(),out,&success);

  if(out)
    *out	<< "\nLet v == const_cast<const DVector&>(v3), begin = v.begin()";
  test_access(&const_cast<const DVector&>(v3),const_cast<const DVector&>(v3).begin()
    ,const_cast<const DVector&>(v3).end(),out,&success);

  if(out)
    *out	<< "\nLet v == vs3, begin = v.begin()";
  test_access(&vs3,vs3.begin(),vs3.end(),out,&success);

  if(out)
    *out	<< "\nLet v == const_cast<const DVectorSlice&>(vs3), begin = v.begin()";
  test_access(&const_cast<const DVectorSlice&>(vs3),const_cast<const DVectorSlice&>(vs3).begin()
    ,const_cast<const DVectorSlice&>(vs3).end(),out,&success);

  if(out)
    *out	<< "\nLet v == v3.rev(), begin = v3.rbegin()";
  test_access(&v3.rev(),v3.rbegin(),v3.rend(),out,&success);

  if(out)
    *out	<< "\nLet v == const_cast<const DVector&>(v3).rev(), begin = const_cast<const DVector&>(v3).rbegin()";
  test_access(&const_cast<const DVector&>(v3).rev(),const_cast<const DVector&>(v3).rbegin()
    ,const_cast<const DVector&>(v3).rend(),out,&success);

  if(out)
    *out	<< "\nLet v == vs3.rev(), begin = vs3.rbegin()";
  test_access(&vs3.rev(),vs3.rbegin(),vs3.rend(),out,&success);

  if(out)
    *out	<< "\nLet v == const_cast<const DVectorSlice&>(vs3).rev()"
          ", begin = const_cast<const DVectorSlice&>(vs3).rbegin()";
  test_access(&const_cast<const DVectorSlice&>(vs3).rev(),const_cast<const DVectorSlice&>(vs3).rbegin()
    ,const_cast<const DVectorSlice&>(vs3).rend(),out,&success);

#ifdef LINALGPACK_CHECK_RANGE

  if(out) *out << "\n*** Test subscriping out of bounds\n";

  if(out) *out << "\nv3(0); (should throw std::out_of_range)\n";
  try{
    v3(0);
    result = false;	
  }
  catch(const std::out_of_range&) {
    result = true;
  }
  if(out) *out << "v3(0) threw std::out_of_range : " << result << std::endl;
  update_success( result, &success );
  
  if(out) *out << "\nv3(n+1); (should throw std::out_of_range)\n";
  try{
    v3(n+1);
    result = false;	
  }
  catch(const std::out_of_range&) {
    result = true;
  }
  if(out) *out << "v3(n+1) threw std::out_of_range : " << result << std::endl;
  update_success( result, &success );

  if(out) *out << "\nvs3(0); (should throw std::out_of_range)\n";
  try{
    vs3(0);
    result = false;	
  }
  catch(const std::out_of_range &) {
    result = true;
  }
  if(out) *out << "vs3(0) threw std::out_of_range : " << result << std::endl;
  update_success( result, &success );
  
  if(out) *out << "\nvs3(n+1); (should throw std::out_of_range)\n";
  try{
    vs3(n+1);
    result = false;	
  }
  catch(const std::out_of_range&) {
    result = true;
  }
  if(out) *out << "vs3(n+1) threw std::out_of_range : " << result << std::endl;
  update_success( result, &success );

#endif

  // ////////////////////////////////
  // Test Subregion Access

  if(out) *out << "\n***\n*** Testing Subregion Access\n***\n";

  // DVector Subregions
  Range1D rng;

  if(out) *out << "\nv = v3, rng = [1,n/2], vs = v(rng)";
  rng = Range1D(1,n/2);
  test_subregion_access( &v3, &v3(rng), rng, out, &success );

  if(out) *out << "\nv = v3, rng = [n/3,2*n/3], vs = v(rng)";
  rng = Range1D(n/3,2*n/3);
  test_subregion_access( &v3, &v3(rng), rng, out, &success );
  
  if(out) *out << "\nv = v3, rng = [n/2,n], vs = v(rng)";
  rng = Range1D(n/2,n);
  test_subregion_access( &v3, &v3(rng), rng, out, &success );

  if(out) *out << "\nv = const_cast<const DVector&>(v3), rng = [n/2,n], vs = v(rng)";
  rng = Range1D(n/2,n);
  test_subregion_access( &v3, &const_cast<const DVector&>(v3)(rng), rng, out, &success );

  if(out) *out << "\nv = v3, rng = [1,n/2], vs = v(1,n/2)";
  rng = Range1D(1,n/2);
  test_subregion_access( &v3, &v3(1,n/2), rng, out, &success );

  if(out) *out << "\nv = const_cast<const DVector&>(v3), rng = [n/2,n], vs = v(n/2,n)";
  rng = Range1D(n/2,n);
  test_subregion_access( &v3, &const_cast<const DVector&>(v3)(n/2,n), rng, out, &success );

  // DVectorSlice Subregions

  if(out) *out << "\nv = vs3, rng = [1,n/2], vs = v(rng)";
  rng = Range1D(1,n/2);
  test_subregion_access( &vs3, &vs3(rng), rng, out, &success );

  if(out) *out << "\nv = vs3, rng = [n/3,2*n/3], vs = v(rng)";
  rng = Range1D(n/3,2*n/3);
  test_subregion_access( &vs3, &vs3(rng), rng, out, &success );
  
  if(out) *out << "\nv = vs3, rng = [n/2,n], vs = v(rng)";
  rng = Range1D(n/2,n);
  test_subregion_access( &vs3, &vs3(rng), rng, out, &success );

  if(out) *out << "\nv = const_cast<const DVectorSlice&>(vs3), rng = [n/2,n], vs = v(rng)";
  rng = Range1D(n/2,n);
  test_subregion_access( &vs3, &const_cast<const DVectorSlice&>(vs3)(rng), rng, out, &success );

  if(out) *out << "\nv = vs3, rng = [1,n/2], vs = v(1,n/2)";
  rng = Range1D(1,n/2);
  test_subregion_access( &vs3, &vs3(1,n/2), rng, out, &success );

  if(out) *out << "\nv = const_cast<const DVectorSlice&>(vs3), rng = [n/2,n], vs = v(n/2,n)";
  rng = Range1D(n/2,n);
  test_subregion_access( &vs3, &const_cast<const DVectorSlice&>(vs3)(n/2,n), rng, out, &success );

  // ///////////////////////
  // Test Assignment

  if(out) *out << "\n***\n*** Testing assignment operators\n***\n";

  // DVector Assignment

  if(out) *out << "\nv1.resize(n); v1 = 0.0;\n";
  v1.resize(n);
  v1 = 0.0;
  result = update_success( comp( v1, 0.0 ), &success );
  if(out)
    *out	<< "v1 =\n" << v1
        << "v1 == 0.0 : " << result << std::endl;
  
  if(out) *out << "\nv1 = 0.0; v1 = vs3;\n";
  v1 = 0.0;
  v1 = vs3;
  result = update_success( comp( v1, vs3 ), &success );
  if(out)
    *out	<< "v1 =\n" << v1
        << "v1 == vs3 : " << result << std::endl;

  if(out) *out << "\nv1 = 0.0; v1 = v3;\n";
  v1 = 0.0;
  v1 = v3;
  result = update_success( comp( v1, v3 ), &success );
  if(out)
    *out	<< "v1 =\n" << v1
        << "v1 == v3 : " << result << std::endl;

  // DVectorSlice Assignment

  if(out) *out << "\nv1.resize(n); v1 = 0.0; vs1.bind(v1());\n";
  v1.resize(n);
  v1 = 0.0;
  vs1.bind(v1());

  if(out) *out << "\nvs1 = alpha1;\n";
  vs1 = alpha1;
  result = update_success( comp( vs1, alpha1 ), &success );
  if(out)
    *out	<< "vs1 =\n" << v1
        << "vs1 == alpha1 : " << result << std::endl;
  
  if(out) *out << "\nvs1 = 0.0; vs1 = vs3;\n";
  vs1 = 0.0;
  vs1 = vs3;
  result = update_success( comp( vs1, vs3 ), &success );
  if(out)
    *out	<< "vs1 =\n" << vs1
        << "vs1 == vs3 : " << result << std::endl;

  // ////////////////////////
  // Test overlap()

  if(out) *out << "\n***\n*** Testing overlap\n***\n";

  EOverLap ovlap;

  // DVector overlap

  if(out) *out << "\n*** DVector overlap\n";

  if(out) *out << "(v1.overlap(v3) -> ";
  ovlap = v1.overlap(v3);
  result = update_success( ovlap == NO_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == NO_OVERLAP : " << result << std::endl;
  
  if(out) *out << "(v3.overlap(v3) -> ";
  ovlap = v3.overlap(v3);
  result = update_success( ovlap == SAME_MEM, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SAME_MEM : " << result << std::endl;

  if(out) *out << "(v3.overlap(v3(1,n-1)) -> ";
  ovlap = v3.overlap(v3(1,n-1));
  result = update_success( ovlap == SOME_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

  if(out) *out << "(v3.overlap(v3(2,n-1)) -> ";
  ovlap = v3.overlap(v3(2,n-1));
  result = update_success( ovlap == SOME_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

  // DVectorSlice overlap

  if(out) *out << "\n*** DVectorSlice overlap\n"
          << "vs1.bind(v3());\n";

  vs1.bind(v3());

  if(out) *out << "(vs1.overlap(v2) -> ";
  ovlap = vs1.overlap(v2);
  result = update_success( ovlap == NO_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == NO_OVERLAP : " << result << std::endl;

  if(out) *out << "(vs1.overlap(vs1) -> ";
  ovlap = vs1.overlap(vs1);
  result = update_success( ovlap == SAME_MEM, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SAME_MEM : " << result << std::endl;

  if(out) *out << "(vs1(1,n/2).overlap(vs1(n/2+1,n)) -> ";
  ovlap = vs1(1,n/2).overlap(vs1(n/2+1,n));
  result = update_success( ovlap == NO_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == NO_OVERLAP : " << result << std::endl;

  if(out) *out << "(vs1(1,n/2).overlap(vs1(n/3,2*n/3)) -> ";
  ovlap = vs1(1,n/2).overlap(vs1(n/3,2*n/3));
  result = update_success( ovlap == SOME_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

  if(out) *out << "(vs1(n/3,2*n/3).overlap(vs1(n/2+1,n)) -> ";
  ovlap = vs1(n/3,2*n/3).overlap(vs1(n/2+1,n));
  result = update_success( ovlap == SOME_OVERLAP, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

  if(out) *out << "(vs1(n/3,2*n/3).overlap(vs1(n/3,2*n/3)) -> ";
  ovlap = vs1(n/3,2*n/3).overlap(vs1(n/3,2*n/3));
  result = update_success( ovlap == SAME_MEM, &success );
  if(out)	*out	<< overlap_str(ovlap) << ") == SAME_MEM : " << result << std::endl;

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
        << "\n*** Congradulations, DVector and DVectorSlice seem to check out. ***\n";
    else
      (*out)
        << "\n*** Oops, all of the tests for DVector and DVectorSlice "
          "where not successful. ***\n";
  }

  return success;
}

