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
//
// 2/10/00: g++ 2.95.2 is giving me some trouble with update_success(...) in template
// functions for some reason?
//

#include <iomanip>
#include <ostream>
#include <vector>
#include <typeinfo>

#include "DenseLinAlgPack_TestDenseLinAlgPack.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_MatVecCompare.hpp"

namespace {

using DenseLinAlgPack::size_type;
using DenseLinAlgPack::value_type;
using DenseLinAlgPack::DVector;
using DenseLinAlgPack::DVectorSlice;
using DenseLinAlgPack::DMatrix;
using DenseLinAlgPack::DMatrixSlice;

// vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs
template<class M_t>
void test_MtV( M_t& M_rhs1, BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2
  , const DVectorSlice& expected_MtV, DVector* tmp1, DVector* tmp2
  , std::ostream* out, bool* success )
{
  using BLAS_Cpp::rows;
  using DenseLinAlgPack::Vt_S;
  using DenseLinAlgPack::Vp_StV;
  using DenseLinAlgPack::Vp_StMtV;
  using DenseLinAlgPack::comp;
  using TestingHelperPack::update_success;
  
  // Check alpha = 1.0, 0.5.  Check beta = 0.0, 0.5, 1.0
  bool result = true;
  value_type
    alpha = 1.0,
    beta = 0.0;
  
  const size_type
    m = M_rhs1.rows(),
    n = M_rhs1.cols();
  
  tmp1->resize( rows(m,n,trans_rhs1) );
  tmp2->resize( tmp1->dim() );

  for( int i = 0; i < 2; ++i, alpha -= 0.5 ) {
    beta = 0.0;
    for( int j = 0; j < 3; ++j, beta += 0.5 ) {
      if(out)
        *out	<< "vs_lhs = " << alpha << " * M_rhs1"
            << ( trans_rhs1 == BLAS_Cpp::trans ? '\'' : ' ' )
            << " * vs_rhs2 + " << beta << " * vs_lhs : ";
      // Initialize vs_lhs
      for( int k = 1; k <= tmp1->dim(); ++k )
        (*tmp1)(k) = (*tmp2)(k) = k;
      // Compute expected
      Vt_S( &(*tmp1)(), beta );
      Vp_StV( &(*tmp1)(), alpha, expected_MtV );
      // Computed
      Vp_StMtV( &(*tmp2)(), alpha, M_rhs1, trans_rhs1, vs_rhs2, beta );
      // Compare
      update_success( result = comp( *tmp1, *tmp2 ), success );
      if(out) *out << result << std::endl;
      if(out && !result) *out << "expected_vs_lhs =\n" << *tmp1
                  << "computed_vs_lhs =\n"<< *tmp2;
    }
  }
}

// Print a '\'' for trans and ' ' for no_trans
inline char trans_char(BLAS_Cpp::Transp _trans) {
  return (_trans == BLAS_Cpp::trans ? '\'' : ' ' );
}

// Test DenseLinAlgPack compatable Level-3 BLAS (Mp_StMtM(...))
template <class M_t>
void test_MtM( M_t& B, BLAS_Cpp::Transp trans_B, const DMatrixSlice& I
  , BLAS_Cpp::Transp trans_I, std::ostream* out
  , DMatrix* C1, DMatrix* C2, bool* success )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  using DenseLinAlgPack::assign;
  using DenseLinAlgPack::Mp_StM;
  using DenseLinAlgPack::Mp_StMtM;
  using DenseLinAlgPack::comp;
  using TestingHelperPack::update_success;
  
  value_type alpha = 1.0, beta = 0.0;
  bool result = true;

  for( int i = 1; i <= 3; ++i ) {
    // Set alpha + (beta)(0.5) = 1.0
    switch(i) {
      case 1:
        alpha = 1.0;
        beta = 0.0;
        break;
      case 2:
        alpha = 0.75;
        beta = 0.5;
        break;
      case 3:
        alpha = 0.5;
        beta = 1.0;
        break;
    }
    // Perform tests
    const char
      B_trans_chr = (trans_B==no_trans?' ':'\''),
      I_trans_chr = (trans_I==no_trans?' ':'\'');
    // expected result
    const size_type
      rows_B = rows( B.rows(), B.cols(), trans_B ),
      cols_B = cols( B.rows(), B.cols(), trans_B );
    C1->resize(rows_B,cols_B);
    (*C1) = 0.0;
    Mp_StM( &(*C1)(), 1.0, B, trans_B );
    // (left)
    if(out) *out << "C = " << alpha << " * B" << B_trans_chr
        << "* I" << I_trans_chr << " + " << beta << " * (0.5) * B" << B_trans_chr
        << " : ";
    C2->resize(rows_B,cols_B);
    (*C2) = 0.0;
    Mp_StM( &(*C2)(), 0.5, B, trans_B );
    size_type n = cols( B.rows(), B.cols(), trans_B );
    Mp_StMtM( &(*C2)(), alpha, B, trans_B, I(1,n,1,n), trans_I, beta );
    update_success( result = comp( (*C1)(), (*C2)() ), success );
    if(out) *out << result << std::endl;
    if( out && !result ) *out << "C1 =\n" << (*C1) << "C2 =\n" << (*C2);
    // (right)
    if(out) *out << "C = " << alpha << " * I" << I_trans_chr
        << "* B" << B_trans_chr << " + " << beta << " * (0.5) * B" << B_trans_chr
        << " : ";
    (*C2) = 0.0;
    Mp_StM( &(*C2)(), 0.5, B, trans_B );
    n = rows( B.rows(), B.cols(), trans_B );
    Mp_StMtM( &(*C2)(), alpha, I(1,n,1,n), trans_I, B, trans_B, beta );
    update_success( result = comp( (*C1)(), (*C2)() ), success );
    if(out) *out << result << std::endl;
    if( out && !result ) *out << "C1 =\n" << (*C1) << "C2 =\n" << (*C2);

  }
}

}	// end namespace

bool DenseLinAlgPack::TestingPack::TestGenMatrixOp(std::ostream* out)
{

  using BLAS_Cpp::trans;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::upper;
  using BLAS_Cpp::lower;
  using BLAS_Cpp::unit;
  using BLAS_Cpp::nonunit;
  using BLAS_Cpp::trans_to_bool;
  using DenseLinAlgPack::comp;
  using DenseLinAlgPack::sqrt_eps;

  bool success = true;
  bool result, result1, result2;

  if(out)
    *out	<< "\n**************************************"
        << "\n*** Testing GenMatrixOp Operations ***"
        << "\n**************************************\n"
        << std::boolalpha << std::setw(10);

  // Arrays for looping through options.

  BLAS_Cpp::Uplo		a_uplo[2]			= { lower,		upper		};
  char				str_uplo[2][10]		= { "lower",	"upper"		};

  BLAS_Cpp::Diag		a_diag[2]			= { unit,		nonunit		};
  char				str_diag[2][10]		= { "unit",		"nonunit"	};

  BLAS_Cpp::Transp	a_trans[2]			= { trans,		no_trans	};
//	char				str_trans[2][10]	= { "trans",	"no_trans"	};

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

  DMatrix gm_lhs(m,n);
  const DMatrixSlice gms_rhs(const_cast<value_type*>(ptr),m*n,m,m,n);

  // ////////////////////////////////////////////////////////////////////////////////
  // Test Element-wise Assignment functions not already tested in TestGenMatrixClass

  if(out) *out << "\n***\n*** Test Element-wise Assignment functions not already "
          "tested in TestGenMatrixClass(...)\n***\n";
  
  // gms_lhs = op(gms_rhs)
  if(out) *out << "\ngms_lhs = op(gms_rhs)\n"
         << "\ngm_lhs.resize(n,m); assign( &gm_lhs(), gms_rhs, trans );\n";
  gm_lhs.resize(n,m);
  assign( &gm_lhs(), gms_rhs, trans );
  update_success( result = comp( gm_lhs(), no_trans, gms_rhs, trans ), &success );
  if(out && !result) *out << "gm_lhs =\n" << gm_lhs();
  if(out) *out << "gm_lhs == gms_rhs' : " << result << std::endl;

  // gms_lhs = op(gms_rhs)
  if(out) *out << "\ngm_lhs = op(gms_rhs)\n"
         << "\ngm_lhs.resize(1,1); assign( &gm_lhs, gms_rhs, trans );\n";
  gm_lhs.resize(1,1);
  assign( &gm_lhs, gms_rhs, trans );
  update_success( result = comp( gm_lhs(), no_trans, gms_rhs, trans ), &success );
  if(out && !result) *out << "gm_lhs =\n" << gm_lhs();
  if(out) *out << "gm_lhs == gms_rhs' : " << result << std::endl;

  // tri_ele_gms_lhs = alpha (elementwise)

  if(out) *out << "\ntri_ele_gms_lhs = alpha (elementwise)\n";

  if(out) *out << "\ngms_lhs.resize(m,m); gm_lhs = 1.0; assign(&nonconst_tri_ele(gm_lhs(),lower),2.0);\n";
  gm_lhs.resize(m,m);
  gm_lhs = 1.0;
  assign(&nonconst_tri_ele(gm_lhs(),lower),2.0);
  update_success( result1 = comp(tri_ele(gm_lhs(),lower),2.0), &success );
  update_success( result2 = comp(tri_ele(gm_lhs(1,m-1,2,m),upper),1.0), &success );
  if(out && (!result1 || !result2) )
    *out << "gm_lhs =\n" << gm_lhs();
  if(out) *out << "tri_ele(gm_lhs(),lower) == 2.0 : " << result1 << std::endl
         << "tri_ele(gm_lhs(1,m-1,2,m),upper) == 1.0 : " << result2 << std::endl;

  if(out) *out << "\nassign(&tri_ele(gm_lhs(1,m-1,2,m),upper),3.0);\n";
  assign(&nonconst_tri_ele(gm_lhs(1,m-1,2,m),upper),3.0);
  update_success( result1 = comp(tri_ele(gm_lhs(),lower),2.0), &success );
  update_success( result2 = comp(tri_ele(gm_lhs(1,m-1,2,m),upper),3.0), &success );
  if(out && (!result1 || !result2) )
    *out << "gm_lhs =\n" << gm_lhs();
  if(out) *out << "tri_ele(gm_lhs(),lower) == 2.0 : " << result1 << std::endl
         << "tri_ele(gm_lhs(1,m-1,2,m),upper) == 3.0 : " << result2 << std::endl;

  // tri_ele_gms_lhs = tri_ele_gms_rhs

  if(out) *out << "\ntri_ele_gms_lhs = tri_ele_gms_rhs\n"
         << "\nassign(&tri_ele(gm_lhs(2,m,1,m-1),lower),tri_ele(gms_rhs(1,m-1,2,m),upper));\n";
  assign(&nonconst_tri_ele(gm_lhs(2,m,1,m-1),lower),tri_ele(gms_rhs(1,m-1,2,m),upper));
  if(out) *out << "assign(&tri_ele(gm_lhs(1,m-1,2,m),upper),tri_ele(gms_rhs(2,m,1,m-1),lower));\n";
  assign(&nonconst_tri_ele(gm_lhs(1,m-1,2,m),upper),tri_ele(gms_rhs(2,m,1,m-1),lower));
  if(out) *out << "gm_lhs.diag() = gms_rhs.diag();\n";
  gm_lhs.diag() = gms_rhs.diag();
  update_success( result1 = comp( gm_lhs(1,m,1,m), no_trans, gms_rhs(1,m,1,m), trans ), &success );
  update_success( result2 = comp( tri_ele(gm_lhs(2,m,1,m-1),lower), tri_ele(gms_rhs(1,m-1,2,m),upper) ), &success );
  if(out && (!result1 || !result2) )
    *out << "gm_lhs =\n" << gm_lhs();
  if(out) *out << "gm_lhs(1,m,1,m) == gms_rhs(1,m,1,m)' : " << result1 << std::endl
         << "tri_ele(gm_lhs(2,m,1,m-1),lower) == tri_ele(gms_rhs(1,m-1,2,m),upper) : " << result2 << std::endl;

  // //////////////////////////////////////////
  // Test Element-wise Algebraic Operations

  if(out) *out << "\n***\n*** Test Element-wise Algebraic Operations\n***\n";

  // gms_lhs *= alpha
  if(out) *out << "\ngms_lhs *= alpha\n"
         << "\ngm_lhs = 1.0; Mt_S( &gm_lhs(), 2.0 );\n";
  gm_lhs = 1.0;
  Mt_S( &gm_lhs(), 2.0 );
  update_success( result = comp( gm_lhs(), 2.0 ), &success );
  if(out && !result)
    *out << "gm_lhs =\n" << gm_lhs();
  if(out) *out << "gm_lhs == 2.0: " << result << std::endl;

  // tri_lhs *= alpha
  {
    if(out) *out << "\ntri_lhs *= alpha\n"
           << "\ngm_lhs = 1.0;\nLet tm1 = tri_ele(gm_lhs(1,m,1,m),BLAS_Cpp::lower), "
            "tm2 = tri_ele(gm_lhs(1,m-1,2,m),BLAS_Cpp::upper)\n";
    gm_lhs = 1.0;
    DMatrixSliceTriEle
      tm1 = nonconst_tri_ele(gm_lhs(1,m,1,m),BLAS_Cpp::lower),
      tm2 = nonconst_tri_ele(gm_lhs(1,m-1,2,m),BLAS_Cpp::upper);
    if(out) *out << "Mt_S( &tm1, 2.0 );\n";
    Mt_S( &tm1, 2.0 );
    if(out) *out << "Mt_S( &tm2, 3.0 );\n";
    Mt_S( &tm2, 3.0 );
    update_success( result1 = comp( tm1, 2.0 ), &success );
    update_success( result2 = comp( tm2, 3.0 ), &success );
    if(out && (!result1 || !result2) )
      *out << "gm_lhs =\n" << gm_lhs();
    if(out) *out << "tm1 == 2.0 : " << result1 << std::endl
           << "tm2 == 3.0 : " << result2 << std::endl;
  }

  // tri_lhs += alpha * tri_rhs
  if(out) *out << "\ntri_lhs += alpha * tri_rhs\n"
            "gm_lhs = 1.0;\n";
  gm_lhs = 1.0;
  if(out) *out << "Mp_StM( &tri_ele(gm_lhs(2,m,1,m-1),lower), 2.0, tri_ele(gm_lhs(1,m-1,2,m),upper) );\n";
  Mp_StM( &nonconst_tri_ele(gm_lhs(2,m,1,m-1),lower), 2.0, tri_ele(gm_lhs(1,m-1,2,m),upper) );
  update_success( result1 = comp( tri_ele(gm_lhs(2,m,1,m-1),lower), 3.0 ), &success );
  update_success( result2 = comp( tri_ele(gm_lhs(1,m,1,m),upper), 1.0 ), &success );
  if(out && (!result1 || !result2) )
    *out << "gm_lhs =\n" << gm_lhs();
  if(out) *out << "tri_ele(gm_lhs(2,m,1,m-1),lower) == 3.0 : " << result1 << std::endl
         << "tri_ele(gm_lhs(1,m,1,m),upper) == 1.0 : " << result2 << std::endl;

  // LinAlgOpPack compatable (Mp_StM(...)).
  if(out) *out << "\n****** LinAlgOpPack compatable (Mp_StM(...))\n";

  // gms_lhs += alpha * op(gms_rhs)
  if(out) *out << "\ngms_lhs += alpha * op(gms_rhs)\n"
            "--------------------------\n"
          "gm_lhs.resize(m,n); gm_lhs = 0.0;"
          "Mp_StM( &gm_lhs(), 0.5, gm_rhs(), no_trans ); "
          "Mp_StM( &gm_lhs(), 0.5, gm_rhs(), no_trans );\n";
  gm_lhs.resize(m,n);
  gm_lhs = 0.0;
  Mp_StM( &gm_lhs(), 0.5, gms_rhs, no_trans );
  Mp_StM( &gm_lhs(), 0.5, gms_rhs, no_trans );
  update_success( result = comp( gm_lhs(), gms_rhs ), &success );
  if(out && !result)
    *out << "gm_lhs =\n" << gm_lhs();
  if(out) *out << "gm_lhs == gm_rhs : " << result << std::endl;

  // gms_lhs += alpha * op(sym_rhs)
  if(out) *out << "\ngms_lhs += alpha * op(sym_rhs)\n"
            "------------------------------\n";
  gm_lhs.resize(m,m);
  {for(int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    for(int i_trans = 0; i_trans < 2; ++i_trans) {
      const BLAS_Cpp::Uplo	_uplo	= a_uplo[i_uplo];
      const BLAS_Cpp::Transp	_trans	= a_trans[i_trans];
      if(out)
        *out
          << "gms_lhs += alpha * sym(gms_rhs(1,m,1,m),"
          << str_uplo[i_uplo] << ")" << (_trans == trans ? '\'' : ' ' ) << " : ";
      gm_lhs = 0.0;
      const DMatrixSliceSym M = sym(gms_rhs(1,m,1,m),_uplo);
      Mp_StM( &gm_lhs(), 0.5, M, _trans );
      Mp_StM( &gm_lhs(), 0.5, M, _trans );
      update_success( result1 = comp( tri_ele(gm_lhs(),lower), tri_ele(M.gms(),_uplo) )
        , &success );
      update_success( result2 = comp( tri_ele(gm_lhs(),upper), tri_ele(M.gms(),_uplo) )
        , &success );
      if(out) *out << ( result1 && result2 ) << std::endl;
      if( out && ( !result1 || !result2 ) )
        *out << "gm_lhs =\n" << gm_lhs();
    }
  }}

  // gms_lhs += alpha * op(tri_rhs)
  if(out) *out << "\ngms_lhs += alpha * op(tri_rhs)\n"
            "------------------------------\n";
  gm_lhs.resize(m,m);
  {for(int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    for(int i_diag = 0; i_diag < 2; ++i_diag) {
      for(int i_trans = 0; i_trans < 2; ++i_trans) {
        const BLAS_Cpp::Uplo	_uplo	= a_uplo[i_uplo];
        const BLAS_Cpp::Diag	_diag	= a_diag[i_diag];
        const BLAS_Cpp::Transp	_trans	= a_trans[i_trans];
        if(out)
          *out
            << "gms_lhs += alpha * tri(gms_rhs(1,m,1,m),"
            << str_uplo[i_uplo] << "," << str_diag[i_diag] << ")"
            << (_trans == trans ? '\'' : ' ' ) << " : ";
        // compute
        gm_lhs = 0.0;
        const DMatrixSliceTri M = tri(gms_rhs(1,m,1,m),_uplo,_diag);
        Mp_StM( &gm_lhs(), 0.5, M, _trans );
        Mp_StM( &gm_lhs(), 0.5, M, _trans );
        // test diagonal
        if( _diag == nonunit )
          result = comp( gm_lhs.diag(), gms_rhs.diag() );
        else
          result = comp( gm_lhs.diag(), 1.0 );
        update_success( result, &success ); 
        // test the rest
        const BLAS_Cpp::Uplo
          as_uplo = (		( _uplo == lower && _trans == no_trans )
                ||	( _uplo == upper && _trans == trans )
                ?	lower : upper											),
          oth_as_uplo = ( as_uplo == lower ? upper : lower );
        const DMatrixSliceTriEle
          M_ele = tri_ele( ( _uplo==lower ? M.gms()(2,m,1,m-1) : M.gms()(1,m-1,2,m) )
                    , _uplo ),
          tri_reg_ele = tri_ele( ( as_uplo==lower ? gm_lhs(2,m,1,m-1) : gm_lhs(1,m-1,2,m) )
                    , as_uplo ),
          oth_tri_reg_ele = tri_ele( ( oth_as_uplo==lower ? gm_lhs(2,m,1,m-1) : gm_lhs(1,m-1,2,m) )
                    , oth_as_uplo );
        update_success( result1 = comp( tri_reg_ele, M_ele ), &success );
        update_success( result2 = comp( oth_tri_reg_ele, 0.0 ), &success );
        if(out) *out << ( result && result1 && result2 ) << std::endl;
        if( out && ( !result || !result1 || !result2 ) )
          *out << "gm_lhs =\n" << gm_lhs();
      }
    }
  }}

  // //////////////////////////////////////////
  // Level-2 BLAS

  if(out) *out << "\n***\n*** Level-2 BLAS \n***\n";

  DVector tmp1, tmp2, vs_rhs(n), expected_MtV;
  {for(int k = 1; k <= n; ++k) vs_rhs(k) = k; }

  // *** Triagnular matrices
  if(out) *out << "\n*** BLAS-2 Equivalent Triangular Matrix operations\n";

  {for(int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    for(int i_diag = 0; i_diag < 2; ++i_diag) {
      for(int i_trans = 0; i_trans < 2; ++i_trans) {
        if(out)
          *out << "\nLet M = tri(gms_rhs(1,m,1,m)," << str_uplo[i_uplo]
             << "," << str_diag[i_diag] << ")"
             << (a_trans[i_trans] == trans ? '\'' : ' ' ) << std::endl;
        DMatrixSliceTri M = tri(gms_rhs(1,m,1,m),a_uplo[i_uplo],a_diag[i_diag]);
        // v_lhs
        tmp1.resize(1);
        tmp2.resize(1);
        V_InvMtV( &tmp1, M, a_trans[i_trans], vs_rhs(1,m) );
        V_MtV( &tmp2, M, a_trans[i_trans], tmp1() );
        update_success( result = comp( tmp2, vs_rhs(1,m)), &success );
        if(out) *out << "(v_lhs)... M * (inv(M)*vs_rhs(1,m)) == vs_rhs(1,m) : "
               << result << std::endl;
        // vs_lhs
        V_InvMtV( &tmp1(), M, a_trans[i_trans], vs_rhs(1,m) );
        V_MtV( &tmp1(), M, a_trans[i_trans], tmp1() );
        update_success( result = comp( tmp1, vs_rhs(1,m)), &success );
        if(out) *out << "(vs_lhs)... M * (inv(M)*vs_rhs(1,m)) == vs_rhs(1,m) : "
               << result << std::endl;
      }
    }
  }}

  // *** LinAlgOpPack complient
  if(out) *out << "\n*** Test DenseLinAlgPack complient (Vp_StMtV(...))\n";

  // vs_lhs = alpha * op(gms_rhs1) * vs_rhs2 + beta * vs_lhs (xGEMV)

  if(out) *out << "\n****** General Matrix MtV\n";

  // no_trans
  if(out) *out << "\nvs_lhs = alpha * gms_rhs1 * vs_rhs2 + beta * vs_lhs (xGEMV)\n";
  expected_MtV.resize(m);
  expected_MtV = 0.0;
  {for(int i=1;i<=m;++i)
    for(int j=1;j<=n;++j)
      expected_MtV(i) += gms_rhs(i,j) * vs_rhs(j);
  }
  test_MtV( gms_rhs, no_trans, vs_rhs(1,n), expected_MtV, &tmp1, &tmp2
    , out, &success );

  // trans
  if(out) *out << "\nvs_lhs = alpha * gms_rhs1' * vs_rhs2 + beta * vs_lhs (xGEMV)\n";
  expected_MtV.resize(n);
  expected_MtV = 0.0;
  {for(int i=1;i<=n;++i)
    for(int j=1;j<=m;++j)
      expected_MtV(i) += gms_rhs(j,i) * vs_rhs(j);
  }
  test_MtV( gms_rhs, trans, vs_rhs(1,m), expected_MtV, &tmp1, &tmp2
    , out, &success );

  // vs_lhs = alpha * op(sym_rhs1) * vs_rhs2 + beta * vs_lhs (xSYMV)

  if(out) *out << "\n****** Symmetric Matrix MtV\n";

  {for(int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    for(int i_trans = 0; i_trans < 2; ++i_trans) {
      if(out)
        *out
          << "\nvs_lhs = alpha * sym(gms_rhs(1,m,1,m),"
          << str_uplo[i_uplo] << ")"
          << (a_trans[i_trans] == trans ? '\'' : ' ' )
          << " * vs_rhs2 + beta * vs_lhs (xSYMV)\n";
      // compute expected MtV
      expected_MtV.resize(m);
      expected_MtV = 0.0;
      switch(a_uplo[i_uplo]) {
        case lower: {
          for(int i=1;i<=m;++i)
            for(int j=1;j<=m;++j)
              expected_MtV(i)
                += (j < i ? gms_rhs(i,j) : gms_rhs(j,i) ) * vs_rhs(j);
          break;
        }
        case upper: {
          for(int i=1;i<=m;++i)
            for(int j=1;j<=m;++j)
              expected_MtV(i)
                += (j > i ? gms_rhs(i,j) : gms_rhs(j,i) ) * vs_rhs(j);
          break;
        }
      }
      test_MtV( sym(gms_rhs(1,m,1,m),a_uplo[i_uplo]), a_trans[i_trans]
        , vs_rhs(1,m), expected_MtV, &tmp1, &tmp2, out, &success );
    }
  }}


  // vs_lhs = alpha * op(tri_rhs1) * vs_rhs2 + beta * vs_lhs.

  if(out) *out << "\n****** Triangular Matrix MtV\n";

  {for(int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    for(int i_diag = 0; i_diag < 2; ++i_diag) {
      for(int i_trans = 0; i_trans < 2; ++i_trans) {
        const BLAS_Cpp::Uplo	_uplo	= a_uplo[i_uplo];
        const BLAS_Cpp::Diag	_diag	= a_diag[i_diag];
        const BLAS_Cpp::Transp	_trans	= a_trans[i_trans];
        if(out)
          *out
            << "\nvs_lhs = alpha * tri(gms_rhs(1,m,1,m),"
            << str_uplo[i_uplo] << "," << str_diag[i_diag] 
            << ")" << (_trans == trans ? '\'' : ' ' )
            << " * vs_rhs2 + beta * vs_lhs\n";
        DMatrixSliceTri M = tri(gms_rhs(1,m,1,m),_uplo,_diag);
        //
        // Compute expected MtV
        //
        // Determine conceptually lower or upper triangular 
        const BLAS_Cpp::Uplo
          as_uplo = (		( _uplo == lower && _trans == no_trans )
                ||	( _uplo == upper && _trans == trans )
                ? lower : upper );
        // Compute expected
        expected_MtV.resize(m);
        expected_MtV = 0.0;
        switch(as_uplo) {
          case lower: {
            for(int i=1;i<=m;++i) {
              for(int j=1;j<i;++j) {
                const int
                  _i = ( _trans == no_trans ? i : j ),
                  _j = ( _trans == no_trans ? j : i );
                expected_MtV(i) += gms_rhs(_i,_j) * vs_rhs(j);
              }
              expected_MtV(i) += (_diag==unit?1.0:gms_rhs(i,i))*vs_rhs(i);
            }
            break;
          }
          case upper: {
            for(int i=1;i<=m;++i) {
              expected_MtV(i) += (_diag==unit?1.0:gms_rhs(i,i))*vs_rhs(i);
              for(int j=i+1;j<=m;++j) {
                const int
                  _i = ( _trans == no_trans ? i : j ),
                  _j = ( _trans == no_trans ? j : i );
                expected_MtV(i) += gms_rhs(_i,_j) * vs_rhs(j);
              }
            }
            break;
          }
        }
        test_MtV( tri(gms_rhs(1,m,1,m),_uplo,_diag), _trans, vs_rhs(1,m)
          , expected_MtV, &tmp1, &tmp2, out, &success );
      }
    }
  }}

  // *** Symmetric Matrix Updates.

  if(out) *out << "\n*** Symmetric Matrix Updates\n";

  // sym_lhs = alpha * vs_rhs * vs_rhs' + sym_lhs (BLAS xSYR).
  if(out) *out << "\n****** sym_lhs = alpha * vs_rhs * vs_rhs' + sym_lhs (BLAS xSYR)\n";

  gm_lhs.resize(m,m);
  gm_lhs = 0.0;
  DMatrix ex_gm_lhs(0.0,m,m);
  value_type alpha = 2.0;
    
  {for( int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    const BLAS_Cpp::Uplo _uplo = a_uplo[i_uplo];
    if(out)
      *out << "\nFor sym_lhs = sym(gms_rhs(1,m,1,m)," << str_uplo[i_uplo]
         << ")\n";
    assign( &nonconst_tri_ele(gm_lhs(),_uplo), tri_ele(gms_rhs(1,m,1,m),_uplo) );
    assign( &nonconst_tri_ele(ex_gm_lhs(),_uplo), tri_ele(gms_rhs(1,m,1,m),_uplo) );
    // Compute expected
    for(int i = 1; i<=m; ++i) {
      for(int j = i; j<=m; ++j) {	// upper triangular
        const int
          _i = (_uplo == upper ? i : j),	// adjust for lower triangular
          _j = (_uplo == upper ? j : i);
        ex_gm_lhs(_i,_j) += alpha * vs_rhs(i) * vs_rhs(j);
      }
    }
    // Compute update
    syr( alpha, vs_rhs(1,m), &nonconst_sym(gm_lhs(),_uplo) );
    // Compare
    update_success( result = comp( tri_ele(gm_lhs(),_uplo), tri_ele(ex_gm_lhs(),_uplo) )
      , &success );
    if( out && !result )
      *out << "gm_lhs =\n" << gm_lhs << "ex_gm_lhs =\n" << ex_gm_lhs();
    if(out) *out << "    gm_lhs == expected_gm_lhs : " << result << std::endl;
  }}

  // sym_lhs = alpha * vs_rhs1 * vs_rhs2' + alpha * vs_rhs2 * vs_rhs1' + sym_lhs (BLAS xSYR2).
  if(out) *out << "\n****** sym_lhs = alpha * vs_rhs1 * vs_rhs2' + alpha * vs_rhs2 * vs_rhs1' + sym_lhs (BLAS xSYR2)\n";

  gm_lhs = 0.0;
  ex_gm_lhs = 0.0;
  alpha = 2.0;
  DVector vs_rhs2 = vs_rhs.rev();
    
  {for( int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    const BLAS_Cpp::Uplo _uplo = a_uplo[i_uplo];
    if(out)
      *out << "\nFor sym_lhs = sym(gms_rhs(1,m,1,m)," << str_uplo[i_uplo]
         << ")\n";
    assign( &nonconst_tri_ele(gm_lhs(),_uplo), tri_ele(gms_rhs(1,m,1,m),_uplo) );
    assign( &nonconst_tri_ele(ex_gm_lhs(),_uplo), tri_ele(gms_rhs(1,m,1,m),_uplo) );
    // Compute expected
    for(int i = 1; i<=m; ++i) {
      for(int j = i; j<=m; ++j) {	// upper triangular
        const int
          _i = (_uplo == upper ? i : j),	// adjust for lower triangular
          _j = (_uplo == upper ? j : i);
        ex_gm_lhs(_i,_j) += alpha * (vs_rhs(i) * vs_rhs2(j) + vs_rhs2(i) * vs_rhs(j));
      }
    }
    // Compute update
    syr2( alpha, vs_rhs(1,m), vs_rhs2(1,m), &nonconst_sym(gm_lhs(),_uplo) );
    // Compare
    update_success( result = comp( tri_ele(gm_lhs(),_uplo), tri_ele(ex_gm_lhs(),_uplo) )
      , &success );
    if( out && !result )
      *out << "gm_lhs =\n" << gm_lhs() << "ex_gm_lhs =\n" << ex_gm_lhs();
    if(out) *out << "    gm_lhs == expected_gm_lhs : " << result << std::endl;
  }}

  // /////////////////////////////////////////
  // Level-3 BLAS

  if(out) *out << "\n***\n*** Level-3 BLAS\n***\n";

  // Triangular Matrices
  if(out) *out << "\n*** BLAS-2 Equivalent Triangular Matrix operations\n";

  DMatrix Tmp1, Tmp2;

  {for(int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    for(int i_diag = 0; i_diag < 2; ++i_diag) {
      for(int i_trans1 = 0; i_trans1 < 2; ++i_trans1) {
        for(int i_trans2 = 0; i_trans2 < 2; ++i_trans2) {
          const BLAS_Cpp::Transp
            _trans1 = a_trans[i_trans1],
            _trans2 = a_trans[i_trans2]; 
          if(out)
            *out << "\nLet M = tri(gms_rhs(1,m,1,m)," << str_uplo[i_uplo]
               << "," << str_diag[i_diag] << ")"
               << trans_char(_trans1) << std::endl;
          DMatrixSliceTri M = tri(gms_rhs(1,m,1,m),a_uplo[i_uplo],a_diag[i_diag]);
          // gm_lhs (left)
          // gm_lhs = alpha * op(tri_rhs1) * op(gms_rhs2) (left) (BLAS xTRMM).
          // gm_lhs = alpha * inv(op(tri_rhs1)) * op(gms_rhs2) (left) (BLAS xTRSM).
          Tmp1.resize(1,1);
          Tmp2.resize(1,1);
          M_StInvMtM( &Tmp1, 2.0, M, _trans1, gms_rhs(1,m,1,m), _trans2 );
          M_StMtM( &Tmp2, 0.5, M, _trans1, Tmp1(), no_trans );
          update_success( result = comp( Tmp2(), no_trans, gms_rhs(1,m,1,m), _trans2 )
            , &success );
          if(out) *out << "(gm_lhs,left)...0.5*M*(2*inv(M)*gms_rhs(1,m,1,m)"
                 << trans_char(_trans2) << ")"
                  " == gms_rhs(1,m,1,m)" << trans_char(_trans2) << " : "
                  << result << std::endl;
          // gms_lhs (left)
          // gms_lhs = alpha * op(tri_rhs1) * op(gms_rhs2) (left) (BLAS xTRMM).
          // gms_lhs = alpha * inv(op(tri_rhs1)) * op(gms_rhs2) (left) (BLAS xTRSM).
          assign( &Tmp2, gms_rhs(1,m,1,m), _trans2 );
          M_StInvMtM( &Tmp1(), 2.0, M, _trans1, Tmp2(), no_trans );
          M_StMtM( &Tmp2, 0.5, M, _trans1, Tmp1(), no_trans );
          update_success( result = comp( Tmp2(), no_trans, gms_rhs(1,m,1,m), _trans2 )
            , &success );
          if(out) {
             *out
               << "(gms_lhs,left)...0.5*M*(2*inv(M)*gms_rhs(1,m,1,m)"
              << trans_char(_trans2) << ")"
                " == gms_rhs(1,m,1,m)" << trans_char(_trans2) << " : "
              << result << std::endl;
            if(!result) {
              *out
                << "\ngms_lhs =\n" << Tmp2
                << "\ngms_rhs(1,m,1,m) =\n" << gms_rhs(1,m,1,m);
            }
          }
          // gm_lhs (right)
          // gm_lhs = alpha * op(gms_rhs1) * op(tri_rhs2) (right) (BLAS xTRMM).
          // gm_lhs = alpha * op(gms_rhs1) * inv(op(tri_rhs2)) (right) (BLAS xTRSM).
          Tmp1.resize(1,1);
          Tmp2.resize(1,1);
          M_StMtInvM( &Tmp1, 2.0, gms_rhs(1,m,1,m), _trans1, M, _trans2 );
          M_StMtM( &Tmp2, 0.5, Tmp1(), no_trans, M, _trans2 );
          update_success( result = comp( Tmp2(), no_trans, gms_rhs(1,m,1,m), _trans1 )
            , &success );
          if(out) *out << "(gm_lhs,right)...5.0*(2.0*gms_rhs(1,m,1,m)"
                 << trans_char(_trans1) << "*inv(M))*M"
                  " == gms_rhs(1,m,1,m)" << trans_char(_trans1) << " : "
                  << result << std::endl;
          // gms_lhs (right)
          // gms_lhs = alpha * op(gms_rhs1) * op(tri_rhs2) (right) (BLAS xTRMM).
          // gms_lhs = alpha * op(gms_rhs1) * inv(op(tri_rhs2)) (right) (BLAS xTRSM).
          assign( &Tmp2(), gms_rhs(1,m,1,m), _trans1 );
          M_StMtInvM( &Tmp2(), 2.0, Tmp2(), no_trans, M, _trans2 );
          M_StMtM( &Tmp2(), 0.5, Tmp2(), no_trans, M, _trans2 );
          update_success( result = comp( Tmp2(), no_trans, gms_rhs(1,m,1,m), _trans1 )
            , &success );
          if(out) *out << "(gms_lhs,right)...5.0*(2.0*gms_rhs(1,m,1,m)"
                 << trans_char(_trans1) << "*inv(M))*M"
                  " == gms_rhs(1,m,1,m)" << trans_char(_trans1) << " : "
                  << result << std::endl;
        }
      }
    }
  }}

  // *** LinAlgOpPack complient
  if(out) *out << "\n*** Test DenseLinAlgPack complient (Mp_StMtM(...))\n";

  DMatrix I(0.0,n,n);
  I.diag() = 1.0;

  // ****** Rectangular Matrices 
  if(out) *out << "\n****** Rectangular MtM (BLAS xGEMV)\n";

  // gms_lhs = alpha * op(gms_rhs1) * op(gms_rhs2) + beta * gms_lhs (BLAS xGEMV).
  {for(int i_trans_B = 0; i_trans_B < 2; ++i_trans_B) {
    for(int i_trans_I = 0; i_trans_I < 2; ++i_trans_I) {
      const BLAS_Cpp::Transp
        _trans_B = a_trans[i_trans_B],
        _trans_I = a_trans[i_trans_I];
      if(out) *out << "\nLet B = gms_rhs\n";
      test_MtM( gms_rhs, _trans_B, I, _trans_I, out, &Tmp1, &Tmp2, &success );
    }
  }}

  // ****** Symmetric Matrices
  if(out) *out << "\n****** Symmetric MtM\n";

  // gms_lhs = alpha * op(sym_rhs1) * op(gms_rhs2) + beta * gms_lhs (left) (BLAS xSYMM).
  // gms_lhs = alpha * op(gms_rhs1) * op(sym_rhs2) + beta * gms_lhs (right) (BLAS xSYMM).
  {for(int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    for(int i_trans_B = 0; i_trans_B < 2; ++i_trans_B) {
      for(int i_trans_I = 0; i_trans_I < 2; ++i_trans_I) {
        const BLAS_Cpp::Uplo
          _uplo = a_uplo[i_uplo];
        const BLAS_Cpp::Transp
          _trans_B = a_trans[i_trans_B],
          _trans_I = a_trans[i_trans_I]; 
        if(out) *out << "\nLet B = sym(gms_rhs(1,m,1,m)," << str_uplo[i_uplo] << ")\n";
        test_MtM( sym(gms_rhs(1,m,1,m),_uplo), _trans_B, I, _trans_I, out
          , &Tmp1, &Tmp2, &success );
      }
    }
  }}

  // ****** Triangular Matrices
  if(out) *out << "\n****** Triangular MtM\n";

  // gms_lhs = alpha * op(tri_rhs1) * op(gms_rhs2) + beta * gms_lhs (left).
  // gms_lhs = alpha * op(gms_rhs1) * op(tri_rhs2) + beta * gms_lhs (right).

  {for(int i_uplo = 0; i_uplo < 2; ++i_uplo) {
    for(int i_diag = 0; i_diag < 2; ++i_diag) {
      for(int i_trans_B = 0; i_trans_B < 2; ++i_trans_B) {
        for(int i_trans_I = 0; i_trans_I < 2; ++i_trans_I) {
          const BLAS_Cpp::Uplo
            _uplo = a_uplo[i_uplo];
          const BLAS_Cpp::Diag
            _diag = a_diag[i_diag];
          const BLAS_Cpp::Transp
            _trans_B = a_trans[i_trans_B],
            _trans_I = a_trans[i_trans_I]; 
          if(out) *out << "\nLet B = tri(gms_rhs(1,m,1,m)," << str_uplo[i_uplo]
                 << "," << str_diag[i_diag] << ")\n";
          test_MtM( tri(gms_rhs(1,m,1,m),_uplo,_diag), _trans_B, I, _trans_I, out
            , &Tmp1, &Tmp2, &success );
        }
      }
    }
  }}

  // *** Symmetric Matrix updating
  if(out) *out << "\n*** Symmetric Matrix updating\n";

  if(out) *out << "\nWarning! Not Tested!\n";

  // sym_lhs = alpha * op(gms_rhs) * op(gms_rhs')  + beta * sym_lhs (BLAS xSYRK).
  // sym_lhs = alpha * op(gms_rhs1) * op(gms_rhs2') + alpha * op(gms_rhs2) * op(gms_rhs1') + beta * sym_lhs (BLAS xSYR2K)

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
        << "\n*** Congradulations, GenMatrixOp operations seem to check out. ***\n";
    else
      (*out)
        << "\n*** Oops, all of the tests for GenMatrixOp operations "
          "where not successful. ***\n";
  }


  return success;
}
