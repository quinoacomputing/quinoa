/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file LVQDTargetTest.cpp
 *  \brief unit tests for LVQDTargetCalculator class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_LVQDTargetCalculator.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_PatchData.hpp"
#include "UnitUtil.hpp"
#include "PatchDataInstances.hpp"

#include <iostream>

using namespace Mesquite;

class LVQDTargetTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(LVQDTargetTest);

    // Test that LVQDTargetCalculator accepts a NULL component
    // and treats the corresponding value as I
  CPPUNIT_TEST (test_LVQD_default_is_I_2D);
  CPPUNIT_TEST (test_LVQD_default_is_I_surface);
  CPPUNIT_TEST (test_LVQD_default_is_I_3D);
    // Test that LVQDTargetCalculator returns the product of
    // its components.
  CPPUNIT_TEST (test_LVQD_product_2D);
  CPPUNIT_TEST (test_LVQD_product_surface);
  CPPUNIT_TEST (test_LVQD_product_3D);

  CPPUNIT_TEST_SUITE_END();

    // Define a few pre-factored matrices to use in tests.
    // Initialized by SetUp() method.
  MsqMatrix<3,3> V3D_Z45, V3D_X90, Q3D_45, D3D_123, I33;
  MsqMatrix<3,2> V2D_Z45, V2D_X90, I32;
  MsqMatrix<2,2> Q2D_45, D2D_21, I22;
    // PatchDatas for use in calls to 2D or 3D versions of functions.
  PatchData pd2D, pd3D;

    // Use LVQDTargetCalculator to calculate the product of the passed values.
  MsqMatrix<3,2> target( const double* L, 
                         const MsqMatrix<3,2>* V, 
                         const MsqMatrix<2,2>* Q, 
                         const MsqMatrix<2,2>* D );
    // Use LVQDTargetCalculator to calculate the product of the passed values.
  MsqMatrix<2,2> target( const double* L, 
                         const MsqMatrix<2,2>* Q, 
                         const MsqMatrix<2,2>* D );
    // Use LVQDTargetCalculator to calculate the product of the passed values.
  MsqMatrix<3,3> target( const double* L, 
                         const MsqMatrix<3,3>* V, 
                         const MsqMatrix<3,3>* Q, 
                         const MsqMatrix<3,3>* D );

public:

  void setUp();

    // Test that LVQDTargetCalculator accepts a NULL component
    // and treats the corresponding value as I
  void test_LVQD_default_is_I_2D();
  void test_LVQD_default_is_I_3D();
  void test_LVQD_default_is_I_surface();
   // Test that LVQDTargetCalculator returns the product of
    // its components.
  void test_LVQD_product_2D();
  void test_LVQD_product_3D();
  void test_LVQD_product_surface();

    // Helper class: return constant values for target matrices.
  class ConstantTarget : public TargetCalculator
  {
    private:
      MsqMatrix<3,3> target3D;
      MsqMatrix<3,2> targetSurf;
      MsqMatrix<2,2> target2D;
      bool have3D, haveSurf, have2D;
      bool flagError;

    public:
      ConstantTarget( MsqMatrix<3,3> val3D, MsqMatrix<3,2> val2D ) 
        : target3D(val3D), targetSurf(val2D), have3D(true), haveSurf(true), have2D(false) 
        {}
      ConstantTarget( MsqMatrix<3,3> val3D, MsqMatrix<2,2> val2D ) 
        : target3D(val3D), target2D(val2D), have3D(true), haveSurf(false), have2D(true) 
        { }
      ConstantTarget( double C, bool surf )
        : target3D(C), targetSurf(C), target2D(C), have3D(true), haveSurf(surf), have2D(!surf) 
        {}
      ConstantTarget( MsqMatrix<3,3> val3D ) 
        : target3D(val3D), target2D(0.0), have3D(true), haveSurf(false), have2D(false) 
        {}
      ConstantTarget( MsqMatrix<3,2> val2D ) 
        : targetSurf(val2D), have3D(false), haveSurf(true), have2D(false) 
        {}
      ConstantTarget( MsqMatrix<2,2> val2D ) 
        : target2D(val2D), have3D(false), haveSurf(false), have2D(true) 
        { }

      virtual bool get_3D_target( PatchData&, size_t, Sample, MsqMatrix<3,3>& result, MsqError& err )
        { CPPUNIT_ASSERT(have3D); result = target3D; return have3D; }

      virtual bool get_surface_target( PatchData&, size_t, Sample, MsqMatrix<3,2>& result, MsqError& err )
        { CPPUNIT_ASSERT(haveSurf); result = targetSurf; return haveSurf; }

      virtual bool get_2D_target( PatchData&, size_t, Sample, MsqMatrix<2,2>& result, MsqError& err )
        { CPPUNIT_ASSERT(have2D); result = target2D; return have2D; }
        
      virtual bool have_surface_orient() const { return haveSurf; }
  };

    // Helper class: return 'invalid' for target matrices, and
    // optionally flag an error.
  class TargetError : public TargetCalculator
  {
      bool flagError;
      bool surfOrient;

    public:
      TargetError( bool flag_error = true, bool orient = false ) 
        : flagError(flag_error), surfOrient(orient) {}

      bool get_3D_target( PatchData&, size_t, Sample, MsqMatrix<3,3>&, MsqError& err)
        { if (flagError) MSQ_SETERR(err)(MsqError::INVALID_MESH); return false; }

      bool get_2D_target( PatchData&, size_t, Sample, MsqMatrix<2,2>&, MsqError& err )
        { if (flagError) MSQ_SETERR(err)(MsqError::INVALID_MESH); return false; }

      bool get_surface_target( PatchData&, size_t, Sample, MsqMatrix<3,2>&, MsqError& err )
        { if (flagError) MSQ_SETERR(err)(MsqError::INVALID_MESH); return false; }
        
      virtual bool have_surface_orient() const { return surfOrient; }
  };
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LVQDTargetTest, "LVQDTargetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LVQDTargetTest, "Unit");


void LVQDTargetTest::setUp()
{
  const double cos45 = MSQ_SQRT_TWO/2.0;
  const double rotation_3D_Z45[9] = { cos45, -cos45, 0,
                                      cos45,  cos45, 0,
                                          0,      0, 1 };
  const double rotation_2D_Z45[6] = { cos45, -cos45,
                                      cos45,  cos45,
                                          0,      0 };
                                          
  const double rotation_3D_X90[9] = { 1,  0,  0,
                                      0,  0, -1,
                                      0,  1,  0 };
  const double rotation_2D_X90[6] = { 1,  0,
                                      0,  0,
                                      0,  1 };
                                      
  const double rc45 = sqrt(cos45);
  const double skew_2D_45[4] = { 1/rc45, rc45,
                                 0,      rc45 };
  const double skew_3D_45[9] = { 1, cos45, cos45,
                                 0, cos45, 1 - cos45,
                                 0,     0, sqrt(MSQ_SQRT_TWO-1) };
                                 
  const double aspect_2D_2x[4] = { MSQ_SQRT_TWO, 0,
                                   0,            MSQ_SQRT_TWO/2 };
  const double r6 = Mesquite::cbrt(1.0/6.0);
  const double aspect_3D_123[9] = { r6,    0,    0,
                                     0, 2*r6,    0,
                                     0,    0, 3*r6 };
                                     
  V3D_Z45 = MsqMatrix<3,3>(rotation_3D_Z45);
  V3D_X90 = MsqMatrix<3,3>(rotation_3D_X90);
  Q3D_45  = MsqMatrix<3,3>(skew_3D_45);
  Q3D_45  *= 1/Mesquite::cbrt(det(Q3D_45));
  D3D_123 = MsqMatrix<3,3>(aspect_3D_123);
  
  V2D_Z45 = MsqMatrix<3,2>(rotation_2D_Z45);
  V2D_X90 = MsqMatrix<3,2>(rotation_2D_X90);
  Q2D_45  = MsqMatrix<2,2>(skew_2D_45);
  D2D_21  = MsqMatrix<2,2>(aspect_2D_2x);
  
  I33 = MsqMatrix<3,3>(1.0);
  I32 = MsqMatrix<3,2>(1.0);
  I22 = MsqMatrix<2,2>(1.0);
  
  MsqError err;
  create_one_tet_patch( pd3D, err );
  ASSERT_NO_ERROR(err);
  create_one_tri_patch( pd2D, err );
  ASSERT_NO_ERROR(err);
};

MsqMatrix<3,2> LVQDTargetTest::target( const double* L, 
                                       const MsqMatrix<3,2>* V, 
                                       const MsqMatrix<2,2>* Q, 
                                       const MsqMatrix<2,2>* D )
{
  ConstantTarget W_size  ( L ? *L : 1.0, true );
  ConstantTarget W_orient( V ? *V : I32 );
  ConstantTarget W_skew  ( Q ? *Q : I22 );
  ConstantTarget W_aspect( D ? *D : I22 );
  LVQDTargetCalculator LVQD( L ? &W_size   : NULL,
                             V ? &W_orient : NULL,
                             Q ? &W_skew   : NULL,
                             D ? &W_aspect : NULL );
  MsqError err;
  MsqMatrix<3,2> W;
  bool v;
  if (!V) {
    MsqMatrix<2,2> W_2D;
    v = LVQD.get_2D_target( pd2D, 0, Sample(0,0), W_2D, err );
    W.set_row( 0, W_2D.row(0) );
    W.set_row( 1, W_2D.row(1) );
    W.set_row( 2, MsqMatrix<1,2>(0.0) );
  }
  else {
    v = LVQD.get_surface_target( pd2D, 0, Sample(0,0), W, err );
  }
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return W;
}

MsqMatrix<2,2> LVQDTargetTest::target( const double* L, 
                                       const MsqMatrix<2,2>* Q, 
                                       const MsqMatrix<2,2>* D )
{
  ConstantTarget W_size  ( L ? *L : 1.0, false );
  ConstantTarget W_skew  ( Q ? *Q : I22 );
  ConstantTarget W_aspect( D ? *D : I22 );
  LVQDTargetCalculator LVQD( L ? &W_size   : NULL,
                             NULL,
                             Q ? &W_skew   : NULL,
                             D ? &W_aspect : NULL );
  MsqError err;
  MsqMatrix<2,2> W;
  bool v = LVQD.get_2D_target( pd2D, 0, Sample(0,0), W, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return W;
}

MsqMatrix<3,3> LVQDTargetTest::target( const double* L, 
                                       const MsqMatrix<3,3>* V, 
                                       const MsqMatrix<3,3>* Q, 
                                       const MsqMatrix<3,3>* D )
{
  ConstantTarget W_size  ( L ? *L : 1.0, true );
  ConstantTarget W_orient( V ? *V : I33 );
  ConstantTarget W_skew  ( Q ? *Q : I33 );
  ConstantTarget W_aspect( D ? *D : I33 );
  LVQDTargetCalculator LVQD( L ? &W_size   : NULL,
                             V ? &W_orient : NULL,
                             Q ? &W_skew   : NULL,
                             D ? &W_aspect : NULL );
  MsqError err;
  MsqMatrix<3,3> W;
  bool v = LVQD.get_3D_target( pd3D, 0, Sample(0,0), W, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(v);
  return W;
}

void LVQDTargetTest::test_LVQD_default_is_I_2D()
{
  double s = 3.2;
  ASSERT_MATRICES_EQUAL( I22, target(0,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<2,2>(s)), target(&s,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( Q2D_45, target(0,&Q2D_45,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( D2D_21, target(0,0,&D2D_21), 1e-8 );
}

void LVQDTargetTest::test_LVQD_default_is_I_surface()
{
  double s = 3.2;
  MsqMatrix<3,2>* null_V = 0;
  ASSERT_MATRICES_EQUAL( I32, target(0,null_V,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,2>(s)), target(&s,null_V,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45, target(0,&V2D_Z45,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V2D_Z45, target(&s,&V2D_Z45,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45*Q2D_45, target(0,&V2D_Z45,&Q2D_45,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45*D2D_21, target(0,&V2D_Z45,0,&D2D_21), 1e-8 );
}

void LVQDTargetTest::test_LVQD_default_is_I_3D()
{
  double s = 2.6;
  MsqMatrix<3,3>* null_V = 0;
  ASSERT_MATRICES_EQUAL( I33, target(0,null_V,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,3>(s)), target(&s,null_V,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45, target(0,&V3D_Z45,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V3D_Z45, target(&s,&V3D_Z45,0,0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45*Q3D_45,  target(0,&V3D_Z45,&Q3D_45, 0), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45*D3D_123, target(0,&V3D_Z45,0,&D3D_123), 1e-8 );
}

void LVQDTargetTest::test_LVQD_product_2D()
{
  double s = 3.2;
  double o = 1.0;
  ASSERT_MATRICES_EQUAL( I22, target(&o,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<2,2>(s)), target(&s,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( Q2D_45, target(&o,&Q2D_45,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( D2D_21, target(&o,&I22,&D2D_21), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*Q2D_45*D2D_21, target(&s,&Q2D_45,&D2D_21), 1e-8 );
}

void LVQDTargetTest::test_LVQD_product_surface()
{
  double s = 3.2;
  double o = 1.0;
  ASSERT_MATRICES_EQUAL( I32, target(&o,&I32,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,2>(s)), target(&s,&I32,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45, target(&o,&V2D_Z45,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V2D_Z45, target(&s,&V2D_Z45,&I22,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45*Q2D_45, target(&o,&V2D_Z45,&Q2D_45,&I22), 1e-8 );
  ASSERT_MATRICES_EQUAL( V2D_Z45*D2D_21, target(&o,&V2D_Z45,&I22,&D2D_21), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V2D_Z45*Q2D_45*D2D_21, target(&s,&V2D_Z45,&Q2D_45,&D2D_21), 1e-8 );
}

void LVQDTargetTest::test_LVQD_product_3D()
{
  double s = 2.6;
  double o = 1.0;
  ASSERT_MATRICES_EQUAL( I33, target(&o,&I33,&I33,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( (MsqMatrix<3,3>(s)), target(&s,&I33,&I33,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45, target(&o,&V3D_Z45,&I33,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V3D_Z45, target(&s,&V3D_Z45,&I33,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45*Q3D_45, target(&o,&V3D_Z45,&Q3D_45,&I33), 1e-8 );
  ASSERT_MATRICES_EQUAL( V3D_Z45*D3D_123, target(&o,&V3D_Z45,&I33,&D3D_123), 1e-8 );
  ASSERT_MATRICES_EQUAL( s*V3D_Z45*Q3D_45*D3D_123, target(&s,&V3D_Z45,&Q3D_45,&D3D_123), 1e-8 );
}
