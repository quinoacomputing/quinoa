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


/** \file IdealTargetTest.cpp
 *  \brief Test the IdealShapeTarget
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_IdealShapeTarget.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_Vector3D.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_MappingFunction.hpp"
#include "Mesquite_Settings.hpp"
#include "Mesquite_IdealElements.hpp"
#include "Mesquite_ElemSampleQM.hpp"
#include "Mesquite_cppunit/extensions/HelperMacros.h"
#include "UnitUtil.hpp"

using namespace Mesquite;

class IdealTargetTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(IdealTargetTest);
  CPPUNIT_TEST (test_tri_corner);
  CPPUNIT_TEST (test_tri_edge);
  CPPUNIT_TEST (test_tri_center);
  CPPUNIT_TEST (test_hex_corner);
  CPPUNIT_TEST (test_hex_edge);
  CPPUNIT_TEST (test_hex_face);
  CPPUNIT_TEST (test_hex_center);
  CPPUNIT_TEST_SUITE_END();

public:

  void test_tri_corner();
  void test_tri_edge();
  void test_tri_center();
  void test_hex_corner();
  void test_hex_edge();
  void test_hex_face();
  void test_hex_center();
  
private:

  void get_calc_target( EntityTopology type, 
                        Sample sample,
                        MsqMatrix<3,3>&, 
                        MsqMatrix<2,2>& );
                        
  void get_ideal_target( EntityTopology type, 
                         Sample sample,
                         MsqMatrix<3,3>&, 
                         MsqMatrix<2,2>& );
                        
  void do_test( EntityTopology type, Sample location );
  
  Settings settings;
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(IdealTargetTest, "IdealTargetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(IdealTargetTest, "Unit");

void IdealTargetTest::test_tri_corner()
{
  do_test( TRIANGLE, Sample(0, 0) );
  do_test( TRIANGLE, Sample(0, 1) );
  do_test( TRIANGLE, Sample(0, 2) );
}

void IdealTargetTest::test_tri_edge()
{
  do_test( TRIANGLE, Sample(1, 0) );
  do_test( TRIANGLE, Sample(1, 1) );
  do_test( TRIANGLE, Sample(1, 2) );
}

void IdealTargetTest::test_tri_center()
{
  do_test( TRIANGLE, Sample(2, 0) );
}

void IdealTargetTest::test_hex_corner()
{
  do_test( HEXAHEDRON, Sample(0, 0) );
  do_test( HEXAHEDRON, Sample(0, 1) );
  do_test( HEXAHEDRON, Sample(0, 6) );
  do_test( HEXAHEDRON, Sample(0, 7) );
}

void IdealTargetTest::test_hex_edge()
{
  do_test( HEXAHEDRON, Sample(1, 1) );
  do_test( HEXAHEDRON, Sample(1, 4) );
  do_test( HEXAHEDRON, Sample(1, 11) );
}

void IdealTargetTest::test_hex_face()
{
  do_test( HEXAHEDRON, Sample(2, 0) );
  do_test( HEXAHEDRON, Sample(2, 1) );
  do_test( HEXAHEDRON, Sample(2, 2) );
  do_test( HEXAHEDRON, Sample(2, 3) );
  do_test( HEXAHEDRON, Sample(2, 4) );
  do_test( HEXAHEDRON, Sample(2, 5) );
}

void IdealTargetTest::test_hex_center()
{
  do_test( HEXAHEDRON, Sample(3, 0) );
}

void IdealTargetTest::get_calc_target( EntityTopology type, 
                                       Sample location,
                                       MsqMatrix<3,3>& w3, 
                                       MsqMatrix<2,2>& w2 )
{
  MsqPrintError err( std::cout );
  const int elem_dim = TopologyInfo::dimension(type);
  
    // create a patch -- actual coords and such don't really matter
  std::vector<double> coords( 24, 0.0 );
  const size_t conn[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  PatchData pd;
  pd.fill( 8, arrptr(coords), 1, type, conn, 0, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  pd.attach_settings( &settings );
  
  IdealShapeTarget tc;
  if (elem_dim == 2)
    tc.get_2D_target( pd, 0, location, w2, err );
  else
    tc.get_3D_target( pd, 0, location, w3, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
}

void IdealTargetTest::get_ideal_target( EntityTopology type,
                                        Sample location,
                                        MsqMatrix<3,3>& w3, 
                                        MsqMatrix<2,2>& w2 )
{
  MsqPrintError err( std::cout );
  const unsigned elem_dim = TopologyInfo::dimension(type);
  
    // get the target matrix for an ideal element
  size_t indices[100];
  size_t num_vtx;
  const Vector3D* coords = unit_element( type );
  Vector3D c[3];
  if (elem_dim == 2) {
    MsqVector<2> derivs[100];
    const MappingFunction2D* func = settings.get_mapping_function_2D( type );
    func->derivatives( location, NodeSet(), indices, derivs, num_vtx, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));

    MsqMatrix<3,2> J;
    for (size_t i = 0; i < num_vtx; ++i) 
      for (unsigned j = 0; j < 2; ++j)
        c[j] += derivs[i][j] * coords[indices[i]];

    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 2; ++j)
        J(i,j) = c[j][i];
    
    w2 = TargetCalculator::skew(J);
  }
  else {
    MsqVector<3> derivs[100];
    const MappingFunction3D* func = settings.get_mapping_function_3D( type );
    func->derivatives( location, NodeSet(), indices, derivs, num_vtx, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));

    for (size_t i = 0; i < num_vtx; ++i) 
      for (unsigned j = 0; j < 3; ++j)
        c[j] += derivs[i][j] * coords[indices[i]];

    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        w3(i,j) = c[j][i];

    w3 = TargetCalculator::skew(w3);
  }
  
}

void IdealTargetTest::do_test( EntityTopology type, Sample location )
{
  MsqMatrix<3,3> w3_calc, w3_exp;
  MsqMatrix<2,2> w2_calc, w2_exp;
  get_calc_target( type, location, w3_calc, w2_calc );
  get_ideal_target( type, location, w3_exp, w2_exp );
  if (TopologyInfo::dimension(type) == 2) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(w2_calc), 1e-6 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(w2_exp), 1e-6 );
    // a rotation of the expected matrix is acceptable.
    MsqMatrix<2,2> R = inverse(w2_calc) * w2_exp;
    ASSERT_MATRICES_EQUAL( transpose(R), inverse(R), 1e-6 );
  }
  else {
    MsqMatrix<3,3> R = inverse(w3_calc) *w3_exp;
    ASSERT_MATRICES_EQUAL( transpose(R), inverse(R), 1e-6 );
  }
}

 
