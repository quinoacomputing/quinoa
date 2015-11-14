/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file HigherOrderTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#define HAVE_HO_HEX
#define TEST_HO_QUAD


#include "Mesquite.hpp"
#include "Mesquite_MsqError.hpp"

#include "Mesquite_ArrayMesh.hpp"
#include "Mesquite_DomainClassifier.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_MeshDomain1D.hpp"

#include "Mesquite_IdealShapeTarget.hpp"
#include "Mesquite_TInverseMeanRatio.hpp"
#include "Mesquite_TShapeSizeOrientNB1.hpp"
#include "Mesquite_TShapeSizeNB3.hpp"
#include "Mesquite_TQualityMetric.hpp"
#include "Mesquite_PMeanPTemplate.hpp"
#include "Mesquite_LPtoPTemplate.hpp"
#include "Mesquite_SteepestDescent.hpp"
#include "Mesquite_FeasibleNewton.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_TerminationCriterion.hpp"

//#include "Mesquite_TriLagrangeShape.hpp"
//#include "Mesquite_TetLagrangeShape.hpp"
//#include "Mesquite_QuadLagrangeShape.hpp"
#ifdef HAVE_HO_HEX
# include "HexLagrangeShape.hpp"
#endif

#include "UnitUtil.hpp"

using namespace Mesquite;

#include <iostream>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

const int MAX_ITERATIONS = 1000;
const double QEL = 1.0; // edge length of ideal quad

class HigherOrderTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(HigherOrderTest);
  CPPUNIT_TEST (test_tri_basic_ideal);
  CPPUNIT_TEST (test_tri_basic_mid_spin);
  CPPUNIT_TEST (test_tri_basic_mid_convex);
  CPPUNIT_TEST (test_tri_basic_peak_up);
  CPPUNIT_TEST (test_tri_basic_peak_down);
  CPPUNIT_TEST (test_tri_basic_peak_over);
#ifdef TEST_HO_QUAD
  CPPUNIT_TEST (test_quad_basic_ideal);
  CPPUNIT_TEST (test_quad_basic_mid_spin);
  CPPUNIT_TEST (test_quad_basic_mid_convex);
  CPPUNIT_TEST (test_quad_basic_left_down);
  CPPUNIT_TEST (test_quad_basic_top_down);
  CPPUNIT_TEST (test_quad_basic_right_up);
  CPPUNIT_TEST (test_quad_basic_left_over);
#endif
  CPPUNIT_TEST (test_tet_basic_ideal);
  CPPUNIT_TEST (test_tet_basic_mid_spin);
  CPPUNIT_TEST (test_tet_basic_mid_convex);
  CPPUNIT_TEST (test_tet_basic_apex_up);
  CPPUNIT_TEST (test_tet_basic_apex_down);
  CPPUNIT_TEST (test_tet_basic_apex_over);
#ifdef HAVE_HO_HEX
  CPPUNIT_TEST (test_hex_basic);
#endif
//  CPPUNIT_TEST (test_tri_open_domain);
//  CPPUNIT_TEST (test_tri_slac);
  CPPUNIT_TEST_SUITE_END();

//  TriLagrangeShape tri_shape;
//  TetLagrangeShape tet_shape;
//  QuadLagrangeShape quad_shape;
#ifdef HAVE_HO_HEX
  HexLagrangeShape hex_shape;
#endif
  InstructionQueue q;
  TShapeSizeNB3 tm;
  IdealShapeTarget tc;
  TQualityMetric metric;
  PMeanPTemplate func;
  SteepestDescent solver;
  TerminationCriterion crit, outer;

public:

  HigherOrderTest() : metric( &tc, &tm ), 
                      func( 1, &metric ),
                      solver( &func )
  {
    MsqError err;
//    q.set_mapping_function( &tri_shape );
//    q.set_mapping_function( &tet_shape );
//    q.set_mapping_function( &quad_shape );
#ifdef HAVE_HO_HEX
    q.set_mapping_function( &hex_shape );
#endif
    q.set_master_quality_improver( &solver, err ); 
    
    q.set_slaved_ho_node_mode( Settings::SLAVE_NONE );
    
    outer.add_iteration_limit( 1 );
    crit.add_absolute_vertex_movement( 1e-6 );
    crit.add_iteration_limit( MAX_ITERATIONS );
    solver.set_outer_termination_criterion( &outer );
    solver.set_inner_termination_criterion( &crit );
  }
  
  bool hit_iteration_limit() const
    { return crit.get_iteration_count() >= MAX_ITERATIONS; }

  void test_tri_basic_ideal();
  void test_tri_basic_mid_spin();
  void test_tri_basic_mid_convex();
  void test_tri_basic_peak_up();
  void test_tri_basic_peak_down();
  void test_tri_basic_peak_over();
  void test_quad_basic_ideal();
  void test_quad_basic_mid_spin();
  void test_quad_basic_mid_convex();
  void test_quad_basic_left_down();
  void test_quad_basic_top_down();
  void test_quad_basic_right_up();
  void test_quad_basic_left_over();
  void test_tet_basic_ideal();
  void test_tet_basic_mid_spin();
  void test_tet_basic_mid_convex();
  void test_tet_basic_apex_up();
  void test_tet_basic_apex_down();
  void test_tet_basic_apex_over();
  void test_hex_basic();
  void test_tri_open_domain();
  void test_tri_slac();
  
  void test_tri_open_domain( double& x1, double& x3, double& x4,
                             double  y2, double  y5, double& y4 );
                           
 
  void basic_tri_test( double& x2, double& y2,
                       double& x3, double& y3,
                       double& x4, double& y4,
                       double& x5, double& y5,
                       MsqError& err );

  void basic_tet_test( Vector3D& p3,
                       Vector3D& p4,
                       Vector3D& p5,
                       Vector3D& p6,
                       Vector3D& p7,
                       Vector3D& p8,
                       Vector3D& p9,
                       MsqError& err );

  void basic_quad_test( Vector3D& p2,
                        Vector3D& p3,
                        Vector3D& p4,
                        Vector3D& p5,
                        Vector3D& p6,
                        Vector3D& p7,
                        MsqError& err );

  void basic_hex_test( Vector3D& p4,
                       Vector3D& p5,
                       Vector3D& p6,
                       Vector3D& p7,
                       Vector3D& p8,
                       Vector3D& p9,
                       Vector3D& p10,
                       Vector3D& p11,
                       Vector3D& p12,
                       Vector3D& p13,
                       Vector3D& p14,
                       Vector3D& p15,
                       Vector3D& p16, 
                       Vector3D& p17,
                       Vector3D& p18,
                       Vector3D& p19,
                       MsqError& err );

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(HigherOrderTest, "HigherOrderTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(HigherOrderTest, "Regression");

static inline double dist( double x1, double y1, double x2, double y2 )
  { return sqrt( (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) ); }

// ideal triangle is defined such that J has unit determinant
const double IDEAL_TRI_SIDE = sqrt(2.0/sqrt(3.0));
const double IDEAL_TRI_HEIGHT = 1/IDEAL_TRI_SIDE;

int tri_mid_edge_nodes_edge_center( double x2, double y2,
                                    double x3, double y3,
                                    double x4, double y4,
                                    double x5, double y5,
                                    double epsilon )
{
  double x0 = 0.0, y0 = 0.0;
  double x1 = IDEAL_TRI_SIDE, y1 = 0.0;
  double x01 = 0.5 * (x0 + x1);
  double x12 = 0.5 * (x1 + x2);
  double x20 = 0.5 * (x2 + x0);
  double y01 = 0.5 * (y0 + y1);
  double y12 = 0.5 * (y1 + y2);
  double y20 = 0.5 * (y2 + y0);
  
  int result = 0;
  if (dist(x3,y3,x01,y01) > epsilon)
    result |= 1;
  if (dist(x4,y4,x12,y12) > epsilon)
    result |= 2;
  if (dist(x5,y5,x20,y20) > epsilon)
    result |= 4;
  
  return result;
}

void HigherOrderTest::test_tri_basic_ideal()
{
  MsqPrintError err(cerr);
  const double eps = 1e-4;
  const double x2eq = 0.5 * IDEAL_TRI_SIDE;
  const double y2eq = IDEAL_TRI_HEIGHT;
  
    // try starting with the optimal result
  double x2 = x2eq;
  double y2 = y2eq;
  double x3 = 0.5*IDEAL_TRI_SIDE;
  double y3 = 0.0;
  double x4 = 0.75*IDEAL_TRI_SIDE;
  double y4 = 0.5*IDEAL_TRI_HEIGHT;
  double x5 = 0.25*IDEAL_TRI_SIDE;
  double y5 = 0.5*IDEAL_TRI_HEIGHT;
//  crit.write_mesh_steps( "test_tri_basic_ideal", TerminationCriterion::VTK );
  basic_tri_test( x2, y2, x3, y3, x4, y4, x5, y5, err );
//  crit.write_mesh_steps( "", TerminationCriterion::NOTYPE );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( x2eq, x2, eps );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( y2eq, y2, eps );
  int midok = tri_mid_edge_nodes_edge_center( x2, y2, x3, y3, x4, y4, x5, y5, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tri_basic_mid_spin()
{
  MsqPrintError err(cerr);
  const double eps = 1e-4;
  const double x2eq = 0.5 * IDEAL_TRI_SIDE;
  const double y2eq = IDEAL_TRI_HEIGHT;
  
    // try moving the mid-edge nodes along the edge away from the edge center
  double x0 = 0.0, y0 = 0.0;
  double x1 = IDEAL_TRI_SIDE, y1 = 0.0;
  double x2 = x2eq;
  double y2 = y2eq;
  const double f = 0.7;
  double x3 =  f * x0 + (1-f) * x1;
  double y3 =  f * y0 + (1-f) * y1;
  double x4 =  f * x1 + (1-f) * x2;
  double y4 =  f * y1 + (1-f) * y2;
  double x5 =  f * x2 + (1-f) * x0;
  double y5 =  f * y2 + (1-f) * y0;
  //crit.write_mesh_steps( "test_tri_basic_mid_spin", TerminationCriterion::VTK );
  basic_tri_test( x2, y2, x3, y3, x4, y4, x5, y5, err );
  //crit.write_mesh_steps( "", TerminationCriterion::NOTYPE );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( x2eq, x2, eps );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( y2eq, y2, eps );
  int midok = tri_mid_edge_nodes_edge_center( x2, y2, x3, y3, x4, y4, x5, y5, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tri_basic_mid_convex()
{
  MsqPrintError err(cerr);
  const double eps = 1e-4;
  const double x2eq = 0.5 * IDEAL_TRI_SIDE;
  const double y2eq = IDEAL_TRI_HEIGHT;
  
    // try equilateral corners with all egdes convex
  double x2 = x2eq;
  double y2 = y2eq;
  double x3 = 0.5*IDEAL_TRI_SIDE;
  double y3 = -0.3;
  double x4 = IDEAL_TRI_SIDE;
  double y4 = 0.5;
  double x5 = 0.0;
  double y5 = 0.5;
  basic_tri_test( x2, y2, x3, y3, x4, y4, x5, y5, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( x2eq, x2, eps );
  //CPPUNIT_ASSERT_DOUBLES_EQUAL( y2eq, y2, eps );
  int midok = tri_mid_edge_nodes_edge_center( x2, y2, x3, y3, x4, y4, x5, y5, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tri_basic_peak_up()
{
  MsqPrintError err(cerr);
  const double eps = 1e-4;
  const double x2eq = 0.5 * IDEAL_TRI_SIDE;
  const double y2eq = IDEAL_TRI_HEIGHT;
  
    // try moving the top vertex up, also move mid-edge nodes proportionally
  double x2 = x2eq;
  double y2 = 2.0 * y2eq;
  double x3 =  0.5 * IDEAL_TRI_SIDE;
  double y3 =  0.0;
  double x4 =  1.5 * x2;
  double y4 =  0.5 * y2;
  double x5 =  0.5 * x2;
  double y5 =  0.5 * y2;
  basic_tri_test( x2, y2, x3, y3, x4, y4, x5, y5, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( x2eq, x2, eps );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( y2eq, y2, eps );
  int midok = tri_mid_edge_nodes_edge_center( x2, y2, x3, y3, x4, y4, x5, y5, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tri_basic_peak_down()
{
  MsqPrintError err(cerr);
  const double eps = 1e-4;
  const double x2eq = 0.5 * IDEAL_TRI_SIDE;
  const double y2eq = IDEAL_TRI_HEIGHT;
  
    // try moving the top vertex down, also move mid-edge nodes proportionally
  double x2 = x2eq;
  double y2 = 0.5 * y2eq;
  double x3 =  0.5 * IDEAL_TRI_SIDE;
  double y3 =  0.0;
  double x4 =  1.5 * x2;
  double y4 =  0.5 * y2;
  double x5 =  0.5 * x2;
  double y5 =  0.5 * y2;
  basic_tri_test( x2, y2, x3, y3, x4, y4, x5, y5, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( x2eq, x2, eps );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( y2eq, y2, eps );
  int midok = tri_mid_edge_nodes_edge_center( x2, y2, x3, y3, x4, y4, x5, y5, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tri_basic_peak_over()
{
  MsqPrintError err(cerr);
  const double eps = 1e-4;
  const double x2eq = 0.5 * IDEAL_TRI_SIDE;
  const double y2eq = IDEAL_TRI_HEIGHT;

    // try moving the top vertex to the right, also move mid-edge nodes proportionally
  double x2 = x2eq + 0.5;
  double y2 = y2eq;
  double x3 =  0.5 * IDEAL_TRI_SIDE;
  double y3 =  0.0;
  double x4 =  1.5 * x2;
  double y4 =  0.5 * y2;
  double x5 =  0.5 * x2;
  double y5 =  0.5 * y2;
  basic_tri_test( x2, y2, x3, y3, x4, y4, x5, y5, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( x2eq, x2, eps );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( y2eq, y2, eps );
  int midok = tri_mid_edge_nodes_edge_center( x2, y2, x3, y3, x4, y4, x5, y5, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

int quad_all_in_xy_plane( const Vector3D& p2,
                          const Vector3D& p3,
                          const Vector3D& p4,
                          const Vector3D& p5,
                          const Vector3D& p6,
                          const Vector3D& p7,
                          double epsilon )
{
  Vector3D list[] = { p2, p3, p4, p5, p6, p7 };
  const int n = sizeof(list)/sizeof(list[0]);
  int result = 0;
  for (int i = 0; i < n; ++i)
    if (fabs(list[i].z()) > epsilon)
      result |= (1 << i);
  return result;
}

int quad_mid_edge_nodes_edge_center( const Vector3D& p2,
                                     const Vector3D& p3,
                                     const Vector3D& p4,
                                     const Vector3D& p5,
                                     const Vector3D& p6,
                                     const Vector3D& p7,
                                     double epsilon )
{
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(QEL, 0.0, 0.0);
  const Vector3D e0 = 0.5*(p0+p1);
  const Vector3D e1 = 0.5*(p1+p2);
  const Vector3D e2 = 0.5*(p2+p3);
  const Vector3D e3 = 0.5*(p3+p0);
  
  int result = 0;
  if ((p4-e0).length() > epsilon)
    result |= 1;
  if ((p5-e1).length() > epsilon)
    result |= 2;
  if ((p6-e2).length() > epsilon)
    result |= 4;
  if ((p7-e3).length() > epsilon)
    result |= 8;
  
  return result;
}

static void get_ideal_quad( Vector3D& p2,
                            Vector3D& p3,
                            Vector3D& p4,
                            Vector3D& p5,
                            Vector3D& p6,
                            Vector3D& p7 )
{
  p2.set( QEL,   QEL,   0 );
  p3.set( 0.0,   QEL,   0 );
  p4.set( QEL/2, 0.0,   0 );
  p5.set( QEL,   QEL/2, 0 );
  p6.set( QEL/2, QEL,   0 );
  p7.set( 0.0,   QEL/2, 0 );
}

void HigherOrderTest::test_quad_basic_ideal()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  Vector3D p2, p3, p4, p5, p6, p7;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(QEL, 0.0, 0.0);
  
    // try starting with the optimal result
  get_ideal_quad( p2, p3, p4, p5, p6, p7 );
  basic_quad_test( p2, p3, p4, p5, p6, p7, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_EQUAL( 0, quad_all_in_xy_plane( p2, p3, p4, p5, p6, p7, eps ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( QEL, QEL, 0 ), p2, eps );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( 0,   QEL, 0 ), p3, eps );
  int midok = quad_mid_edge_nodes_edge_center( p2, p3, p4, p5, p6, p7, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_quad_basic_mid_spin()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  Vector3D p2, p3, p4, p5, p6, p7;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(QEL, 0.0, 0.0);
  
    // try moving the mid-edge nodes along the edge away from the edge center
  p2.set(QEL,QEL,0);
  p3.set(0  ,QEL,0);
  double f = 0.4;
  p4 = f*p0 + (1-f)*p1;
  p5 = f*p1 + (1-f)*p2;
  p6 = f*p2 + (1-f)*p3;
  p7 = f*p3 + (1-f)*p0;
//  crit.write_mesh_steps( "quad_basic_mid_spin", TerminationCriterion::VTK );
  basic_quad_test( p2, p3, p4, p5, p6, p7, err );
//  crit.write_mesh_steps( "", TerminationCriterion::NOTYPE );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_EQUAL( 0, quad_all_in_xy_plane( p2, p3, p4, p5, p6, p7, eps ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( QEL, QEL, 0 ), p2, eps );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( 0,   QEL, 0 ), p3, eps );
  int midok = quad_mid_edge_nodes_edge_center( p2, p3, p4, p5, p6, p7, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_quad_basic_mid_convex()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  Vector3D p2, p3, p4, p5, p6, p7;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(QEL, 0.0, 0.0);
  
    // try square corners with all egdes convex
  get_ideal_quad( p2, p3, p4, p5, p6, p7 );
  p4 += Vector3D( 0.0,     -0.2*QEL, 0);
  p5 += Vector3D( 0.2*QEL,  0.0,     0);
  p6 += Vector3D( 0.0,      0.2*QEL, 0);
  p7 += Vector3D(-0.2*QEL,  0.0,     0);
//  crit.write_mesh_steps( "quad_basic_mid_convex", TerminationCriterion::VTK );
  basic_quad_test( p2, p3, p4, p5, p6, p7, err );
//  crit.write_mesh_steps( "", TerminationCriterion::NOTYPE );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_EQUAL( 0, quad_all_in_xy_plane( p2, p3, p4, p5, p6, p7, eps ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( QEL, QEL, 0 ), p2, eps );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(   0, QEL, 0 ), p3, eps );
  int midok = quad_mid_edge_nodes_edge_center( p2, p3, p4, p5, p6, p7, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_quad_basic_left_down()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  Vector3D p2, p3, p4, p5, p6, p7;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(QEL, 0.0, 0.0);
  
    // try moving the top left vertex down, also move mid-edge nodes proportionally
  get_ideal_quad( p2, p3, p4, p5, p6, p7 );
  p3 -= Vector3D(0.0,0.0,QEL/2);
  p6 = 0.5*(p2+p3);
  p7 = 0.5*(p0+p3);
  basic_quad_test( p2, p3, p4, p5, p6, p7, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_EQUAL( 0, quad_all_in_xy_plane( p2, p3, p4, p5, p6, p7, eps ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( QEL, QEL, 0 ), p2, eps );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(   0, QEL, 0 ), p3, eps );
  int midok = quad_mid_edge_nodes_edge_center( p2, p3, p4, p5, p6, p7, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_quad_basic_top_down()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  Vector3D p2, p3, p4, p5, p6, p7;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(QEL, 0.0, 0.0);
  
    // try moving the top two vertices down, also move mid-edge nodes proportionally
  get_ideal_quad( p2, p3, p4, p5, p6, p7 );
  p2 -= Vector3D(0.0,QEL/2,0.0);
  p3 -= Vector3D(0.0,QEL/2,0.0);
  p5 = 0.5*(p1+p2);
  p6 = 0.5*(p2+p3);
  p7 = 0.5*(p0+p3);
  basic_quad_test( p2, p3, p4, p5, p6, p7, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_EQUAL( 0, quad_all_in_xy_plane( p2, p3, p4, p5, p6, p7, eps ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( QEL, QEL, 0 ), p2, eps );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(   0, QEL, 0 ), p3, eps );
  int midok = quad_mid_edge_nodes_edge_center( p2, p3, p4, p5, p6, p7, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_quad_basic_right_up()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  Vector3D p2, p3, p4, p5, p6, p7;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(QEL, 0.0, 0.0);
  
    // try moving the top right vertex up, also move mid-edge nodes proportionally
  get_ideal_quad( p2, p3, p4, p5, p6, p7 );
  p2 += Vector3D(0.0,QEL*2,0.0);
  p5 = 0.5*(p1+p2);
  p6 = 0.5*(p2+p3);
  basic_quad_test( p2, p3, p4, p5, p6, p7, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_EQUAL( 0, quad_all_in_xy_plane( p2, p3, p4, p5, p6, p7, eps ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( QEL, QEL, 0 ), p2, eps );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(   0, QEL, 0 ), p3, eps );
  int midok = quad_mid_edge_nodes_edge_center( p2, p3, p4, p5, p6, p7, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_quad_basic_left_over()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  Vector3D p2, p3, p4, p5, p6, p7;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(QEL, 0.0, 0.0);
  
    // try moving the top left vertex to the right, also move mid-edge nodes proportionally
  get_ideal_quad( p2, p3, p4, p5, p6, p7 );
  p3 -= Vector3D(QEL/2, 0.0, 0.0);
  p6 = 0.5*(p2+p3);
  p7 = 0.5*(p0+p3);
  basic_quad_test( p2, p3, p4, p5, p6, p7, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_EQUAL( 0, quad_all_in_xy_plane( p2, p3, p4, p5, p6, p7, eps ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( QEL, QEL, 0 ), p2, eps );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(   0, QEL, 0 ), p3, eps );
  int midok = quad_mid_edge_nodes_edge_center( p2, p3, p4, p5, p6, p7, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}


// ideal tet is defined to have unit lambda (cube root of det == 1)
// const double IDEAL_TET_SIDE = sqrt(2.0) * pow(3.0, 1.0/3.0); // this was for unit area
const double IDEAL_TET_SIDE = 2.0/pow(33,1.0/6.0);
const double IDEAL_TET_BASE = sqrt(3.0) * 0.5 * IDEAL_TET_SIDE;
const double IDEAL_TET_HEIGHT = sqrt(2.0/3.0) * IDEAL_TET_SIDE;

int tet_mid_edge_nodes_edge_center( const Vector3D& p3,
                                    const Vector3D& p4,
                                    const Vector3D& p5,
                                    const Vector3D& p6,
                                    const Vector3D& p7,
                                    const Vector3D& p8,
                                    const Vector3D& p9,
                                    double epsilon )
{
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  const Vector3D e0 = 0.5*(p0+p1);
  const Vector3D e1 = 0.5*(p1+p2);
  const Vector3D e2 = 0.5*(p2+p0);
  const Vector3D e3 = 0.5*(p0+p3);
  const Vector3D e4 = 0.5*(p1+p3);
  const Vector3D e5 = 0.5*(p2+p3);
  
  int result = 0;
  if ((p4-e0).length() > epsilon)
    result |= 1;
  if ((p5-e1).length() > epsilon)
    result |= 2;
  if ((p6-e2).length() > epsilon)
    result |= 4;
  if ((p7-e3).length() > epsilon)
    result |= 8;
  if ((p8-e4).length() > epsilon)
    result |= 16;
  if ((p9-e5).length() > epsilon)
    result |= 32;
  
  return result;
}

inline static void get_ideal_tet( Vector3D& p3,
                                  Vector3D& p4,
                                  Vector3D& p5,
                                  Vector3D& p6,
                                  Vector3D& p7,
                                  Vector3D& p8,
                                  Vector3D& p9 )
{
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  p3.set( 0.5 * IDEAL_TET_SIDE, IDEAL_TET_BASE/3.0, IDEAL_TET_HEIGHT );
  p4.set(0.5*(p0+p1));
  p5.set(0.5*(p1+p2));
  p6.set(0.5*(p2+p0));
  p7.set(0.5*(p0+p3));
  p8.set(0.5*(p1+p3));
  p9.set(0.5*(p2+p3));
}

void HigherOrderTest::test_tet_basic_ideal()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  const Vector3D p3eq( 0.5 * IDEAL_TET_SIDE, IDEAL_TET_BASE/3.0, IDEAL_TET_HEIGHT );
  Vector3D p3, p4, p5, p6, p7, p8, p9;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  
    // try starting with the optimal result
  get_ideal_tet( p3, p4, p5, p6, p7, p8, p9 );
  basic_tet_test( p3, p4, p5, p6, p7, p8, p9, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_VECTORS_EQUAL( p3eq, p3, eps );
  int midok = tet_mid_edge_nodes_edge_center( p3, p4, p5, p6, p7, p8, p9, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tet_basic_mid_spin()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  const Vector3D p3eq( 0.5 * IDEAL_TET_SIDE, IDEAL_TET_BASE/3.0, IDEAL_TET_HEIGHT );
  Vector3D p3, p4, p5, p6, p7, p8, p9;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  
    // try moving the mid-edge nodes along the edge away from the edge center
  p3 = p3eq;
  p4 = 0.6*p0 + 0.4*p1;
  p5 = 0.6*p1 + 0.4*p2;
  p6 = 0.6*p2 + 0.4*p0;
  p7 = 0.6*p0 + 0.4*p3;
  p8 = 0.4*p1 + 0.6*p3;
  p9 = 0.3*p2 + 0.7*p3;
  basic_tet_test( p3, p4, p5, p6, p7, p8, p9, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_VECTORS_EQUAL( p3eq, p3, eps );
  int midok = tet_mid_edge_nodes_edge_center( p3, p4, p5, p6, p7, p8, p9, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tet_basic_mid_convex()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  const Vector3D p3eq( 0.5 * IDEAL_TET_SIDE, IDEAL_TET_BASE/3.0, IDEAL_TET_HEIGHT );
  Vector3D p3, p4, p5, p6, p7, p8, p9;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  
    // try equilateral corners with all egdes convex
  get_ideal_tet( p3, p4, p5, p6, p7, p8, p9 );
  p4 += Vector3D( 0.0, -0.2, -0.2);
  p5 += Vector3D( 0.2,  0.2, -0.2);
  p6 += Vector3D(-0.2,  0.2, -0.2);
  p7 += Vector3D(-0.2, -0.2,  0.2);
  p8 += Vector3D( 0.2, -0.2,  0.2);
  p9 += Vector3D( 0.0,  0.2,  0.2);
  basic_tet_test( p3, p4, p5, p6, p7, p8, p9, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_VECTORS_EQUAL( p3eq, p3, eps );
  int midok = tet_mid_edge_nodes_edge_center( p3, p4, p5, p6, p7, p8, p9, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tet_basic_apex_down()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  const Vector3D p3eq( 0.5 * IDEAL_TET_SIDE, IDEAL_TET_BASE/3.0, IDEAL_TET_HEIGHT );
  Vector3D p3, p4, p5, p6, p7, p8, p9;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  
    // try moving the top vertex down, also move mid-edge nodes proportionally
  get_ideal_tet( p3, p4, p5, p6, p7, p8, p9 );
  p3 -= Vector3D(0.0,0.0,0.5);
  p7 -= Vector3D(0.0,0.0,0.25);
  p8 -= Vector3D(0.0,0.0,0.25);
  p9 -= Vector3D(0.0,0.0,0.25);
  basic_tet_test( p3, p4, p5, p6, p7, p8, p9, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_VECTORS_EQUAL( p3eq, p3, eps );
  int midok = tet_mid_edge_nodes_edge_center( p3, p4, p5, p6, p7, p8, p9, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tet_basic_apex_up()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  const Vector3D p3eq( 0.5 * IDEAL_TET_SIDE, IDEAL_TET_BASE/3.0, IDEAL_TET_HEIGHT );
  Vector3D p3, p4, p5, p6, p7, p8, p9;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  
    // try moving the top vertex up, also move mid-edge nodes proportionally
  get_ideal_tet( p3, p4, p5, p6, p7, p8, p9 );
  p3 += Vector3D(0.0,0.0,3.0);
  p7 += Vector3D(0.0,0.0,1.5);
  p8 += Vector3D(0.0,0.0,1.5);
  p9 += Vector3D(0.0,0.0,1.5);
  basic_tet_test( p3, p4, p5, p6, p7, p8, p9, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_VECTORS_EQUAL( p3eq, p3, eps );
  int midok = tet_mid_edge_nodes_edge_center( p3, p4, p5, p6, p7, p8, p9, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_tet_basic_apex_over()
{
  MsqPrintError err(cerr);
  const double eps = 5e-2;
  const Vector3D p3eq( 0.5 * IDEAL_TET_SIDE, IDEAL_TET_BASE/3.0, IDEAL_TET_HEIGHT );
  Vector3D p3, p4, p5, p6, p7, p8, p9;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  
    // try moving the top vertex to the right, also move mid-edge nodes proportionally
  get_ideal_tet( p3, p4, p5, p6, p7, p8, p9 );
  p3 -= Vector3D(0.3, 0.0, 0.0);
  p7 = 0.5*(p0 + p3);
  p8 = 0.5*(p1 + p3);
  p9 = 0.5*(p2 + p3);
  basic_tet_test( p3, p4, p5, p6, p7, p8, p9, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  CPPUNIT_ASSERT_VECTORS_EQUAL( p3eq, p3, eps );
  int midok = tet_mid_edge_nodes_edge_center( p3, p4, p5, p6, p7, p8, p9, eps );
  CPPUNIT_ASSERT_EQUAL( 0, midok );
}

void HigherOrderTest::test_hex_basic()
{
  MsqError err;
  const double P = 0.25; // fraction to perturb higher-order nodes: 0.5->optimal
  const double Q = 1-P;
  const double EPSILON = 1e-4;
  const unsigned long num_vtx = 27;

  // Define a 27-node hex with corner nodes fixed and higher-order
  // nodes perturbed a bit.  Smooth to see if HO nodes are moved
  // to center of edge/face/volume.
  double coords[3*num_vtx] = {
    // bottom corners
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0,
    // top corners
    0, 0, 1,
    1, 0, 1,
    1, 1, 1,
    0, 1, 1,
    // bottom mid-edge
    P, 0, 0,
    1, Q, 0,
    Q, 1, 0,
    0, P, 0,
    // vertical mid-edge
    0, 0, P,
    1, 0, P,
    1, 1, Q,
    0, 1, Q,
    // top mid-edge
    Q, 0, 1,
    1, P, 1,
    P, 1, 1,
    0, Q, 1,
    // vertical faces
    P, 0, Q,
    1, Q, P,
    Q, 1, P,
    0, P, Q,
    // bottom and top mid-face
    P, Q, 0,
    Q, P, 1,
    // mid-element
    P, Q, P };
  unsigned long conn[num_vtx];
  for (unsigned i = 0; i < num_vtx; i++)
    conn[i] = i;
  int fixed[num_vtx];
  std::fill( fixed, fixed+8, 1 );
  std::fill( fixed+8, fixed+num_vtx, 0 );
  const EntityTopology type = HEXAHEDRON;
  unsigned long offsets[] = { 0, num_vtx };
  ArrayMesh one_hex( 3, num_vtx, coords, fixed, 1, &type, conn, offsets );
  
    // smooth
  q.run_instructions( &one_hex, err ); ASSERT_NO_ERROR(err);
 
    // test that higher-order nodes were moved to the expected locations
    // bottom mid-edge
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 0.0, 0), Vector3D(coords+3* 8), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1.0, 0.5, 0), Vector3D(coords+3* 9), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 1.0, 0), Vector3D(coords+3*10), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.0, 0.5, 0), Vector3D(coords+3*11), EPSILON );
    // vertical mid-edge
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0, 0, 0.5), Vector3D(coords+3*12), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1, 0, 0.5), Vector3D(coords+3*13), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1, 1, 0.5), Vector3D(coords+3*14), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0, 1, 0.5), Vector3D(coords+3*15), EPSILON );
    // top mid-edge
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 0.0, 1), Vector3D(coords+3*16), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1.0, 0.5, 1), Vector3D(coords+3*17), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 1.0, 1), Vector3D(coords+3*18), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.0, 0.5, 1), Vector3D(coords+3*19), EPSILON );
    // vertical faces
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 0.0, 0.5), Vector3D(coords+3*20), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1.0, 0.5, 0.5), Vector3D(coords+3*21), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 1.0, 0.5), Vector3D(coords+3*22), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.0, 0.5, 0.5), Vector3D(coords+3*23), EPSILON );
    // bottom and top mid-face
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 0.5, 0), Vector3D(coords+3*24), EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 0.5, 1), Vector3D(coords+3*25), EPSILON );
    // mid-element
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0.5, 0.5, 0.5), Vector3D(coords+3*26), EPSILON );
}

/*
    Jason,
 
        Here is a very simple test we should run.
 
       A single quadratic triangle.   IMR metric with W = equilateral ideal    No geometry.
 
      Nodes:
 
         (x,y)            Type               F/S/F            
       -----------------------------------------------------------
        (0,0)           Corner               Fixed
        (1,0)           Corner               Fixed
        (x2,y2)         Corner                Free
        (x3,y3)         Mid                   Free
        (x4,y4)         Mid                   Free
        (x5,y5)         Mid                   Free
       --------------------------------------------------------
 
       Start with sample points at each of the node.
       The optimal element should be equilateral.
       Try different initial values for x2,y2,x3, etc.
       Try only including some of the sample points to see if still get equilateral.
       Probably should add this one to the testSuite.
 
   -- Pat
*/
void HigherOrderTest::basic_tri_test( double& x2, double& y2,
                                      double& x3, double& y3,
                                      double& x4, double& y4,
                                      double& x5, double& y5,
                                      MsqError& err )
{
    // Construct Mesh instance
  const int DIM = 3;
  const int NVTX = 6;
  const int NELEM = 1;
  double coords[DIM*NVTX] = { 0., 0., 0,
                              IDEAL_TRI_SIDE, 0., 0,
                              x2, y2, 0,
                              x3, y3, 0,
                              x4, y4, 0,
                              x5, y5, 0 };
  const int fixed[NVTX] = { 1, 1, 0, 0, 0, 0 };
  const unsigned long conn[NVTX*NELEM] = { 0, 1, 2, 3, 4, 5, };
  ArrayMesh mesh( DIM, NVTX, coords, fixed,
                  NELEM, TRIANGLE, conn, false, NVTX );
  PlanarDomain xy( PlanarDomain::XY );
  
    // Solve
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &xy);
  q.run_instructions( &mesh_and_domain, err ); MSQ_ERRRTN(err);

    // Pass back modified coordinates
  x2 = coords[2*DIM+0]; y2 = coords[2*DIM+1];
  x3 = coords[3*DIM+0]; y3 = coords[3*DIM+1];
  x4 = coords[4*DIM+0]; y4 = coords[4*DIM+1];
  x5 = coords[5*DIM+0]; y5 = coords[5*DIM+1];
}

void HigherOrderTest::basic_tet_test( Vector3D& p3,
                                      Vector3D& p4,
                                      Vector3D& p5,
                                      Vector3D& p6,
                                      Vector3D& p7,
                                      Vector3D& p8,
                                      Vector3D& p9,
                                      MsqError& err )
{
    // Construct Mesh instance
  const int DIM = 3;
  const int NVTX = 10;
  const int NELEM = 1;
  const Vector3D p0(0.0, 0.0, 0.0);
  const Vector3D p1(IDEAL_TET_SIDE, 0.0, 0.0);
  const Vector3D p2(0.5*IDEAL_TET_SIDE, IDEAL_TET_BASE, 0.0 );
  double coords[DIM*NVTX] = { p0.x(), p0.y(), p0.z(),
                              p1.x(), p1.y(), p1.z(),
                              p2.x(), p2.y(), p2.z(),
                              p3.x(), p3.y(), p3.z(),
                              p4.x(), p4.y(), p4.z(),
                              p5.x(), p5.y(), p5.z(),
                              p6.x(), p6.y(), p6.z(),
                              p7.x(), p7.y(), p7.z(),
                              p8.x(), p8.y(), p8.z(),
                              p9.x(), p9.y(), p9.z() };
                              
  const int fixed[NVTX] = { 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 };
  const unsigned long conn[NVTX*NELEM] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  ArrayMesh mesh( DIM, NVTX, coords, fixed,
                  NELEM, TETRAHEDRON, conn, false, NVTX );

  
    // Solve
  q.run_instructions( &mesh, err ); MSQ_ERRRTN(err);

    // Pass back modified coordinates
  p3.set( coords + 3*DIM );
  p4.set( coords + 4*DIM );
  p5.set( coords + 5*DIM );
  p6.set( coords + 6*DIM );
  p7.set( coords + 7*DIM );
  p8.set( coords + 8*DIM );
  p9.set( coords + 9*DIM );
}

void HigherOrderTest::basic_quad_test( Vector3D& p2,
                                       Vector3D& p3,
                                       Vector3D& p4,
                                       Vector3D& p5,
                                       Vector3D& p6,
                                       Vector3D& p7,
                                       MsqError& err )
{
    // Construct Mesh instance
  const int DIM = 3;
  const int NVTX = 8;
  const int NELEM = 1;
  double coords[DIM*NVTX] = {    0.0,    0.0,    0.0,
                                 QEL,    0.0,    0.0,
                              p2.x(), p2.y(), p2.z(),
                              p3.x(), p3.y(), p3.z(),
                              p4.x(), p4.y(), p4.z(),
                              p5.x(), p5.y(), p5.z(),
                              p6.x(), p6.y(), p6.z(),
                              p7.x(), p7.y(), p7.z() };
                              
  const int fixed[NVTX] = { 1, 1, 0, 0, 0, 0, 0, 0 };
  const unsigned long conn[NVTX*NELEM] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  ArrayMesh mesh( DIM, NVTX, coords, fixed,
                  NELEM, QUADRILATERAL, conn, false, NVTX );
  PlanarDomain xy( PlanarDomain::XY );

  
    // Solve
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &xy);
  q.run_instructions( &mesh_and_domain, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());

    // Pass back modified coordinates
  p2.set( coords + 2*DIM );
  p3.set( coords + 3*DIM );
  p4.set( coords + 4*DIM );
  p5.set( coords + 5*DIM );
  p6.set( coords + 6*DIM );
  p7.set( coords + 7*DIM );
}

void HigherOrderTest::basic_hex_test( Vector3D& p4,
                                      Vector3D& p5,
                                      Vector3D& p6,
                                      Vector3D& p7,
                                      Vector3D& p8,
                                      Vector3D& p9,
                                      Vector3D& p10,
                                      Vector3D& p11,
                                      Vector3D& p12,
                                      Vector3D& p13,
                                      Vector3D& p14,
                                      Vector3D& p15,
                                      Vector3D& p16, 
                                      Vector3D& p17,
                                      Vector3D& p18,
                                      Vector3D& p19,
                                      MsqError& err )
{
    // Construct Mesh instance
  const int DIM = 3;
  const int NVTX = 20;
  const int NELEM = 1;
  double coords[DIM*NVTX] = { 0.0, 0.0, 0,
                              2.0, 0.0, 0,
                              2.0, 2.0, 0,
                              0.0, 2.0, 0,
                               p4.x(), p4.y(), p4.z(),
                               p5.x(), p5.y(), p5.z(),
                               p6.x(), p6.y(), p6.z(),
                               p7.x(), p7.y(), p7.z(),
                               p8.x(), p8.y(), p8.z(),
                               p9.x(), p9.y(), p9.z(),
                              p10.x(), p10.y(), p10.z(),
                              p11.x(), p11.y(), p11.z(),
                              p12.x(), p12.y(), p12.z(),
                              p13.x(), p13.y(), p13.z(),
                              p14.x(), p14.y(), p14.z(),
                              p15.x(), p15.y(), p15.z(),
                              p16.x(), p16.y(), p16.z(),
                              p17.x(), p17.y(), p17.z(),
                              p18.x(), p18.y(), p18.z(),
                              p19.x(), p19.y(), p19.z() };
  const int fixed[NVTX] = { 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  const unsigned long conn[NVTX*NELEM] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                          10,11,12,13,14,15,16,17,18,19 };
  ArrayMesh mesh( DIM, NVTX, coords, fixed,
                  NELEM, HEXAHEDRON, conn, false, NVTX );
  
    // Solve
  q.run_instructions( &mesh, err ); MSQ_ERRRTN(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());

    // Pass back modified coordinates
  p4.set( coords + 4*DIM );
  p5.set( coords + 5*DIM );
  p6.set( coords + 6*DIM );
  p7.set( coords + 7*DIM );
  p8.set( coords + 8*DIM );
  p9.set( coords + 9*DIM );
  p10.set( coords + 10*DIM );
  p11.set( coords + 11*DIM );
  p12.set( coords + 12*DIM );
  p13.set( coords + 13*DIM );
  p14.set( coords + 14*DIM );
  p15.set( coords + 15*DIM );
  p16.set( coords + 16*DIM );
  p17.set( coords + 17*DIM );
  p18.set( coords + 18*DIM );
  p19.set( coords + 19*DIM );
}

/*       The mesh consists of a single triangle with 3 corner
         and 3 mid-side nodes.  Two sides of the triangle are
         on the boundary, which consists of two curves (maybe
         just straight lines), so one side is on each curve.
         The corner vertex at the intersection of the two curves is
         fixed, the two mid-side nodes on the boundary are
         free (with snap-to), one of the other two corner vertices
         is fixed, the other free (with snap-to).  The remaining
         mid-node vertex is free, being in the interior.  To close
         the domain, I leave up to you...

      ----2*
      |   |
      |   |
      | --5       4
     y2 | |
      |y5 |
      | | |
      ----0*------3-----------1-----
          |<--x3->|           |
          |<-------x1-------->|
          
    * Fixed vertices
          
*/
void HigherOrderTest::test_tri_open_domain( double& x1, double& x3, double& x4,
                                            double  y2, double  y5, double& y4 )
{
  MsqPrintError err(cerr);

    // Validate input
  CPPUNIT_ASSERT(x3 < x1);
  CPPUNIT_ASSERT(y5 < y2);
  CPPUNIT_ASSERT(x3 > 0.0);
  CPPUNIT_ASSERT(y5 > 0.0);

    // Construct Mesh instance
  const int DIM = 3;
  const int NVTX = 6;
  const int NELEM = 1;
  double coords[DIM*NVTX] = { 0, 0, 0,
                             x1, 0, 0,
                              0,y2, 0,
                             x3, 0, 0,
                             x4,y4, 0,
                              0,y5, 0 };
  const int fixed[NVTX] = { 1, 0, 1, 0, 0, 0 };
  const unsigned long conn[NVTX*NELEM] = { 0, 1, 2, 3, 4, 5, };
  ArrayMesh mesh( DIM, NVTX, coords, fixed,
                  NELEM, TRIANGLE, conn, false, NVTX );

    // Construct Geometry
  LineDomain xaxis( Vector3D(0,0,0), Vector3D(1,0,0) );
  LineDomain yaxis( Vector3D(0,0,0), Vector3D(0,1,0) );
  PlanarDomain xyplane( PlanarDomain::XY );
  DomainClassifier::DomainSet sets[] = { &xaxis, &yaxis, &xyplane };
  
  vector<Mesh::VertexHandle> verts;
  vector<Mesh::ElementHandle> elems;
  mesh.get_all_vertices( verts, err ); ASSERT_NO_ERROR(err);
  mesh.get_all_elements( elems, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, elems.size() );
  
    // Associate mesh with geometry
  sets[0].vertices.push_back( verts[1] );
  sets[0].vertices.push_back( verts[3] );
  sets[1].vertices.push_back( verts[2] );
  sets[1].vertices.push_back( verts[5] );
  sets[2].elements.push_back( elems[0] );
  sets[2].vertices.push_back( verts[4] );
  DomainClassifier geom;
  DomainClassifier::classify_by_handle( geom, &mesh, sets, 3, err );
  ASSERT_NO_ERROR(err);
  
    // Solve
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &geom);
  q.run_instructions( &mesh_and_domain, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
 
    // Pass back modified coordinate values
  x1 = coords[DIM*1+0];
  x3 = coords[DIM*3+0];
  x4 = coords[DIM*4+0];
  y4 = coords[DIM*4+1];
}

void HigherOrderTest::test_tri_open_domain( )
{
  double x1 = 2;
  double y2 = 2;
  double x3 = 1;
  double y5 = 1;
  double x4 = 1;
  double y4 = 1;
  test_tri_open_domain( x1, x3, x4, y2, y5, y4 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, x1, 1e-3 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, y2, 1e-3 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, x3, 1e-3 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, x4, 1e-3 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, y4, 1e-3 );
}
  
/*
   Jason,
 
        Here is an interesting test we should run that's slightly more like the SLAC problem.
 
       A single quadratic triangle.   IMR metric with W = equilateral ideal    No geometry.
 
      Nodes:
 
         (x,y)            Type               F/S/F            
       -----------------------------------------------------------
        (0,0)           Corner               Fixed
        (1,0)           Corner               Fixed
        (0.5,10)       Corner               Free
        (0.50,0.49)     Mid                 Fixed
        (0.75,0.5)       Mid                Free
        (0.25,0.5)       Mid                Free
       --------------------------------------------------------
 
       Start with sample points at each of the node.
 
       The Jacobian at the fixed mid-side node is nearly zero
         and the goal is to see if this can be improved by moving the 3 free nodes.
 
         The node at (0.5,1.0) should move downward\u2026. Moving it upward makes
           the mid-node Jacobian negative, and the IMR barrier should prevent that.
 
       Try including only some of the sample points to see if all are needed.
 
   -- Pat
*/
void HigherOrderTest::test_tri_slac()
{
  MsqPrintError err(cerr);

    // Construct Mesh instance
  const int DIM = 3;
  const int NVTX = 6;
  const int NELEM = 1;
  double coords[DIM*NVTX] = { 0.00, 0.00, 0,
                              1.00, 0.00, 0,
                              0.50,10.0,  0,
                              0.50, 0.30, 0,
                              0.75, 4.50, 0,
                              0.25, 4.50, 0 };
  const int fixed[NVTX] = { 1, 1, 0, 1, 0, 0 };
  const unsigned long conn[NVTX*NELEM] = { 0, 1, 2, 3, 4, 5, };
  ArrayMesh mesh( DIM, NVTX, coords, fixed,
                  NELEM, TRIANGLE, conn, false, NVTX );
  PlanarDomain xy( PlanarDomain::XY );
  
    // Solve
  q.run_instructions( &mesh, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!hit_iteration_limit());
  
  const Vector3D v0(coords+0*DIM);
  const Vector3D v1(coords+1*DIM);
  const Vector3D v2(coords+2*DIM);
  const Vector3D v3(coords+3*DIM);
  const Vector3D v4(coords+4*DIM);
  const Vector3D v5(coords+5*DIM);
  
    // Expect vertex 2 to have moved downwards signficantly
  CPPUNIT_ASSERT( v2.y() < 2.0 );
  CPPUNIT_ASSERT( v2.y() > 0.5 );
    // Expect it to have stayed in the center along the X axis
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, v2.x(), 1e-3 );
    // Expect free mid-edge nodes to be "between" adjacent corners
  CPPUNIT_ASSERT( v4.x() > 0.5 );
  CPPUNIT_ASSERT( v4.x() < 1.0 );
  CPPUNIT_ASSERT( v5.x() > 0.0 );
  CPPUNIT_ASSERT( v5.x() < 0.5 );
    // Expect two free mid-edge nodes to be slightly above the
    // line between adjacent corners (slightly convex edge)
    // because the fixed mid-edge node on the bottom is moved upward
  Vector3D m12 = v2 - v1;
  Vector3D m02 = v2 - v0;
  double t1 = (m12 % (v4 - v0)) / (m12 % m12);
  double t2 = (m02 % (v5 - v0)) / (m02 % m02);
  Vector3D l4 = v0 * t1 * m12; // closed point to v4 on line v1,v2
  Vector3D l5 = v0 * t2 * m02; // closed point to v5 on line v0,v2
    // Check not to far from line
  CPPUNIT_ASSERT( (l4 - v4).length() < 0.1 );
  CPPUNIT_ASSERT( (l5 - v5).length() < 0.1 );
    // Check that edges are convex
  CPPUNIT_ASSERT( l4.x() < v4.x() );
  CPPUNIT_ASSERT( l5.x() > v5.x() );
}

