/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */

#include "QualityMetricTester.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_QualityMetric.hpp"
#include "Mesquite_IdealElements.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_ElemSampleQM.hpp"
#include "Mesquite_TopologyInfo.hpp"
#include "Mesquite_EdgeQM.hpp"

#include <cppunit/extensions/HelperMacros.h>

static std::vector<EntityTopology> types_in_group( QualityMetricTester::ElemTypeGroup group )
{
  std::vector<EntityTopology> types;
  switch (group) {
    default:
    case QualityMetricTester::ALL:
      types.push_back( POLYHEDRON );
      types.push_back( POLYGON );
    case QualityMetricTester::ALL_FE:
      types.push_back( SEPTAHEDRON );
    case QualityMetricTester::ALL_FE_EXCEPT_SEPTAHEDRON:
      types.push_back( PYRAMID );
      types.push_back( PRISM );
    case QualityMetricTester::NON_MIXED_FE:
      types.push_back( HEXAHEDRON );
      types.push_back( QUADRILATERAL );
    case QualityMetricTester::SIMPLICIES:
      types.push_back(TETRAHEDRON);
      types.push_back(TRIANGLE);
      break;
      
    case QualityMetricTester::THREE_D:
      types.push_back( POLYHEDRON );
    case QualityMetricTester::THREE_D_FE:
      types.push_back( SEPTAHEDRON );
    case QualityMetricTester::THREE_D_FE_EXCEPT_SEPTAHEDRON:
      types.push_back( PYRAMID );
      types.push_back( PRISM );
    case QualityMetricTester::THREE_D_NON_MIXED_FE:
      types.push_back( HEXAHEDRON );
      types.push_back(TETRAHEDRON);
      break;
      
    case QualityMetricTester::TWO_D:
      types.push_back( POLYGON );
    case QualityMetricTester::TWO_D_FE:
      types.push_back( TRIANGLE );
      types.push_back( QUADRILATERAL );
      break;
  }
  std::reverse( types.begin(), types.end() );
  return types;
}

QualityMetricTester::QualityMetricTester( 
                     const EntityTopology* supported_elem_types,
                     size_t len,
                     const Settings* settings )
  : degenHexPyramid(false), 
    types( len ), 
    mSettings( settings ), 
    geomPlane( Vector3D(0,0,1), Vector3D(0,0,0) )
{
  std::copy( supported_elem_types, supported_elem_types+len, types.begin() );
}

QualityMetricTester::QualityMetricTester( 
                     ElemTypeGroup group,
                     const Settings* settings )
  : degenHexPyramid(false), 
    types( types_in_group(group) ),
    mSettings( settings ), 
    geomPlane( Vector3D(0,0,1), Vector3D(0,0,0) )
{
}

void QualityMetricTester::get_ideal_tris( PatchData& pd, bool unit_area )
{
  //      6 ------- 5     .
  //       /\     /\      .
  //      /  \   /  \     .
  //     /    \ /    \    .
  //  1 <------X------> 4 .
  //     \    /0\    /    .
  //      \  /   \  /     .
  //       \/     \/      .
  //      2 ------- 3     .
  static const double y3 = MSQ_SQRT_THREE_DIV_TWO;
  double ideal_tri_verts[] = { 0.0, 0.0, 0.0,
                              -1.0, 0.0, 0.0,
                              -0.5, -y3, 0.0,
                               0.5, -y3, 0.0,
                               1.0, 0.0, 0.0,
                               0.5,  y3, 0.0,
                              -0.5,  y3, 0.0 };
  const size_t ideal_tri_elems[] = { 0, 1, 2, 
                                     0, 2, 3, 
                                     0, 4, 5, 
                                     0, 4, 5, 
                                     0, 5, 6, 
                                     0, 6, 1 };
  const bool fixed[] = { false, true, true, true, true, true, true };
  
  if (unit_area) {
    const double unit_tri_scale = 2.0/std::sqrt(6.0);
    for (size_t i = 0; i < sizeof(ideal_tri_verts)/sizeof(double); ++i)
      ideal_tri_verts[i] *= unit_tri_scale;
  }
  
  pd.set_domain( &geomPlane );
  
  MsqPrintError err( std::cout );
  pd.fill( 7, ideal_tri_verts, 6, TRIANGLE, ideal_tri_elems, fixed, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
}
  
void QualityMetricTester::get_ideal_quads( PatchData& pd )
{
  //            7
  //   8 +------+------+ 6
  //     |      |      |
  //     |      |      |
  //     |      |0     |
  //   1 +------+------+ 5
  //     |      |      |
  //     |      |      |
  //     |      |      |
  //   2 +------+------+ 4
  //            3
  static const double ideal_quad_verts[] = { 0.0, 0.0, 0.0,
                                            -1.0, 0.0, 0.0,
                                            -1.0,-1.0, 0.0,
                                             0.0,-1.0, 0.0,
                                             1.0,-1.0, 0.0,
                                             1.0, 0.0, 0.0,
                                             1.0, 1.0, 0.0,
                                             0.0, 1.0, 0.0,
                                            -1.0, 1.0, 0.0 };
  static const size_t ideal_quad_elems[] = { 0, 1, 2, 3,
                                             0, 3, 4, 5, 
                                             0, 5, 6, 7, 
                                             0, 7, 8, 1 };
  const bool fixed[] = { false, true, true, true, true, true, true, true, true };
  
  pd.set_domain( &geomPlane );
  
  MsqPrintError err( std::cout );
  pd.fill( 9, ideal_quad_verts, 4, QUADRILATERAL, ideal_quad_elems, fixed, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
}

  
void QualityMetricTester::get_ideal_hexes( PatchData& pd )
{
  static const double ideal_hex_verts[] = {  0.0, 0.0, 0.0,
                                            -1.0, 0.0, 0.0,
                                            -1.0,-1.0, 0.0,
                                             0.0,-1.0, 0.0,
                                             1.0,-1.0, 0.0,
                                             1.0, 0.0, 0.0,
                                             1.0, 1.0, 0.0,
                                             0.0, 1.0, 0.0,
                                            -1.0, 1.0, 0.0, 
                                             0.0, 0.0,-1.0,
                                            -1.0, 0.0,-1.0,
                                            -1.0,-1.0,-1.0,
                                             0.0,-1.0,-1.0,
                                             1.0,-1.0,-1.0,
                                             1.0, 0.0,-1.0,
                                             1.0, 1.0,-1.0,
                                             0.0, 1.0,-1.0,
                                            -1.0, 1.0,-1.0, 
                                             0.0, 0.0, 1.0,
                                            -1.0, 0.0, 1.0,
                                            -1.0,-1.0, 1.0,
                                             0.0,-1.0, 1.0,
                                             1.0,-1.0, 1.0,
                                             1.0, 0.0, 1.0,
                                             1.0, 1.0, 1.0,
                                             0.0, 1.0, 1.0,
                                            -1.0, 1.0, 1.0 };
  static const size_t ideal_hex_elems[] = {  9, 10, 11, 12, 0, 1, 2, 3,
                                             9, 12, 13, 14, 0, 3, 4, 5,
                                             9, 14, 15, 16, 0, 5, 6, 7,
                                             9, 16, 17, 10, 0, 7, 8, 1,
                                             0, 1, 2, 3, 18, 19, 20, 21,
                                             0, 3, 4, 5, 18, 21, 22, 23,
                                             0, 5, 6, 7, 18, 23, 24, 25, 
                                             0, 7, 8, 1, 18, 25, 26, 19 };
  bool fixed[27];
  fixed[0] = false;
  for (size_t f = 1; f < 27; ++f)
    fixed[f] = true;
  
  MsqPrintError err( std::cout );
  pd.fill( 27, ideal_hex_verts, 8, HEXAHEDRON, ideal_hex_elems, fixed, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
}

// free_vertex_index:
//  >= 0 : all vertice except specified one are fixed
//    -1 : all vertices are free
//    -2 : first vertex is fixed, all others are free
void QualityMetricTester::get_ideal_element( EntityTopology type, 
                                             bool unit_area,
                                             PatchData& pd,
                                             int free_vertex_index )
{
  if (TopologyInfo::dimension(type) == 2)
    pd.set_domain( &geomPlane );
  
  const size_t n = TopologyInfo::corners(type);
  const Vector3D* coords = unit_area ? unit_element(type,degenHexPyramid) : unit_edge_element(type,degenHexPyramid);
  CPPUNIT_ASSERT(coords != 0);
  
  CPPUNIT_ASSERT(sizeof(double)*3 == sizeof(Vector3D));
  const double* elem_verts = reinterpret_cast<const double*>(coords);
  
  bool* fixed = new bool[n];
  std::vector<size_t> conn(n);
  for (size_t i = 0; i < n; ++i) {
    conn[i] = i;
    fixed[i] = (free_vertex_index >= 0);
  }
  if (free_vertex_index == -2)
    fixed[0] = true;
  else if (free_vertex_index >= 0)
    fixed[free_vertex_index] = false;

  MsqPrintError err( std::cout );
  pd.fill( n, elem_verts, 1, type, arrptr(conn), fixed, err );
  delete [] fixed;
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
}
void QualityMetricTester::get_ideal_element( EntityTopology type, 
                                             bool unit_area,
                                             PatchData& pd,
                                             bool first_vertex_fixed )
{
  get_ideal_element( type, unit_area, pd, first_vertex_fixed ? -2 : -1 );
}

/** Begin with an ideal element, and move the first vertex one half
  * of the distnace to the centroid */
void QualityMetricTester::get_nonideal_element( EntityTopology type, 
                                                PatchData& pd, 
                                                int free_vtx_type )
{
  MsqError err;
  get_ideal_element( type, false, pd, free_vtx_type );
  const size_t* verts = pd.element_by_index(0).get_vertex_index_array();
  const MsqVertex* coords = pd.get_vertex_array(err);
  size_t n = pd.element_by_index(0).vertex_count();
  
  Vector3D sum(0,0,0);
  for (unsigned i = 0; i < n; ++i)
    sum += coords[verts[i]];
  
  sum /= n;
  sum += coords[verts[0]];
  sum *= 0.5;
  pd.set_vertex_coordinates( sum, verts[0], err );
}

void QualityMetricTester::get_nonideal_element( EntityTopology type, 
                                                PatchData& pd, 
                                                bool first_vertex_fixed )
{
  get_nonideal_element( type, pd, first_vertex_fixed ? -2 : -1 );
}

/** Collapse first and second vertices of an ideal element */
void QualityMetricTester::get_degenerate_element( EntityTopology type, PatchData& pd )
{
  MsqError err;
  get_ideal_element( type, false, pd );
  const size_t* verts = pd.element_by_index(0).get_vertex_index_array();
  const MsqVertex* coords = pd.get_vertex_array(err);
  pd.set_vertex_coordinates( coords[verts[1]], verts[0], err );
}

/** Create element with zero area/volume
  * For quads and hexes, the results are also inverted elements.
  */
void QualityMetricTester::get_zero_element( EntityTopology type, PatchData& pd )
{
  MsqError err;
  get_ideal_element( type, false, pd );
  MsqMeshEntity& elem = pd.element_by_index(0);
  const size_t* verts = elem.get_vertex_index_array();
  const MsqVertex* coords = pd.get_vertex_array(err);
  unsigned i;
  Vector3D sum(0,0,0);
 
  switch (type) {
    case TRIANGLE:
      pd.set_vertex_coordinates( 0.5*(coords[verts[0]]+coords[verts[1]]), verts[2], err );
      break;
    case QUADRILATERAL:
      pd.set_vertex_coordinates( coords[verts[0]], verts[1], err );
      pd.set_vertex_coordinates( coords[verts[3]], verts[2], err );
      break;
    case TETRAHEDRON:
      for (i = 0; i < 3; ++i)
        sum += coords[verts[i]];
      pd.set_vertex_coordinates( sum/3.0, verts[3], err );
      break;
    case PYRAMID:
      for (i = 0; i < 4; ++i)
        sum += coords[verts[i]];
      pd.set_vertex_coordinates( 0.25*sum, verts[4], err );
      break;
    case PRISM:
      pd.set_vertex_coordinates( 0.5*(coords[verts[0]]+coords[verts[1]]), verts[2], err );
      pd.set_vertex_coordinates( 0.5*(coords[verts[3]]+coords[verts[4]]), verts[5], err );
      break;
    case HEXAHEDRON:
      pd.set_vertex_coordinates( coords[verts[0]], verts[4], err );
      pd.set_vertex_coordinates( coords[verts[1]], verts[5], err );
      pd.set_vertex_coordinates( coords[verts[2]], verts[6], err );
      pd.set_vertex_coordinates( coords[verts[3]], verts[7], err );
      break;
    default:
      CPPUNIT_ASSERT(false);
      break;
  }
}

/** Create inverted elements.
  * For tri and tet elements, reflect one vertex about the opposite side.
  * For all other elements, introduce a concave vertex in one of the 
  * quadrilateral faces.
  */
void QualityMetricTester::get_inverted_element( EntityTopology type, PatchData& pd )
{
  MsqError err;
  get_ideal_element( type, false, pd );
  MsqMeshEntity& elem = pd.element_by_index(0);
  const size_t* verts = elem.get_vertex_index_array();
  const MsqVertex* coords = pd.get_vertex_array(err);
  unsigned i;
  Vector3D sum(0,0,0);
  
//  switch (type) {
//    case TRIANGLE:
//      coords[verts[2]] = coords[verts[0]] + coords[verts[1]] - coords[verts[2]];
//      break;
//    case QUADRILATERAL:
//    case PYRAMID:
//    case HEXAHEDRON:
//      coords[verts[2]] += 0.75 * (coords[verts[0]] - coords[verts[2]]);
//      break;
//    case TETRAHEDRON:
//      for (i = 0; i < 3; ++i)
//        sum += coords[verts[i]];
//      coords[verts[3]] += 0.5*sum - 1.5*coords[verts[3]];
//      break;
//    case PRISM:
//      coords[verts[4]] += 0.75 * (coords[verts[0]] - coords[verts[4]]);
//      break;
//    default:
//      CPPUNIT_ASSERT(false);
//      break;
//  }
  if (type == TRIANGLE)
      pd.set_vertex_coordinates( coords[verts[0]] + coords[verts[1]] - coords[verts[2]], verts[2], err );
  else if (type == QUADRILATERAL || type == PYRAMID || type == HEXAHEDRON)
      pd.set_vertex_coordinates( 0.75 * (coords[verts[0]] - coords[verts[2]]), verts[2], err );
  else if (type == TETRAHEDRON) {
      for (i = 0; i < 3; ++i)
        sum += coords[verts[i]];
      pd.set_vertex_coordinates( 0.5*sum - 1.5*coords[verts[3]], verts[3], err );
  }
  else if (type == PRISM) 
      pd.set_vertex_coordinates( 0.75 * (coords[verts[0]] - coords[verts[4]]), verts[4], err );
  else
      CPPUNIT_ASSERT(false);
  
  if (TopologyInfo::dimension(type) == 3 || pd.domain_set()) {
    int inverted, total;
    elem.check_element_orientation(pd,inverted,total,err);
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(total);
  }
}

void QualityMetricTester::test_type_is_supported( EntityTopology type, QualityMetric* qm )
{
  MsqPrintError err(std::cout);
  bool rval;
    // create patch containing only elements of sepcified type
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  get_ideal_element( type, false, pd );
    // get list of evaluation locations in patch
  std::vector<size_t> handles;
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!handles.empty());
    // evaluate the metric at one location
  double value = -HUGE_VAL;
  rval = qm->evaluate( pd, handles[0], value, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(rval);
    // if metric returned input value for 'value' try again
    // with a different initial value to verify that the metric
    // is setting 'value' to something.
  if (value == -HUGE_VAL) {
    value = 0.0;
    rval = qm->evaluate( pd, handles[0], value, err );
    CPPUNIT_ASSERT(!err && rval);
    CPPUNIT_ASSERT( value != 0.0 );
  }
  std::vector<size_t> indices;
  rval = qm->evaluate_with_indices( pd, handles[0], value, indices, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err)); CPPUNIT_ASSERT(rval);
  std::vector<Vector3D> grad;
  rval = qm->evaluate_with_gradient( pd, handles[0], value, indices, grad, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err)); CPPUNIT_ASSERT(rval);
  std::vector<Matrix3D> hess;
  rval = qm->evaluate_with_Hessian( pd, handles[0], value, indices, grad, hess, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err)); CPPUNIT_ASSERT(rval);
}

void QualityMetricTester::test_type_is_not_supported( EntityTopology type, QualityMetric* qm )
{
  MsqError err;
  bool rval;
    // create patch containing only elements of sepcified type
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  get_ideal_element( type, false, pd );
    // get list of evaluation locations in patch
  std::vector<size_t> handles;
  qm->get_evaluations( pd, handles, false, err );
  if (err.error_code() == MsqError::UNSUPPORTED_ELEMENT)
    return;
  CPPUNIT_ASSERT(!err);
  if (handles.empty())
    return;
    // evaluate the metric at one location
  double value;
  rval = qm->evaluate( pd, handles[0], value, err );
  CPPUNIT_ASSERT_EQUAL(MsqError::UNSUPPORTED_ELEMENT, err.error_code());
  std::vector<size_t> indices;
  rval = qm->evaluate_with_indices( pd, handles[0], value, indices, err );
  CPPUNIT_ASSERT_EQUAL(MsqError::UNSUPPORTED_ELEMENT, err.error_code());
  std::vector<Vector3D> grad;
  rval = qm->evaluate_with_gradient( pd, handles[0], value, indices, grad, err );
  CPPUNIT_ASSERT_EQUAL(MsqError::UNSUPPORTED_ELEMENT, err.error_code());
  std::vector<Matrix3D> hess;
  rval = qm->evaluate_with_Hessian( pd, handles[0], value, indices, grad, hess, err );
  CPPUNIT_ASSERT_EQUAL(MsqError::UNSUPPORTED_ELEMENT, err.error_code());
}

void QualityMetricTester::test_supported_element_types( QualityMetric* qm )
{
  for (int i = TRIANGLE; i < MIXED; ++i) {
    if (i == POLYGON || i == POLYHEDRON || i == SEPTAHEDRON)
      continue;
    else if (type_is_supported( (EntityTopology)i ))
      test_type_is_supported( (EntityTopology)i, qm );
    else
      test_type_is_not_supported( (EntityTopology)i, qm );
  }
}

static void test_evaluate( PatchData& pd, QualityMetric* qm, double value )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles;
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!handles.empty());
  for (unsigned i = 0; i < handles.size(); ++i) {
    double qmval;
    bool rval = qm->evaluate( pd, handles[i], qmval, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT(rval);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( value, qmval, 1e-6 );
  }
}

void QualityMetricTester::test_evaluate_unit_element( QualityMetric* qm, EntityTopology type, double value )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  get_ideal_element( type, true, pd );
  test_evaluate( pd, qm, value );
}

void QualityMetricTester::test_evaluate_unit_edge_element( QualityMetric* qm, EntityTopology type, double value )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  get_ideal_element( type, false, pd );
  test_evaluate( pd, qm, value );
}

void QualityMetricTester::test_evaluate_unit_tris_about_vertex( QualityMetric* qm, double value )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  get_ideal_tris( pd, true );
  test_evaluate( pd, qm, value );
}

void QualityMetricTester::test_evaluate_unit_quads_about_vertex( QualityMetric* qm, double value )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  get_ideal_quads( pd );
  test_evaluate( pd, qm, value );
}

void QualityMetricTester::test_evaluate_unit_hexes_about_vertex( QualityMetric* qm, double value )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  get_ideal_hexes( pd );
  test_evaluate( pd, qm, value );
}

void QualityMetricTester::test_evaluate_unit_edge_tris_about_vertex( QualityMetric* qm, double value )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  get_ideal_tris( pd, false );
  test_evaluate( pd, qm, value );
}

void QualityMetricTester::test_get_element_evaluations( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  unsigned count = 0;
  std::vector<size_t> handles;

  if (type_is_supported(HEXAHEDRON)) {
    get_ideal_hexes( pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
    std::sort( handles.begin(), handles.end() );
    handles.erase( std::unique( handles.begin(), handles.end() ), handles.end() );
    CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
    ++count;
  }

  if (type_is_supported(QUADRILATERAL)) {
    get_ideal_quads( pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
    std::sort( handles.begin(), handles.end() );
    handles.erase( std::unique( handles.begin(), handles.end() ), handles.end() );
    CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
    ++count;
  }

  if (type_is_supported(TRIANGLE)) {
    get_ideal_tris( pd, false );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
    std::sort( handles.begin(), handles.end() );
    handles.erase( std::unique( handles.begin(), handles.end() ), handles.end() );
    CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
    ++count;
  }
  
  CPPUNIT_ASSERT(count > 0);
}

void QualityMetricTester::test_get_vertex_evaluations( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  unsigned count = 0;
  std::vector<size_t> handles;

  if (type_is_supported(HEXAHEDRON)) {
    get_ideal_hexes( pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
    ++count;
  }

  if (type_is_supported(TRIANGLE)) {
    get_ideal_tris( pd, false );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
    ++count;
  }

  if (type_is_supported(QUADRILATERAL)) {
    get_ideal_quads( pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
    ++count;
  }
  
  CPPUNIT_ASSERT(count > 0);
}

void QualityMetricTester::test_get_sample_evaluations( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  unsigned count = 0;
  std::vector<size_t> handles;
  
  if (type_is_supported(HEXAHEDRON)) {
    get_ideal_hexes( pd );
    size_t expected_evals = 0;
    for (size_t i= 0; i < pd.num_elements(); ++i)
      expected_evals += pd.get_samples(i).num_nodes();

    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( expected_evals, handles.size() );

    std::sort( handles.begin(), handles.end() );
    handles.erase( std::unique( handles.begin(), handles.end() ), handles.end() );
    CPPUNIT_ASSERT_EQUAL(expected_evals, handles.size() );

    ++count;
  }

  if (type_is_supported(TRIANGLE)) {
    get_ideal_tris( pd, false );
    size_t expected_evals = 0;
    for (size_t i= 0; i < pd.num_elements(); ++i)
      expected_evals += pd.get_samples(i).num_nodes();

    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( expected_evals, handles.size() );

    std::sort( handles.begin(), handles.end() );
    handles.erase( std::unique( handles.begin(), handles.end() ), handles.end() );
    CPPUNIT_ASSERT_EQUAL(expected_evals, handles.size() );

    ++count;
  }

  if (type_is_supported(QUADRILATERAL)) {
    get_ideal_quads( pd );
    size_t expected_evals = 0;
    for (size_t i= 0; i < pd.num_elements(); ++i)
      expected_evals += pd.get_samples(i).num_nodes();

    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT_EQUAL( expected_evals, handles.size() );

    std::sort( handles.begin(), handles.end() );
    handles.erase( std::unique( handles.begin(), handles.end() ), handles.end() );
    CPPUNIT_ASSERT_EQUAL(expected_evals, handles.size() );

    ++count;
  }
  
  CPPUNIT_ASSERT(count > 0);
}
  void get_ideal_element( EntityTopology type, 
                          bool unit_area, 
                          PatchData& pd,
                          int free_vertex_index );

void QualityMetricTester::test_get_edge_evaluations( EdgeQM* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles;

  CPPUNIT_ASSERT(!types.empty());
  for (size_t i = 0; i < types.size(); ++i) {
    get_ideal_element( types[i], true, pd );
    handles.clear();
    qm->get_evaluations( pd, handles, false, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL( TopologyInfo::edges(types[i]), (unsigned)handles.size() );
  }
}

void QualityMetricTester::test_get_in_element_evaluations( ElemSampleQM* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles1, handles2;
  
  CPPUNIT_ASSERT(!types.empty());
  
  get_ideal_element( types[0], false, pd );
  CPPUNIT_ASSERT(pd.num_elements() == 1);
  qm->get_evaluations( pd, handles1, false, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!handles1.empty());
  qm->get_element_evaluations( pd, 0, handles2, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!handles2.empty());
  std::sort( handles1.begin(), handles1.end() );
  std::sort( handles2.begin(), handles2.end() );
  CPPUNIT_ASSERT( handles1 == handles2 );
  
  get_ideal_element( types.back(), false, pd );
  CPPUNIT_ASSERT(pd.num_elements() == 1);
  qm->get_evaluations( pd, handles1, false, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!handles1.empty());
  qm->get_element_evaluations( pd, 0, handles2, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!handles2.empty());
  std::sort( handles1.begin(), handles1.end() );
  std::sort( handles2.begin(), handles2.end() );
  CPPUNIT_ASSERT( handles1 == handles2 );
}

void QualityMetricTester::test_get_indices_fixed( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd1, pd2;
  if (mSettings) {
    pd1.attach_settings(mSettings);  
    pd2.attach_settings(mSettings);
  }
  std::vector<size_t> handles1, handles2, indices1, indices2;
  double qm_val1, qm_val2;
  
  CPPUNIT_ASSERT(!types.empty());
    // get element with no fixed vertices
  get_ideal_element( types.back(), false, pd1, false );
  qm->get_evaluations( pd1, handles1, false, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!handles1.empty());
    // call evaluate with indices
  size_t handle = handles1.back();
  qm->evaluate_with_indices( pd1, handle, qm_val1, indices1, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!indices1.empty());
    // get element with one fixed vertex
  get_ideal_element( types.back(), false, pd2, true );
  qm->get_evaluations( pd2, handles2, false, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!handles2.empty());
    // If vertex-based metric, need to find corresponding handle
    // in second patch.  The vertex order changed as a result of
    // one of the vertices beging fixed in the second patch.
  if (qm->get_metric_type() == QualityMetric::VERTEX_BASED )
  {
    handle = handles1.back();
    const size_t* conn1 = pd1.element_by_index(0).get_vertex_index_array();
    const size_t* conn2 = pd2.element_by_index(0).get_vertex_index_array();
    const size_t len = TopologyInfo::corners( types.back() );
    const size_t* ptr = std::find( conn1, conn1+len, handle );
    handle = conn2[ptr - conn1];
  }
    // call evaluate with indices
  qm->evaluate_with_indices( pd2, handle, qm_val2, indices2, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(!indices2.empty());
    // make sure we got the same QM value 
    // (shoudn't be affected by fixed/non-fixed)
  CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    // should have gotten back one less index (for element-based metrics)
  if (handles1.size() == 1) {
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size()+1 );
  }
    // indices2 shouldn't contain any fixed vertices
  std::sort( indices2.begin(), indices2.end() );
  CPPUNIT_ASSERT( indices2.back() < pd2.num_free_vertices() );
    // indices2 shouldn/t contain any duplicates
  CPPUNIT_ASSERT( std::unique(indices2.begin(), indices2.end()) == indices2.end() );
}

void QualityMetricTester::test_get_element_indices( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices, verts;
  double qm_val;
  
  CPPUNIT_ASSERT(qm->get_metric_type() == QualityMetric::ELEMENT_BASED);
  CPPUNIT_ASSERT( !types.empty() );
  
  for (size_t i = 0; i < types.size(); ++i) {
      // construct patch w/ one fixed vertex
    get_ideal_element( types[i], false, pd, true );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
      // get vertex list
    verts.clear();
    pd.element_by_index(0).get_vertex_indices( verts );
      // remove fixed vertex
    for (size_t j = 0; j <verts.size(); ++j)
      if (verts[j] >= pd.num_free_vertices()) {
        verts.erase( verts.begin() + j);
        break;
      }
      // evaluate metric
    qm->evaluate_with_indices( pd, handles[0], qm_val, indices, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      // index list should match vertex list (element metric
      // should depend on all free vertices in element).
    CPPUNIT_ASSERT_EQUAL( verts.size(), indices.size() );
    std::sort( verts.begin(), verts.end() );
    std::sort( indices.begin(), indices.end() );
    CPPUNIT_ASSERT( verts == indices );
  }
}  

void QualityMetricTester::test_get_vertex_indices( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> indices, verts;
  double qm_val;
  
  CPPUNIT_ASSERT(qm->get_metric_type() == QualityMetric::VERTEX_BASED);
  CPPUNIT_ASSERT( !types.empty() );
  
  for (size_t i = 0; i < types.size(); ++i) {
      // construct patch w/ one fixed vertex
    get_ideal_element( types[i], false, pd, true );
      // get adjacent vertices
    verts.clear();
    pd.get_adjacent_vertex_indices( 0, verts, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      // remove fixed vertex
    for (size_t j = 0; j <verts.size(); ++j)
      if (verts[j] >= pd.num_free_vertices()) {
        verts.erase( verts.begin() + j);
        break;
      }
    verts.push_back( 0 );
      // evaluate metric
    qm->evaluate_with_indices( pd, 0, qm_val, indices, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      // vertex metric should depend on all adjacent free vertices
    CPPUNIT_ASSERT_EQUAL( verts.size(), indices.size() );
    std::sort( verts.begin(), verts.end() );
    std::sort( indices.begin(), indices.end() );
    CPPUNIT_ASSERT( verts == indices );
  }
}

void QualityMetricTester::test_get_edge_indices( EdgeQM* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd, pd2;
  if (mSettings) {
    pd.attach_settings(mSettings);
    pd2.attach_settings(mSettings);
  }
  std::vector<size_t> handles, indices, indices2, verts;
  std::vector<size_t>::iterator it;
  double qm_val;
  
  CPPUNIT_ASSERT(qm->get_metric_type() == QualityMetric::VERTEX_BASED);
  CPPUNIT_ASSERT( !types.empty() );
  
  for (size_t i = 0; i < types.size(); ++i) {
      // construct patch w/ no free vertices
    get_ideal_element( types[i], false, pd, false );
    
    qm->get_evaluations( pd, handles, false, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL( (size_t)TopologyInfo::edges(types[i]), handles.size() );

    for (size_t j = 0; j < handles.size(); ++j) {
        // evaluate metric
      qm->evaluate_with_indices( pd, handles[j], qm_val, indices, err );
      ASSERT_NO_ERROR(err);
        // evaluation at each edge should depend on at least the two end vertices
      CPPUNIT_ASSERT( indices.size() >= (size_t)2 );
    }
  }
  
  for (size_t i = 0; i < types.size(); ++i) {
      // construct patch w/ one free vertex   
    get_ideal_element( types[i], false, pd, true );
    const size_t fixed_vertex = pd.num_free_vertices();
    
    qm->get_evaluations( pd, handles, false, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL( (size_t)TopologyInfo::edges(types[i]), handles.size() );

    for (size_t j = 0; j < handles.size(); ++j) {
        // evaluate metric
      qm->evaluate_with_indices( pd, handles[j], qm_val, indices, err );
      ASSERT_NO_ERROR(err);
        // evaluation at each edge should depend on at least the two end vertices
      CPPUNIT_ASSERT( !indices.empty() );
      
        // indices should never contain the index of a fixed vertex
      it = std::find(indices.begin(), indices.end(), fixed_vertex );
      CPPUNIT_ASSERT( it == indices.end() );
    }
  }
}  

void QualityMetricTester::test_get_sample_indices( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices, verts;
  double qm_val;
  
  CPPUNIT_ASSERT(qm->get_metric_type() == QualityMetric::ELEMENT_BASED);
  CPPUNIT_ASSERT( !types.empty() );
  
  for (size_t i = 0; i < types.size(); ++i) {
      // get evaluation locations
    get_ideal_element( types[i], false, pd, true );
    const size_t count = pd.get_samples(0).num_nodes();
    CPPUNIT_ASSERT( count > 0 );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      // should have one evaluation for each sample point
    CPPUNIT_ASSERT_EQUAL( count, handles.size() );
      // get vertex list
    verts.clear();
    pd.element_by_index(0).get_vertex_indices( verts );
      // remove fixed vertex
    for (size_t j = 0; j <verts.size(); ++j)
      if (verts[j] >= pd.num_free_vertices()) {
        verts.erase( verts.begin() + j);
        break;
      }
      // evaluate metric
    qm->evaluate_with_indices( pd, handles[0], qm_val, indices, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      // only one fixed vertex, and presumably every possible evaluation
      // depends on more than one vertex, so should be at least one
      // other, non-fixed, vertex.
    CPPUNIT_ASSERT( indices.size() > 0 );
      // make sure no duplicates
    std::sort( indices.begin(), indices.end() );
    CPPUNIT_ASSERT( std::unique(indices.begin(), indices.end()) == indices.end() );
      // check that all indices are in vertex list
    CPPUNIT_ASSERT( indices.size() <= verts.size() );
    for (std::vector<size_t>::iterator j = indices.begin(); j != indices.end(); ++j)
      CPPUNIT_ASSERT( std::find( verts.begin(), verts.end(), *j ) != verts.end() );
  }
}  

void QualityMetricTester::compare_eval_and_eval_with_indices( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices;
  double qm_val1, qm_val2;
  bool rval;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_ideal_element( types[i], false, pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    for (size_t j = 0; j < handles.size(); ++j) {
      rval = qm->evaluate( pd, handles[j], qm_val1, err );
      CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      CPPUNIT_ASSERT( rval );
      rval = qm->evaluate_with_indices( pd, handles[j], qm_val2, indices, err );
      CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
      CPPUNIT_ASSERT( rval );
      
      CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    }
  }
}

void QualityMetricTester::compare_eval_with_indices_and_eval_with_gradient( QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->evaluate_with_indices( pd, handles[j], qm_val1, indices1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_gradient( pd, handles[j], qm_val2, indices2, grad, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );

    std::sort( indices1.begin(), indices1.end() );
    std::sort( indices2.begin(), indices2.end() );
    CPPUNIT_ASSERT( indices1 == indices2 );
    CPPUNIT_ASSERT( !indices1.empty() );
  }
}

void QualityMetricTester::compare_eval_with_indices_and_eval_with_gradient( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_ideal_element( types[i], false, pd, true );
    compare_eval_with_indices_and_eval_with_gradient( qm, pd );
  }
}

void QualityMetricTester::compare_eval_with_indices_and_eval_with_hessian( 
                             QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad;
  std::vector<Matrix3D> Hess;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->evaluate_with_indices( pd, handles[j], qm_val1, indices1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_Hessian( pd, handles[j], qm_val2, indices2, grad, Hess, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );

    std::sort( indices1.begin(), indices1.end() );
    std::sort( indices2.begin(), indices2.end() );
    CPPUNIT_ASSERT( indices1 == indices2 );
    CPPUNIT_ASSERT( !indices1.empty() );
  }
}

void QualityMetricTester::compare_eval_with_indices_and_eval_with_hessian( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_nonideal_element( types[i], pd, true );
    compare_eval_with_indices_and_eval_with_hessian( qm, pd );
  }
}

void QualityMetricTester::compare_eval_with_indices_and_eval_with_diagonal( 
                             QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad;
  std::vector<SymMatrix3D> Hess;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->evaluate_with_indices( pd, handles[j], qm_val1, indices1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_Hessian_diagonal( pd, handles[j], qm_val2, indices2, grad, Hess, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );

    std::sort( indices1.begin(), indices1.end() );
    std::sort( indices2.begin(), indices2.end() );
    CPPUNIT_ASSERT( indices1 == indices2 );
    CPPUNIT_ASSERT( !indices1.empty() );
  }
}

void QualityMetricTester::compare_eval_with_indices_and_eval_with_diagonal( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_nonideal_element( types[i], pd, true );
    compare_eval_with_indices_and_eval_with_diagonal( qm, pd );
  }
}

void QualityMetricTester::compare_eval_with_grad_and_eval_with_hessian( 
                               QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad1, grad2;
  std::vector<Matrix3D> Hess;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->evaluate_with_gradient( pd, handles[j], qm_val1, indices1, grad1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_Hessian( pd, handles[j], qm_val2, indices2, grad2, Hess, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );

    for (size_t k = 0; k < indices1.size(); ++k) {
      std::vector<size_t>::iterator it =
        std::find( indices2.begin(), indices2.end(), indices1[k] );
      CPPUNIT_ASSERT( it != indices2.end() );
      size_t m = it - indices2.begin();
      CPPUNIT_ASSERT_VECTORS_EQUAL( grad1[k], grad2[m], 1e-6 );
    }
  }
}

void QualityMetricTester::compare_eval_with_grad_and_eval_with_hessian( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
 
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_nonideal_element( types[i], pd, true );
    compare_eval_with_grad_and_eval_with_hessian( qm, pd );
  }
}

void QualityMetricTester::compare_eval_with_grad_and_eval_with_diagonal( 
                               QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad1, grad2;
  std::vector<SymMatrix3D> Hess;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->evaluate_with_gradient( pd, handles[j], qm_val1, indices1, grad1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_Hessian_diagonal( pd, handles[j], qm_val2, indices2, grad2, Hess, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );

    for (size_t k = 0; k < indices1.size(); ++k) {
      std::vector<size_t>::iterator it =
        std::find( indices2.begin(), indices2.end(), indices1[k] );
      CPPUNIT_ASSERT( it != indices2.end() );
      size_t m = it - indices2.begin();
      CPPUNIT_ASSERT_VECTORS_EQUAL( grad1[k], grad2[m], 1e-6 );
    }
  }
}

void QualityMetricTester::compare_eval_with_grad_and_eval_with_diagonal( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
 
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_nonideal_element( types[i], pd, true );
    compare_eval_with_grad_and_eval_with_diagonal( qm, pd );
  }
}

void QualityMetricTester::compare_eval_with_diag_and_eval_with_hessian( 
                               QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad1, grad2;
  std::vector<SymMatrix3D> hess1;
  std::vector<Matrix3D> hess2;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->evaluate_with_Hessian_diagonal( pd, handles[j], qm_val1, indices1, grad1, hess1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_Hessian( pd, handles[j], qm_val2, indices2, grad2, hess2, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );
    CPPUNIT_ASSERT_EQUAL( grad1.size(), grad2.size() );
    CPPUNIT_ASSERT_EQUAL( hess2.size(), hess1.size() * (hess1.size()+1) / 2 );
    unsigned h2step = indices2.size();
    std::vector<Matrix3D>::const_iterator h2i = hess2.begin();

    for (size_t k = 0; k < indices2.size(); ++k) {
      std::vector<size_t>::iterator it =
        std::find( indices1.begin(), indices1.end(), indices2[k] );
      CPPUNIT_ASSERT( it != indices1.end() );
      size_t m = it - indices1.begin();
      CPPUNIT_ASSERT_VECTORS_EQUAL( grad1[m], grad2[k], 5e-5 );
      CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(hess1[m]), *h2i, 5e-5 ); 
      h2i += h2step--;
    }
  }
}

void QualityMetricTester::compare_eval_with_diag_and_eval_with_hessian( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
 
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_nonideal_element( types[i], pd, true );
    compare_eval_with_diag_and_eval_with_hessian( qm, pd );
  }
}

void QualityMetricTester::compare_analytical_and_numerical_gradients( 
                                  QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad1, grad2;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->QualityMetric::evaluate_with_gradient( pd, handles[j], qm_val1, indices1, grad1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_gradient( pd, handles[j], qm_val2, indices2, grad2, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );
    CPPUNIT_ASSERT( !indices1.empty() );

    std::vector<size_t>::iterator it1, it2;
    for (it1 = indices1.begin(); it1 != indices1.end(); ++it1) {
      it2 = std::find( indices2.begin(), indices2.end(), *it1 );
      CPPUNIT_ASSERT( it2 != indices2.end() );

      size_t idx1 = it1 - indices1.begin();
      size_t idx2 = it2 - indices2.begin();
      CPPUNIT_ASSERT_VECTORS_EQUAL( grad1[idx1], grad2[idx2], 0.01 );
    }
  }
}

void QualityMetricTester::compare_analytical_and_numerical_gradients( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_nonideal_element( types[i], pd, true );
    compare_analytical_and_numerical_gradients( qm, pd );
  }
}


void QualityMetricTester::compare_analytical_and_numerical_hessians( 
                            QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad1, grad2;
  std::vector<Matrix3D> Hess1, Hess2;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->QualityMetric::evaluate_with_Hessian( pd, handles[j], qm_val1, indices1, grad1, Hess1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_Hessian( pd, handles[j], qm_val2, indices2, grad2, Hess2, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );
    CPPUNIT_ASSERT( !indices1.empty() );

    std::vector<size_t>::iterator it;
    unsigned h = 0;
    for (unsigned r = 0; r < indices1.size(); ++r) {
      it = std::find( indices2.begin(), indices2.end(), indices1[r] );
      CPPUNIT_ASSERT( it != indices2.end() );
      unsigned r2 = it - indices2.begin();

      for (unsigned c = r; c < indices1.size(); ++c, ++h) {
        it = std::find( indices2.begin(), indices2.end(), indices1[c] );
        CPPUNIT_ASSERT( it != indices2.end() );
        unsigned c2 = it - indices2.begin();

        unsigned h2;
        if (r2 <= c2) 
          h2 = indices2.size()*r - r*(r+1)/2 + c;
        else
          h2 = indices2.size()*c - c*(c+1)/2 + r;

        //if (!utest_mat_equal(Hess1[h],Hess2[h2],0.001))
        //  assert(false);
        CPPUNIT_ASSERT_MATRICES_EQUAL( Hess1[h], Hess2[h2], 0.001 );
      }
    }
  }
}

void QualityMetricTester::compare_analytical_and_numerical_hessians( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_nonideal_element( types[i], pd, true );
    compare_analytical_and_numerical_hessians( qm, pd );
  }
}

void QualityMetricTester::compare_analytical_and_numerical_diagonals( 
                            QualityMetric* qm, PatchData& pd )
{
  MsqPrintError err( std::cout );
  std::vector<size_t> handles, indices1, indices2;
  std::vector<Vector3D> grad1, grad2;
  std::vector<Matrix3D> Hess1;
  std::vector<SymMatrix3D> Hess2;
  double qm_val1, qm_val2;
  bool rval;
  
  qm->get_evaluations( pd, handles, false, err );
  CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
  CPPUNIT_ASSERT( !handles.empty() );
  for (size_t j = 0; j < handles.size(); ++j) {
    rval = qm->QualityMetric::evaluate_with_Hessian( pd, handles[j], qm_val1, indices1, grad1, Hess1, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );
    rval = qm->evaluate_with_Hessian_diagonal( pd, handles[j], qm_val2, indices2, grad2, Hess2, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( rval );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( qm_val1, qm_val2, 1e-6 );
    CPPUNIT_ASSERT_EQUAL( indices1.size(), indices2.size() );
    CPPUNIT_ASSERT( !indices1.empty() );
    CPPUNIT_ASSERT_EQUAL( indices1.size() * (indices1.size()+1) / 2, Hess1.size() );
    CPPUNIT_ASSERT_EQUAL( indices2.size(), Hess2.size() );

    size_t h = 0;
    std::vector<size_t>::iterator it;
    for (unsigned r = 0; r < indices1.size(); ++r) {
      it = std::find( indices2.begin(), indices2.end(), indices1[r] );
      CPPUNIT_ASSERT( it != indices2.end() );
      unsigned r2 = it - indices2.begin();
      //if (!utest_mat_equal(Hess1[h],Hess2[r2],0.001))
      //  assert(false);
      CPPUNIT_ASSERT_MATRICES_EQUAL( Hess1[h], Hess2[r2], 0.001 );
      h += indices1.size() - r;
    }
  }
}

void QualityMetricTester::compare_analytical_and_numerical_diagonals( QualityMetric* qm )
{
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  MsqError err;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_nonideal_element( types[i], pd, true );
    compare_analytical_and_numerical_diagonals( qm, pd );
  }
}

class PatchTranslate : public QualityMetricTester::PatchXform
{
    Vector3D offset;
  public:
    inline PatchTranslate( Vector3D s ) : offset(s) {}
    virtual ~PatchTranslate();
    virtual void xform( PatchData& pd, PlanarDomain* dom ) ;
    virtual void xform_grad( std::vector<Vector3D>& );
};


PatchTranslate::~PatchTranslate() {}
 void PatchTranslate::xform( PatchData& pd, PlanarDomain* dom )
{
    // If patch has associated planar geometry, translate that too
  MeshDomain* geom = pd.get_domain();
  if (geom) {
    PlanarDomain* plane = dynamic_cast<PlanarDomain*>(geom);
    CPPUNIT_ASSERT(plane);
    CPPUNIT_ASSERT(dom);

    Vector3D norm = plane->get_normal();
    Vector3D point = plane->get_origin() + offset;
    dom->set_plane( norm, point );
    pd.set_domain( dom );
  }

    // Translate mesh vertices
  MsqError err;
  for (size_t i = 0; i < pd.num_nodes(); ++i)
    pd.move_vertex( offset, i, err );
}
void PatchTranslate::xform_grad( std::vector<Vector3D>& ) {}


class PatchScale : public QualityMetricTester::PatchXform
{
    double scaleFactor;
  public:
    inline PatchScale( double s ) : scaleFactor(s) {}
    virtual ~PatchScale();
    virtual void xform( PatchData& pd, PlanarDomain* dom ) ;
    virtual void xform_grad( std::vector<Vector3D>& );
};

PatchScale::~PatchScale() {}
void PatchScale::xform( PatchData& pd, PlanarDomain* dom ) 
{
    // If patch has associated planar geometry, scale distance from origin
  MeshDomain* geom = pd.get_domain();
  if (geom) {
    PlanarDomain* plane = dynamic_cast<PlanarDomain*>(geom);
    CPPUNIT_ASSERT(plane);
    CPPUNIT_ASSERT(dom);

    Vector3D norm = plane->get_normal();
    Vector3D point = scaleFactor * plane->get_origin();
    dom->set_plane( norm, point );
    pd.set_domain( dom );
  }

    // Scale about origin
  MsqError err;
  for (size_t i = 0; i < pd.num_nodes(); ++i)
    pd.set_vertex_coordinates( pd.vertex_by_index(i) * scaleFactor, i, err );
}
void PatchScale::xform_grad( std::vector<Vector3D>& ) {}

class PatchRotate : public QualityMetricTester::PatchXform
{
    Matrix3D rotation;
  public:
    inline PatchRotate( Matrix3D s ) : rotation(s) {}
    virtual ~PatchRotate();
    virtual void xform( PatchData& pd, PlanarDomain* dom ) ;
    virtual void xform_grad( std::vector<Vector3D>& grad );
};

PatchRotate::~PatchRotate() {}
void PatchRotate::xform( PatchData& pd, PlanarDomain* dom ) 
{
    // If patch has associated planar geometry, rotate that also
  MeshDomain* geom = pd.get_domain();
  if (geom) {
    PlanarDomain* plane = dynamic_cast<PlanarDomain*>(geom);
    CPPUNIT_ASSERT(plane);
    CPPUNIT_ASSERT(dom);

    Vector3D norm = rotation * plane->get_normal();
    Vector3D point = rotation * plane->get_origin();
    dom->set_plane( norm, point );
    pd.set_domain( dom );
  }

    // Scale about origin
  MsqError err;
  for (size_t i = 0; i < pd.num_nodes(); ++i)
    pd.set_vertex_coordinates( rotation * pd.vertex_by_index(i), i, err );
}
void PatchRotate::xform_grad( std::vector<Vector3D>& grad )
{
  for (unsigned i = 0; i < grad.size(); ++i)
    grad[i] = rotation * grad[i];
}

void QualityMetricTester::test_transform_invariant( QualityMetric* qm,
                                QualityMetricTester::PatchXform& xform,
                                bool untangle )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles;
  double qmval;
  bool rval;
  PlanarDomain dom(Vector3D(0,0,1),Vector3D(0,0,0));
  size_t i, j;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (i = 0; i < types.size(); ++i) {
    if (untangle)
      get_inverted_element( types[i], pd );
    else
      get_nonideal_element( types[i], pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    
      // get value(s) for ideal element
    std::vector<double> values(handles.size());
    for (j = 0; j < handles.size(); ++j) {
      rval = qm->evaluate( pd, handles[j], values[j], err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
    }
    
      // compare values for transformed patch
    xform.xform( pd, &dom );
    for (j = 0; j < handles.size(); ++j) {
      rval = qm->evaluate( pd, handles[j], qmval, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( values[j], qmval, 1e-3 );
    }
  }
}

void QualityMetricTester::test_grad_transform_invariant( QualityMetric* qm, 
                                    QualityMetricTester::PatchXform& xform,
                                    bool untangle )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices2;
  std::vector<Vector3D> grad2;
  double qmval;
  bool rval;
  PlanarDomain dom(Vector3D(0,0,1),Vector3D(0,0,0));
  size_t i, j, k;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (i = 0; i < types.size(); ++i) {
    if (untangle)
      get_inverted_element( types[i], pd );
    else
      get_nonideal_element( types[i], pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    
      // get value(s) for orignal patch
    std::vector< std::vector<Vector3D> > grad1(handles.size());
    std::vector< std::vector<size_t> > indices1(handles.size());
    for (j = 0; j < handles.size(); ++j) {
      rval = qm->evaluate_with_gradient( pd, handles[j], qmval, indices1[j], grad1[j], err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
      CPPUNIT_ASSERT( !grad1[j].empty() );
    }
    
      // transform patch
    xform.xform( pd, &dom );
    
      // compare values
    for (j = 0; j < handles.size(); ++j) {
      rval = qm->evaluate_with_gradient( pd, handles[j], qmval, indices2, grad2, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
    
      CPPUNIT_ASSERT_EQUAL( indices1[j].size(), indices2.size() );
      CPPUNIT_ASSERT_EQUAL( grad1[j].size(), grad2.size() );
      CPPUNIT_ASSERT( indices1[j] == indices2 );
    
      xform.xform_grad( grad1[j] );
      for (k = 0; k < grad2.size(); ++k)
        CPPUNIT_ASSERT_VECTORS_EQUAL( grad1[j][k], grad2[k], 1e-3 );
    }
  }
}

void QualityMetricTester::test_hessian_transform_invariant( QualityMetric* qm, 
                                       QualityMetricTester::PatchXform& xform,
                                       bool untangle )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices2;
  std::vector<Vector3D> grad1, grad2;
  std::vector<Matrix3D> hess2;
  double qmval;
  bool rval;
  PlanarDomain dom(Vector3D(0,0,1),Vector3D(0,0,0));
  size_t i, j, k;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (i = 0; i < types.size(); ++i) {
    if (untangle)
      get_inverted_element( types[i], pd );
    else
      get_nonideal_element( types[i], pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    
    std::vector< std::vector<size_t> > indices1(handles.size());
    std::vector< std::vector<Matrix3D> > hess1(handles.size());
    for (j = 0; j < handles.size(); ++j) {
      rval = qm->evaluate_with_Hessian( pd, handles[j], qmval, indices1[j], grad1, hess1[j], err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
      CPPUNIT_ASSERT(!hess1[j].empty());
    }
    
    xform.xform( pd, &dom );
    
    for (j = 0; j < handles.size(); ++j) {
      rval = qm->evaluate_with_Hessian( pd, handles[j], qmval, indices2, grad2, hess2, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
    
      CPPUNIT_ASSERT_EQUAL( indices1[j].size(), indices2.size() );
      CPPUNIT_ASSERT_EQUAL( hess1[j].size(), hess2.size() );
      CPPUNIT_ASSERT( indices1[j] == indices2 );
    
      for (k = 0; k < hess2.size(); ++k)
        CPPUNIT_ASSERT_MATRICES_EQUAL( hess1[j][k], hess2[k], 1e-3 );
    }
  }
}

void QualityMetricTester::test_location_invariant( QualityMetric* qm, bool untangle )
{
  PatchTranslate xform( Vector3D(1,1,-1) );
  test_transform_invariant( qm, xform, untangle );
}

void QualityMetricTester::test_scale_invariant( QualityMetric* qm, bool untangle )
{
  PatchScale scale1( 0.5 ), scale2( 2.0 );
  test_transform_invariant( qm, scale1, untangle );
  test_transform_invariant( qm, scale2, untangle );
}

void QualityMetricTester::test_orient_invariant( QualityMetric* qm, bool untangle )
{
    // rotate 45 degrees about x-axis
  const double f = std::sqrt(2.0)/2;
  Matrix3D rotation( 1, 0, 0,
                     0, f, f, 
                     0,-f, f );
  PatchRotate xform(rotation);
  test_transform_invariant( qm, xform, untangle );
}

void QualityMetricTester::test_grad_location_invariant( QualityMetric* qm, bool untangle )
{
  PatchTranslate xform( Vector3D(1,1,-1) );
  test_grad_transform_invariant( qm, xform, untangle );
}

void QualityMetricTester::test_hessian_location_invariant( QualityMetric* qm, bool untangle )
{
  PatchTranslate xform( Vector3D(1,1,-1) );
  test_hessian_transform_invariant( qm, xform, untangle );
}

void QualityMetricTester::test_grad_orient_invariant( QualityMetric* qm, bool untangle )
{
    // rotate 45 degrees about x-axis
  const double f = std::sqrt(2.0)/2;
  Matrix3D rotation( 1, 0, 0,
                     0, f, f, 
                     0,-f, f );
  PatchRotate xform(rotation);
  test_grad_transform_invariant( qm, xform, untangle );
}

void QualityMetricTester::test_measures_transform( QualityMetric* qm,
                               QualityMetricTester::PatchXform& xform,
                               bool unit_area )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles;
  double qmval1, qmval2;
  bool rval;
  PlanarDomain dom(Vector3D(0,0,1),Vector3D(0,0,0));
  bool decreasing = false;
  if (qm->get_negate_flag() != 1) {
    CPPUNIT_ASSERT_EQUAL( qm->get_negate_flag(), -1 );
    decreasing = true;
  }
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    get_ideal_element( types[i], unit_area, pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    
    rval = qm->evaluate( pd, handles[0], qmval1, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT(rval);
    
    xform.xform( pd, &dom );
    rval = qm->evaluate( pd, handles[0], qmval2, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT(rval);
    if (decreasing)
      CPPUNIT_ASSERT( qmval2 < qmval1 );
    else
      CPPUNIT_ASSERT( qmval2 > qmval1 );
  }
}

void QualityMetricTester::test_measures_size( QualityMetric* qm, bool unit_area )
{
  PatchScale s1(0.5), s2(2.0);
  test_measures_transform( qm, s1, unit_area );
  test_measures_transform( qm, s2, unit_area );
}

void QualityMetricTester::test_measures_in_plane_orientation( QualityMetric* qm )
{
    // rotate 45 degrees about z-axis
  const double f = std::sqrt(2.0)/2;
  Matrix3D rotation( f,-f, 0,
                     f, f, 0, 
                     0, 0, 1 );
  PatchRotate xform( rotation );
  test_measures_transform( qm, xform, false );
}

void QualityMetricTester::test_measures_out_of_plane_orientation( QualityMetric* qm )
{
    // rotate 45 degrees about z-axis
  const double f = std::sqrt(2.0)/2;
  Matrix3D rotation( 1, 0, 0,
                     0, f, f, 
                     0,-f, f );
  PatchRotate xform( rotation );
  test_measures_transform( qm, xform, false );
}

typedef void (QualityMetricTester::*patch_func_t)( EntityTopology, PatchData& );
static void test_bad_element( QualityMetricTester* instance,
                              QualityMetric* qm, 
                              patch_func_t type_func, 
                              bool should_succeed,
                              const std::vector<EntityTopology>& types,
                              const Settings* mSettings )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles;
  double qmval;
  bool rval;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
    (instance->*type_func)( types[i], pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    
    rval = true;
    for (size_t j = 0; j < handles.size(); ++j) {
      rval = rval && qm->evaluate( pd, handles[j], qmval, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    }
    CPPUNIT_ASSERT_EQUAL( should_succeed, rval );
  }
}

void QualityMetricTester::test_evaluate_inverted_element( QualityMetric* qm, bool should_succeed )
  { test_bad_element( this, qm, &QualityMetricTester::get_inverted_element, should_succeed, types, mSettings ); }
void QualityMetricTester::test_evaluate_degenerate_element( QualityMetric* qm, bool should_succeed )
  { test_bad_element( this, qm, &QualityMetricTester::get_degenerate_element, should_succeed, types, mSettings ); }
void QualityMetricTester::test_evaluate_zero_element( QualityMetric* qm, bool should_succeed )
  { test_bad_element( this, qm, &QualityMetricTester::get_zero_element, should_succeed, types, mSettings ); }

void QualityMetricTester::test_measures_quality( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles;
  double qmval1, qmval2;
  bool rval;
  bool decreasing = false;
  if (qm->get_negate_flag() != 1) {
    CPPUNIT_ASSERT_EQUAL( qm->get_negate_flag(), -1 );
    decreasing = true;
  }
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
      // get a patch
    get_ideal_element( types[i], false, pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
      // evaluate for ideal element
    rval = qm->evaluate( pd, handles[0], qmval1, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT( rval );
      // calculate element centroid
    const size_t* verts = pd.element_by_index(0).get_vertex_index_array();
    const MsqVertex* coords = pd.get_vertex_array(err);
    size_t n = pd.element_by_index(0).vertex_count();
    Vector3D cent(0,0,0);
    Vector3D coords0 = coords[verts[0]];
    for (unsigned j = 0; j < n; ++j)
      cent += coords[verts[j]];
    cent /= n;
      // evaluate for non-ideal element (decreased area)
    pd.set_vertex_coordinates( 0.5*(cent + coords0), verts[0], err );
    rval = qm->evaluate( pd, handles[0], qmval2, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT( rval );
      // check values
    if (decreasing) {
      CPPUNIT_ASSERT( qmval2 < qmval1 );
    }
    else {
      CPPUNIT_ASSERT( qmval2 > qmval1 );
    }
      // evaluate for non-ideal element (increased area)
    pd.set_vertex_coordinates( 3*coords0 - 2*cent, verts[0], err );
    rval = qm->evaluate( pd, handles[0], qmval2, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT( rval );
      // check values
    if (decreasing) {
      CPPUNIT_ASSERT( qmval2 < qmval1 );
    }
    else {
      CPPUNIT_ASSERT( qmval2 > qmval1 );
    }
  }
}

void QualityMetricTester::test_domain_deviation_quality( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices;
  size_t i;
  double qmval1, qmval2;
  bool rval;
  bool decreasing = false;
  if (qm->get_negate_flag() != 1) {
    CPPUNIT_ASSERT_EQUAL( qm->get_negate_flag(), -1 );
    decreasing = true;
  }
  
  CPPUNIT_ASSERT( type_is_supported(TRIANGLE)||type_is_supported(QUADRILATERAL) );
  
  if (type_is_supported(TRIANGLE)) {
      // get a patch
    get_ideal_tris( pd, false );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
      // find an evaluation that depends on the free vertex
    for (i = 0; i < handles.size(); ++i) {
      rval = qm->evaluate_with_indices( pd, handles[i], qmval1, indices, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
      if (!indices.empty())
        break;
    }
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() ); // one free vertex in patch
    CPPUNIT_ASSERT_EQUAL( (size_t)0, indices.front() ); // one free vertex in patch
    const size_t handle = handles[i];
      // rotate patch about origin such that elements are no longer in the plane
    const Matrix3D rot( M_PI/4.0, Vector3D( 1, 1, 0 ) );
    for (i = 0; i < pd.num_nodes(); ++i)
      pd.set_vertex_coordinates( rot * pd.vertex_by_index(i), i, err );
      // evaluate for rotated patch
    rval = qm->evaluate( pd, handle, qmval2, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT(rval);
      // check quality
    if (decreasing) {
      CPPUNIT_ASSERT( qmval1 > qmval2 );
    }
    else {
      CPPUNIT_ASSERT( qmval1 < qmval2 );
    }
  }
  
  if (type_is_supported(QUADRILATERAL)) {
      // get a patch
    get_ideal_quads( pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
      // find an evaluation that depends on the free vertex
    for (i = 0; i < handles.size(); ++i) {
      rval = qm->evaluate_with_indices( pd, handles[i], qmval1, indices, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
      if (!indices.empty())
        break;
    }
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() ); // one free vertex in patch
    CPPUNIT_ASSERT_EQUAL( (size_t)0, indices.front() ); // one free vertex in patch
    const size_t handle = handles[i];
      // rotate patch about origin such that elements are no longer in the plane
    const Matrix3D rot( M_PI/3.0, Vector3D( 1, 0, 0 ) );
    for (i = 0; i < pd.num_nodes(); ++i)
      pd.set_vertex_coordinates( rot * pd.vertex_by_index(i), i, err );
      // evaluate for rotated patch
    rval = qm->evaluate( pd, handle, qmval2, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT(rval);
      // check quality
    if (decreasing) {
      CPPUNIT_ASSERT( qmval1 > qmval2 );
    }
    else {
      CPPUNIT_ASSERT( qmval1 < qmval2 );
    }
  }
}
 

void QualityMetricTester::test_domain_deviation_gradient( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices;
  std::vector<Vector3D> grad;
  size_t i;
  double qmval;
  bool rval;
  
  CPPUNIT_ASSERT( type_is_supported(TRIANGLE)||type_is_supported(QUADRILATERAL) );
  
  if (type_is_supported(TRIANGLE)) {
      // get a patch
    get_ideal_tris( pd, false );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
      // move the free vertex out of the plane
    pd.move_vertex( Vector3D( 0, 0, 0.2 ), 0, err );
      // find an evaluation that depends on the free vertex
    for (i = 0; i < handles.size(); ++i) {
      rval = qm->evaluate_with_gradient( pd, handles[i], qmval, indices, grad, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
      if (!grad.empty())
        break;
    }
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() ); // one free vertex in patch
    CPPUNIT_ASSERT_EQUAL( (size_t)0, indices.front() ); // one free vertex in patch
    CPPUNIT_ASSERT_EQUAL( (size_t)1, grad.size() );
    
    Vector3D v = qm->get_negate_flag() * Vector3D( 0, 0, 1 );
    CPPUNIT_ASSERT( v % grad[0] > DBL_EPSILON );
  }
  
  if (type_is_supported(QUADRILATERAL)) {
      // get a patch
    get_ideal_quads( pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
      // move the free vertex out of the plane
    pd.move_vertex( Vector3D( 0, 0, 0.2 ), 0, err );
      // find an evaluation that depends on the free vertex
    for (i = 0; i < handles.size(); ++i) {
      rval = qm->evaluate_with_gradient( pd, handles[i], qmval, indices, grad, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT(rval);
      if (!grad.empty())
        break;
    }
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() ); // one free vertex in patch
    CPPUNIT_ASSERT_EQUAL( (size_t)0, indices.front() ); // one free vertex in patch
    CPPUNIT_ASSERT_EQUAL( (size_t)1, grad.size() );
    
    Vector3D v = qm->get_negate_flag() * Vector3D( 0, 0, 1 );
    CPPUNIT_ASSERT( v % grad[0] > DBL_EPSILON );
  }

}

void QualityMetricTester::test_gradient_reflects_quality( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices;
  std::vector<Vector3D> grad;
  double qmval;
  bool rval;
  bool decreasing = false;
  if (qm->get_negate_flag() != 1) {
    CPPUNIT_ASSERT_EQUAL( qm->get_negate_flag(), -1 );
    decreasing = true;
  }
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
      // get a patch
    get_ideal_element( types[i], false, pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
      // make element non-ideal
    const size_t* verts = pd.element_by_index(0).get_vertex_index_array();
    const MsqVertex* coords = pd.get_vertex_array(err);
    size_t n = pd.element_by_index(0).vertex_count();
    Vector3D sum(0,0,0);
    for (unsigned j = 0; j < n; ++j)
      sum += coords[verts[j]];
    const Vector3D init_pos = coords[verts[0]];
    const Vector3D new_pos = 0.5 * coords[verts[0]] + sum/(2*n);
    pd.set_vertex_coordinates( new_pos, verts[0], err );
      // evaluate for non-ideal element
    rval = qm->evaluate_with_gradient( pd, handles[0], qmval, indices, grad, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT( rval );
      // find gradient for vertex 0
    std::vector<size_t>::iterator it = std::find( indices.begin(), indices.end(), verts[0] );
    CPPUNIT_ASSERT( it != indices.end() );
    Vector3D g = grad[it - indices.begin()];
      // check gradient direction
    Vector3D v = init_pos - coords[verts[0]];
    double dotprod = v % g;
    if (decreasing)
      CPPUNIT_ASSERT( dotprod > 1e-6 );
    else
      CPPUNIT_ASSERT( dotprod < -1e-6 );
  }
}

void QualityMetricTester::test_measures_vertex_quality( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  double qmval1, qmval2;
  bool rval;
  bool increasing = true;
  if (qm->get_negate_flag() != 1) {
    CPPUNIT_ASSERT_EQUAL( qm->get_negate_flag(), -1 );
    increasing = true;
  }
  
  unsigned count = 0;
  for (size_t i = 0; i < 3; ++i) {
      // get a patch
    switch (i) { 
      case 0: 
        if (!type_is_supported(TRIANGLE)) continue;
        get_ideal_tris( pd, false ); break;
      case 1: 
        if (!type_is_supported(QUADRILATERAL)) continue;
        get_ideal_quads( pd ); break;
      case 2: 
        if (!type_is_supported(HEXAHEDRON)) continue;
        get_ideal_hexes( pd ); break;
    }
    ++count;
    
      // evaluate for ideal element
    rval = qm->evaluate( pd, 0, qmval1, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT( rval );
      // move center vertex
    pd.move_vertex( Vector3D(0.3,0.3,0.3), 0, err );
      // evaluate for non-ideal element
    rval = qm->evaluate( pd, 0, qmval2, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT( rval );
      // check values
    if (increasing) {
      CPPUNIT_ASSERT( qmval2 > qmval1 );
    }
    else {
      CPPUNIT_ASSERT( qmval2 < qmval1 );
    }
  }
    // make sure we tested something
  CPPUNIT_ASSERT( count > 0 );
}

void QualityMetricTester::test_vertex_gradient_reflects_quality( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> indices;
  std::vector<Vector3D> grad;
  double qmval;
  bool rval;
  bool increasing = true;
  if (qm->get_negate_flag() != 1) {
    CPPUNIT_ASSERT_EQUAL( qm->get_negate_flag(), -1 );
    increasing = false;
  }
  
  unsigned count = 0;
  for (size_t i = 0; i < 3; ++i) {
      // get a patch
    switch (i) { 
      case 0: 
        if (!type_is_supported(TRIANGLE)) continue;
        get_ideal_tris( pd, false ); break;
      case 1: 
        if (!type_is_supported(QUADRILATERAL)) continue;
        get_ideal_quads( pd ); break;
      case 2: 
        if (!type_is_supported(HEXAHEDRON)) continue;
        get_ideal_hexes( pd ); break;
    }
    
      // move center vertex
    pd.move_vertex( Vector3D(0.3,0.3,0.3), 0, err );

      // evaluate for ideal non-element
    rval = qm->evaluate_with_gradient( pd, 0, qmval, indices, grad, err );
    CPPUNIT_ASSERT(!MSQ_CHKERR(err));
    CPPUNIT_ASSERT( rval );
      // find gradient for vertex 0
    if (indices.empty())
      continue;
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() ); // only 1 free vertex
    CPPUNIT_ASSERT_EQUAL( (size_t)0, indices[0] ); // and that vertex is vertex 0
    ++count;
    Vector3D g = grad[0];
      // check gradient direction
    if (increasing)
      CPPUNIT_ASSERT( g[2] > 1e-6 );
    else
      CPPUNIT_ASSERT( g[2] < -1e-6 );
  }
    // make sure we tested something
  CPPUNIT_ASSERT( count > 0 );
}

void QualityMetricTester::test_ideal_element_zero_gradient( QualityMetric* qm, bool unit_area )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices;
  std::vector<Vector3D> grad;
  double qmval;
  bool rval;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
      // get a patch
    get_ideal_element( types[i], unit_area, pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    for (size_t j = 0; j < handles.size(); ++j) {
        // evaluate 
      rval = qm->evaluate_with_gradient( pd, handles[j], qmval, indices, grad, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT( rval );
      CPPUNIT_ASSERT_EQUAL( indices.size(), grad.size() );
        // check that all gradients are zero
      for (size_t k = 0; k < grad.size(); ++k) {
        CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,0), grad[k], 1e-4 );
      } // for(k < grad.size())
    } // for(j < handles.size())
  } // for(i < types.size())
}

void QualityMetricTester::test_ideal_element_zero_vertex_gradient( 
                                  QualityMetric* qm, bool unit_area )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices;
  std::vector<Vector3D> grad;
  double qmval;
  bool rval;
  
  unsigned count = 0;
  for (size_t i = 0; i < 3; ++i) {
      // get a patch
    switch (i) { 
      case 0: 
        if (!type_is_supported(TRIANGLE)) continue;
        get_ideal_tris( pd, unit_area ); break;
      case 1: 
        if (!type_is_supported(QUADRILATERAL)) continue;
        get_ideal_quads( pd ); break;
      case 2: 
        if (!type_is_supported(HEXAHEDRON)) continue;
        get_ideal_hexes( pd ); break;
    }
     
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    for (size_t j = 0; j < handles.size(); ++j) {
        // evaluate 
      rval = qm->evaluate_with_gradient( pd, handles[j], qmval, indices, grad, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT( rval );
      CPPUNIT_ASSERT_EQUAL( indices.size(), grad.size() );
      if (grad.empty())
        continue;
      ++count;
      CPPUNIT_ASSERT_EQUAL( (size_t)1, grad.size() ); // only 1 free vertex in patch
        // check that all gradients are zero
      CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,0), grad[0], 1e-4 );
    } // for(j < handles.size())
  } // for(i < 3)
    // make sure we tested something
  CPPUNIT_ASSERT( count > 0 );
}

static double** allocate_matrix( unsigned n )
{
  unsigned num_ptr = n;
  if (num_ptr%2) ++num_ptr;
  
  void* storage = malloc( n*n*sizeof(double) + num_ptr*sizeof(double*) );
  double** ptrs = (double**)storage;
  double* data = (double*)(ptrs+num_ptr);
  for (unsigned i = 0; i < n; ++i)
    ptrs[i] = data + i*n;
  return ptrs;
}

static double value( const Matrix3D* blocks, unsigned n, unsigned i, unsigned j )
{
  if (i > j)
    std::swap(i,j);
  
  unsigned mi = i / 3;
  unsigned mj = j / 3;
  unsigned idx = n*mi - mi*(mi+1)/2 + mj;
  return blocks[idx][i%3][j%3];
}
/*
static bool positive_definite( const Matrix3D* blocks, unsigned n )
{
    // Do Cholesky-Banachiewicz decompositon of the matrix.
    // If this results in any imaginary values for diagonal terms,
    // then the matrix is not positive-definite.
  const int N = 3*n;
  double** mat = allocate_matrix(N);
  int i, j, k;
  double s;
  for (i = 0; i < N; ++i)
  {
    for (j = 0; j < i; ++j)
    {
      s = value(blocks,n,i,j);
      for (k = 0; k < j - 1; ++k)
        s -= mat[i][k]*mat[j][k];
      mat[i][j] = s / mat[j][j];
    }
    s = value(blocks,n,i,i);
    for (k = 0; k < i - 1; ++k)
      s -= mat[i][k]*mat[i][k];
    if (s < 0.0)
    {
      free( mat );
      return false;
    }
    mat[i][i] = sqrt(s);
  }
  return true;
}
*/

    /** Test that Hessian is positive-definite for ideal elements */
void QualityMetricTester::test_ideal_element_positive_definite_Hessian( 
                                       QualityMetric* qm, bool unit_area )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices;
  std::vector<Vector3D> grad;
  std::vector<Matrix3D> Hess;
  double qmval;
  bool rval;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
      // get a patch
    get_ideal_element( types[i], unit_area, pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    for (size_t j = 0; j < handles.size(); ++j) {
        // evaluate 
      rval = qm->evaluate_with_Hessian( pd, handles[j], qmval, indices, grad, Hess, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT( rval );
      const size_t n = indices.size();
      CPPUNIT_ASSERT_EQUAL( n, grad.size() );
      CPPUNIT_ASSERT_EQUAL( n*(n+1)/2, Hess.size() );
        // if necessary, adjust Hessian for negate_flag();
      if (qm->get_negate_flag() == -1)
        for (size_t k = 0; k < Hess.size(); ++k)
          Hess[k] *= -1.0;
        // check that all diagonal blocks are positive definite
      size_t h_idx = 0;
      for (size_t k = 0; k < n; h_idx+=(n-k), ++k) {
        CPPUNIT_ASSERT( Hess[h_idx].positive_definite() );
      }
        // check that the entire Hessian is positive definite
      //CPPUNIT_ASSERT( positive_definite( arrptr(Hess), n ) );
    } // for(j < handles.size())
  } // for(i < types.size())
}

    /** Test that diagonal bocks of Hessian are symmetric */
void QualityMetricTester::test_symmetric_Hessian_diagonal_blocks( QualityMetric* qm )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (mSettings)
    pd.attach_settings(mSettings);  
  std::vector<size_t> handles, indices;
  std::vector<Vector3D> grad;
  std::vector<Matrix3D> Hess;
  double qmval;
  bool rval;
  
  CPPUNIT_ASSERT( !types.empty() );
  for (size_t i = 0; i < types.size(); ++i) {
      // get a patch
    get_ideal_element( types[i], false, pd );
    qm->get_evaluations( pd, handles, false, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) );
    CPPUNIT_ASSERT( !handles.empty() );
    for (size_t j = 0; j < handles.size(); ++j) {
        // evaluate 
      rval = qm->evaluate_with_Hessian( pd, handles[j], qmval, indices, grad, Hess, err );
      CPPUNIT_ASSERT(!MSQ_CHKERR(err));
      CPPUNIT_ASSERT( rval );
      const size_t n = indices.size();
      CPPUNIT_ASSERT_EQUAL( n, grad.size() );
      CPPUNIT_ASSERT_EQUAL( n*(n+1)/2, Hess.size() );
        // check that diagonal Hessian blocks are symmetric
      size_t h_idx = 0;
      for (size_t k = 0; k < n; h_idx+=(n-k), ++k) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL( Hess[h_idx][0][1], Hess[h_idx][1][0], 5e-3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL( Hess[h_idx][0][2], Hess[h_idx][2][0], 5e-3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL( Hess[h_idx][1][2], Hess[h_idx][2][1], 5e-3 );
      } // for(k < Hess.size())
    } // for(j < handles.size())
  } // for(i < types.size())
}

void QualityMetricTester::test_gradient_with_fixed_vertex( QualityMetric* qm,
                                                           const Settings* settings )
{
  CPPUNIT_ASSERT(!types.empty());
  for (unsigned i = 0; i < types.size(); ++i)
    test_gradient_with_fixed_vertex( types[i], qm, settings );
}

void QualityMetricTester::test_hessian_with_fixed_vertex( QualityMetric* qm,
                                                          const Settings* settings )
{
  CPPUNIT_ASSERT(!types.empty());
  for (unsigned i = 0; i < types.size(); ++i)
    test_hessian_with_fixed_vertex( types[i], qm, settings );
}

void QualityMetricTester::test_diagonal_with_fixed_vertex( QualityMetric* qm, 
                                                           const Settings* settings )
{
  CPPUNIT_ASSERT(!types.empty());
  for (unsigned i = 0; i < types.size(); ++i)
    test_diagonal_with_fixed_vertex( types[i], qm, settings );
}

void QualityMetricTester::test_gradient_with_fixed_vertex( EntityTopology type, 
                                                           QualityMetric* qm,
                                                           const Settings* settings )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (settings)
    pd.attach_settings(settings);
  else if (mSettings)
    pd.attach_settings(mSettings);
  
  std::vector<size_t> handles, indices, indices_fixed, conn;
  std::vector<Vector3D> grad, grad_fixed;
  double val, val_fixed;
  bool rval;
  const int n = TopologyInfo::corners( type );
  
  get_nonideal_element( type, pd );
  pd.element_by_index(0).get_vertex_indices( conn );
  
  qm->get_evaluations( pd, handles, false, err );
  ASSERT_NO_ERROR(err);
    // For sample-based metrics, it is difficult to set up such that
    // both surface and volume elements have only one sample point.
    // Make things easier by skipping element types with no sample points.
  if (handles.empty()) 
    return;
  CPPUNIT_ASSERT_EQUAL_MESSAGE( "test only valid for element-based metrics", 
                                (size_t)1, handles.size() );
  size_t h = handles[0];
  
  rval = qm->evaluate_with_gradient( pd, h, val, indices, grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( (size_t)n, indices.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)n, grad.size() );
  
  for (int i = 0; i < n; ++i) {
    get_nonideal_element( type, pd, i );

    qm->get_evaluations( pd, handles, false, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL_MESSAGE( "test only valid for element-based metrics", 
                                  (size_t)1, handles.size() );
    h = handles[0];
    
    rval = qm->evaluate_with_gradient( pd, h, val_fixed, indices_fixed, grad_fixed, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(rval);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices_fixed.size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, grad_fixed.size() );
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( val, val_fixed, 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( grad[conn[i]], grad_fixed.front(), 1e-5 );
  }
}


void QualityMetricTester::test_hessian_with_fixed_vertex( EntityTopology type, 
                                                          QualityMetric* qm,
                                                          const Settings* settings )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (settings) 
    pd.attach_settings(settings);
  else if (mSettings)
    pd.attach_settings(mSettings);
  std::vector<size_t> handles, indices, indices_fixed, conn;
  std::vector<Vector3D> grad, grad_fixed;
  std::vector<Matrix3D> hess, hess_fixed;
  double val, val_fixed;
  bool rval;
  const int n = TopologyInfo::corners( type );
  
  get_nonideal_element( type, pd );
  pd.element_by_index(0).get_vertex_indices( conn );
  
  qm->get_evaluations( pd, handles, false, err );
  ASSERT_NO_ERROR(err);
    // For sample-based metrics, it is difficult to set up such that
    // both surface and volume elements have only one sample point.
    // Make things easier by skipping element types with no sample points.
  if (handles.empty()) 
    return;
  CPPUNIT_ASSERT_EQUAL_MESSAGE( "test only valid for element-based metrics", 
                                (size_t)1, handles.size() );
  size_t h = handles[0];
  
  rval = qm->evaluate_with_Hessian( pd, h, val, indices, grad, hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( (size_t)n, indices.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)n, grad.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)(n*(n+1)/2), hess.size() );
  
  for (int i = 0; i < n; ++i) {
    get_nonideal_element( type, pd, i );

    qm->get_evaluations( pd, handles, false, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL_MESSAGE( "test only valid for element-based metrics", 
                                  (size_t)1, handles.size() );
    h = handles[0];
    
    rval = qm->evaluate_with_Hessian( pd, h, val_fixed, indices_fixed, grad_fixed, hess_fixed, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(rval);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices_fixed.size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, grad_fixed.size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, hess_fixed.size() );
    
    const size_t j = conn[i];
    CPPUNIT_ASSERT_DOUBLES_EQUAL( val, val_fixed, 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( grad[j], grad_fixed.front(), 1e-5 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( hess[n*j - j*(j-1)/2], hess_fixed.front(), 1e-4 );
  }
}

void QualityMetricTester::test_diagonal_with_fixed_vertex( EntityTopology type, 
                                                           QualityMetric* qm,
                                                           const Settings* settings )
{
  MsqPrintError err( std::cout );
  PatchData pd;
  if (settings) 
    pd.attach_settings(settings);
  else if (mSettings)
    pd.attach_settings(mSettings);
  std::vector<size_t> handles, indices, indices_fixed, conn;
  std::vector<Vector3D> grad, grad_fixed;
  std::vector<SymMatrix3D> hess, hess_fixed;
  double val, val_fixed;
  bool rval;
  const int n = TopologyInfo::corners( type );
  
  get_nonideal_element( type, pd );
  pd.element_by_index(0).get_vertex_indices( conn );
  
  qm->get_evaluations( pd, handles, false, err );
  ASSERT_NO_ERROR(err);
    // For sample-based metrics, it is difficult to set up such that
    // both surface and volume elements have only one sample point.
    // Make things easier by skipping element types with no sample points.
  if (handles.empty()) 
    return;
  CPPUNIT_ASSERT_EQUAL_MESSAGE( "test only valid for element-based metrics", 
                                (size_t)1, handles.size() );
  size_t h = handles[0];
  
  rval = qm->evaluate_with_Hessian_diagonal( pd, h, val, indices, grad, hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( (size_t)n, indices.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)n, grad.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)n, hess.size() );
  
  for (int i = 0; i < n; ++i) {
    get_nonideal_element( type, pd, i );

    qm->get_evaluations( pd, handles, false, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL_MESSAGE( "test only valid for element-based metrics", 
                                  (size_t)1, handles.size() );
    h = handles[0];
    
    rval = qm->evaluate_with_Hessian_diagonal( pd, h, val_fixed, indices_fixed, grad_fixed, hess_fixed, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(rval);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices_fixed.size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, grad_fixed.size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, hess_fixed.size() );
    
    const size_t j = conn[i];
    CPPUNIT_ASSERT_DOUBLES_EQUAL( val, val_fixed, 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( grad[j], grad_fixed.front(), 1e-5 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( hess[j], hess_fixed.front(), 1e-4 );
  }
}
