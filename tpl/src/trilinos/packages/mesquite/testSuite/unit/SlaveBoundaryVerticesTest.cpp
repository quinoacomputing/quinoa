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


/** \file SlaveBoundaryVerticesTest.cpp
 *  \brief Test SlaveBoundaryVertices class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "UnitUtil.hpp"
#include "meshfiles.h"
#include "Mesquite_Settings.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_SlaveBoundaryVertices.hpp"

#include "Mesquite_MsqVertex.hpp"
#include "Mesquite_DomainClassifier.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_MeshDomain1D.hpp"
#include "Mesquite_MeshImpl.hpp"

#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace Mesquite;

class SlaveBoundaryVerticesTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(SlaveBoundaryVerticesTest);
  CPPUNIT_TEST(test_fail_if_slaves_not_calculated);
  CPPUNIT_TEST(test_one_depth_from_fixed);
  CPPUNIT_TEST(test_two_depth_from_fixed);
  CPPUNIT_TEST(test_zero_depth_from_surface);
  CPPUNIT_TEST(test_one_depth_from_surface);
  CPPUNIT_TEST(test_two_depth_from_surface);
  CPPUNIT_TEST(test_zero_depth_from_curve);
  CPPUNIT_TEST(test_one_depth_from_curve);
  CPPUNIT_TEST(test_two_depth_from_curve);
  CPPUNIT_TEST_SUITE_END();

  void test_slaved_common( unsigned depth, unsigned boundary ); 

  void make_mesh( MeshImpl& mesh,
                  DomainClassifier& domain,
                  const int intervals );

public:

  void test_fail_if_slaves_not_calculated();
  void test_one_depth_from_fixed()          { test_slaved_common( 1, 4 ); }
  void test_two_depth_from_fixed()          { test_slaved_common( 2, 4 ); }
  void test_zero_depth_from_surface()       { test_slaved_common( 0, 2 ); }
  void test_one_depth_from_surface()        { test_slaved_common( 1, 2 ); }
  void test_two_depth_from_surface()        { test_slaved_common( 2, 2 ); }
  void test_zero_depth_from_curve()         { test_slaved_common( 0, 1 ); }
  void test_one_depth_from_curve()          { test_slaved_common( 1, 1 ); }
  void test_two_depth_from_curve()          { test_slaved_common( 2, 1 ); }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SlaveBoundaryVerticesTest, "SlaveBoundaryVerticesTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SlaveBoundaryVerticesTest, "Unit");

void SlaveBoundaryVerticesTest::test_fail_if_slaves_not_calculated()
{
  MsqError err;

  Settings settings;
  settings.set_slaved_ho_node_mode( Settings::SLAVE_ALL );
  SlaveBoundaryVertices tool( 1 );
  
  MeshImpl mesh;
  DomainClassifier domain;
  make_mesh( mesh, domain, 2 );
    
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &domain);
  tool.loop_over_mesh( &mesh_and_domain, &settings, err );
  CPPUNIT_ASSERT(err);
  err.clear();
}

void SlaveBoundaryVerticesTest::make_mesh( MeshImpl& mesh,
                                           DomainClassifier& domain,
                                           const int intervals )
{
  MsqPrintError err(std::cerr);
  const char input_file[] = MESH_FILES_DIR "3D/vtk/quadratic/6x6x6-hex20.vtk";
  
  const Vector3D min( -3, -3, -3 );
  const Vector3D max(  3,  3,  3 );
  const Vector3D space = (max - min) * 1.0/(intervals);
  
  mesh.clear();
  mesh.read_vtk( input_file, err );
  ASSERT_NO_ERROR(err);
  
  const Vector3D corners[8] = { Vector3D(min[0],min[1],min[2]),
                                Vector3D(max[0],min[1],min[2]),
                                Vector3D(max[0],max[1],min[2]),
                                Vector3D(min[0],max[1],min[2]),
                                Vector3D(min[0],min[1],max[2]),
                                Vector3D(max[0],min[1],max[2]),
                                Vector3D(max[0],max[1],max[2]),
                                Vector3D(min[0],max[1],max[2]) };
  
  MeshDomain* subdomains[26] = { 
                    new PlanarDomain( PlanarDomain::XZ, min[1] ),
                    new PlanarDomain( PlanarDomain::YZ, max[0] ),
                    new PlanarDomain( PlanarDomain::XZ, max[1] ),
                    new PlanarDomain( PlanarDomain::YZ, min[0] ),
                    new PlanarDomain( PlanarDomain::XY, min[2] ),
                    new PlanarDomain( PlanarDomain::XY, max[2] ),
                    new LineDomain( corners[0], corners[1] - corners[0] ),
                    new LineDomain( corners[1], corners[2] - corners[1] ),
                    new LineDomain( corners[2], corners[3] - corners[2] ),
                    new LineDomain( corners[3], corners[0] - corners[3] ),
                    new LineDomain( corners[0], corners[4] - corners[0] ),
                    new LineDomain( corners[1], corners[5] - corners[1] ),
                    new LineDomain( corners[2], corners[6] - corners[0] ),
                    new LineDomain( corners[3], corners[7] - corners[1] ),
                    new LineDomain( corners[4], corners[5] - corners[4] ),
                    new LineDomain( corners[5], corners[6] - corners[5] ),
                    new LineDomain( corners[6], corners[7] - corners[6] ),
                    new LineDomain( corners[7], corners[4] - corners[7] ),
                    new PointDomain( corners[0] ),
                    new PointDomain( corners[1] ),
                    new PointDomain( corners[2] ),
                    new PointDomain( corners[3] ),
                    new PointDomain( corners[4] ),
                    new PointDomain( corners[5] ),
                    new PointDomain( corners[6] ),
                    new PointDomain( corners[7] ) };
  const int subdims[26] = { 2, 2, 2, 2, 2, 2,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                            0, 0, 0, 0, 0, 0, 0, 0 };
  DomainClassifier::classify_skin_geometrically( domain, &mesh, 1e-6,
                                                 subdomains, subdims, 26,
                                                 err );
  domain.delete_sub_domains(true);
  ASSERT_NO_ERROR(err);
}

void SlaveBoundaryVerticesTest::test_slaved_common( unsigned depth, unsigned boundary )
{
  MeshImpl mesh;
  DomainClassifier domain;
  make_mesh( mesh, domain, 2*depth+2 );

  MsqPrintError err(std::cerr);
  std::vector< std::vector<Mesh::VertexHandle> > depths(depth+1);
  std::set<Mesh::VertexHandle> non_slave;
  std::set<Mesh::VertexHandle>::iterator p;

    // find boundary vertices
  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!verts.empty());
  if (boundary >= 4) {
    std::vector<bool> flags;
    mesh.vertices_get_fixed_flag( arrptr(verts), flags, verts.size(), err );
    ASSERT_NO_ERROR(err);
    for (size_t i = 0; i < verts.size(); ++i)
      if (flags[i]) {
        depths[0].push_back( verts[i] );
        non_slave.insert( verts[i] );
      }
  }
  else {
    std::vector<unsigned short> dim(verts.size());
    domain.domain_DoF( arrptr(verts), arrptr(dim), verts.size(), err );
    ASSERT_NO_ERROR(err);
    for (size_t i = 0; i < verts.size(); ++i)
      if (dim[i] <= boundary) {
        depths[0].push_back( verts[i] );
        non_slave.insert( verts[i] );
      }
  }
  
    // check that our input is usable for this test
  CPPUNIT_ASSERT( !verts.empty() );
  
    // find all vertices up to specified depth
  for (unsigned d = 0; d < depth; ++d) {
    for (size_t i = 0; i < depths[d].size(); ++i) {
      std::vector<Mesh::ElementHandle> adj;
      std::vector<size_t> junk;
      mesh.vertices_get_attached_elements( &depths[d][i], 1, adj, junk, err );
      ASSERT_NO_ERROR(err);
      for(size_t j = 0; j < adj.size(); ++j) {
        junk.clear();
        std::vector<Mesh::VertexHandle> conn;
        mesh.elements_get_attached_vertices( &adj[j], 1, conn, junk, err );
        ASSERT_NO_ERROR(err);
        for (size_t k = 0; k < conn.size(); ++k) {
          p = non_slave.find(conn[k]);
          if (p == non_slave.end()) {
            non_slave.insert( p, conn[k] );
            depths[d+1].push_back( conn[k] );
          }
        }
      }
    }
  }
  
    // Check that our input is usable for this test:
    // Should have some vertices that are not within the specified depth of 
    // the boundary.
  CPPUNIT_ASSERT( non_slave.size() < verts.size() );
  
    // Now build a map of all higher-order nodes in the mesh
  std::set<Mesh::VertexHandle> higher_order;
  std::vector<Mesh::ElementHandle> elems;
  mesh.get_all_elements( elems, err ); 
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!elems.empty());
  std::vector<EntityTopology> types(elems.size());
  mesh.elements_get_topologies( arrptr(elems), arrptr(types), elems.size(), err );
  ASSERT_NO_ERROR(err);
  for (size_t i = 0; i < elems.size(); ++i) {
    std::vector<Mesh::VertexHandle> conn;
    std::vector<size_t> junk;
    mesh.elements_get_attached_vertices( &elems[i], 1, conn, junk, err );
    ASSERT_NO_ERROR(err);
    for (size_t j = TopologyInfo::corners( types[i] ); j < conn.size(); ++j)
      higher_order.insert( conn[j] );
  }
  
    // Check that our input is usable for this test:
    // Should have some higher-order vertices
  CPPUNIT_ASSERT( !higher_order.empty() );
  
    // Now build a map of all fixed vertices
  std::set<Mesh::VertexHandle> fixed_vertices;
  std::vector<bool> fixed;
  mesh.vertices_get_fixed_flag( arrptr(verts), fixed, verts.size(), err );
  ASSERT_NO_ERROR(err);
  for (size_t i = 0; i < verts.size(); ++i)
    if (fixed[i])
      fixed_vertices.insert( verts[i] );

    // Now actually run the tool
  Settings settings;
  settings.set_slaved_ho_node_mode( Settings::SLAVE_CALCULATED );
  SlaveBoundaryVertices tool( depth, boundary );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &domain);
  tool.loop_over_mesh( &mesh_and_domain, &settings, err );
  ASSERT_NO_ERROR(err);
  
    // Now verify the results
  std::vector<unsigned char> bytes( verts.size() );
  mesh.vertices_get_byte( arrptr(verts), arrptr(bytes), verts.size(), err );
  ASSERT_NO_ERROR(err);
  for (size_t i = 0; i < verts.size(); ++i) {
    bool in_non_slave = (non_slave.find( verts[i] ) != non_slave.end());
    bool in_fixed = (fixed_vertices.find( verts[i] ) != fixed_vertices.end());
    bool in_higher_order = (higher_order.find( verts[i] ) != higher_order.end());
    if (bytes[i] & MsqVertex::MSQ_DEPENDENT) { // if slave node
        // must not be within 'depth' of boundary
      CPPUNIT_ASSERT( !in_non_slave );
        // must be a higher-order vertex
      CPPUNIT_ASSERT( in_higher_order );
        // must not be fixed
      CPPUNIT_ASSERT( !in_fixed );
    }
    else {
        // there are three reasons that a vertex isn't slaved
      bool in_non_slave = (non_slave.find( verts[i] ) != non_slave.end());
      bool in_fixed = (fixed_vertices.find( verts[i] ) != fixed_vertices.end());
      bool in_higher_order = (higher_order.find( verts[i] ) != higher_order.end());
      CPPUNIT_ASSERT( in_fixed || !in_higher_order || in_non_slave );
    }
  }
}

