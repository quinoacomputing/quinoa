/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
/*!
  \file   MeshInterface.hpp
  \brief  

  \author Thomas Leurent
  \date   2003-04-17
*/

#include "meshfiles.h"

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "MeshImpl.hpp"
#include "Vector3D.hpp"

#include "cppunit/extensions/HelperMacros.h"

#include <iostream>
#include <algorithm>


using std::cout;
using std::cerr;
using std::endl;
using Mesquite::arrptr;

class MeshInterfaceTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(MeshInterfaceTest);
  CPPUNIT_TEST (test_get_geometric_dimension);
  CPPUNIT_TEST (test_vertices);
  //CPPUNIT_TEST (test_vertices_are_on_boundary);
  /* Not implemented yet
  CPPUNIT_TEST (test_vertex_is_fixed);
  */
  CPPUNIT_TEST (test_vertex_byte);
  CPPUNIT_TEST (test_vertex_get_attached_elements);
  CPPUNIT_TEST (test_elements);
  CPPUNIT_TEST (test_elements_get_attached_vertices);
  CPPUNIT_TEST (test_elements_get_topology);
  CPPUNIT_TEST(test_element_get_attached_vertex_indices);
  CPPUNIT_TEST_SUITE_END();

private:
  Mesquite::MeshImpl *mMesh; 
  std::vector<Mesquite::Mesh::VertexHandle> mConnectivity;
  std::vector<Mesquite::Mesh::VertexHandle> mVertices;
  std::vector<Mesquite::Mesh::ElementHandle> mElements;
  std::vector<size_t> mOffsets;
  Mesquite::MsqError mErr;
  
public:

  MeshInterfaceTest() 
    : mMesh(0) 
    {}

   /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
      // Read a VTK file -- 1 triangle flanked by 1 quad on each side (1 tri + 3 quads)
    mMesh = new Mesquite::MeshImpl;
    mMesh->read_vtk(MESH_FILES_DIR "2D/vtk/mixed/untangled/hybrid_3quad_1tri.vtk", mErr);
    CPPUNIT_ASSERT(!mErr);

      // Get mesh data
    mMesh->get_all_elements( mElements, mErr );
    CPPUNIT_ASSERT(!mErr);
    mMesh->elements_get_attached_vertices( arrptr(mElements),
                                           mElements.size(),
                                           mConnectivity,
                                           mOffsets,
                                           mErr );
    CPPUNIT_ASSERT(!mErr);
    
      // Construct list of vertices w/out duplicates from
      // connectivity list.
    std::vector<Mesquite::Mesh::VertexHandle>::iterator new_end;
    mVertices = mConnectivity;
    std::sort( mVertices.begin(), mVertices.end() );
    new_end = std::unique( mVertices.begin(), mVertices.end() );
    mVertices.resize( new_end - mVertices.begin() );
  }
  
    // Automatically called by CppUnit after each test function.
  void tearDown()
  {

    delete mMesh;
    if(mErr) cout << mErr << endl;
    mErr.clear();
  }
  
public:
  
  void test_get_geometric_dimension()
  {
    int d = mMesh->get_geometric_dimension(mErr);
    CPPUNIT_ASSERT_EQUAL(d,3);
  }

  
  void test_vertices()
  {
    size_t nbVert = mVertices.size();
    CPPUNIT_ASSERT_EQUAL(9,(int)nbVert);

    Mesquite::MsqVertex correct_coords[9], coords[9];
    correct_coords[0].set(1,0,0);
    correct_coords[1].set(0,1.732,0);
    correct_coords[2].set(-1,0,0);
    correct_coords[3].set(-1,-2,0);
    correct_coords[4].set(1,-2,0);
    correct_coords[5].set(2.732,1,0);
    correct_coords[6].set(1.732,2.732,0);
    correct_coords[7].set(-1.732,2.732,0);
    correct_coords[8].set(-2.732,1,0);

    mMesh->vertices_get_coordinates(arrptr(mVertices), coords, nbVert, mErr);
    CPPUNIT_ASSERT(!mErr);
    for (size_t i=0; i<nbVert; ++i) {
      for (int j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(coords[i][j], correct_coords[i][j], .01);
    }

    coords[3].set(2.,3.,4.);
    mMesh->vertex_set_coordinates(mVertices[3], coords[3], mErr);
    CPPUNIT_ASSERT(!mErr);
    Mesquite::MsqVertex coords_2;
    mMesh->vertices_get_coordinates(&mVertices[3], &coords_2, 1, mErr);
    CPPUNIT_ASSERT(!mErr);
    for (int j=0; j<3; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(coords[3][j], coords_2[j], 1e-6);
  }

//  void test_vertices_are_on_boundary()
//  {
//    bool on_bnd[9];
//    mMesh->vertices_are_on_boundary(mVertices, on_bnd, nbVert, mErr);
//    CPPUNIT_ASSERT(!mErr);
//    bool correct_boundary[9] = {false, false, false, true, true, true, true, true, true};
//    for (size_t i=0; i<nbVert; ++i) {
//      CPPUNIT_ASSERT(on_bnd[i] == correct_boundary[i]);
//    }
//  }

  void test_vertex_is_fixed()
  {
    size_t nbVert = mVertices.size();
    bool correct_fixed[9] = {false, false, false, true, true, true, true, true, true};
    std::vector<bool> fixed;
    mMesh->vertices_get_fixed_flag( arrptr(mVertices), fixed, 9, mErr );
    CPPUNIT_ASSERT(!mErr);
    CPPUNIT_ASSERT_EQUAL( (size_t)8, fixed.size() );
    for (size_t i=0; i<nbVert; ++i) {
      CPPUNIT_ASSERT_EQUAL(correct_fixed[i], (bool)fixed[i]);
    }
  }

  void test_vertex_byte()
  {
    size_t nbVert = mVertices.size();
    size_t i;
	unsigned char* bytes = new unsigned char[nbVert];
    mMesh->vertices_get_byte(arrptr(mVertices), bytes, nbVert, mErr); 
    CPPUNIT_ASSERT(!mErr);

    // Asserts all vertex bytes are initialised to 0. 
    for (i=0; i<nbVert; ++i)
      CPPUNIT_ASSERT(bytes[i] == 0);

    // Test various vertex byte read / write routines.
    bytes[3] |= 4;
    mMesh->vertices_set_byte(arrptr(mVertices), bytes, nbVert, mErr); 
    CPPUNIT_ASSERT(!mErr);
    mMesh->vertex_set_byte(mVertices[5], 8, mErr); 
    CPPUNIT_ASSERT(!mErr);
    unsigned char byte;
    mMesh->vertex_get_byte(mVertices[3], &byte, mErr);
    CPPUNIT_ASSERT(!mErr);
    CPPUNIT_ASSERT(byte == 4);
    mMesh->vertices_get_byte(arrptr(mVertices), bytes, nbVert, mErr);
    CPPUNIT_ASSERT(!mErr);
    for (i=0; i<nbVert; ++i) {
      if (i==3)
        CPPUNIT_ASSERT(bytes[i] == 4);
      else if (i==5)
        CPPUNIT_ASSERT(bytes[i] == 8);
      else
        CPPUNIT_ASSERT(bytes[i] == 0);
    }

    delete [] bytes;
  }

  
  void test_vertex_get_attached_elements()
  {
	  size_t i;
    const size_t nbVert = mVertices.size();

    std::vector<Mesquite::Mesh::ElementHandle> elements;
    std::vector<size_t> offsets;
    mMesh->vertices_get_attached_elements( arrptr(mVertices),
                                           mVertices.size(),
                                           elements,
                                           offsets,
                                           mErr );
    CPPUNIT_ASSERT(!mErr);
    CPPUNIT_ASSERT_EQUAL(offsets.size(), mVertices.size() + 1);

    // checks we have 6 vertices contained in 1 element only
    // and 3 vertices contained in 3 elements. 
    int n1=0;
    int n3=0;
    for (i = 1; i <= mVertices.size(); ++i)
    {
      const size_t nev = offsets[i] - offsets[i-1];
      if (nev==1)
        ++n1;
      else if (nev==3)
        ++n3;
      else // failure. 
        CPPUNIT_ASSERT(false);
    }
    CPPUNIT_ASSERT(n1==6);
    CPPUNIT_ASSERT(n3==3);

    // gets the index of a vertex in a corner
    int one_corner_vertex_index=0;
    for (i = 0; i < nbVert; ++i)
    {
      const size_t nev = offsets[i+1] - offsets[i];
      if (1 == nev)
        break;
    }
    CPPUNIT_ASSERT( i < nbVert );
    one_corner_vertex_index = i;

    // retrieve the attached element.
    // This is a poor test.  We allow an element handle of zero,
    // and just testing that the function returned something isn't
    // that useful. - J.K.
    //Mesquite::Mesh::ElementHandle elem=0;
    //mMesh->vertex_get_attached_elements(mVertices[one_corner_vertex_index], &elem, 1, mErr);
    //CPPUNIT_ASSERT(!mErr);
    //CPPUNIT_ASSERT(elem!=0);

  }
  
  void test_elements()
  {
    CPPUNIT_ASSERT_EQUAL(4,(int)mElements.size());

    
  }


  void test_elements_get_attached_vertices()
  {
    const size_t nbElem = mElements.size();
     
    // checks we have 3 elements containing 4 vertices 
    // and 1 element containing 3 vertices. 
    int n3=0;
    int n4=0;
	size_t i;
    for (i=0; i<nbElem; ++i) {
      size_t nve = mOffsets[i+1] - mOffsets[i];
      if (nve==3)
        ++n3;
      else if (nve==4)
        ++n4;
      else // failure. 
        CPPUNIT_ASSERT(false);
    }
    CPPUNIT_ASSERT(n3==1); // 1 triangle
    CPPUNIT_ASSERT(n4==3); // 3 quads

    
    // Make sure CSR data is valid
    std::map<Mesquite::Mesh::VertexHandle,int> vtx_repeated_occurence;
    for (i = 0; i < mVertices.size(); ++i)
      vtx_repeated_occurence[mVertices[i]] = 0;
    for (i = 0; i < mConnectivity.size(); ++i)
      ++vtx_repeated_occurence[mConnectivity[i]];
    for (i=0; i<9; ++i) {
      CPPUNIT_ASSERT( vtx_repeated_occurence[mVertices[i]] <= 3 );
    }

    // Makes sure CSR offsets are valid
    CPPUNIT_ASSERT( mOffsets[0] == 0 );
    CPPUNIT_ASSERT( mOffsets[1] >=3 && mOffsets[1]<=12 );
    CPPUNIT_ASSERT( mOffsets[2] >=3 && mOffsets[2]<=12 );
    CPPUNIT_ASSERT( mOffsets[3] >=3 && mOffsets[3]<=12 );
    CPPUNIT_ASSERT( mOffsets[4] == 15 );
  }

  
  void test_elements_get_topology()
  {
    const size_t nbElem = mElements.size();
    int nb_quads=0;
    int nb_tri=0;
    Mesquite::EntityTopology* topos = new Mesquite::EntityTopology[nbElem];
    mMesh->elements_get_topologies(arrptr(mElements), topos, nbElem, mErr);
    CPPUNIT_ASSERT(!mErr);
    for (size_t i=0; i<nbElem; ++i) {
      switch (topos[i]) {
      case Mesquite::TRIANGLE:
        ++nb_tri;
        break;
      case Mesquite::QUADRILATERAL:
        ++nb_quads;
        break;
      default:
        CPPUNIT_FAIL("Topology should be quad or Hex only.");
      }
    }
    CPPUNIT_ASSERT_EQUAL(1,nb_tri);
    CPPUNIT_ASSERT_EQUAL(3,nb_quads);
    delete []topos;
  }



  void test_element_get_attached_vertex_indices()
  {
    // Find the index of the triangle
    Mesquite::EntityTopology topo=Mesquite::MIXED;
    int tri_index = -1;
    while (topo != Mesquite::TRIANGLE) {
      ++tri_index;
      CPPUNIT_ASSERT((unsigned)tri_index < mElements.size());
      Mesquite::Mesh::ElementHandle handle = mElements[tri_index];
      mMesh->elements_get_topologies(&handle, &topo, 1, mErr);
      CPPUNIT_ASSERT(!mErr);
    }

    // creates list with correct vertices coordinates for the triangle
    std::vector<Mesquite::Vector3D> correct_coords;
    correct_coords.push_back(Mesquite::Vector3D(1.,0.,0.));
    correct_coords.push_back(Mesquite::Vector3D(0.,1.732050807,0.));
    correct_coords.push_back(Mesquite::Vector3D(-1.,0.,0.));

    // Creates same list from the mesh implementation
    std::vector<Mesquite::MsqVertex> tri_coords(3);
    mMesh->vertices_get_coordinates(&mConnectivity[mOffsets[tri_index]],
                                    arrptr(tri_coords), 3, mErr );
    CPPUNIT_ASSERT(!mErr);

    // Makes sure both list contain the same elements (not necessarily in the same order).
    std::vector<Mesquite::Vector3D>::iterator correct_iter;
    std::vector<Mesquite::MsqVertex>::iterator tri_iter;
    for (tri_iter = tri_coords.begin(); tri_iter != tri_coords.end(); ++tri_iter)
    {
      for (correct_iter = correct_coords.begin(); 
           correct_iter != correct_coords.end(); 
           ++correct_iter)
      {
        if (Mesquite::Vector3D::distance_between(*tri_iter, *correct_iter) < 10e-4)   
          break;
      }
      
      // check if a match was found
      CPPUNIT_ASSERT( correct_iter != correct_coords.end() );
      
      // remove match from list
      correct_coords.erase( correct_iter );
    }
    CPPUNIT_ASSERT(correct_coords.empty());
  }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshInterfaceTest, "MeshInterfaceTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshInterfaceTest, "Unit");
