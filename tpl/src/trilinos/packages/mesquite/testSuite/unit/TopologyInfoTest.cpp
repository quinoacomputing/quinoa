/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#include <algorithm>
#include <iostream>
using std::cout;
#include "cppunit/extensions/HelperMacros.h"

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "TopologyInfo.hpp"

using namespace Mesquite;

class TopologyInfoTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(TopologyInfoTest);

   CPPUNIT_TEST (tri);
   CPPUNIT_TEST (tri3);
   CPPUNIT_TEST (tri4);
   CPPUNIT_TEST (tri6);
   CPPUNIT_TEST (tri7);
   
   CPPUNIT_TEST (quad);
   CPPUNIT_TEST (quad4);
   CPPUNIT_TEST (quad5);
   CPPUNIT_TEST (quad8);
   CPPUNIT_TEST (quad9);
   
   CPPUNIT_TEST (tet);
   CPPUNIT_TEST (tet4);
   CPPUNIT_TEST (tet5);
   CPPUNIT_TEST (tet8);
   CPPUNIT_TEST (tet9);
   CPPUNIT_TEST (tet10);
   CPPUNIT_TEST (tet11);
   CPPUNIT_TEST (tet14);
   CPPUNIT_TEST (tet15);

   CPPUNIT_TEST (hex);
   CPPUNIT_TEST (hex8);
   CPPUNIT_TEST (hex9);
   CPPUNIT_TEST (hex14);
   CPPUNIT_TEST (hex15);
   CPPUNIT_TEST (hex20);
   CPPUNIT_TEST (hex21);
   CPPUNIT_TEST (hex26);
   CPPUNIT_TEST (hex27);
 
   CPPUNIT_TEST (pyramid);
   CPPUNIT_TEST (pyramid5);
   CPPUNIT_TEST (pyramid13);
   
   CPPUNIT_TEST (wedge);
   CPPUNIT_TEST (wedge6);
   CPPUNIT_TEST (wedge15);

   CPPUNIT_TEST (polygon);
   CPPUNIT_TEST (polyhedron);
   
   CPPUNIT_TEST (bad_type);
   
   CPPUNIT_TEST (tri_adj_vert);
   CPPUNIT_TEST (quad_adj_vert);
   CPPUNIT_TEST (tet_adj_vert);
   CPPUNIT_TEST (hex_adj_vert);
   CPPUNIT_TEST (pyr_adj_vert);
   CPPUNIT_TEST (wdg_adj_vert);
   
   CPPUNIT_TEST (tri_rev_adj_vert);
   CPPUNIT_TEST (quad_rev_adj_vert);
   CPPUNIT_TEST (tet_rev_adj_vert);
   CPPUNIT_TEST (hex_rev_adj_vert);
   CPPUNIT_TEST (pyr_rev_adj_vert);
   CPPUNIT_TEST (wdg_rev_adj_vert);
   
   CPPUNIT_TEST (higher_order_from_side);
   CPPUNIT_TEST (side_from_higher_order);
   CPPUNIT_TEST (find_edge);
   CPPUNIT_TEST (find_face);
   CPPUNIT_TEST (find_side);
   CPPUNIT_TEST (compare_sides);

   CPPUNIT_TEST_SUITE_END();

public:

  void setUp() {}
  
  void tearDown() {}
  
  TopologyInfoTest() {}
  
  bool compare_edge( const unsigned* a, const unsigned* b );
  
  bool compare_face( unsigned len, const unsigned* a, const unsigned* b );
  
  bool compare_vol( unsigned len, const unsigned* a, const unsigned* b );
  
  void test_face_elem( EntityTopology topo, 
                       unsigned num_nodes,
                       unsigned num_sides );

  
  void test_vol_elem( EntityTopology topo, 
                      unsigned num_nodes,
                      unsigned num_verts,
                      unsigned num_edges,
                      unsigned num_faces );
  
  void test_poly( EntityTopology topo );
  
  void tri();
  void tri3();
  void tri4();
  void tri6();
  void tri7();
  
  void quad();
  void quad4();
  void quad5();
  void quad8();
  void quad9();
  
  void tet();
  void tet4();
  void tet5();
  void tet8();
  void tet9();
  void tet10();
  void tet11();
  void tet14();
  void tet15();

  void hex();
  void hex8();
  void hex9();
  void hex14();
  void hex15();
  void hex20();
  void hex21();
  void hex26();
  void hex27();
  
  void pyramid();
  void pyramid5();
  void pyramid13();
  
  void wedge();
  void wedge6();
  void wedge15();
  
  void polygon();
  void polyhedron();

  void bad_type();
  
  void tri_adj_vert();
  void quad_adj_vert();
  void tet_adj_vert();
  void hex_adj_vert();
  void pyr_adj_vert();
  void wdg_adj_vert();
  void test_adj( Mesquite::EntityTopology, const unsigned expected[][5] );
  
  void tri_rev_adj_vert()  { test_rev_adj( Mesquite::TRIANGLE      ); }
  void quad_rev_adj_vert() { test_rev_adj( Mesquite::QUADRILATERAL ); }
  void tet_rev_adj_vert()  { test_rev_adj( Mesquite::TETRAHEDRON   ); }
  void hex_rev_adj_vert()  { test_rev_adj( Mesquite::HEXAHEDRON    ); }
  void pyr_rev_adj_vert()  { test_rev_adj( Mesquite::PYRAMID       ); }
  void wdg_rev_adj_vert()  { test_rev_adj( Mesquite::PRISM         ); }
  void test_rev_adj( Mesquite::EntityTopology type );
   
  void higher_order_from_side();
  void side_from_higher_order();
  void find_edge();
  void find_face();
  void find_side();
  void compare_sides();
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TopologyInfoTest, "TopologyInfoTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TopologyInfoTest, "Unit");


  
bool TopologyInfoTest::compare_edge( const unsigned* a, const unsigned* b )
  { return (a[0] == b[0] && a[1] == b[1]) ||
           (a[0] == b[1] && a[1] == b[0]); }

bool TopologyInfoTest::compare_face( unsigned len, const unsigned* a, const unsigned* b )
{
  unsigned i, j;
  for (i = 0; i < len; ++i)
  {
    for (j = 0; j < len; ++j)
    {
      if (a[j] != b[(i+j)%len])
        break;
    }
    if (j == len)
      return true;
  }
  return false;
}

bool TopologyInfoTest::compare_vol( unsigned len, const unsigned* a, const unsigned* b )
{
  for (unsigned i = 0; i < len; ++i)
    if (a[i] != b[i])
      return false;
  return true;
}

void TopologyInfoTest::test_face_elem( EntityTopology topo, 
                     unsigned num_nodes,
                     unsigned num_sides )
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  unsigned index = 0;
  TopologyInfo::higher_order( topo, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!vol);

  unsigned nodes = num_sides + edge * num_sides + face;
  CPPUNIT_ASSERT( num_nodes == nodes );

  unsigned side, dim;
  for (index = 0; index < num_sides; ++index)
  {
    TopologyInfo::side_number( topo, num_nodes, index, dim, side, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(dim == 0 && side == index);
  }

  if (edge)
  {
    for (unsigned s = 0; s < num_sides; ++s)
    {
      TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(dim == 1 && side == s);
    }
  }
  if (face)
  {
    TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(dim == 2 && side == 0);
  }
}


void TopologyInfoTest::test_vol_elem( EntityTopology topo, 
                    unsigned num_nodes,
                    unsigned num_verts,
                    unsigned num_edges,
                    unsigned num_faces )
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  unsigned index = 0;
  TopologyInfo::higher_order( topo, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);

  unsigned nodes = num_verts + edge * num_edges + face * num_faces + vol;
  CPPUNIT_ASSERT( num_nodes == nodes );

  unsigned side, dim;
  for (index = 0; index < num_verts; ++index)
  {
    TopologyInfo::side_number( topo, num_nodes, index, dim, side, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(dim == 0 && side == index);
  }

  if (edge)
  {
    for (unsigned s = 0; s < num_edges; ++s)
    {
      TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(dim == 1 && side == s);
    }
  }
  if (face)
  {
    for (unsigned s = 0; s < num_faces; ++s)
    {
      TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(dim == 2 && side == s);
    }
  }
  if (vol)
  {
    TopologyInfo::side_number( topo, num_nodes, index++, dim, side, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(dim == 3 && side == 0);
  }
}





void TopologyInfoTest::tri()
{
  MsqPrintError err(cout);

  CPPUNIT_ASSERT (2 == TopologyInfo::dimension( TRIANGLE ));
  CPPUNIT_ASSERT (3 == TopologyInfo::adjacent( TRIANGLE, 1 ));
  CPPUNIT_ASSERT (3 == TopologyInfo::adjacent( TRIANGLE, 0 ));
  CPPUNIT_ASSERT (1 == TopologyInfo::adjacent( TRIANGLE, 2 ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( TRIANGLE, 3 ));

  CPPUNIT_ASSERT (3 == TopologyInfo::sides( TRIANGLE ));
  CPPUNIT_ASSERT (3 == TopologyInfo::corners( TRIANGLE ));
  CPPUNIT_ASSERT (3 == TopologyInfo::edges( TRIANGLE ));
  CPPUNIT_ASSERT (1 == TopologyInfo::faces( TRIANGLE ));

  const unsigned num_edges = 3;
  const unsigned* side;
  const unsigned edges[num_edges][2] = { {0, 1}, {1, 2}, {2, 0} };
  unsigned count;
  const unsigned face[] = { 0, 1, 2 };

  for (unsigned i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( TRIANGLE, i, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( TRIANGLE, 1, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  }

  side = TopologyInfo::side_vertices( TRIANGLE, 2, 0, count, err ); 
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(3 == count);
  CPPUNIT_ASSERT( compare_face( 3, side, face ) );
}

void TopologyInfoTest::tri3()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 3;

  TopologyInfo::higher_order( TRIANGLE, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);

  test_face_elem( TRIANGLE, num_nodes, 3 );
}

void TopologyInfoTest::tri4()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 4;

  TopologyInfo::higher_order( TRIANGLE, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( TRIANGLE, num_nodes, 3 );
}

void TopologyInfoTest::tri6()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 6;

  TopologyInfo::higher_order( TRIANGLE, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( TRIANGLE, num_nodes, 3 );
}

void TopologyInfoTest::tri7()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 7;

  TopologyInfo::higher_order( TRIANGLE, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( TRIANGLE, num_nodes, 3 );
}


void TopologyInfoTest::quad()
{
  MsqPrintError err(cout);

  CPPUNIT_ASSERT (2 == TopologyInfo::dimension( QUADRILATERAL ));
  CPPUNIT_ASSERT (4 == TopologyInfo::adjacent( QUADRILATERAL, 1 ));
  CPPUNIT_ASSERT (4 == TopologyInfo::adjacent( QUADRILATERAL, 0 ));
  CPPUNIT_ASSERT (1 == TopologyInfo::adjacent( QUADRILATERAL, 2 ));
  CPPUNIT_ASSERT (0 == TopologyInfo::adjacent( QUADRILATERAL, 3 ));

  CPPUNIT_ASSERT (4 == TopologyInfo::sides( QUADRILATERAL ));
  CPPUNIT_ASSERT (4 == TopologyInfo::corners( QUADRILATERAL ));
  CPPUNIT_ASSERT (4 == TopologyInfo::edges( QUADRILATERAL ));
  CPPUNIT_ASSERT (1 == TopologyInfo::faces( QUADRILATERAL ));

  const unsigned num_edges = 4;
  const unsigned* side;
  const unsigned edges[num_edges][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
  unsigned count;
  const unsigned face[] = { 0, 1, 2, 3 };

  for (unsigned i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( QUADRILATERAL, i, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( QUADRILATERAL, 1, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  }

  side = TopologyInfo::side_vertices( QUADRILATERAL, 2, 0, count, err ); 
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( 4 == count );
  CPPUNIT_ASSERT( compare_face( 4, side, face ) );
}

void TopologyInfoTest::quad4()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 4;

  TopologyInfo::higher_order( QUADRILATERAL, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( QUADRILATERAL, num_nodes, 4 );
}

void TopologyInfoTest::quad5()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 5;

  TopologyInfo::higher_order( QUADRILATERAL, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( QUADRILATERAL, num_nodes, 4 );
}

void TopologyInfoTest::quad8()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 8;

  TopologyInfo::higher_order( QUADRILATERAL, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( QUADRILATERAL, num_nodes, 4 );
}

void TopologyInfoTest::quad9()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 9;

  TopologyInfo::higher_order( QUADRILATERAL, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_face_elem( QUADRILATERAL, num_nodes, 4 );
}


void TopologyInfoTest::tet()
{
  MsqPrintError err(cout);

  const unsigned num_verts = 4;
  const unsigned num_edges = 6;
  const unsigned num_faces = 4;
  CPPUNIT_ASSERT (3 == TopologyInfo::dimension( TETRAHEDRON ));
  CPPUNIT_ASSERT (1 == TopologyInfo::adjacent( TETRAHEDRON, 3 ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::adjacent( TETRAHEDRON, 2 ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::adjacent( TETRAHEDRON, 1 ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::adjacent( TETRAHEDRON, 0 ));

  CPPUNIT_ASSERT (num_faces == TopologyInfo::sides( TETRAHEDRON ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::corners( TETRAHEDRON ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::edges( TETRAHEDRON ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::faces( TETRAHEDRON ));

  const unsigned* side;
  unsigned i, count;
  const unsigned vert_per_face = 3;
  unsigned edges[num_edges][2] = { {0, 1}, {1, 2}, {2, 0},
                                   {0, 3}, {1, 3}, {2, 3} };
  unsigned faces[num_faces][vert_per_face] = { { 0, 1, 3 }, { 1, 2, 3 }, 
                                               { 2, 0, 3 }, { 2, 1, 0 } };

  for (i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( TETRAHEDRON, i, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( TETRAHEDRON, 1, i, count, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  } 

  for (i = 0; i < num_faces; ++i)
  {
    side = TopologyInfo::face_vertices( TETRAHEDRON, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == vert_per_face);
    CPPUNIT_ASSERT( compare_face( vert_per_face, side, faces[i] ) );

    side = TopologyInfo::side_vertices( TETRAHEDRON, 2, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == vert_per_face);
    CPPUNIT_ASSERT( compare_face( vert_per_face, side, faces[i] ) );
  } 

  side = TopologyInfo::side_vertices( TETRAHEDRON, 3, 0, count, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(count == num_verts);
  for (i = 0; i < num_verts; ++i)
    CPPUNIT_ASSERT(side[i] == i);
}

void TopologyInfoTest::tet4()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 4;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet5()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 5;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet8()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 8;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet9()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 9;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet10()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 10;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet11()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 11;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet14()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 14;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}

void TopologyInfoTest::tet15()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 15;

  TopologyInfo::higher_order( TETRAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( TETRAHEDRON, num_nodes, 4, 6, 4 );
}


void TopologyInfoTest::hex()
{
  MsqPrintError err(cout);

  const unsigned num_verts = 8;
  const unsigned num_edges = 12;
  const unsigned num_faces = 6;
  CPPUNIT_ASSERT (3 == TopologyInfo::dimension( HEXAHEDRON ));
  CPPUNIT_ASSERT (1 == TopologyInfo::adjacent( HEXAHEDRON, 3 ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::adjacent( HEXAHEDRON, 2 ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::adjacent( HEXAHEDRON, 1 ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::adjacent( HEXAHEDRON, 0 ));

  CPPUNIT_ASSERT (num_faces == TopologyInfo::sides( HEXAHEDRON ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::corners( HEXAHEDRON ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::edges( HEXAHEDRON ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::faces( HEXAHEDRON ));

  const unsigned* side;
  unsigned i, count;
  const unsigned vert_per_face = 4;
  unsigned edges[num_edges][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
                                   { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
                                   { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } };
  unsigned faces[num_faces][vert_per_face] = { { 0, 1, 5, 4 }, 
                                               { 1, 2, 6, 5 },
                                               { 2, 3, 7, 6 },
                                               { 3, 0, 4, 7 },
                                               { 3, 2, 1, 0 },
                                               { 4, 5, 6, 7 } };

  for (i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( HEXAHEDRON, i, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( HEXAHEDRON, 1, i, count, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  } 

  for (i = 0; i < num_faces; ++i)
  {
    side = TopologyInfo::face_vertices( HEXAHEDRON, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == vert_per_face);
    CPPUNIT_ASSERT( compare_face( vert_per_face, side, faces[i] ) );

    side = TopologyInfo::side_vertices( HEXAHEDRON, 2, i, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == vert_per_face);
    CPPUNIT_ASSERT( compare_face( vert_per_face, side, faces[i] ) );
  } 

  side = TopologyInfo::side_vertices( HEXAHEDRON, 3, 0, count, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(count == num_verts);
  for (i = 0; i < num_verts; ++i)
    CPPUNIT_ASSERT(side[i] == i);
}

void TopologyInfoTest:: hex8()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 8;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void TopologyInfoTest:: hex9()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 9;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void TopologyInfoTest:: hex14()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 14;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void  TopologyInfoTest::hex15()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 15;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void TopologyInfoTest:: hex20()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 20;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void  TopologyInfoTest::hex21()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 21;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void  TopologyInfoTest::hex26()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 26;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void  TopologyInfoTest::hex27()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 27;

  TopologyInfo::higher_order( HEXAHEDRON, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT( face);
  CPPUNIT_ASSERT( vol);

  test_vol_elem( HEXAHEDRON, num_nodes, 8, 12, 6 );
}

void TopologyInfoTest::pyramid()
{
  MsqPrintError err(cout);

  const unsigned num_verts = 5;
  const unsigned num_edges = 8;
  const unsigned num_faces = 5;
  CPPUNIT_ASSERT (3 == TopologyInfo::dimension( PYRAMID ));
  CPPUNIT_ASSERT (1 == TopologyInfo::adjacent( PYRAMID, 3 ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::adjacent( PYRAMID, 2 ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::adjacent( PYRAMID, 1 ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::adjacent( PYRAMID, 0 ));

  CPPUNIT_ASSERT (num_faces == TopologyInfo::sides( PYRAMID ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::corners( PYRAMID ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::edges( PYRAMID ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::faces( PYRAMID ));

  const unsigned* side;
  unsigned i, count;
  const unsigned num_tri_faces = 4;
  const unsigned num_quad_faces = 1;
  const bool tri_before_quad = true;
  CPPUNIT_ASSERT( num_tri_faces + num_quad_faces == num_faces );
  unsigned edges[num_edges][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
                                   { 0, 4 }, { 1, 4 }, { 2, 4 }, { 3, 4 } };
  unsigned tris[num_tri_faces][3] = { { 0, 1, 4 }, { 1, 2, 4 }, 
                                       { 2, 3, 4 }, { 3, 0, 4 } };
  unsigned quads[num_quad_faces][4] = { { 3, 2, 1, 0 } };

  for (i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( PYRAMID, i, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( PYRAMID, 1, i, count, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  } 

  const unsigned tri_off = num_quad_faces * !tri_before_quad;
  for (i = 0; i < num_tri_faces; ++i)
  {
    side = TopologyInfo::face_vertices( PYRAMID, i+tri_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 3);
    CPPUNIT_ASSERT( compare_face( 3, side, tris[i] ) );

    side = TopologyInfo::side_vertices( PYRAMID, 2, i+tri_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 3);
    CPPUNIT_ASSERT( compare_face( 3, side, tris[i] ) );
  } 

  const unsigned quad_off = num_tri_faces * tri_before_quad;
  for (i = 0; i < num_quad_faces; ++i)
  {
    side = TopologyInfo::face_vertices( PYRAMID, i+quad_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 4);
    CPPUNIT_ASSERT( compare_face( 4, side, quads[i] ) );

    side = TopologyInfo::side_vertices( PYRAMID, 2, i+quad_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 4);
    CPPUNIT_ASSERT( compare_face( 4, side, quads[i] ) );
  } 

  side = TopologyInfo::side_vertices( PYRAMID, 3, 0, count, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(count == num_verts);
  for (i = 0; i < num_verts; ++i)
    CPPUNIT_ASSERT(side[i] == i);
}

void TopologyInfoTest::pyramid5()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 5;

  TopologyInfo::higher_order( PYRAMID, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( PYRAMID, num_nodes, 5, 8, 5 );
}

void TopologyInfoTest::pyramid13()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 13;

  TopologyInfo::higher_order( PYRAMID, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( PYRAMID, num_nodes, 5, 8, 5 );
}

void TopologyInfoTest::wedge()
{
  MsqPrintError err(cout);

  const unsigned num_verts = 6;
  const unsigned num_edges = 9;
  const unsigned num_faces = 5;
  CPPUNIT_ASSERT (3 == TopologyInfo::dimension( PRISM ));
  CPPUNIT_ASSERT (1 == TopologyInfo::adjacent( PRISM, 3 ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::adjacent( PRISM, 2 ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::adjacent( PRISM, 1 ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::adjacent( PRISM, 0 ));

  CPPUNIT_ASSERT (num_faces == TopologyInfo::sides( PRISM ));
  CPPUNIT_ASSERT (num_verts == TopologyInfo::corners( PRISM ));
  CPPUNIT_ASSERT (num_edges == TopologyInfo::edges( PRISM ));
  CPPUNIT_ASSERT (num_faces == TopologyInfo::faces( PRISM ));

  const unsigned* side;
  unsigned i, count;
  const unsigned num_tri_faces = 2;
  const unsigned num_quad_faces = 3;
  const bool tri_before_quad = false;
  CPPUNIT_ASSERT( num_tri_faces + num_quad_faces == num_faces );
  unsigned edges[num_edges][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, 
                                   { 0, 3 }, { 1, 4 }, { 2, 5 },
                                   { 3, 4 }, { 4, 5 }, { 5, 3 } };
  unsigned tris[num_tri_faces][3] = { { 2, 1, 0 }, { 3, 4, 5 } };
  unsigned quads[num_quad_faces][4] = { { 0, 1, 4, 3 },
                                        { 1, 2, 5, 4 },
                                        { 2, 0, 3, 5 } };

  for (i = 0; i < num_edges; ++i)
  {
    side = TopologyInfo::edge_vertices( PRISM, i, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );

    side = TopologyInfo::side_vertices( PRISM, 1, i, count, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(2 == count);
    CPPUNIT_ASSERT( compare_edge( side, edges[i] ) );
  } 

  const unsigned tri_off = num_quad_faces * !tri_before_quad;
  for (i = 0; i < num_tri_faces; ++i)
  {
    side = TopologyInfo::face_vertices( PRISM, i+tri_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 3);
    CPPUNIT_ASSERT( compare_face( 3, side, tris[i] ) );

    side = TopologyInfo::side_vertices( PRISM, 2, i+tri_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 3);
    CPPUNIT_ASSERT( compare_face( 3, side, tris[i] ) );
  } 

  const unsigned quad_off = num_tri_faces * tri_before_quad;
  for (i = 0; i < num_quad_faces; ++i)
  {
    side = TopologyInfo::face_vertices( PRISM, i+quad_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 4);
    CPPUNIT_ASSERT( compare_face( 4, side, quads[i] ) );

    side = TopologyInfo::side_vertices( PRISM, 2, i+quad_off, count, err ); 
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(count == 4);
    CPPUNIT_ASSERT( compare_face( 4, side, quads[i] ) );
  } 

  side = TopologyInfo::side_vertices( PRISM, 3, 0, count, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(count == num_verts);
  for (i = 0; i < num_verts; ++i)
    CPPUNIT_ASSERT(side[i] == i);
}

void TopologyInfoTest::wedge6()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 6;

  TopologyInfo::higher_order( PRISM, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(!edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( PRISM, num_nodes, 6, 9, 5 );
}

void TopologyInfoTest::wedge15()
{
  MsqPrintError err(cout);
  bool edge, face, vol;
  const int num_nodes = 15;

  TopologyInfo::higher_order( PRISM, num_nodes, edge, face, vol , err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( edge);
  CPPUNIT_ASSERT(!face);
  CPPUNIT_ASSERT(!vol);

  test_vol_elem( PRISM, num_nodes, 6, 9, 5 );
}

void TopologyInfoTest::test_poly( EntityTopology topo )
{
  CPPUNIT_ASSERT( TopologyInfo::adjacent(topo, 1) == 0 );
  CPPUNIT_ASSERT( TopologyInfo::sides(topo) == 0 );
  CPPUNIT_ASSERT( TopologyInfo::corners(topo) == 0 );
  CPPUNIT_ASSERT( TopologyInfo::edges(topo) == 0 );
  CPPUNIT_ASSERT( TopologyInfo::faces(topo) == 0 );
}


void TopologyInfoTest::polygon()
{
  CPPUNIT_ASSERT( TopologyInfo::dimension( POLYGON ) == 2 );
  test_poly( POLYGON );
}

void TopologyInfoTest::polyhedron()
{
  CPPUNIT_ASSERT( TopologyInfo::dimension( POLYHEDRON ) == 3 );
  test_poly( POLYHEDRON );
}


void TopologyInfoTest::bad_type()
{
  Mesquite::MsqError err;
  EntityTopology bad_types[] = { (EntityTopology)0,
                                 (EntityTopology)1,
                                 MIXED,
                                 (EntityTopology)(MIXED + 1)
                               };

  for (unsigned i = 0; i < (sizeof(bad_types)/sizeof(EntityTopology)); ++i)
  {
    CPPUNIT_ASSERT( TopologyInfo::dimension(bad_types[i]) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::adjacent(bad_types[i], 1) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::sides(bad_types[i]) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::corners(bad_types[i]) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::edges(bad_types[i]) == 0 );
    CPPUNIT_ASSERT( TopologyInfo::faces(bad_types[i]) == 0 );

    bool a,b, c;
    TopologyInfo::higher_order( bad_types[i], 20, a, b, c, err );
    CPPUNIT_ASSERT(err);

    const unsigned* ptr;
    ptr = TopologyInfo::edge_vertices( bad_types[i], 0, err );
    CPPUNIT_ASSERT(err);

    unsigned count;
    ptr = TopologyInfo::face_vertices( bad_types[i], 0, count, err );
    CPPUNIT_ASSERT(err);

    for (unsigned j = 0; j < 4; ++j)
    {
      ptr = TopologyInfo::side_vertices( bad_types[i], j, 0, count, err );
      CPPUNIT_ASSERT(err);
    }

    unsigned dim, side;
    for (unsigned idx = 0; idx < 20; ++idx)
    {
      TopologyInfo::side_number( bad_types[i], idx, idx/2, dim, side, err );
      CPPUNIT_ASSERT(err);
    }
  }
}


void TopologyInfoTest::tri_adj_vert()
{
  unsigned data[][5] = { { 2, 1, 2, 0, 0 },
                         { 2, 2, 0, 0, 0 },
                         { 2, 0, 1, 0, 0 } };
  test_adj( Mesquite::TRIANGLE, data );
}

void TopologyInfoTest::quad_adj_vert()
{
  unsigned data[][5] = { { 2, 1, 3, 0, 0 },
                         { 2, 2, 0, 0, 0 },
                         { 2, 3, 1, 0, 0 },
                         { 2, 0, 2, 0, 0 } };
  test_adj( Mesquite::QUADRILATERAL, data );
}

void TopologyInfoTest::tet_adj_vert()
{
  unsigned data[][5] = { { 3, 1, 2, 3, 0 },
                         { 3, 2, 0, 3, 0 },
                         { 3, 0, 1, 3, 0 },
                         { 3, 2, 1, 0, 0 } };
  test_adj( Mesquite::TETRAHEDRON, data );
}

void TopologyInfoTest::hex_adj_vert()
{
  unsigned data[][5] = { { 3, 1, 3, 4, 0 },
                         { 3, 2, 0, 5, 0 },
                         { 3, 3, 1, 6, 0 },
                         { 3, 0, 2, 7, 0 },
                         { 3, 7, 5, 0, 0 },
                         { 3, 4, 6, 1, 0 },
                         { 3, 5, 7, 2, 0 },
                         { 3, 6, 4, 3, 0 } };
  test_adj( Mesquite::HEXAHEDRON, data );
}

void TopologyInfoTest::pyr_adj_vert()
{
  unsigned data[][5] = { { 3, 1, 3, 4, 0 },
                         { 3, 2, 0, 4, 0 },
                         { 3, 3, 1, 4, 0 },
                         { 3, 0, 2, 4, 0 },
                         { 4, 3, 2, 1, 0 } };
  test_adj( Mesquite::PYRAMID, data );
}

void TopologyInfoTest::wdg_adj_vert()
{
  unsigned data[][5] = { { 3, 1, 2, 3, 0 },
                         { 3, 2, 0, 4, 0 },
                         { 3, 0, 1, 5, 0 },
                         { 3, 5, 4, 0, 0 },
                         { 3, 3, 5, 1, 0 },
                         { 3, 4, 3, 2, 0 } };
  test_adj( Mesquite::PRISM, data );
}

void TopologyInfoTest::test_adj( Mesquite::EntityTopology type,
                                 const unsigned data[][5] )
{
    // Get num vertices from type
  unsigned n = TopologyInfo::corners( type );
  CPPUNIT_ASSERT( n > 0 );
  
  // The first index into data is the vertex.
  // Each column of "data", indexed by vertex, is a 
  // vector containing the number of adjacent vertices
  // followed by the list if adjacent vertex indices.
  
    // for each vertex
  for (unsigned i = 0; i < n; ++i)
  {
      // Get the data corresponding to this vertex of the element
    unsigned const * corner = data[i];
    unsigned expected_count = corner[0];
      // Query TopologyInfo for the same data
    unsigned actual_count;
    unsigned const* actual_adj = TopologyInfo::adjacent_vertices( type, i, actual_count );
      // Check result is valid and counts match 
    CPPUNIT_ASSERT( actual_adj != NULL );
    CPPUNIT_ASSERT( expected_count == actual_count );
    
      // For 3-D elements, returned vertices are expected to be oriented 
      // such that  a face bounded by the vertices in the counter-clockwise 
      // order will have a normal pointing away from the input vertex.
      // So the vertices must be in a certain order, but may begin with
      // any of the adjacent vertices.
      
      // Find the location in the result list at which the first
      // vertex in the expected list occurs.
    unsigned j;
    for (j = 0; j < actual_count; ++j)
      if (corner[1] == actual_adj[j])
        break;
      // Asssert that the first expected vertex was somewhere in 
      // the result list.
    CPPUNIT_ASSERT( j < actual_count );
      // Compare the remaining vertices, enforcing the order.
    for (unsigned k = 1; k < actual_count; ++k)
      CPPUNIT_ASSERT( corner[k+1] == actual_adj[(k+j)%actual_count] );
  }
}

void TopologyInfoTest::test_rev_adj( Mesquite::EntityTopology type )
{
    // Get num vertices from type
  unsigned n = TopologyInfo::corners( type );
  CPPUNIT_ASSERT( n > 0 );
  
    // for each vertex
  for (unsigned i = 0; i < n; ++i)
  {
      // get adjacent vertex list
    unsigned num_adj_idx;
    unsigned const* adj_idx = TopologyInfo::adjacent_vertices( type, i, num_adj_idx );
    CPPUNIT_ASSERT( adj_idx != NULL );
    CPPUNIT_ASSERT( num_adj_idx > 1 );
    
      // make sure adjacent vertex list is unique
    std::vector<unsigned> adj_idx_copy( num_adj_idx );
    std::copy( adj_idx, adj_idx+num_adj_idx, adj_idx_copy.begin() );
    std::sort( adj_idx_copy.begin(), adj_idx_copy.end() );
    std::vector<unsigned>::iterator iter;
    iter = std::unique( adj_idx_copy.begin(), adj_idx_copy.end() );
    CPPUNIT_ASSERT( iter == adj_idx_copy.end() );
    
      // Get reverse mapping indices
    unsigned num_rev_idx;
    unsigned const* rev_idx
     = TopologyInfo::reverse_vertex_adjacency_offsets( type, i, num_rev_idx );
    CPPUNIT_ASSERT( rev_idx != NULL );
    CPPUNIT_ASSERT( num_rev_idx == num_adj_idx );
    
      // for each adjacent vertex, test reverse mapping 
    for (unsigned j = 0; j < num_adj_idx; ++j)
    {
      unsigned num_adj_adj_idx;
      unsigned const* adj_adj_idx 
        = TopologyInfo::adjacent_vertices( type, adj_idx[j], num_adj_adj_idx );
      
      CPPUNIT_ASSERT( rev_idx[j] < num_adj_adj_idx );
      CPPUNIT_ASSERT( adj_adj_idx[rev_idx[j]] == i );
    }
  }
}

      
struct ho_result { unsigned dim; unsigned num; unsigned idx; };   

const ho_result HEX27[] = {
    { 1,  0,  8 },
    { 1,  1,  9 },
    { 1,  2, 10 },
    { 1,  3, 11 },
    { 1,  4, 12 },
    { 1,  5, 13 },
    { 1,  6, 14 },
    { 1,  7, 15 },
    { 1,  8, 16 },
    { 1,  9, 17 },
    { 1, 10, 18 },
    { 1, 11, 19 },
    { 2,  0, 20 },
    { 2,  1, 21 },
    { 2,  2, 22 },
    { 2,  3, 23 },
    { 2,  4, 24 },
    { 2,  5, 25 },
    { 3,  0, 26 } };

const ho_result* const HEX20 = HEX27;
const ho_result* const HEX26 = HEX27;

const ho_result HEX15[] = {
  { 2, 0,  8 },
  { 2, 1,  9 },
  { 2, 2, 10 },
  { 2, 3, 11 },
  { 2, 4, 12 },
  { 2, 5, 13 },
  { 3, 0, 14 } };

const ho_result* const HEX14 = HEX15;

const ho_result HEX9[] = { { 3, 0, 8 } };

void TopologyInfoTest::higher_order_from_side()
{
  MsqPrintError err(std::cerr);
    // HEX-27
  for (unsigned i = 0; i < 19; ++i) {
    unsigned result = TopologyInfo::higher_order_from_side( 
                         HEXAHEDRON, 
                         27, 
                         HEX27[i].dim,
                         HEX27[i].num,
                         err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX27[i].idx, result );
  }
    // HEX-26
  for (unsigned i = 0; i < 18; ++i) {
    unsigned result = TopologyInfo::higher_order_from_side( 
                         HEXAHEDRON, 
                         26, 
                         HEX26[i].dim,
                         HEX26[i].num,
                         err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX26[i].idx, result );
  }
    // HEX-20
  for (unsigned i = 0; i < 12; ++i) {
    unsigned result = TopologyInfo::higher_order_from_side( 
                         HEXAHEDRON, 
                         20, 
                         HEX20[i].dim,
                         HEX20[i].num,
                         err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX20[i].idx, result );
  }
    // HEX-15
  for (unsigned i = 0; i < 7; ++i) {
    unsigned result = TopologyInfo::higher_order_from_side( 
                         HEXAHEDRON, 
                         15, 
                         HEX15[i].dim,
                         HEX15[i].num,
                         err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX15[i].idx, result );
  }
    // HEX-14
  for (unsigned i = 0; i < 6; ++i) {
    unsigned result = TopologyInfo::higher_order_from_side( 
                         HEXAHEDRON, 
                         14, 
                         HEX14[i].dim,
                         HEX14[i].num,
                         err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX14[i].idx, result );
  }
    // HEX-9
  for (unsigned i = 0; i < 1; ++i) {
    unsigned result = TopologyInfo::higher_order_from_side( 
                         HEXAHEDRON, 
                         9, 
                         HEX9[i].dim,
                         HEX9[i].num,
                         err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX9[i].idx, result );
  }
}

void TopologyInfoTest::side_from_higher_order()
{
  MsqPrintError err(std::cerr);
  unsigned dim, num;
    // HEX-27
  for (unsigned i = 0; i < 19; ++i) {
    TopologyInfo::side_from_higher_order( HEXAHEDRON, 27, HEX27[i].idx, dim, num, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX27[i].dim, dim );
    CPPUNIT_ASSERT_EQUAL( HEX27[i].num, num );
  }
    // HEX-26
  for (unsigned i = 0; i < 18; ++i) {
    TopologyInfo::side_from_higher_order( HEXAHEDRON, 26, HEX26[i].idx, dim, num, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX26[i].dim, dim );
    CPPUNIT_ASSERT_EQUAL( HEX26[i].num, num );
  }
    // HEX-20
  for (unsigned i = 0; i < 12; ++i) {
    TopologyInfo::side_from_higher_order( HEXAHEDRON, 20, HEX20[i].idx, dim, num, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX20[i].dim, dim );
    CPPUNIT_ASSERT_EQUAL( HEX20[i].num, num );
  }
    // HEX-15
  for (unsigned i = 0; i < 7; ++i) {
    TopologyInfo::side_from_higher_order( HEXAHEDRON, 15, HEX15[i].idx, dim, num, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX15[i].dim, dim );
    CPPUNIT_ASSERT_EQUAL( HEX15[i].num, num );
  }
    // HEX-14
  for (unsigned i = 0; i < 6; ++i) {
    TopologyInfo::side_from_higher_order( HEXAHEDRON, 14, HEX14[i].idx, dim, num, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX14[i].dim, dim );
    CPPUNIT_ASSERT_EQUAL( HEX14[i].num, num );
  }
    // HEX-9
  for (unsigned i = 0; i < 1; ++i) {
    TopologyInfo::side_from_higher_order( HEXAHEDRON, 9, HEX9[i].idx, dim, num, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( HEX9[i].dim, dim );
    CPPUNIT_ASSERT_EQUAL( HEX9[i].num, num );
  }
}




void TopologyInfoTest::find_edge()
{
  MsqPrintError err(std::cerr);
  bool reversed;
  
  EntityTopology types[] = { TRIANGLE, QUADRILATERAL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID };
  const unsigned num_types = sizeof(types)/sizeof(types[0]);
  unsigned idx;
  for (unsigned t = 0; t < num_types; ++t) {
    const EntityTopology type = types[t];
    
    for (unsigned e = 0; e < TopologyInfo::edges(type); ++e) {
      
      const unsigned *const v = TopologyInfo::edge_vertices( type, e, err );
      CPPUNIT_ASSERT(!err);
      idx = TopologyInfo::find_edge( type, v, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( e, idx );
      CPPUNIT_ASSERT( !reversed );
      
      const unsigned switched[2] = { v[1], v[0] };
      idx = TopologyInfo::find_edge( type, switched, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( e, idx );
      CPPUNIT_ASSERT( reversed );
    }
  }
}

static void shift( const unsigned* in, unsigned* out, unsigned size, unsigned offset )
{
  for (unsigned i = 0; i < size; ++i)
    out[i] = in[(i+offset)%size];
}

static void reverse( const unsigned* in, unsigned* out, unsigned size, unsigned offset )
{
  for (unsigned i = 0; i < size; ++i)
    out[i] = in[(offset + size - i - 1)%size];
}

void TopologyInfoTest::find_face()
{
  MsqPrintError err(std::cerr);
  bool reversed;
  unsigned switched[4], idx;
  
  EntityTopology types[] = { TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID };
  const unsigned num_types = sizeof(types)/sizeof(types[0]);
  for (unsigned t = 0; t < num_types; ++t) {
    const EntityTopology type = types[t];
    
    for (unsigned f = 0; f < TopologyInfo::faces(type); ++f) {
      
      unsigned n;
      const unsigned *const v = TopologyInfo::face_vertices( type, f, n, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(n == 3 || n == 4);
      
      idx = TopologyInfo::find_face( type, v, n, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, idx );
      CPPUNIT_ASSERT( !reversed );
      
      // advance by 1 and try again
      shift( v, switched, n, 1 );
      idx = TopologyInfo::find_face( type, switched, n, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, idx );
      CPPUNIT_ASSERT( !reversed );
      
      // advance by 2 and try again
      shift( v, switched, n, 2 );
      idx = TopologyInfo::find_face( type, switched, n, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, idx );
      CPPUNIT_ASSERT( !reversed );
      
      // advance by 3 and try again
      shift( v, switched, n, 3 );
      idx = TopologyInfo::find_face( type, switched, n, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, idx );
      CPPUNIT_ASSERT( !reversed );
      
      // reverse and try again
      reverse( v, switched, n, 0 );
      idx = TopologyInfo::find_face( type, switched, n, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, idx );
      CPPUNIT_ASSERT( reversed );
      
      // reverse, advance by 1 and try again
      reverse( v, switched, n, 1 );
      idx = TopologyInfo::find_face( type, switched, n, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, idx );
      CPPUNIT_ASSERT( reversed );
      
      // reverse, advance by 2 and try again
      reverse( v, switched, n, 2 );
      idx = TopologyInfo::find_face( type, switched, n, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, idx );
      CPPUNIT_ASSERT( reversed );
      
      // reverse, advance by 3 and try again
      reverse( v, switched, n, 3 );
      idx = TopologyInfo::find_face( type, switched, n, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, idx );
      CPPUNIT_ASSERT( reversed );
    }
  }
}


void TopologyInfoTest::find_side()
{
  MsqPrintError err(std::cerr);
  bool reversed;
  unsigned i, dim, switched[4];
  
  EntityTopology types[] = { TRIANGLE, QUADRILATERAL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID };
  const unsigned num_types = sizeof(types)/sizeof(types[0]);
  for (unsigned t = 0; t < num_types; ++t) {
    const EntityTopology type = types[t];
    
    for (unsigned e = 0; e < TopologyInfo::edges(type); ++e) {
      
      const unsigned *const v = TopologyInfo::edge_vertices( type, e, err );
      CPPUNIT_ASSERT(!err);
      TopologyInfo::find_side( type, v, 2, dim, i, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( e, i );
      CPPUNIT_ASSERT_EQUAL( 1u, dim );
      CPPUNIT_ASSERT( !reversed );
      
      switched[0] = v[1]; switched[1] = v[0];
      TopologyInfo::find_side( type, switched, 2, dim, i, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( e, i );
      CPPUNIT_ASSERT_EQUAL( 1u, dim );
      CPPUNIT_ASSERT( reversed );
    }
  }

  for (unsigned t = 2; t < num_types; ++t) {
    const EntityTopology type = types[t];
    
    for (unsigned f = 0; f < TopologyInfo::faces(type); ++f) {
      
      unsigned n;
      const unsigned *const v = TopologyInfo::face_vertices( type, f, n, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(n == 3 || n == 4);
      
      TopologyInfo::find_side( type, v, n, dim, i, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, i );
      CPPUNIT_ASSERT_EQUAL( 2u, dim );
      CPPUNIT_ASSERT( !reversed );

      // reverse and try again
      reverse( v, switched, n, 0 );
      TopologyInfo::find_side( type, switched, n, dim, i, reversed, err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( f, i );
      CPPUNIT_ASSERT_EQUAL( 2u, dim );
      CPPUNIT_ASSERT( reversed );
    }
  }
}

void TopologyInfoTest::compare_sides()
{
  // define two hexes sharing a face:
//    6-------7-------8
//   /|      /|      /|
//  0-------1-------2 |
//  | |     | |     | |
//  | 9-----|-10----|-11
//  |/      |/      |/
//  3-------4-------5
  
  const size_t hex1[] = { 3, 4, 10, 9, 0, 1, 7, 6 };
  const size_t hex2[] = { 4, 5, 11, 10, 1, 2, 8, 7 };
    // shared edges: { hex1_edge, hex2_edge }
  const unsigned edges[][2] = { { 1, 3 },
                                { 5, 4 },
                                { 6, 7 },
                                { 9, 11} };
  const unsigned num_edges = sizeof(edges)/sizeof(edges[0]);
    // shared faces: { hex1_face, hex2_face }
  const unsigned faces[][2] = { { 1, 3 } };
  const unsigned num_faces = sizeof(faces)/sizeof(faces[0]);
  
  MsqPrintError err(std::cerr);
  
    // try every possible edge combination
  for (unsigned e1 = 0; e1 < 12; ++e1) {
    unsigned match;
    for (match = 0; match < num_edges; ++match)
      if (edges[match][0] == e1)
        break;
     
    for (unsigned e2 = 0; e2 < 12; ++e2) {
      const bool expected = (match < num_edges) && (edges[match][1] == e2);
    
      const bool result = TopologyInfo::compare_sides( hex1, HEXAHEDRON, e1,
                                                       hex2, HEXAHEDRON, e2,
                                                       1, err );
      
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( expected, result );
    }
  }
  
    // try every possible face combination
  for (unsigned f1 = 0; f1 < 6; ++f1) {
    unsigned match;
    for (match = 0; match < num_faces; ++match)
      if (faces[match][0] == f1)
        break;
     
    for (unsigned f2 = 0; f2 < 6; ++f2) {
      const bool expected = (match < num_faces) && (faces[match][1] == f2);
    
      const bool result = TopologyInfo::compare_sides( hex1, HEXAHEDRON, f1,
                                                       hex2, HEXAHEDRON, f2,
                                                       2, err );
      
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_EQUAL( expected, result );
    }
  }
}
