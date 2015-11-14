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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 12-Nov-02 at 18:05:56
//  LAST-MOD:  9-Jun-04 at 14:43:39 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file PatchDataInstances.hpp

This header file contains some functions to instantiates particular PatchData Objects.
Those objects can be used in unit tests.
Patches must be allocated and dealocated by the caller. 

\author Thomas Leurent
\author Michael Brewer

*/
// DESCRIP-END.
//

#ifndef PatchDataInstances_hpp
#define PatchDataInstances_hpp

#include "Mesquite_MsqVertex.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_IdealElements.hpp"
#include "Mesquite_TopologyInfo.hpp"

#include <math.h>
#include <iostream>

#include "cppunit/extensions/HelperMacros.h"

namespace MESQUITE_NS
{
  //! must be called in sync with create_...._patch_with_domain
  inline void destroy_patch_with_domain(PatchData &pd)
  {
    MeshDomain* domain = pd.get_domain();
    pd.set_domain( 0 );
    delete domain;
    
    //Mesh* mesh = pd.get_mesh();
    //pd.set_mesh( 0 );
    //delete mesh;
  }
  
  inline void move_vertex( PatchData& pd, 
                           const Vector3D& position,
                           const Vector3D& delta,
                           MsqError& err )
  {
    const MsqVertex* array = pd.get_vertex_array( err ); 
    if (err) return;
    
    int idx = 0, cnt = 0;
    for (size_t i = 0; i < pd.num_nodes(); ++i)
      if ((array[i] - position).length_squared() < DBL_EPSILON) 
        { idx = i; ++cnt; }
    
    CPPUNIT_ASSERT_EQUAL( cnt, 1 );
    
    pd.move_vertex( delta, idx, err );
  }
        
  
  /*! creates a patch containing one ideal hexahedra
  */
   inline void create_one_hex_patch(PatchData &one_hex_patch, MsqError &err)
   {
     double coords[] = { 1.0, 1.0, 1.0,
                         2.0, 1.0, 1.0,
                         2.0, 2.0, 1.0,
                         1.0, 2.0, 1.0,
                         1.0, 1.0, 2.0,
                         2.0, 1.0, 2.0,
                         2.0, 2.0, 2.0,
                         1.0, 2.0, 2.0 };
     
     size_t indices[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
     
     one_hex_patch.fill( 8, coords, 1, HEXAHEDRON, indices, 0, err );
   }
   

   //! creates a Patch containing an ideal tetrahedra
   inline void create_one_tet_patch(PatchData &one_tet_patch, MsqError &err)
   {
     double coords[] = { 1.0, 1.0, 1.0,
                         2.0, 1.0, 1.0,
                         1.5, 1+sqrt(3.0)/2.0, 1.0,
                         1.5, 1+sqrt(3.0)/6.0, 1+sqrt(2.0)/sqrt(3.0) };

     size_t indices[4] = { 0, 1, 2, 3 };

     one_tet_patch.fill( 4, coords, 1, TETRAHEDRON, indices, 0, err );
   }
   
   //! create patch containing one ideal pyramid
   inline void create_one_pyr_patch( PatchData& one_pyr_patch, MsqError& err )
   {
     /* Equilateral triangles
     double coords[] = { 1, -1, 0,
                         1,  1, 0,
                        -1,  1, 0,
                        -1, -1, 0,
                         0,  0, sqrt(2) };
     */
     /* Unit height */
     double coords[] = { 1, -1, 0,
                         1,  1, 0,
                        -1,  1, 0,
                        -1, -1, 0,
                         0,  0, 2 };
                         
     size_t indices[5] = { 0, 1, 2, 3, 4 };
     
     one_pyr_patch.fill( 5, coords, 1, PYRAMID, indices, 0, err );
   } 
   
   //! create patch containing one ideal wedge
   inline void create_one_wdg_patch( PatchData& one_wdg_patch, MsqError& err )
   {
     double hgt = 0.5 * MSQ_SQRT_THREE;
     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         0.5, hgt, 0.0,
                         0.0, 0.0, 1.0,
                         1.0, 0.0, 1.0,
                         0.5, hgt, 1.0 };
                         
     size_t indices[6] = { 0, 1, 2, 3, 4, 5 };
     
     one_wdg_patch.fill( 6, coords, 1, PRISM, indices, 0, err );
   } 

      //! creates a Patch containing an ideal tetrahedra, inverted
   inline void create_one_inverted_tet_patch(PatchData &one_tet_patch,
                                             MsqError &err)
   {
     double coords[] = { 1, 1, 1,
                         2, 1, 1,
                         1.5, 1+sqrt(3.0)/2.0, 1,
                         1.5, 1+sqrt(3.0)/6.0, 1-sqrt(2.0)/sqrt(3.0), };
                         
     size_t indices[4] = { 0, 1, 2, 3 };

     one_tet_patch.fill( 4, coords, 1, TETRAHEDRON, indices, 0, err );
   }
   
   //! creates a Patch containing an ideal quadrilateral
   inline void create_one_quad_patch(PatchData &one_qua_patch, MsqError &err)
   {
     double coords[] = { 1, 1, 1,
                         2, 1, 1,
                         2, 2, 1,
                         1, 2 , 1 };
     
     size_t indices[4] = { 0, 1, 2, 3 };
     
     one_qua_patch.fill( 4, coords, 1, QUADRILATERAL, indices, 0, err );
   }
   
   
     /*! \fn create_one_tri_patch(PatchData &one_tri_patch, MsqError &err)
            2
           / \      creates a Patch containing an ideal triangle
          /   \
         0-----1
         This Patch also has the normal information. 
     */
   inline void create_one_tri_patch(PatchData &one_tri_patch, MsqError &err)
   {
       /* ************** Creates normal info ******************* */
     Vector3D pnt(0,0,0);
     Vector3D s_norm(0,0,3);
     one_tri_patch.set_domain( new PlanarDomain( s_norm, pnt ) );

       /* *********************FILL tri************************* */
     double coords[] = { 1, 1, 1,
                         2, 1, 1,
                         1.5, 1+sqrt(3.0)/2.0, 1 };
     
     size_t indices[3] = { 0, 1, 2 };
     one_tri_patch.fill( 3, coords, 1, TRIANGLE, indices, 0, err );
   }
     
  
   /*! \fn create_two_tri_patch(PatchData &one_tri_patch, MsqError &err)
            2
           / \      creates a Patch containing two ideal triangles
          / 0 \
         0-----1
          \ 1 /
           \ /
            3
         This Patch also has the normal information. 
   */
   inline void create_two_tri_patch(PatchData &pd, MsqError &err)
   {
       /* ************** Creates normal info ******************* */
     Vector3D pnt(0,0,1);
     Vector3D s_norm(0,0,3);
     pd.set_domain( new PlanarDomain(s_norm, pnt) );

       // **********************FILL tri*************************

     double coords[] = { 1, 1, 1,
                         2, 1, 1,
                         1.5, 1+sqrt(3.0)/2.0, 1,
                         1.5, 1-sqrt(3.0)/2.0, 1 };

     size_t indices[] = { 0, 1, 2, 0, 3, 1 };
     
     pd.fill( 4, coords, 2, TRIANGLE, indices, 0, err );
  }
   

   /*! \fn create_four_quads_patch(PatchData &four_quads, MsqError &err)
     our 2D set up: 4 quads, center vertex outcentered by (0,-0.5)
      7____6____5
      |    |    |
      | 2  |  3 |
      8-_  |  _-4       vertex 1 is at (0,0)
      |  -_0_-  |       vertex 5 is at (2,2)
      | 0  |  1 |
      1----2----3
   */
   inline void create_four_quads_patch(PatchData &four_quads, MsqError &err) 
   {
     double coords[] = { 1, .5, 0,
                         0, 0, 0,
                         1, 0, 0,
                         2, 0, 0,
                         2, 1, 0,
                         2, 2, 0,
                         1, 2, 0,
                         0, 2, 0,
                         0, 1, 0 };

     size_t indices[] = { 1, 2, 0, 8, 
                          2, 3, 4, 0,
                          8, 0, 6, 7,
                          0, 4, 5, 6 };
                        
     
     four_quads.fill( 9, coords, 4, QUADRILATERAL, indices, 0, err );
   }
   

   /*! \fn create_six_quads_patch(PatchData &four_quads, MsqError &err)
     our 2D set up: 6 quads, 1 center vertex outcentered by (0,-0.5), the other centered
      7____6____5___11
      |    |    |    |
      | 2  |  3 | 5  |
      8-_  |  _-4---10       vertex 1 is at (0,0)
      |  -_0_-  |    |       vertex 11 is at (3,2)
      | 0  |  1 | 4  |
      1----2----3----9

      use destroy_patch_with_domain() in sync.
   */
   inline void create_six_quads_patch_with_domain(PatchData &pd, MsqError &err) 
   {
     // associates domain
     Vector3D pnt(0,0,0);
     Vector3D s_norm(0,0,3);
     pd.set_domain( new PlanarDomain(s_norm, pnt) );

     double coords[] = { 1,.5, 0,
                         0, 0, 0,
                         1, 0, 0,
                         2, 0, 0,
                         2, 1, 0,
                         2, 2, 0,
                         1, 2, 0,
                         0, 2, 0,
                         0, 1, 0,
                         3, 0, 0,
                         3, 1, 0,
                         3, 2, 0 };

     size_t indices[] = { 1,  2,  0, 8,
                          2,  3,  4, 0,
                          8,  0,  6, 7,
                          0,  4,  5, 6,
                          3,  9, 10, 4,
                          4, 10, 11, 5 };
     
     bool fixed[] = { false, true, true, true, 
                      false, true, true, true, 
                      true, true, true, true };
     
     pd.fill( 12, coords, 6, QUADRILATERAL, indices, fixed, err );
   }
   

   /*! \fn create_six_quads_patch_inverted_with_domain(PatchData &four_quads, MsqError &err)
     our 2D set up: 6 quads, 1 center vertex outcentered by (0,-0.5), the other centered
      7____6____5___11
      |    |    |    |
      | 2  |  3 | 5  |
      8    |    4---10       vertex 1 is at (0,0)
      |\       /|    |       vertex 11 is at (3,2)
      |    |    | 4  |
      1----2----3----9
         \  /
          0      
      use destroy_patch_with_domain() in sync.
   */
   inline void create_six_quads_patch_inverted_with_domain(PatchData &pd, MsqError &err) 
   {
     create_six_quads_patch_with_domain(pd,err); MSQ_CHKERR(err);

     Vector3D displacement(0,-1.5,0);
     
     pd.move_vertex( displacement, 0, err );
   }
   

   /*! \fn create_twelve_hex_patch(PatchData &pd, MsqError &err)
     3D set up: 12 quads, one center vertex outcentered by (0,-0.5),
     the other centered. Vertex 1 is at (0,0,-1). Vertex 35 is at (3,2,1).
     
      7____6____5___11     19___18____17__23     31___30___29___35
      |    |    |    |      |    |    |    |      |    |    |    |
      | 2  |  3 | 5  |      |    |    |    |      | 8  |  9 | 11 |
      8----0----4---10     20-_  |  _16---22     32---24---28---34       
      |    |    |    |      |  -12_-  |    |      |    |    |    |       
      | 0  |  1 | 4  |      |    |    |    |      | 6  |  7 | 10 |
      1----2----3----9     13---14---15---21     25---26---27---33
   */
   inline void create_twelve_hex_patch(PatchData &pd, MsqError &err) 
   {
     double coords[] = { 1, 1, -1,
                         0, 0, -1,
                         1, 0, -1,
                         2, 0, -1,
                         2, 1, -1,
                         2, 2, -1,
                         1, 2, -1,
                         0, 2, -1,
                         0, 1, -1,
                         3, 0, -1,
                         3, 1, -1,
                         3, 2, -1,

                         1,.5, 0,
                         0, 0, 0,
                         1, 0, 0,
                         2, 0, 0,
                         2, 1, 0,
                         2, 2, 0,
                         1, 2, 0,
                         0, 2, 0,
                         0, 1, 0,
                         3, 0, 0,
                         3, 1, 0,
                         3, 2, 0,

                         1, 1, 1,
                         0, 0, 1,
                         1, 0, 1,
                         2, 0, 1,
                         2, 1, 1,
                         2, 2, 1,
                         1, 2, 1,
                         0, 2, 1,
                         0, 1, 1,
                         3, 0, 1,
                         3, 1, 1,
                         3, 2, 1 };
     
     size_t connectivity[] = { 1, 2, 0, 8, 13, 14, 12, 20, // 0
                               2, 3, 4, 0, 14, 15, 16, 12, // 1
                               8, 0, 6, 7, 20, 12, 18, 19, // 2
                               0, 4, 5, 6, 12, 16, 17, 18, // 3
                               3, 9, 10, 4, 15, 21, 22, 16, // 4
                               4, 10, 11, 5, 16, 22, 23, 17, // 5
                               13, 14, 12, 20, 25, 26, 24, 32, // 6
                               14, 15, 16, 12, 26, 27, 28, 24, // 7
                               20, 12, 18, 19, 32, 24, 30, 31, // 8
                               12, 16, 17, 18, 24, 28, 29, 30, // 9
                               15, 21, 22, 16, 27, 33, 34, 28, // 10
                               16, 22, 23, 17, 28, 34, 35, 29 }; // 11
    
     bool fixed[] = { true,  true,  true,  true,  true,  true,  true,  true,
                      true,  true,  true,  true,  false, true,  true,  true,
                      false, true,  true,  true,  true,  true,  true,  true,
                      true,  true,  true,  true,  true,  true,  true,  true,
                      true,  true,  true,  true };


     pd.fill( 36, coords, 12, HEXAHEDRON, connectivity, fixed, err );
   }
   
   inline void create_twelve_hex_patch_inverted(PatchData &pd, MsqError &err)
   {
     create_twelve_hex_patch(pd,err); MSQ_CHKERR(err); 
     move_vertex( pd, Vector3D(2,1,0), Vector3D(0,0,1.5), err ); MSQ_CHKERR(err);
   }
     
     
   /* Patch used in several quality metric tests.
      Our triangular patch is made of two tris.  tri_1 is a perfect
      equilateral (the ideal for most metrics).  tri_2 is an arbitrary
      triangle.
      Memory allocated in this function must be deallocated with
      destroy_patch_with_domain().
   */
   inline void create_qm_two_tri_patch_with_domain(PatchData &triPatch, MsqError &err) 
   {
     Vector3D pnt(0,0,0);
     Vector3D s_norm(0,0,3);
     triPatch.set_domain( new PlanarDomain(s_norm, pnt) );

     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         0.5, sqrt(3.0)/2.0, 0.0,
                         2.0, -4.0, 2.0 };
     

     const size_t conn[] = { 0, 1, 2, 0, 3, 1 };
     
     triPatch.fill( 4, coords, 2, TRIANGLE, conn, 0, err );
   }
   
     /* Patch used in several quality metric tests.
       Our quad patch is made of two quads.  quad_1 is a perfect
       square (the ideal for most metrics).  quad_2 is an arbitrary
       quad.
       Memory allocated in this function must be deallocated with
       destroy_patch_with_domain().
     */
   inline void create_qm_two_quad_patch_with_domain(PatchData &quadPatch, MsqError &err)
   {
     Vector3D pnt(0,0,0);
     Vector3D s_norm(0,0,3);
     quadPatch.set_domain( new PlanarDomain(s_norm, pnt) );

     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         1.0, 1.0, 0.0,
                         0.0, 1.0, 0.0,
                         2.0, -1.0, .5,
                         1.5, 1.0, 1.0 };
     
     const size_t conn[] = { 0, 1, 2, 3, 1, 4, 5, 2 };
     
     quadPatch.fill( 6, coords, 2, QUADRILATERAL, conn, 0, err );
   }
  
     /* Patch used in several quality metric tests.
        Our tet patch is made of two tets.  tet_1 is a perfect
        equilateral (the ideal for most metrics).  tet_2 is an arbitrary
        tet.
     */
   inline void create_qm_two_tet_patch(PatchData &tetPatch, MsqError &err)
   {
     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         0.5, sqrt(3.0)/2.0, 0.0,
                         0.5, sqrt(3.0)/6.0, sqrt(2.0)/sqrt(3.0),
                         2.0, 3.0, -.5 };
     

     const size_t conn[] = { 0, 1, 2, 3, 1, 4, 2, 3 };
     
     tetPatch.fill( 5, coords, 2, TETRAHEDRON, conn, 0, err );
   }
  
     /* Patch used in several quality metric tests.
        Our pyr patch is made of two pyramids.  The first is a perfect
        pyramid (the ideal for most metrics).  The second is an arbitrary
        pyramid.
     */
   inline void create_qm_two_pyr_patch(PatchData &pyrPatch, MsqError &err)
   {
     /* Equilateral triangles
     double coords[] = { 1, -1, 0,
                         1,  1, 0,
                        -1,  1, 0,
                        -1, -1, 0,
                         0,  0, sqrt(2) };
     */
     /* Unit height */
     double coords[] = { 
                         /* Equilateral triangles */
                    /*   1, -1, 0,
                         1,  1, 0,
                        -1,  1, 0,
                        -1, -1, 0,
                         0,  0, sqrt(2)  */                   
                         /* Unit height */
                         1, -1, 0,
                         1,  1, 0,
                        -1,  1, 0,
                        -1, -1, 0,
                         0,  0, 2,
                         /* Apex for a squashed pyramid */
                         0,  0, -1
                         };
     

     const size_t conn[] = { 0, 1, 2, 3, 4, 
                             3, 2, 1, 0, 5 };
     
     pyrPatch.fill( 6, coords, 2, PYRAMID, conn, 0, err );
   }
  
     /* Patch used in several quality metric tests.
        Our prism patch is made of two prisms.  The first is a perfect
        prism (the ideal for most metrics).  The second is an arbitrary
        wedge.
     */
   inline void create_qm_two_wdg_patch(PatchData &wdgPatch, MsqError &err)
   {
     double hgt = 0.5 * MSQ_SQRT_THREE;
     double coords[] = {  // ideal prism vertices
                         0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         0.5, hgt, 0.0,
                         0.0, 0.0, 1.0,
                         1.0, 0.0, 1.0,
                         0.5, hgt, 1.0,
                          // top vertices for stretched wedge
                         0.5,-3.0, 0.0,
                         0.5,-4.0, 1.0 };

     const size_t conn[] = { 0, 1, 2, 3, 4, 5,
                             1, 0, 6, 4, 3, 7 };
     
     wdgPatch.fill( 8, coords, 2, PRISM, conn, 0, err );
   }

   /* Patch used in seveal quality metric tests.
      Our hex patch is made of two hexes.  hex_1 is a perfect
      unit cube (the ideal for most metrics).  hex_2 is an arbitrary
      hex.
   */
   inline void create_qm_two_hex_patch(PatchData &hexPatch, MsqError &err)
   {
     double coords[] = { 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         1.0, 1.0, 0.0,
                         0.0, 1.0, 0.0,
                         0.0, 0.0, 1.0,
                         1.0, 0.0, 1.0,
                         1.0, 1.0, 1.0,
                         0.0, 1.0, 1.0,
                         2.0, 0.0, 0.0,
                         2.0, 1.0, 0.0,
                         2.0,-1.0, 1.0,
                         3.0, 2.0, 1.0 };
     
     const size_t conn[] = { 0, 1, 2, 3, 4, 5, 6, 7,
                             1, 8, 9, 2, 5, 10, 11, 6 };
                             
     hexPatch.fill( 12, coords, 2, HEXAHEDRON, conn, 0, err ); 
   }
   
   // Create patch containing one ideal element, optionally higher-order.
   // For 2D elements, will attach appropriate planar domain.
   inline void create_ideal_element_patch( PatchData& pd,
                                           EntityTopology type, 
                                           size_t num_nodes,
                                           MsqError& err )
   {
      static PlanarDomain zplane(PlanarDomain::XY);
      static Settings settings;
      settings.set_slaved_ho_node_mode( Settings::SLAVE_NONE );
      pd.attach_settings( &settings );
      
   
        // build list of vertex coordinates
      const Vector3D* corners = unit_edge_element( type );
      std::vector<Vector3D> coords( corners, corners+TopologyInfo::corners(type) );
      bool mids[4] = {false};
      TopologyInfo::higher_order( type, num_nodes, mids[1], mids[2], mids[3], err );
      MSQ_ERRRTN(err);
      std::vector<size_t> conn(coords.size());
      for (unsigned i = 0; i < coords.size(); ++i)
        conn[i] = i;
  
      for (unsigned dim = 1; dim <= TopologyInfo::dimension(type); ++dim) {
        if (!mids[dim])
          continue;
        
        int num_side;
        if (dim == TopologyInfo::dimension(type))
          num_side = 1;
        else
          num_side = TopologyInfo::adjacent( type, dim );
        
        for (int s = 0; s < num_side; ++s) {
          unsigned idx = TopologyInfo::higher_order_from_side( type, num_nodes, dim, s, err );
          MSQ_ERRRTN(err);
          conn.push_back(idx);
        
          unsigned n;
          const unsigned* side = TopologyInfo::side_vertices( type, dim, s, n, err );
          MSQ_ERRRTN(err);
          Vector3D avg = coords[side[0]];
          for (unsigned v = 1; v < n; ++v)
            avg += coords[side[v]];
          avg *= 1.0/n;
          coords.push_back(avg);
        }
      }
      
      bool* fixed = new bool[coords.size()];
      std::fill( fixed, fixed+coords.size(), false );
      pd.fill( coords.size(), coords[0].to_array(), 1, &type, 
               &num_nodes, &conn[0], fixed, err );
      delete [] fixed;
      MSQ_ERRRTN(err);
      
      if (TopologyInfo::dimension(type) == 2)
        pd.set_domain( &zplane );
   }
   
} // namespace

#endif // PatchDataInstances_hpp
