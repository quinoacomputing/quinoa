/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>

#include <util/TPI.h>
#include <util/ParallelComm.hpp>

#include <mesh/FieldTraits.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/Comm.hpp>

#include <element/Hexahedron_Topologies.hpp>

using namespace phdmesh ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef Field<double,Cartesian> CoordinateField ;

void test_simple_mesh( ParallelMachine pm , std::istream & )
{
  typedef Hexahedron<> Hex ;

  static const char method[] = "test_simple_mesh" ;

  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );

  //--------------------------------------------------------------------
  // Define a mesh mesh_meta_data: the parts and fields.

  MetaData S ;

  // Get some of the predefined parts for later use...
  Part * const owns_part = & S.locally_owned_part();
  Part * const univ_part = & S.universal_part();

  // Declare a part for the element block and side set,
  // these are automatically a subset of the universal part.

  Part * const elem_part = & S.declare_part( std::string("element_block") );
  Part * const face_part = & S.declare_part( std::string("side_set") );

  // Nodal coordinate field dimensioned to 3 everywhere in the mesh

  CoordinateField & node_coordinates =
    S.declare_field<CoordinateField>( std::string("coordinates") );

  S.put_field( node_coordinates , Node , *univ_part , 3 );
  
  // Done defining the mesh_meta_data, commit it.

  S.commit();

  //--------------------------------------------------------------------
  // Create mesh bulk data conformal to the mesh_meta_data.

  const unsigned kernel_capacity = 100 ;

  BulkData M( S , pm , kernel_capacity );

  // Define a trivial mesh, stack of hex elements
  // with one hex element per processor ordered by
  // processor rank.
  // Attach the +X face, use ExodusII element-node ordering
  // and element-face orientation.
  // Node and element identifiers must not be zero,
  // an identifier of zero is reserved for 'undefined'.

  const entity_key_type node_key[8] = {
    // Base of this processor's hex
    entity_key( Node , p_rank * 4 + 1 ) ,
    entity_key( Node , p_rank * 4 + 2 ) ,
    entity_key( Node , p_rank * 4 + 3 ) ,
    entity_key( Node , p_rank * 4 + 4 ) ,

    // Top of this processor's hex
    entity_key( Node , p_rank * 4 + 5 ) ,
    entity_key( Node , p_rank * 4 + 6 ) ,
    entity_key( Node , p_rank * 4 + 7 ) ,
    entity_key( Node , p_rank * 4 + 8 ) };

  const entity_key_type elem_key = entity_key( Element , p_rank + 1 );
  const entity_key_type face_key = entity_key( Face ,    p_rank + 1 );

  // Part membership for the elements and nodes:
  // 'owns_part'  Assume this processor owns everything it declares,
  //              will resolve parallel sharing later.

  std::vector<Part*> add_parts ;
  add_parts.push_back( owns_part );
  add_parts.push_back( elem_part );

  // Declare node and element entities:

  Entity & elem = M.declare_entity( elem_key , add_parts );

  for ( unsigned i = 0 ; i < 8 ; ++i ) {
    Entity & node = M.declare_entity( node_key[i] , add_parts );

    // Declare element <-> node relations
    // These are required to have unique identifiers
    // by providing the 'method' argument.
    // If non-unique then an exception is thrown that includes
    // the text contained in the 'method' string.

    M.declare_relation( elem , node , i );
  }

  // Declare the face entity:

  add_parts.push_back( face_part );
  
  Entity & face = M.declare_entity( face_key , add_parts );

  // Declare element <-> face relation

  M.declare_relation( elem , face , 0 );

  // Declare face <-> node relations

  StaticAssert< Hex::side<0,0>::node == 0 >::ok();
  StaticAssert< Hex::side<0,1>::node == 1 >::ok();
  StaticAssert< Hex::side<0,2>::node == 5 >::ok();
  StaticAssert< Hex::side<0,3>::node == 4 >::ok();

  const unsigned elem_face_node[4] = {
    Hex::side<0,0>::node ,
    Hex::side<0,1>::node ,
    Hex::side<0,2>::node ,
    Hex::side<0,3>::node };

  PairIterRelation elem_node = elem.relations( Node );

  for ( unsigned i = 0 ; i < 4 ; ++i ) {
    Entity & node = * elem_node[ elem_face_node[i] ].entity();
    M.declare_relation( face , node , i );

    // Update the nodes on the face to also be members of the face part.

    M.change_entity_parts( node , add_parts );
  }

  // Set node coordinates:

  Entity & node_0 = * elem_node[0].entity();
  Entity & node_1 = * elem_node[1].entity();
  Entity & node_2 = * elem_node[2].entity();
  Entity & node_3 = * elem_node[3].entity();
  Entity & node_4 = * elem_node[4].entity();
  Entity & node_5 = * elem_node[5].entity();
  Entity & node_6 = * elem_node[6].entity();
  Entity & node_7 = * elem_node[7].entity();

  field_data_valid( node_coordinates , node_0 , method );
  field_data_valid( node_coordinates , node_1 , method );
  field_data_valid( node_coordinates , node_2 , method );
  field_data_valid( node_coordinates , node_3 , method );
  field_data_valid( node_coordinates , node_4 , method );
  field_data_valid( node_coordinates , node_5 , method );
  field_data_valid( node_coordinates , node_6 , method );
  field_data_valid( node_coordinates , node_7 , method );

  double * const node_0_coord = field_data( node_coordinates , node_0 );
  double * const node_1_coord = field_data( node_coordinates , node_1 );
  double * const node_2_coord = field_data( node_coordinates , node_2 );
  double * const node_3_coord = field_data( node_coordinates , node_3 );
  double * const node_4_coord = field_data( node_coordinates , node_4 );
  double * const node_5_coord = field_data( node_coordinates , node_5 );
  double * const node_6_coord = field_data( node_coordinates , node_6 );
  double * const node_7_coord = field_data( node_coordinates , node_7 );

  node_0_coord[0] = 0 ; node_0_coord[1] = 0 ; node_0_coord[2] = p_rank ;
  node_1_coord[0] = 1 ; node_1_coord[1] = 0 ; node_1_coord[2] = p_rank ;
  node_2_coord[0] = 1 ; node_2_coord[1] = 1 ; node_2_coord[2] = p_rank ;
  node_3_coord[0] = 0 ; node_3_coord[1] = 1 ; node_3_coord[2] = p_rank ;
  node_4_coord[0] = 0 ; node_4_coord[1] = 0 ; node_4_coord[2] = p_rank + 1 ;
  node_5_coord[0] = 1 ; node_5_coord[1] = 0 ; node_5_coord[2] = p_rank + 1 ;
  node_6_coord[0] = 1 ; node_6_coord[1] = 1 ; node_6_coord[2] = p_rank + 1 ;
  node_7_coord[0] = 0 ; node_7_coord[1] = 1 ; node_7_coord[2] = p_rank + 1 ;

  // Determine proper parallel sharing and ownership
  comm_mesh_discover_sharing( M );

  // Generate the parallel ghosting 'aura'
  comm_mesh_regenerate_aura( M );

  // Verify that the parallel sharing and aura were generated properly
  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << method
                << " FAILED: is not parallel consistent"
                << std::endl ;
    }
    return ;
  }

  // Get the global counts and identifier stats
  {
    entity_id_type counts[ EntityTypeEnd ];
    entity_id_type max_id[ EntityTypeEnd ];

    comm_mesh_stats( M , counts , max_id );

    if ( p_rank == 0 ) {
      std::cout << method
                << " Stats for { node , edge , face , element , other }"
                << std::endl ;
      std::cout << "  Global Counts = {" 
                << " " << counts[0]
                << " " << counts[1]
                << " " << counts[2]
                << " " << counts[3]
                << " " << counts[4]
                << " " << counts[5]
                << " }" << std::endl ;
      std::cout << "  Global MaxId  = {" 
                << " " << max_id[0]
                << " " << max_id[1]
                << " " << max_id[2]
                << " " << max_id[3]
                << " " << max_id[4]
                << " " << max_id[5]
                << " }" << std::endl ;
    }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      parallel_machine_barrier( pm );
      if ( p_rank == p ) {
        count_entities( M , S.locally_used_part() , counts );

        std::cout << "  P" << p_rank << " Uses  Counts = {" 
                  << " " << counts[0]
                  << " " << counts[1]
                  << " " << counts[2]
                  << " " << counts[3]
                  << " " << counts[4]
                  << " " << counts[5]
                  << " }" << std::endl ;

        count_entities( M , S.locally_owned_part() , counts );

        std::cout << "  P" << p_rank << " Owns  Counts = {" 
                  << " " << counts[0]
                  << " " << counts[1]
                  << " " << counts[2]
                  << " " << counts[3]
                  << " " << counts[4]
                  << " " << counts[5]
                  << " }" << std::endl ;

        std::cout.flush();
      }
      parallel_machine_barrier( pm );
    }
  }
  parallel_machine_barrier( pm );

  //------------------------------

  if ( p_rank == 0 ) {
    std::cout << method << " successful" << std::endl ;
    std::cout.flush();
  }
}


