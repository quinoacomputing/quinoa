/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <stk_percept/fixtures/PyramidFixture.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_io/IossBridge.hpp>

//----------------------------------------------------------------------

namespace stk_classic{
  namespace percept {

    typedef shards::Hexahedron<8>          Hex8;
    typedef shards::Wedge<6>               Wedge6;
    typedef shards::Tetrahedron<4>         Tet4;
    typedef shards::Pyramid<5>             Pyramid5;

    typedef shards::ShellQuadrilateral<4>  ShellQuad4;
    typedef shards::ShellTriangle<3>       ShellTriangle3;

    typedef shards::Quadrilateral<4>  Quad4;
    typedef shards::Triangle<3>       Triangle3;

    PyramidFixture::PyramidFixture( stk_classic::ParallelMachine comm, bool doCommit, bool do_sidesets ) :
      m_spatial_dimension(3)
      , m_metaData(m_spatial_dimension, stk_classic::mesh::fem::entity_rank_names(m_spatial_dimension) )
      , m_bulkData( stk_classic::mesh::fem::FEMMetaData::get_meta_data(m_metaData) , comm )
      , m_block_pyramid(    m_metaData.declare_part< Pyramid5 >( "block_4" ))
      , m_sideset_quad(0), m_sideset_quad_subset(0)
      , m_sideset_tri(0), m_sideset_tri_subset(0)
      , m_elem_rank( m_metaData.element_rank() )
      , m_coordinates_field( m_metaData.declare_field< VectorFieldType >( "coordinates" ))
      , m_centroid_field(    m_metaData.declare_field< VectorFieldType >( "centroid" ))
      , m_temperature_field( m_metaData.declare_field< ScalarFieldType >( "temperature" ))
      , m_volume_field( m_metaData.declare_field< ScalarFieldType >( "volume" ))
      , m_element_node_coordinates_field( m_metaData.declare_field< ElementNodePointerFieldType >( "elem_node_coord" ))
    {
      // Define where fields exist on the mesh:
      stk_classic::mesh::Part & universal = m_metaData.universal_part();

      if (do_sidesets)
        {
          m_sideset_quad_subset = &m_metaData.declare_part(std::string("surface_pyramid5_quad2d2_1"), m_metaData.face_rank());
          m_sideset_quad =        &m_metaData.declare_part(std::string("surface_1"), m_metaData.face_rank());
          stk_classic::mesh::fem::set_cell_topology< Quad4 >(*m_sideset_quad_subset);
          stk_classic::io::put_io_part_attribute(*m_sideset_quad_subset);
          stk_classic::io::put_io_part_attribute(*m_sideset_quad);
          m_metaData.declare_part_subset(*m_sideset_quad, *m_sideset_quad_subset);

          m_sideset_tri_subset = &m_metaData.declare_part(std::string("surface_pyramid5_tri2d2_1"), m_metaData.face_rank());
          m_sideset_tri =        &m_metaData.declare_part(std::string("surface_2"), m_metaData.face_rank());
          stk_classic::mesh::fem::set_cell_topology< Triangle3 >(*m_sideset_tri_subset);
          stk_classic::io::put_io_part_attribute(*m_sideset_tri_subset);
          stk_classic::io::put_io_part_attribute(*m_sideset_tri);
          m_metaData.declare_part_subset(*m_sideset_tri, *m_sideset_tri_subset);
        }
      put_field( m_coordinates_field , stk_classic::mesh::fem::FEMMetaData::NODE_RANK , universal );
      put_field( m_centroid_field , m_elem_rank , universal );
      put_field( m_temperature_field, stk_classic::mesh::fem::FEMMetaData::NODE_RANK, universal );

      put_field( m_volume_field, m_elem_rank, m_block_pyramid );

      // Define the field-relation such that the values of the
      // 'element_node_coordinates_field' are pointers to the
      // element's nodal 'coordinates_field'.
      // I.e., let:
      //   double *const* elem_node_coord =
      //     field_data( m_element_node_coordinates_field , element );
      // then
      //     elem_node_coord[n][0..2] is the coordinates of element node 'n'
      //     that are attached to that node.

      m_metaData.declare_field_relation(
                                        m_element_node_coordinates_field ,
                                        stk_classic::mesh::fem::get_element_node_stencil(3) ,
                                        m_coordinates_field
                                        );

      // Define element node coordinate field for all element parts
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_pyramid, Pyramid5::node_count );
      stk_classic::io::put_io_part_attribute(  m_block_pyramid );

      if (doCommit)
        m_metaData.commit();
    }

    PyramidFixture::~PyramidFixture()
    { }

    //------------------------------------------------------------------------------
    // Use case specific mesh generation data:

    enum { SpatialDim = 3 };
    enum { node_count = 7 };
    enum { number_pyramid = 2 };
    enum { number_quad = 2 };
    enum { number_tri = 6 };

    namespace {

      // Hard coded node coordinate data for all the nodes in the entire mesh
      static const double node_coord_data[ node_count ][ SpatialDim ] = {

        { 0 , 0 , -1 } , { 1 , 0 , -1 } ,  { 2 , 0 , -1 } , 

        { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , 

        { 1 , 1 , -2 } 
      };

      // Hard coded pyramid node ids for all the pyramid nodes in the entire mesh
      static const stk_classic::mesh::EntityId pyramid_node_ids[number_pyramid][ Pyramid5::node_count ] = {
        { 1 , 4 , 5 , 2 , 7 } ,
        { 2 , 5 , 6 , 3 , 7 } };

      // Hard coded quad node ids for all the quad nodes in the entire mesh
      static const stk_classic::mesh::EntityId quad_node_ids[number_quad][ Quad4::node_count ] = {
        { 1 , 4 , 5 , 2  } ,
        { 2 , 5 , 6 , 3  } };

      // Hard coded tri node ids for all the tri nodes in the entire mesh
      static const stk_classic::mesh::EntityId tri_node_ids[number_tri][ Triangle3::node_count ] = {
        { 1 , 7 , 2 },
        { 4, 7, 1 },
        { 5, 7, 4 },
        { 6, 7, 5 },
        { 3, 7, 6 },
        { 2, 7, 3 }
        };

    }

    //------------------------------------------------------------------------------

    void PyramidFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          stk_classic::mesh::EntityId curr_elem_id = 1;

          // For each element topology declare elements

          stk_classic::mesh::Entity *pyramids[2];
          for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++curr_elem_id ) {
            pyramids[i] = &stk_classic::mesh::fem::declare_element( m_bulkData, m_block_pyramid, curr_elem_id, pyramid_node_ids[i] );
          }

          if (m_sideset_quad)
            {
              for ( unsigned i = 0 ; i < number_quad ; ++i , ++curr_elem_id ) {
                stk_classic::mesh::fem::declare_element_side( m_bulkData, 
                                                      curr_elem_id, //side_id,
                                                      *pyramids[i], // element,
                                                      4,            //j_side, // local_side_ord,
                                                      m_sideset_quad_subset);
              }
            }

          if (m_sideset_tri)
            {
              unsigned j_side=0;
              for ( unsigned i = 0 ; i < 3 ; ++i , ++curr_elem_id ) {
                if (i == 2) ++j_side;
                stk_classic::mesh::fem::declare_element_side( m_bulkData, 
                                                      curr_elem_id, //side_id,
                                                      *pyramids[0], // element,
                                                      j_side,            //j_side, // local_side_ord,
                                                      m_sideset_tri_subset);
                ++j_side;
              }
              j_side=1;
              for ( unsigned i = 0 ; i < 3 ; ++i , ++curr_elem_id ) {
                stk_classic::mesh::fem::declare_element_side( m_bulkData, 
                                                      curr_elem_id, //side_id,
                                                      *pyramids[1], // element,
                                                      j_side,            //j_side, // local_side_ord,
                                                      m_sideset_tri_subset);
                ++j_side;
              }
            }

          // For all nodes assign nodal coordinates
          for ( unsigned i = 0 ; i < node_count ; ++i ) {
            stk_classic::mesh::Entity * const node = m_bulkData.get_entity( stk_classic::mesh::fem::FEMMetaData::NODE_RANK , i + 1 );
            double * const coord = field_data( m_coordinates_field , *node );
            coord[0] = node_coord_data[i][0] ;
            coord[1] = node_coord_data[i][1] ;
            coord[2] = node_coord_data[i][2] ;
          }

        }
      m_bulkData.modification_end();

    }

    // Verify mesh for 6 different parts
    bool verifyMesh( const PyramidFixture & mesh )
    {
      bool result = true;

      const stk_classic::mesh::BulkData & bulkData = mesh.m_bulkData ;
      //const VectorFieldType & node_coord = mesh.m_coordinates_field ;
      //const ElementNodePointerFieldType & elem_node_coord  =  mesh.m_element_node_coordinates_field ;

      std::vector<stk_classic::mesh::Bucket *> element_buckets = bulkData.buckets( mesh.m_elem_rank );

      // Create a pair containing Part and matching node_count

      typedef std::pair<stk_classic::mesh::Part*, unsigned> PartNodeCountPair;
      std::vector<PartNodeCountPair> part_and_node_counts;
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_pyramid, Pyramid5::node_count));

      // Verify that entities in each part are set up correctly.
      // Use a PartVector iterator for parts_to_check and call
      // verify_elem_node_coord_by_part in UseCase_Common.cpp for
      // each part in turn.
#if 0
      for( std::vector<PartNodeCountPair>::const_iterator  i = part_and_node_counts.begin() ; i != part_and_node_counts.end() ; ++i )
        {
          result = result &&
            verify_elem_node_coord_by_part(
                                           *(i->first),
                                           element_buckets,
                                           elem_node_coord,
                                           node_coord,
                                           i->second
                                           );
        }
#endif

      // Check that all the nodes were allocated.
      for ( unsigned i = 0 ; i < node_count ; ++i ) {
        stk_classic::mesh::Entity * const node = bulkData.get_entity( stk_classic::mesh::fem::FEMMetaData::NODE_RANK , i + 1 );
        if ( node == NULL ) {
          std::cerr << "Error!  Invalid null pointer for node returned from "
                    << "bulkData.get_entity( stk_classic::mesh::fem::FEMMetaData::NODE_RANK, " << i+1 << " ) " << std::endl;
          result = false;
        }
      }

      return result;
    }

  } //namespace percept
} //namespace stk_classic
