#include <iostream>
#include <set>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "UseCase_ChangeOwner.hpp"

namespace {

typedef shards::Quadrilateral<4> Quad4;

const stk_classic::mesh::EntityRank NODE_RANK = stk_classic::mesh::fem::FEMMetaData::NODE_RANK;

const unsigned spatial_dimension = 2;

}

Grid2D_Fixture::Grid2D_Fixture( stk_classic::ParallelMachine comm )
  : m_fem_meta_data( spatial_dimension ),
    m_bulk_data( stk_classic::mesh::fem::FEMMetaData::get_meta_data(m_fem_meta_data) , comm , 100 ),
    m_quad_part( stk_classic::mesh::fem::declare_part<Quad4>( m_fem_meta_data, "quad")),
    m_coord_field( m_fem_meta_data.declare_field< VectorField >( "coordinates" ) ),
    m_elem_rank( m_fem_meta_data.element_rank() ),
    m_node_rank( m_fem_meta_data.node_rank() )
{
  stk_classic::mesh::put_field( m_coord_field , m_node_rank , m_fem_meta_data.universal_part() );

  m_fem_meta_data.commit();
}

bool Grid2D_Fixture::test_change_owner( unsigned nx , unsigned ny )
{
  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();

  int result = 1 ;

  m_bulk_data.modification_begin();

  // First of all work out the node ids and declare element elem
  if ( p_rank == 0 ) {
    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk_classic::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk_classic::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk_classic::mesh::fem::declare_element( m_bulk_data , m_quad_part , elem , nodes );
      }
    }
  }

  // Only P0 has any nodes or elements
  if ( p_rank == 0 ) {
    result = ! m_bulk_data.buckets( m_node_rank ).empty() &&
      ! m_bulk_data.buckets( m_elem_rank ).empty();
  }
  else {
    result = m_bulk_data.buckets( m_node_rank ).empty() &&
             m_bulk_data.buckets( m_elem_rank ).empty();
  }

  m_bulk_data.modification_end();

  if ( 1 < p_size ) {
    std::vector< stk_classic::mesh::EntityProc > change ;

    if ( p_rank == 0 ) {
      if ( p_size == 3 ) {
        const unsigned nnx = nx + 1 ;
        const unsigned nny = ny + 1 ;
        for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
            stk_classic::mesh::EntityId id = 1 + ix + iy * nnx ;
            unsigned proc = ix < nx/2 ? 1 : 2;
            stk_classic::mesh::EntityProc tmp( m_bulk_data.get_entity( m_node_rank , id ) , proc );
            change.push_back( tmp );
          }
        }
        for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
            stk_classic::mesh::EntityId id = 1 + ix + iy * nx ;
            unsigned proc = ix < nx/2 ? 1 : 2;
            stk_classic::mesh::EntityProc tmp( m_bulk_data.get_entity( m_elem_rank , id ) , proc );
            change.push_back( tmp );
          }
        }
      }
      else
      {
        const unsigned nnx = nx + 1 ;
        const unsigned nny = ny + 1 ;
        for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
            stk_classic::mesh::EntityId id = 1 + ix + iy * nnx ;
            stk_classic::mesh::EntityProc tmp( m_bulk_data.get_entity( m_node_rank , id ) , 1 );
            change.push_back( tmp );
          }
        }
        for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
            stk_classic::mesh::EntityId id = 1 + ix + iy * nx ;
            stk_classic::mesh::EntityProc tmp( m_bulk_data.get_entity( m_elem_rank , id ) , 1 );
            change.push_back( tmp );
          }
        }
      }
    }

    m_bulk_data.modification_begin();
    m_bulk_data.change_entity_owner( change );
    m_bulk_data.modification_end();

    change.clear();

    if ( p_rank == 0 ) {
      const unsigned nnx = nx + 1 ;
      const unsigned nny = ny + 1 ;
      for ( unsigned iy = 0 ; iy < nny / 2 ; ++iy ) {
        for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
          stk_classic::mesh::EntityId id = 1 + ix + iy * nnx ;
          stk_classic::mesh::EntityProc tmp( m_bulk_data.get_entity( m_node_rank , id ) , 1 );
          change.push_back( tmp );
        }
      }
      for ( unsigned iy = 0 ; iy < ny / 2 ; ++iy ) {
        for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
          stk_classic::mesh::EntityId id = 1 + ix + iy * nx ;
          stk_classic::mesh::EntityProc tmp( m_bulk_data.get_entity( m_elem_rank , id ) , 1 );
          change.push_back( tmp );
        }
      }
    }

    m_bulk_data.modification_begin();
    m_bulk_data.change_entity_owner( change );
    m_bulk_data.modification_end();

    if ( p_size == 3 ) {
      change.clear();

      if ( p_rank == 2 ) {
        const unsigned nnx = nx + 1 ;
        const unsigned nny = ny + 1 ;
        for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
          for ( unsigned ix = nx / 2 ; ix < nnx ; ++ix ) {
            stk_classic::mesh::EntityId id = 1 + ix + iy * nnx ;
            unsigned proc = 1;
            stk_classic::mesh::EntityProc tmp( m_bulk_data.get_entity( m_node_rank , id ) , proc );
            change.push_back( tmp );
          }
        }
        for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
          for ( unsigned ix = nx / 2 ; ix < nx ; ++ix ) {
            stk_classic::mesh::EntityId id = 1 + ix + iy * nx ;
            unsigned proc = 1;
            stk_classic::mesh::EntityProc tmp( m_bulk_data.get_entity( m_elem_rank , id ) , proc );
            change.push_back( tmp );
          }
        }
      }

      m_bulk_data.modification_begin();
      m_bulk_data.change_entity_owner( change );
      m_bulk_data.modification_end();
    }

    // Only P1 has any nodes or elements
    if ( p_rank == 1 ) {
      result = ! m_bulk_data.buckets( m_node_rank ).empty() &&
               ! m_bulk_data.buckets( m_elem_rank ).empty();
    }
    else {
      result = m_bulk_data.buckets( m_node_rank ).empty() &&
               m_bulk_data.buckets( m_elem_rank ).empty();
    }
  }

  stk_classic::all_reduce( m_bulk_data.parallel() , stk_classic::ReduceMin<1>( & result ) );

  return result ;
}

//----------------------------------------------------------------------------

bool test_change_owner_with_constraint( stk_classic::ParallelMachine pm )
{
  bool success = true ;

  const unsigned p_rank = stk_classic::parallel_machine_rank( pm );
  const unsigned p_size = stk_classic::parallel_machine_size( pm );

  if ( p_size != 2 ) { return success ; }

  std::vector<std::string> rank_names = stk_classic::mesh::fem::entity_rank_names(spatial_dimension);
  const stk_classic::mesh::EntityRank constraint_rank = rank_names.size();
  rank_names.push_back("Constraint");
  stk_classic::mesh::fem::FEMMetaData fem_meta_data( spatial_dimension, rank_names );
  const stk_classic::mesh::EntityRank element_rank = fem_meta_data.element_rank();

  VectorField * coordinates_field =
    & put_field(
        fem_meta_data.declare_field<VectorField>("coordinates"),
        NODE_RANK,
        fem_meta_data.universal_part() ,
        3
        );

  stk_classic::mesh::Part & owned_part = fem_meta_data.locally_owned_part();
  stk_classic::mesh::Part & quad_part  = stk_classic::mesh::fem::declare_part<Quad4>( fem_meta_data, "quad");

  fem_meta_data.commit();

  stk_classic::mesh::BulkData bulk_data( stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta_data),
                                 pm,
                                 100 );
  bulk_data.modification_begin();

  unsigned nx = 3;
  unsigned ny = 3;

  if ( p_rank==0 )
  {
    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk_classic::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk_classic::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk_classic::mesh::fem::declare_element( bulk_data , quad_part , elem , nodes );
      }
    }

    for ( unsigned iy = 0 ; iy < ny+1 ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx+1 ; ++ix ) {
        stk_classic::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk_classic::mesh::Entity * n = bulk_data.get_entity( NODE_RANK, nid );
        double * const coord = stk_classic::mesh::field_data( *coordinates_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  if ( p_size>1 )
  {
    std::vector<stk_classic::mesh::EntityProc> ep;

    if ( p_rank==0 )
    {
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  3 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  4 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  7 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  8 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 3 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 6 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 1 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();

    bulk_data.modification_begin();

    if ( p_rank==1 )
    {
      // create constraint

      stk_classic::mesh::Entity * n10 = bulk_data.get_entity( NODE_RANK, 10 );
      stk_classic::mesh::Entity * n11 = bulk_data.get_entity( NODE_RANK, 11 );
      stk_classic::mesh::Entity * n12 = bulk_data.get_entity( NODE_RANK, 12 );

      stk_classic::mesh::PartVector add;
      add.push_back( &owned_part );
      const stk_classic::mesh::EntityId c_entity_id = 1;
      stk_classic::mesh::Entity & c = bulk_data.declare_entity( constraint_rank, c_entity_id, add );
      bulk_data.declare_relation( c , *n10 , 0 );
      bulk_data.declare_relation( c , *n11 , 1 );
      bulk_data.declare_relation( c , *n12 , 2 );
    }

    bulk_data.modification_end();

    stk_classic::mesh::Entity * n10 = bulk_data.get_entity( NODE_RANK, 10 );

    if ( p_rank==0 or p_rank==1 )
    {
      ThrowErrorMsgIf( !stk_classic::mesh::in_shared( *n10 ), "NODE[10] not shared" );
    }

    bulk_data.modification_begin();

    if ( p_rank==1 )
    {
      // destroy constraint

      stk_classic::mesh::Entity * c1 = bulk_data.get_entity( constraint_rank, 1 );

      ThrowErrorMsgIf( !bulk_data.destroy_entity( c1 ),
                       "failed to destroy constraint" );
    }

    bulk_data.modification_end();

    if ( p_rank==0 or p_rank==1 )
    {
      ThrowErrorMsgIf( stk_classic::mesh::in_shared( *n10 ), "NODE[10] shared" );
    }
  }

  return success ;
}

//----------------------------------------------------------------------------

bool test_change_owner_2( stk_classic::ParallelMachine pm )
{
  const unsigned p_rank = stk_classic::parallel_machine_rank( pm );
  const unsigned p_size = stk_classic::parallel_machine_size( pm );

  if ( p_size != 3 ) { return true ; }

  stk_classic::mesh::fem::FEMMetaData fem_meta_data( spatial_dimension );
  const stk_classic::mesh::EntityRank element_rank = fem_meta_data.element_rank();

  VectorField * coordinates_field =
    & put_field(
        fem_meta_data.declare_field<VectorField>("coordinates"),
        NODE_RANK,
        fem_meta_data.universal_part() ,
        3
        );

  stk_classic::mesh::Part & quad_part  = stk_classic::mesh::fem::declare_part<Quad4>( fem_meta_data, "quad");

  fem_meta_data.commit();

  stk_classic::mesh::BulkData bulk_data( stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta_data),
                                 pm,
                                 100 );
  bulk_data.modification_begin();

  unsigned nx = 3;
  unsigned ny = 3;

  if ( p_rank==0 )
  {
    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk_classic::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk_classic::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk_classic::mesh::fem::declare_element( bulk_data , quad_part , elem , nodes );
      }
    }

    for ( unsigned iy = 0 ; iy < ny+1 ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx+1 ; ++ix ) {
        stk_classic::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk_classic::mesh::Entity * n = bulk_data.get_entity( NODE_RANK, nid );
        double * const coord = stk_classic::mesh::field_data( *coordinates_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  if ( p_size==3 )
  {
    std::vector<stk_classic::mesh::EntityProc> ep;

    if ( p_rank==0 )
    {
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  9 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 13 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 14 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 7 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 8 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  3 ), 2 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  4 ), 2 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  7 ), 2 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  8 ), 2 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 2 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 2 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 3 ), 2 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 6 ), 2 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();

    ep.clear();

    if ( p_rank==1 )
    {
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  9 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 13 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 14 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 0 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 7 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 8 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 0 ) );
    }

    if ( p_rank==2 )
    {
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  3 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  4 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  7 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  8 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 0 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 3 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 6 ), 0 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();
  }

  return true ;
}

//----------------------------------------------------------------------------

bool test_change_owner_3( stk_classic::ParallelMachine pm )
{
  const int p_rank = stk_classic::parallel_machine_rank( pm );
  const unsigned p_size = stk_classic::parallel_machine_size( pm );

  stk_classic::mesh::fem::FEMMetaData fem_meta_data( spatial_dimension );
  const stk_classic::mesh::EntityRank element_rank = fem_meta_data.element_rank();

  VectorField * coordinates_field =
    & put_field(
        fem_meta_data.declare_field<VectorField>("coordinates"),
        NODE_RANK,
        fem_meta_data.universal_part() ,
        3
        );

  stk_classic::mesh::Part & quad_part  = stk_classic::mesh::fem::declare_part<Quad4>( fem_meta_data, "quad");

  fem_meta_data.commit();

  stk_classic::mesh::BulkData bulk_data( stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta_data),
                                 pm,
                                 100 );
  bulk_data.modification_begin();

  unsigned nx = 3;
  unsigned ny = 3;

  if ( p_rank==0 )
  {
    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk_classic::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk_classic::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk_classic::mesh::fem::declare_element( bulk_data , quad_part , elem , nodes );
      }
    }

    for ( unsigned iy = 0 ; iy < ny+1 ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx+1 ; ++ix ) {
        stk_classic::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk_classic::mesh::Entity * n = bulk_data.get_entity( NODE_RANK, nid );
        double * const coord = stk_classic::mesh::field_data( *coordinates_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  if ( p_size>1 ) {

    std::vector<stk_classic::mesh::EntityProc> ep;

    if ( p_rank==0 )
    {
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  2 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  3 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  4 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  6 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  7 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  8 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 10 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 14 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 2 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 5 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 8 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 3 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 6 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 1 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();

    // output to debug

    ep.clear();

    if ( p_rank==0 )
    {
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  1 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  5 ), 1 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 1 ), 1 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 4 ), 1 ) );
    }
    else if ( p_rank==1 )
    {
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 10 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 0 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 14 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 0 ) );

      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 8 ), 0 ) );
      ep.push_back( stk_classic::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 0 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();
  }

  return true ;
}


