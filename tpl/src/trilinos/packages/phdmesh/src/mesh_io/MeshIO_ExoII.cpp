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

#include <string.h>
#include <stdio.h>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include <strings.h>
#include <util/ParallelReduce.hpp>
#include <util/ParallelComm.hpp>

#include <mesh/MetaData.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/Comm.hpp>

#include <mesh_io/ExoII.hpp>

#include <element/Declarations.hpp>
#include <element/Hexahedron_Topologies.hpp>
#include <element/Tetrahedron_Topologies.hpp>
#include <element/Pyramid_Topologies.hpp>
#include <element/Triangle_Topologies.hpp>
#include <element/Quadrilateral_Topologies.hpp>
#include <element/Wedge_Topologies.hpp>

//----------------------------------------------------------------------

namespace phdmesh {
namespace exodus {

const GlobalLocalIndex & GlobalLocalIndex::tag()
{ static const GlobalLocalIndex self ; return self ; }

const char * GlobalLocalIndex::name() const
{ static const char n[] = "GlobalLocalIndex" ; return n ; }

std::string GlobalLocalIndex::to_string( unsigned size , unsigned index ) const
{
  static const char g[] = "global" ;
  static const char l[] = "local" ;

  if ( 2 != size || size <= index ) {
    std::ostringstream msg ;
    msg << tag().name();
    msg << " ERROR Size = " << size ;
    msg << " Index = " << index ;
    throw std::runtime_error( msg.str() );
  }

  return std::string( index ? g : l );
}

unsigned GlobalLocalIndex::to_index( unsigned size , const std::string & arg )
  const
{
  static const char g[] = "global" ;
  static const char l[] = "local" ;

  const unsigned index = ! strcasecmp( arg.c_str() , g ) ? 0 : (
                         ! strcasecmp( arg.c_str() , l ) ? 1 : 2 );

  if ( 2 == index ) {
    std::ostringstream msg ;
    msg << tag().name();
    msg << " ERROR size = " << size ;
    msg << " label = " << arg ;
    throw std::runtime_error( msg.str() );
  }

  return index ;
}

const ElementAttributes & ElementAttributes::tag()
{ static const ElementAttributes self ; return self ; }

const char * ElementAttributes::name() const
{ static const char n[] = "ElementAttributes" ; return n ; }


}
}


//----------------------------------------------------------------------

#if defined( HAVE_MPI )

#include <mpi.h>

namespace {

void allgather( phdmesh::ParallelMachine c ,
                unsigned * local , unsigned * global , int n )
{
  MPI_Allgather( local , n , MPI_UNSIGNED , global , n , MPI_UNSIGNED , c );
}

}

#else

namespace {

void allgather( phdmesh::ParallelMachine ,
                unsigned * local , unsigned * global , int n )
{
  for ( int i = 0 ; i < n ; ++i ) { global[i] = local[i] ; }
}

}

#endif

//----------------------------------------------------------------------

namespace phdmesh {
namespace exodus {

void verify_node_coordinate_field( const FieldBase & f )
{
  static const char method[] = "phdmesh::exodus::FileSchema::FileSchema" ;

  // Verify the node coordinate field

  const bool bad_type = ! f.type_is<double>();
  const bool bad_rank = 1 != f.rank();

  if ( bad_type || bad_rank ) {
    std::ostringstream msg ;
    msg << method << " FAILED with bad node coordinates field: " << f ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

class FilePart {
public:
  Part               & m_part ;
  const int            m_identifier ;
  const EntityType     m_entity_type ;
  const CellTopology * m_topology ;
  const int            m_number_attr ;

  FilePart( Part      & arg_part ,
            int         arg_id ,
            EntityType  arg_entity_type ,
            const CellTopology * arg_topology ,
            int         arg_number_attributes )
    : m_part( arg_part ),
      m_identifier( arg_id ),
      m_entity_type( arg_entity_type ),
      m_topology( arg_topology ),
      m_number_attr( arg_number_attributes )
    {}
};

struct FilePartLess {
  bool operator()( const FilePart * lhs , int rhs ) const
    { return lhs->m_identifier < rhs ; }
  bool operator()( int lhs, const FilePart * rhs ) const
    { return lhs < rhs->m_identifier; }
  bool operator()( const FilePart * lhs, const FilePart * rhs ) const
    { return lhs->m_identifier < rhs->m_identifier; }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

// Exodus indices: [ global , local ]
const FileSchema::IndexField & exo_index( MetaData & arg_mesh_meta_data )
{
  typedef FileSchema::IndexField IndexField ;

  static const char name[] = "exodusII_local_index" ;

  const Part & upart = arg_mesh_meta_data.universal_part();

  IndexField & f = arg_mesh_meta_data.declare_field< IndexField >( std::string(name) );

  arg_mesh_meta_data.put_field( f , Node ,    upart , 2 );
  arg_mesh_meta_data.put_field( f , Edge ,    upart , 2 );
  arg_mesh_meta_data.put_field( f , Face ,    upart , 2 );
  arg_mesh_meta_data.put_field( f , Element , upart , 2 );

  return f ;
}

}

FileSchema::FileSchema(
  MetaData                      & arg_mesh_meta_data ,
  const FieldBase                   & arg_node_coordinates ,
  const FileSchema::AttributeField  & arg_elem_attributes ,
  const unsigned                      arg_writer_rank )
  : m_schema( arg_mesh_meta_data ),
    m_io_rank( arg_writer_rank ),
    m_dimension( arg_node_coordinates.max_size( Node ) ),
    m_field_node_coord( arg_node_coordinates ),
    m_field_elem_attr(  arg_elem_attributes ),
    m_field_index( exo_index( arg_mesh_meta_data ) )
{
  verify_node_coordinate_field( arg_node_coordinates );
}

FileSchema::~FileSchema() {}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void map_to_exodus( const CellTopology & top ,
                    const unsigned space_dim ,
                    std::string & element_type ,
                    int & number_attributes )
{
  number_attributes = 0 ;

  switch( top.key ) {
  case Node_Traits::key :
    break ;

  case Line<2>::key :
    if ( 1 < space_dim ) {
      element_type.assign( "BEAM_2" );
      number_attributes = space_dim == 2 ? 3 : (
                          space_dim == 3 ? 7 : 0 );
    }
    break ;

  case Line<3>::key :
    if ( 1 < space_dim ) {
      element_type.assign( "BEAM_3" );
      number_attributes = space_dim == 2 ? 3 : (
                          space_dim == 3 ? 7 : 0 );
    }
    break ;

  case ShellLine<2>::key :
    element_type.assign( "SHELL_2" );
    number_attributes = 1 ;
    break ;

  case ShellLine<3>::key :
    element_type.assign( "SHELL_3" );
    number_attributes = 1 ;
    break ;

  case Triangle<3>::key :
    element_type.assign( "TRIANGLE_3" );
    break ;

  case Triangle<6>::key :
    element_type.assign( "TRIANGLE_6" );
    break ;

  case ShellTriangle<3>::key :
    element_type.assign( "TRIANGLE_3" );
    number_attributes = 1 ;
    break ;

  case ShellTriangle<6>::key :
    element_type.assign( "TRIANGLE_6" );
    number_attributes = 1 ;
    break ;

  case Quadrilateral<4>::key :
    element_type.assign( "QUAD_4" );
    break ;

  case Quadrilateral<8>::key :
    element_type.assign( "QUAD_8" );
    break ;

  case Quadrilateral<9>::key :
    element_type.assign( "QUAD_9" );
    break ;

  case ShellQuadrilateral<4>::key :
    element_type.assign( "SHELL_4" );
    number_attributes = 1 ;
    break ;

  case ShellQuadrilateral<8>::key :
    element_type.assign( "SHELL_8" );
    number_attributes = 1 ;
    break ;

  case ShellQuadrilateral<9>::key :
    element_type.assign( "SHELL_9" );
    number_attributes = 1 ;
    break ;

  case Tetrahedron< 4>::key :
    element_type.assign( "TETRA_4" );
    break ;

  case Tetrahedron<10>::key :
    element_type.assign( "TETRA_10" );
    break ;

  case Pyramid< 5>::key :
    element_type.assign( "PYRAMID_5" );
    break ;

  case Pyramid<13>::key :
    element_type.assign( "PYRAMID_13" );
    break ;

  case Pyramid<14>::key :
    element_type.assign( "PYRAMID_14" );
    break ;

  case Wedge< 6>::key :
    element_type.assign( "WEDGE_6" );
    break ;

  case Wedge<15>::key :
    element_type.assign( "WEDGE_15" );
    break ;

  case Wedge<18>::key :
    element_type.assign( "WEDGE_18" );
    break ;

  case Hexahedron< 8>::key :
    element_type.assign( "HEX_8" );
    break ;

  case Hexahedron<20>::key :
    element_type.assign( "HEX_20" );
    break ;

  case Hexahedron<27>::key :
    element_type.assign( "HEX_27" );
    break ;
  }
}

const FilePart * internal_declare_part(
  MetaData & arg_mesh_meta_data ,
  Part         & arg_part ,
  int            arg_id ,
  EntityType     arg_type ,
  const CellTopology * arg_topology ,
  int            arg_space_dim ,
  int            arg_number_attr ,
  const FileSchema::AttributeField  & arg_attributes )
{
  typedef FileSchema::AttributeField AttributeField ;

  static const char method[] = "phdmesh::exodus::FileSchema::declare_part" ;

  std::string element_type ;
  int number_attributes = 0 ;

  if ( arg_topology ) {
    map_to_exodus( * arg_topology , arg_space_dim ,
                   element_type , number_attributes );
  }

  const FilePart * file_part = arg_part.attribute<FilePart>();

  if ( ! arg_number_attr ) { arg_number_attr = number_attributes ; }

  if ( ! file_part ) {

    if ( arg_number_attr ) {
      arg_mesh_meta_data.put_field(
        const_cast<AttributeField &>(arg_attributes),
        Element, arg_part, arg_number_attr );
    }

    file_part = new FilePart( arg_part, arg_id, arg_type,
                              arg_topology, arg_number_attr );

    arg_mesh_meta_data.declare_attribute_with_delete<FilePart>(
      arg_part , file_part );
  }

  if ( file_part == NULL ||
       & file_part->m_part         != & arg_part ||
         file_part->m_identifier   != arg_id ||
         file_part->m_entity_type  != arg_type ||
         file_part->m_topology     != arg_topology ||
         file_part->m_number_attr  != arg_number_attr ) {
    std::ostringstream msg ;
    msg << method << " FAILED redeclaration" ;
    throw std::runtime_error( msg.str() );
  }

  return file_part ;
}

}

void FileSchema::declare_part( Part & arg_part , int arg_id )
{
  std::vector< const FilePart * >::iterator
    i = std::lower_bound( m_parts[Element].begin() ,
                          m_parts[Element].end() ,
                          arg_id , FilePartLess() );

  if ( i != m_parts[ Element ].end() && (*i)->m_identifier == arg_id ) {
    std::ostringstream msg ;
    msg << "phdmesh::exodus::FileSchema::declare_part( " ;
    msg << arg_part.name() ;
    msg << " , " ;
    msg << arg_id ;
    msg << " ) FAILED, already called declare_part( " ;
    msg << (*i)->m_part.name() ;
    msg << " , " ;
    msg << arg_id ;
    msg << " )" ;
    throw std::runtime_error( msg.str() );
  }

  const CellTopology * const top = get_cell_topology( arg_part );

  if ( NULL == top ) {
    std::ostringstream msg ;
    msg << "phdmesh::exodus::FileSchema::declare_part( " ;
    msg << arg_part.name();
    msg << " , " ;
    msg << arg_id ;
    msg << " ) FAILED, Part does not have a CellTopology" ;
    throw std::runtime_error( msg.str() );
  }

  const FilePart * const fp =
    internal_declare_part( m_schema , arg_part , arg_id ,
                           Element , top , m_dimension ,
                           0 , m_field_elem_attr );

  m_parts[ Element ].insert( i , fp );
}

void FileSchema::declare_part(
  Part     & arg_part ,
  int        arg_id ,
  EntityType arg_type )
{
  std::vector< const FilePart * >::iterator
    i = std::lower_bound( m_parts[ arg_type ].begin() ,
                          m_parts[ arg_type ].end() ,
                          arg_id , FilePartLess() );

  if ( i != m_parts[ arg_type ].end() && (*i)->m_identifier == arg_id ) {
    std::ostringstream msg ;
    msg << "phdmesh::exodus::FileSchema::declare_part( " ;
    msg << arg_part.name() ;
    msg << " , " ;
    msg << arg_id ;
    msg << " , " ;
    msg << entity_type_name( arg_type );
    msg << " ) FAILED, already called declare_part( " ;
    msg << (*i)->m_part.name() ;
    msg << " , " ;
    msg << arg_id ;
    msg << " , " ;
    msg << entity_type_name( arg_type );
    msg << " )" ;
    throw std::runtime_error( msg.str() );
  }

  const FilePart * const fp =
     internal_declare_part( m_schema , arg_part , arg_id ,
                            arg_type , NULL , m_dimension ,
                            0 , m_field_elem_attr );

  m_parts[ arg_type ].insert( i , fp );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

struct BlockEntityProc {
  entity_id_type entity_id ;
  unsigned       block_id ;
  unsigned       processor ;
};

struct LessIndices {

  bool operator()( const std::pair<int,const Entity*> & L ,
                   const std::pair<int,const Entity*> & R ) const
  {
    return L.first != R.first ? L.first < R.first :
           L.second->identifier() < R.second->identifier() ;
  }

  bool operator ()( const BlockEntityProc & L ,
                    const BlockEntityProc & R ) const
  {
     return L.block_id  != R.block_id ?  L.block_id  < R.block_id : (
            L.entity_id != R.entity_id ? L.entity_id < R.entity_id :
                                         L.processor < R.processor );
  }
};

void assign_contiguous_indices(
  BulkData & M ,
  unsigned entity_type ,
  const std::vector< std::pair<int,const Entity*> > & entities ,
  const FileSchema::IndexField & f )
{
  static const char method[] = "phdmesh::exodus::assign_contiguous_indices" ;

  const unsigned u_zero = 0 ;

  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();
  const ParallelMachine p_comm = M.parallel();

  std::vector< std::pair<int,const Entity*> >::const_iterator i ;

  unsigned unique_count = 0 ;

  // Get the total count, split evenly among processors, and then map
 
  entity_id_type end_id = 1 ;
  {
    const EntitySet & es = M.entities( entity_type );
    if ( ! es.empty() ) {
      end_id += (-- es.end() )->identifier();
    }
  }

  all_reduce( p_comm , Max<1>( & end_id ) );

  const double p_map = ((double) p_size) / ((double) end_id) ;

  CommAll all( p_comm );

  for ( i = entities.begin() ; i != entities.end() ; ++i ) {
    const unsigned p = (unsigned)( p_map * i->second->identifier() );
    all.send_buffer( p ).skip<entity_id_type>(2);
  }

  all.allocate_buffers( p_size / 4 );

  for ( i = entities.begin() ; i != entities.end() ; ++i ) {
    entity_id_type data[2] ;
    data[0] = i->first ;
    data[1] = i->second->identifier();

    const unsigned p = (unsigned)( p_map * data[1] );
    all.send_buffer( p ).pack<entity_id_type>( data , 2);
  }

  all.communicate();

  std::vector<BlockEntityProc> received ;

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {
      entity_id_type data[2] ;
      buf.unpack<entity_id_type>( data , 2 );
      BlockEntityProc tmp ;
      tmp.block_id  = data[0] ;
      tmp.entity_id = data[1] ;
      tmp.processor = p ;
      received.push_back( tmp );
    }
  }

  std::sort( received.begin() , received.end() , LessIndices() );

  unique_count = 0 ;

  {
    entity_id_type id = 0 ;

    for ( std::vector<BlockEntityProc>::iterator
          j = received.begin() ; j != received.end() ; ++j ) {
      if ( id != j->entity_id ) { id = j->entity_id ; ++unique_count ; }
    }
  }

  std::vector<unsigned> global_count( p_size , u_zero );

  allgather( p_comm , & unique_count , & global_count[0] , 1 );

  entity_id_type count = 0 ;

  for ( unsigned p = 0 ; p < p_rank ; ++p ) { count += global_count[p] ; }

  entity_id_type global_index = count ;

  for ( unsigned p = p_rank ; p < p_size ; ++p ) { count += global_count[p] ; }

  all.swap_send_recv();
  all.reset_buffers();

  {
    entity_id_type id = 0 ;

    for ( std::vector<BlockEntityProc>::iterator
          j = received.begin() ; j != received.end() ; ++j ) {

      if ( id != j->entity_id ) { id = j->entity_id ; ++global_index ; }

      entity_id_type data[2] ;
      data[0] = global_index ;
      data[1] = id ;
      all.send_buffer( j->processor ).pack<entity_id_type>( data , 2 );
    }
  }

  all.communicate();

  entity_id_type local_index = 0 ;

  for ( i = entities.begin() ; i != entities.end() ; ++i ) {
    const entity_id_type id = i->second->identifier();

    const unsigned p = (unsigned)( p_map * id );
    CommBuffer & buf = all.recv_buffer( p );
    entity_id_type tmp[2] ;
    buf.unpack<entity_id_type>( tmp , 2 );

    int * const entity_index = field_data( f , * i->second );

    entity_index[0] = tmp[0] ; // global index
    entity_index[1] = ++local_index ; // local index

    if ( id != tmp[1] ) {
      std::ostringstream msg ;
      msg << method << " FAILED" ;
      throw std::runtime_error(msg.str());
    }
  }
}

void assign_contiguous_indices( BulkData & M , unsigned entity_type ,
                                const FileSchema::IndexField & f )
{
  const EntitySet & es = M.entities( entity_type );

  std::vector< std::pair<int,const Entity*> > entities ;

  entities.reserve( es.size() );

  std::pair<int,const Entity*> tmp ; tmp.first = 0 ;

  for ( EntitySet::const_iterator j = es.begin() ; j != es.end() ; ++j ) {
    tmp.second = & *j ;
    entities.push_back( tmp );
  }

  assign_contiguous_indices( M , entity_type , entities , f );
}

}


void FileSchema::assign_indices( BulkData & arg_mesh ) const
{
  static const char method[] = "phdmesh::FileSchema::assign_indices" ;

  BulkData & M = arg_mesh ;
  const MetaData & S = m_schema ;

  S.assert_same_mesh_meta_data( method , M.mesh_meta_data() );

  //--------------------------------------------------------------------

  assign_contiguous_indices( M , Node , m_field_index );
  assign_contiguous_indices( M , Edge , m_field_index );
  assign_contiguous_indices( M , Face , m_field_index );

  //--------------------------------------------------------------------
  // Assign element indices by block and then element id order.
  // Count the number of owned elements per block.

  std::vector< std::pair<int,const Entity*> > entities ;

  entities.reserve( M.entities( Element ).size() );

  const std::vector<const FilePart*> & elem_parts = m_parts[ Element ];

  const unsigned num_elem_blk = elem_parts.size();

  const KernelSet & ks = M.kernels( Element );

  for ( KernelSet::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {

    const Kernel & k = *ik ;

    bool found = false ;

    for ( unsigned i = 0 ; i < num_elem_blk ; ++i ) {

      if ( k.has_superset( elem_parts[i]->m_part ) ) {

        if ( found ) {
          std::ostringstream msg ;
          msg << method
              << " FAILED: element contained in multiple blocks" ;
          throw std::runtime_error( msg.str() );
        }

        found = true ;

        std::pair<int,const Entity*> tmp ;
        tmp.first = elem_parts[i]->m_identifier ;

        for ( unsigned j = 0 ; j < k.size() ; ++j ) {
          tmp.second = k[j] ;
          entities.push_back( tmp );
        }
      }
    }
  }

  std::sort( entities.begin() , entities.end() , LessIndices() );

  // Assign element indices in element block order and then element id order

  assign_contiguous_indices( M , Element , entities , m_field_index );
}

}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#if defined( HAVE_SNL_EXODUSII )

#include <exodusII.h>
#include <ne_nemesisI.h>

enum { Maximum_Name_Length = MAX_STR_LENGTH };

#if defined( HAVE_MPI )

namespace {

void broadcast( phdmesh::ParallelMachine c , int p_root , int * value , int n )
{ MPI_Bcast( value , n , MPI_INT , p_root , c ); }

void broadcast( phdmesh::ParallelMachine c , int p_root , char * value , int n )
{ MPI_Bcast( value , n , MPI_BYTE , p_root , c ); }

void scatter( phdmesh::ParallelMachine c ,
              const unsigned p_root ,
              const unsigned * const send_size ,
              const unsigned recv_size ,
              void * const data )
{
  const unsigned p_rank = phdmesh::parallel_machine_rank( c );
  const unsigned p_size = phdmesh::parallel_machine_size( c );

  std::vector<int> count ;
  std::vector<int> displ ;

  if ( p_rank == p_root ) {
    count.resize( p_size );
    displ.resize( p_size );
    unsigned n = 0 ;
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      count[p] = p_root != p ? send_size[p] : 0 ;
      displ[p] = n ;
      n += send_size[p] ;
    }
  }

  int * const p_count = p_rank == p_root ? & count[0] : NULL ;
  int * const p_displ = p_rank == p_root ? & displ[0] : NULL ;
  int p_recv = p_rank == p_root ? 0 : recv_size ;

  MPI_Scatterv( data , p_count , p_displ , MPI_BYTE ,
                data , p_recv ,            MPI_BYTE ,
                p_root , c );

  // Now the local copy, if any

  if ( p_rank == p_root && send_size[ p_root ] && displ[ p_root ] ) {
    unsigned char * const dst = reinterpret_cast<unsigned char *>( data );
    unsigned char * const src = dst + displ[ p_root ] ;
    memmove( dst , src , send_size[ p_root ] );
  }
}

}

#else

namespace {

void broadcast( phdmesh::ParallelMachine , int , char * , int ) {}
void broadcast( phdmesh::ParallelMachine , int , int * , int ) {}

void scatter( phdmesh::ParallelMachine ,
              const unsigned ,
              const unsigned * const ,
              const unsigned ,
              void * const ) {}

}

#endif


namespace phdmesh {
namespace exodus {

//----------------------------------------------------------------------

FileOutput::~FileOutput()
{
  const BulkData       & M  = m_mesh ;
  const FileSchema & FS = m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = FS.m_io_rank ;

  if ( p_rank == p_write ) { ex_close( m_exo_id ); }
}

//----------------------------------------------------------------------

namespace {

template<class Op>
void chunk( const std::vector<const Entity*> & entities ,
            const FileSchema::IndexField & field_index ,
            const int index_begin ,
            const int index_end ,
            Op & op )
{
  const int maximum = op.max();

  const std::vector<const Entity *>::const_iterator i_end = entities.end();
        std::vector<const Entity *>::const_iterator i = entities.begin();

  int index = index_begin ;

  while ( index < index_end ) {

    const int index_end_chunk = index_end < index + maximum ?
                                index_end : index + maximum ;

    const std::vector<const Entity *>::const_iterator i_beg = i ;

    while ( i != i_end &&
            field_data( field_index , **i )[0] < index_end_chunk ) { ++i ; }

    op( index , index_end_chunk , i_beg , i );

    index = index_end_chunk ;
  }
}

//----------------------------------------------------------------------

struct WriteNodeIndexCoord {
  // Two buffers: communication and file data
  enum { size_per_node = 2 * ( sizeof(int) + 3 * sizeof(double) ) };

  FileOutput & output ;
  unsigned     maximum ;
  bool         writer ;
  int          exo_error ;
  const char * exo_func ;
  std::vector<int>    identifiers ;
  std::vector<double> coord_x_y_z ;

  unsigned max() const { return maximum ; }

  void operator()( const int , const int ,
                   const std::vector<const Entity *>::const_iterator ,
                   const std::vector<const Entity *>::const_iterator );

  WriteNodeIndexCoord( FileOutput & arg_output , unsigned );
};

WriteNodeIndexCoord::WriteNodeIndexCoord(
  FileOutput & arg_output , unsigned max_buffer )
  : output( arg_output ),
    maximum( max_buffer / size_per_node ),
    writer( false ),
    exo_error( 0 ),
    exo_func( NULL ),
    identifiers(),
    coord_x_y_z()
{
  const BulkData       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = SF.m_io_rank ;

  writer = p_rank == p_write ;
}

void WriteNodeIndexCoord::operator()(
  const int index_begin ,
  const int index_end ,
  const std::vector<const Entity *>::const_iterator send_begin ,
  const std::vector<const Entity *>::const_iterator send_end )
{
  const BulkData & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned send_size =
     ( sizeof(int) * 2 + sizeof(double) * 3 ) * ( send_end - send_begin );

  CommGather gather( p_comm , p_write , send_size );

  {
    CommBuffer & buf = gather.send_buffer();

    for ( std::vector<const Entity *>::const_iterator
          i = send_begin ; i != send_end ; ++i ) {
      const Entity & node = **i ;
      int ids[2] ;
      ids[0] = * field_data( SF.m_field_index , node );
      ids[1] =   node.identifier();
      const double * const coord =
        (double *) field_data( SF.m_field_node_coord , node );

      buf.pack<int>( ids , 2 );
      buf.pack<double>( coord , 3 );
    }
  }

  gather.communicate();

  if ( writer ) {

    const int exo_id = output.exo_id() ;
    const int number = index_end - index_begin ;
    const unsigned n = number ;

    if ( identifiers.size() < n ) { identifiers.resize( n ); }
    if ( coord_x_y_z.size() < n ) { coord_x_y_z.resize( n * 3 ); }

    double * const x = & coord_x_y_z[0] ;
    double * const y = x + n ;
    double * const z = y + n ;

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = gather.recv_buffer(p);
      while ( buf.remaining() ) {
        int ids[2] ;    buf.unpack<int>( ids , 2 );
        double xyz[3] ; buf.unpack<double>( xyz , 3 );

        const unsigned offset = ids[0] - index_begin ;
        identifiers[ offset ] = ids[1] ;
        x[ offset ] = xyz[0] ;
        y[ offset ] = xyz[1] ;
        z[ offset ] = xyz[2] ;
      }
    }

    if ( ! exo_error ) {
      exo_func = "ne_put_n_node_num_map" ;
      exo_error =
        ne_put_n_node_num_map( exo_id, index_begin, number, & identifiers[0] );
    }

    if ( ! exo_error ) {
      exo_func = "ne_put_n_coord" ;
      exo_error = ne_put_n_coord(exo_id,index_begin,number,x,y,z);
    }
  }
}

//----------------------------------------------------------------------

struct WriteElemIndexRelation {
  FileOutput & output ;
  const FilePart & part ;
  int          index_part ;
  unsigned     maximum ;
  bool         writer ;
  int          exo_error ;
  const char * exo_func ;

  std::vector<int> send_data ;
  std::vector<int> identifiers ;
  std::vector<int> relationivity ;

  WriteElemIndexRelation( FileOutput & ,
                         const FilePart & ,
                         int ,
                         unsigned max_buffer );

  unsigned max() const { return maximum ; }

  void operator()( const int , const int ,
                   const std::vector<const Entity *>::const_iterator ,
                   const std::vector<const Entity *>::const_iterator );
};

WriteElemIndexRelation::WriteElemIndexRelation(
  FileOutput & arg_output ,
  const FilePart   & arg_part ,
  int          arg_index_part ,
  unsigned max_buffer )
  : output( arg_output ), part( arg_part ),
    index_part( arg_index_part ), maximum( 0 ),
    writer( false ), exo_error( 0 ), exo_func( NULL ),
    send_data(), identifiers(), relationivity()
{
  const int i_zero = 0 ;
  const BulkData       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned send_data_num = part.m_topology->node_count + 2 ;

  // Two buffers: communication and file data

  maximum = max_buffer / ( 2 * send_data_num * sizeof(int) );

  send_data.resize( send_data_num , i_zero );

  writer = p_rank == p_write ;
}

void WriteElemIndexRelation::operator()(
  const int index_begin ,
  const int index_end ,
  const std::vector<const Entity *>::const_iterator send_begin ,
  const std::vector<const Entity *>::const_iterator send_end )
{
  const BulkData       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned block_id           = part.m_identifier ;
  const unsigned num_nodes_per_elem = part.m_topology->node_count ;
  const unsigned send_data_num      = num_nodes_per_elem + 2 ;

  const unsigned send_size =
    sizeof(int) * send_data_num * ( send_end - send_begin );

  CommGather gather( p_comm , p_write , send_size );

  {
    CommBuffer & buf = gather.send_buffer();

    for ( std::vector<const Entity *>::const_iterator
          i = send_begin ; i != send_end ; ++i ) {
      const Entity & elem = **i ;

      send_data[0] = * field_data( SF.m_field_index , elem );
      send_data[1] =   elem.identifier();

      PairIterRelation con = elem.relations( Node );

      if ( con.size() < num_nodes_per_elem ) {
        std::ostringstream msg ;
        msg << "P" << M.parallel_rank();
        msg << ": phdmesh::exodus::WriteElemIndexRelation FAILED, " ;
        msg << "expected " << num_nodes_per_elem ;
        msg << " but given " << con.size();
        throw std::runtime_error( msg.str() );
      }

      for ( unsigned j = 0 ; j < num_nodes_per_elem ; ++j , ++con ) {
        Entity & node = * con->entity();
        send_data[j+2] = * field_data( SF.m_field_index , node );
      }

      buf.pack<int>( & send_data[0] , send_data_num );
    }
  }

  gather.communicate();

  if ( writer ) {

    const int exo_id = output.exo_id() ;
    const int number = index_end - index_begin ;
    const unsigned n = number ;

    if ( identifiers.size() < n ) { identifiers.resize( n ); }
    if ( relationivity.size() < n * part.m_topology->node_count ) {
         relationivity.resize( n * part.m_topology->node_count );
    }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {

      CommBuffer & buf = gather.recv_buffer(p);

      while ( buf.remaining() ) {

        buf.unpack<int>( & send_data[0] , send_data_num );

        int index = send_data[0] ;

        unsigned offset = index - index_begin ;

        identifiers[offset] = send_data[1] ;

        int * const node_id = & relationivity[ offset * num_nodes_per_elem ];

        for ( unsigned j = 0 ; j < num_nodes_per_elem ; ++j ) {
          node_id[j] = send_data[j+2] ;
        }
      }
    }

    if ( ! exo_error ) {
      exo_func = "ne_put_n_elem_num_map" ;
      exo_error =
        ne_put_n_elem_num_map(exo_id,index_begin,number, & identifiers[0] );
    }

    if ( ! exo_error ) {
      const int index_begin_part = 1 + index_begin - index_part ;

      exo_func = "ne_put_n_elem_conn" ;
      exo_error =
        ne_put_n_elem_conn( exo_id , block_id ,
                            index_begin_part , number , & relationivity[0] );
    }
  }
}

unsigned count( const KernelSet & ks , const PartSet & ps )
{
  unsigned n = 0 ;

  KernelSet::const_iterator i ;
  for ( i = ks.begin() ; i != ks.end() ; ++i ) {
    if ( i->has_superset( ps ) ) {
      n += i->size();
    }
  }
  return n ;
}
 
std::string variable_name( EntityType        r ,
                           const Part      & p , 
                           const FieldBase & f ,
                           const unsigned    k )
{
  const unsigned f_num_dim = f.rank();

  std::string name( f.name() );

  if ( f_num_dim ) {
    unsigned dims[   MaximumFieldDimension ];
    unsigned indices[ MaximumFieldDimension ];

    const FieldBase::Restriction & dim = f.restriction( r , p );

    array_stride_to_natural_dimensions( f_num_dim , dim.stride , dims );
    array_stride_to_natural_indices( f_num_dim , dim.stride , k , indices );

    for ( unsigned i = 0 ; i < f_num_dim ; ++i ) {

      // Natural, as opposed to Fortran, ordering of dimensions:
      const ArrayDimTag & tag = * f.dimension_tags()[ ( f_num_dim - 1 ) - i ];

      std::string tmp = tag.to_string(dims[i],indices[i]);

      if ( tmp.size() ) {
        name.append( "_" );
        name.append( tmp );
      }
    }
  }

  return name ;
}


void variable_add( EntityType entity_type ,
                   const FieldBase                     & field ,
                   const std::vector<const FilePart *> & parts ,
                   std::vector< std::string >          & var_names ,
                   std::vector< FieldIO >              & spec )
{
  const unsigned size_begin = var_names.size();
  const unsigned field_num_dim = field.rank();

  for ( unsigned j = 0 ; j < parts.size() ; ++j ) {
    Part & p = parts[j]->m_part ;

    const FieldBase::Restriction & d = field.restriction( entity_type , p );

    if ( d.stride[0] ) { // Exists
      const unsigned n = array_stride_size( field_num_dim , d.stride );

      for ( unsigned k = 0 ; k < n ; ++k ) {
        const std::string name( variable_name(entity_type, p, field, k ) );
        var_names.push_back( name );
      }
    }
  }

  std::vector<std::string>::iterator i_beg = var_names.begin();

  std::advance( i_beg , size_begin );

  std::sort( i_beg , var_names.end() );

  var_names.erase( std::unique( i_beg , var_names.end() ) ,
                                        var_names.end() );
  
  FieldIO tmp ;

  tmp.m_field     = & field ;
  tmp.m_part      = NULL ;
  tmp.m_offset    = 0 ;
  tmp.m_var_index = 0 ;

  for ( unsigned j = 0 ; j < parts.size() ; ++j ) {
    tmp.m_part = parts[j] ;

    Part & p = parts[j]->m_part ;

    const FieldBase::Restriction & d = field.restriction( entity_type , p );
    const unsigned n = array_stride_size( field_num_dim , d.stride );

    for ( unsigned k = 0 ; k < n ; ++k ) {

      const std::string name( variable_name( entity_type , p , field , k ) );

      const std::vector<std::string>::iterator i =
        std::lower_bound( i_beg , var_names.end() , name );

      tmp.m_offset = k ;
      tmp.m_var_index = 1 + std::distance( var_names.begin() , i );

      spec.push_back( tmp );
    }
  }
}


struct EntityLess {
  bool operator()( const Entity * const lhs , const Entity * const rhs )
    { return lhs->identifier() < rhs->identifier(); }
};

void get_entities(
  std::vector< const Entity * > & tmp ,
  const BulkData & mesh ,
  EntityType type ,
  Part & part )
{
  const KernelSet & ks = mesh.kernels( type );

  unsigned count = 0 ;

  for ( KernelSet::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {
    if ( ik->has_superset( part ) ) { count += ik->size(); }
  }

  {
    const Entity * const en = NULL ;
    tmp.assign( count , en );
  }

  count = 0 ;
  for ( KernelSet::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {
    if ( ik->has_superset( part ) ) {
      for ( unsigned i = 0 ; i < ik->size() ; ++i , ++count ) {
        tmp[ count ] = (*ik)[i] ;
      }
    }
  }

  std::sort( tmp.begin() , tmp.end() , EntityLess() );
}

void get_entities(
  std::vector< const Entity * > & tmp ,
  const BulkData & mesh ,
  EntityType type ,
  Part & part1 ,
  Part & part2 )
{
  const KernelSet & ks = mesh.kernels( type );

  unsigned count = 0 ;

  for ( KernelSet::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {
    if ( ik->has_superset( part1 ) &&
         ik->has_superset( part2 ) ) { count += ik->size(); }
  }

  {
    const Entity * const en = NULL ;
    tmp.assign( count , en );
  }

  count = 0 ;
  for ( KernelSet::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {
    if ( ik->has_superset( part1 ) &&
         ik->has_superset( part2 ) ) {
      for ( unsigned i = 0 ; i < ik->size() ; ++i , ++count ) {
        tmp[ count ] = (*ik)[i] ;
      }
    }
  }

  std::sort( tmp.begin() , tmp.end() , EntityLess() );
}

}

//----------------------------------------------------------------------

FileOutput::FileOutput(
  const FileSchema  & arg_schema ,
  const BulkData & arg_mesh ,
  const std::string & arg_file_path ,
  const std::string & arg_title ,
  const bool          arg_storage_double ,
  const std::vector<const FieldBase * > & arg_fields ,
  const int * const arg_processor )
  : m_schema( arg_schema ),
    m_mesh( arg_mesh ),
    m_exo_id( 0 ),
    m_counter( 0 ),
    m_max_buffer( 0x0200000 )
{
  static const char method[] = "phdmesh::exodus::FileOutput::FileOutput" ;

  const int i_zero = 0 ;

  const BulkData       & M  = arg_mesh ;
  const MetaData     & SM = M.mesh_meta_data();
  const FileSchema & FS = arg_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = FS.m_io_rank ;

  const bool writer = p_write == p_rank ;

  const Part & universal_part = SM.universal_part();
        Part & owns_part      = SM.locally_owned_part();

  const std::vector<const FilePart *> & node_parts = FS.parts( Node );
  const std::vector<const FilePart *> & edge_parts = FS.parts( Edge );
  const std::vector<const FilePart *> & face_parts = FS.parts( Face );
  const std::vector<const FilePart *> & elem_parts = FS.parts( Element );

  //--------------------------------------------------------------------

  std::vector< std::string > name_node_var ;
  std::vector< std::string > name_elem_var ;
  std::vector<int> exist_elem_var ; 

  if ( arg_processor ) {
    const std::string name_processor("processor");

    FieldIO tmp ;
    tmp.m_field = NULL ;
    tmp.m_part  = NULL ;
    tmp.m_offset = 0 ;
    tmp.m_var_index = 1 ;

    if ( arg_processor[ Node ] ) {
      name_node_var.push_back( name_processor );
      m_field_node_universal.push_back( tmp );
    }

    if ( arg_processor[ Element ] ) {
      name_elem_var.push_back( name_processor );
      m_field_elem.push_back( tmp );
    }
  }

  // Universal nodal variables:

  for ( std::vector< const FieldBase * >::const_iterator
        i = arg_fields.begin() ; i != arg_fields.end() ; ++i ) {

    const FieldBase & f = **i ;
    const unsigned f_num_dim = f.rank();
    const FieldBase::Restriction & d = f.restriction( Node , universal_part );

    if ( d.stride[0] ) {
      FieldIO tmp ;

      tmp.m_field     = & f ;
      tmp.m_part      = NULL ;
      tmp.m_offset    = 0 ;
      tmp.m_var_index = 0 ;

      if ( ! f_num_dim ) {
        // Scalar
        name_node_var.push_back( f.name() );
        tmp.m_offset    = 0 ;
        tmp.m_var_index = name_node_var.size();
        m_field_node_universal.push_back( tmp );
      }
      else {
        const unsigned n = array_stride_size( f_num_dim , d.stride );

        for ( unsigned k = 0 ; k < n ; ++k ) {
          name_node_var.push_back( variable_name(Node, universal_part, f, k) );
          tmp.m_offset    = k ;
          tmp.m_var_index = name_node_var.size();
          m_field_node_universal.push_back( tmp );
        }
      }
    }
  }

  // Element variables:

  for ( std::vector< const FieldBase * >::const_iterator
        i = arg_fields.begin() ; i != arg_fields.end() ; ++i ) {

    const FieldBase & f = **i ;

    bool is_element_var = false ;

    for ( std::vector<FieldBase::Restriction>::const_iterator
          j =  f.restrictions().begin() ;
          j != f.restrictions().end() && ! is_element_var ; ++j ) {

      is_element_var = entity_type( j->key ) == Element ;
    }
    if ( is_element_var ) {
      variable_add( Element , f , elem_parts , name_elem_var , m_field_elem );
    }
  }

  exist_elem_var.resize( elem_parts.size() * name_elem_var.size() , i_zero);

  for ( unsigned j = 0 ; j < elem_parts.size() ; ++j ) {

    int * const exist = & exist_elem_var[ j * name_elem_var.size() ];

    for ( unsigned k = 0 ; k < m_field_elem.size() ; ++k ) {
      if ( m_field_elem[k].m_part == NULL ||
           m_field_elem[k].m_part == elem_parts[k] ) { 
        const int offset = m_field_elem[k].m_var_index - 1 ;
        exist[ offset ] = 1 ;
      }
    }
  }

  //--------------------------------------------------------------------
  // Sizes and counts

  const int num_dim       = FS.m_dimension ;
  const int num_node_sets = node_parts.size();
  const int num_edge_sets = edge_parts.size();
  const int num_face_sets = face_parts.size();
  const int num_side_sets = num_edge_sets + num_face_sets ;
  const int num_elem_blk  = elem_parts.size();

  const unsigned num_parts = num_node_sets + num_side_sets + num_elem_blk ;

  m_global_counts.resize( 4 + num_parts , i_zero );

  std::vector<int> part_count_local(  4 + num_parts , i_zero );

  {
    PartSet ps ;
    insert( ps , owns_part );

    part_count_local[ 0 ] = count( M.kernels( Node ) , ps );
    part_count_local[ 1 ] = count( M.kernels( Edge ) , ps );
    part_count_local[ 2 ] = count( M.kernels( Face ) , ps );
    part_count_local[ 3 ] = count( M.kernels( Element ) , ps );
  }

  unsigned n = 4 ;
  for ( int i = 0 ; i < num_node_sets ; ++i , ++n ) {
    Part & part = node_parts[i]->m_part ;
    PartSet ps ;
    insert( ps , part );
    insert( ps , owns_part );
    part_count_local[ n ] = count( M.kernels( Node ) , ps );
  }

  for ( int i = 0 ; i < num_edge_sets ; ++i , ++n ) {
    Part & part = edge_parts[i]->m_part ;
    PartSet ps ;
    insert( ps , part );
    insert( ps , owns_part );
    part_count_local[ n ] = count( M.kernels( Edge ) , ps );
  }

  for ( int i = 0 ; i < num_face_sets ; ++i , ++n ) {
    Part & part = face_parts[i]->m_part ;
    PartSet ps ;
    insert( ps , part );
    insert( ps , owns_part );
    part_count_local[ n ] = count( M.kernels( Face ) , ps );
  }

  for ( int i = 0 ; i < num_elem_blk ; ++i , ++n ) {
    Part & part = elem_parts[i]->m_part ;
    PartSet ps ;
    insert( ps , part );
    insert( ps , owns_part );
    part_count_local[ n ] = count( M.kernels( Element ) , ps );
  }

  all_reduce_sum( p_comm , & part_count_local[0] , & m_global_counts[0] , n );

  const int num_nodes_global = m_global_counts[0] ;
  const int num_elems_global = m_global_counts[3] ;

  const int * const global_num_elem_in_block =
    & m_global_counts[ 4 + num_node_sets + num_edge_sets + num_face_sets ];

  //--------------------------------------------------------------------

  int exo_error = 0 ;
  const char * exo_func = NULL ;

  if ( writer ) {

    // Create file:
    {
      int comp_ws = sizeof(double);
      int io_ws   = arg_storage_double ? sizeof(double) : sizeof(float) ;

      exo_func = "ex_create" ;
      m_exo_id = ex_create( arg_file_path.c_str() ,
                            EX_CLOBBER , & comp_ws , & io_ws );

      if ( m_exo_id < 0 ) { exo_error = -1 ; }
    }

    // Put title and sizes:
    if ( ! exo_error ) {
      exo_func = "ex_put_init" ;
      exo_error = ex_put_init( m_exo_id , arg_title.c_str() ,
                               num_dim , num_nodes_global , num_elems_global ,
                               num_elem_blk , num_node_sets , num_side_sets );
    }

    // Put the model nodal coordinate names:
    if ( ! exo_error ) {
      char coord_name_0[ Maximum_Name_Length ] ;
      char coord_name_1[ Maximum_Name_Length ] ;
      char coord_name_2[ Maximum_Name_Length ] ;

      char * coord_names[3] = { coord_name_0 , coord_name_1 , coord_name_2 };

      const ArrayDimTag & dim_tag = *FS.m_field_node_coord.dimension_tags()[0];

      const std::string & coord_name_base = FS.m_field_node_coord.name();

      for ( int i = 0 ; i < num_dim ; ++i ) {
        strcpy( coord_names[i] , coord_name_base.c_str() );
        strcat( coord_names[i] , "_" );
        strcat( coord_names[i] , dim_tag.to_string( num_dim , i ).c_str() );
      }

      exo_func = "ex_put_coord_names" ;
      exo_error = ex_put_coord_names( m_exo_id , coord_names );
    }

    // Put element blocks' description:
    if ( ! exo_error ) {
      char * const char_ptr_null = NULL ;

      std::vector<char> tmp( num_elem_blk * Maximum_Name_Length );

      std::vector<int>    elem_blk_id( num_elem_blk , 0 );
      std::vector<char *> elem_blk_name( num_elem_blk , char_ptr_null );
      std::vector<char *> elem_type( num_elem_blk , char_ptr_null );
      std::vector<int>    num_nodes_per_elem( num_elem_blk , 0 );
      std::vector<int>    num_attr( num_elem_blk , 0 );

      for ( int i = 0 ; i < num_elem_blk ; ++i ) {
        const FilePart   & fp = * elem_parts[i] ;
        elem_blk_id[i]        = fp.m_identifier ;
        num_nodes_per_elem[i] = fp.m_topology->node_count ;
        num_attr[i]           = fp.m_number_attr ;
        elem_type[i]          = & tmp[ i * Maximum_Name_Length ] ;
        elem_blk_name[i]      = const_cast<char*>( fp.m_part.name().c_str() );

        std::string tmp_elem_type ;
        int tmp_num_attr = 0 ;

        map_to_exodus( * fp.m_topology , num_dim ,
                       tmp_elem_type , tmp_num_attr );

        strcpy( elem_type[i] , tmp_elem_type.c_str() );
      }

      exo_func = "ex_put_concat_elem_block" ;
      exo_error = ex_put_concat_elem_block( m_exo_id ,
                                            & elem_blk_id[0] ,
                                            & elem_type[0] ,
                                              global_num_elem_in_block ,
                                            & num_nodes_per_elem[0] ,
                                            & num_attr[0] ,
                                            1 /* have identifier maps */ );

      if ( ! exo_error ) {
        exo_func = "ex_put_names(element)" ;
        exo_error = ex_put_names( m_exo_id, EX_ELEM_BLOCK, & elem_blk_name[0]);
      }
    }

    // Put results variable description:

    if ( ! exo_error ) {

      int num_var_global = 0 ;
      int num_var_node = name_node_var.size();
      int num_var_elem = name_elem_var.size();
      int num_var_side_set = 0 ;
      int num_var_node_set = 0 ;

      exo_func = "ex_put_all_var_param" ;
      exo_error = ex_put_all_var_param( m_exo_id ,
                                        num_var_global ,
                                        num_var_node ,
                                        num_var_elem ,
                                        & exist_elem_var[0] ,
                                        num_var_node_set ,
                                        NULL ,
                                        num_var_side_set ,
                                        NULL );

      if ( ! exo_error && num_var_node ) {
        std::vector<char *> ptr_var_names( name_node_var.size() );

        for ( unsigned i = 0 ; i < name_node_var.size() ; ++i ) {
          ptr_var_names[i] = const_cast<char*>( name_node_var[i].c_str() );
        }

        exo_func = "ex_put_var_names(node)" ;
        exo_error = ex_put_var_names( m_exo_id , "n" ,
                                      num_var_node ,
                                      & ptr_var_names[0] );
      }

      if ( ! exo_error && num_var_elem ) {
        std::vector<char *> ptr_var_names( name_elem_var.size() );

        for ( unsigned i = 0 ; i < name_elem_var.size() ; ++i ) {
          ptr_var_names[i] = const_cast<char*>( name_elem_var[i].c_str() );
        }

        exo_func = "ex_put_var_names(element)" ;
        exo_error = ex_put_var_names( m_exo_id , "e" ,
                                      num_var_elem ,
                                      & ptr_var_names[0] );
      }
    }
  }

  broadcast( p_comm , p_write , & exo_error , 1 );

  //--------------------------------------------------------------------
  // Write nodal identifiers and genesis coordinates.

  if ( ! exo_error ) {
    WriteNodeIndexCoord work( *this , m_max_buffer );

    const int index = 1 ;
    const int index_end = index + num_nodes_global ;

    std::vector< const Entity * > tmp ;
    get_entities( tmp , M , Node , owns_part );

    chunk( tmp , FS.m_field_index , index , index_end , work );

    exo_error = work.exo_error ;

    broadcast( p_comm , p_write , & exo_error , 1 );

    exo_func  = work.exo_func ;
  }

  //--------------------------------------------------------------------
  // Write element block identifiers and relations,
  // element indices are partitioned by element block.

  if ( ! exo_error ) {

    int index = 1 ;

    for ( int j = 0 ; j < num_elem_blk ; ++j ) {

      const int index_end = index + global_num_elem_in_block[j] ;

      WriteElemIndexRelation work(*this, * elem_parts[j], index, m_max_buffer);

      std::vector<const Entity *> tmp ;
      get_entities( tmp , M , Element , owns_part , elem_parts[j]->m_part );

      chunk( tmp , FS.m_field_index , index , index_end , work );

      index = index_end ;

      if ( ! exo_error ) {
        exo_error = work.exo_error ;
        exo_func  = work.exo_func ;
      }
    }

    broadcast( p_comm , p_write , & exo_error , 1 );
  }

  if ( exo_error ) {
    std::ostringstream msg ;
    if ( writer ) {
      msg << method << " FAILED calling " << exo_func << " with " << exo_error ;
    }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void pack_entity_proc(
  const std::vector<const Entity *>::const_iterator send_begin ,
  const std::vector<const Entity *>::const_iterator send_end ,
  const FileSchema::IndexField      & f_index ,
  CommBuffer & buf )
{
  double item_buffer[2] ;

  for ( std::vector<const Entity *>::const_iterator
        i = send_begin ; i != send_end ; ++i ) {
    const Entity & e = **i ;

    item_buffer[0] = * field_data( f_index , e );
    item_buffer[1] =   e.owner_rank();

    buf.pack<double>( item_buffer , 2 );
  }
}

// Work-around for shoddy pathscale compiler:

void pack_data_2( const double * const b , CommBuffer & buf )
{ buf.pack<double>( b , 2 ); }

template<typename T>
void pack_entity_data(
  const std::vector<const Entity *>::const_iterator ib ,
  const std::vector<const Entity *>::const_iterator ie ,
  const FileSchema::IndexField      & f_index ,
  const FieldBase                   & f_data ,
  const unsigned                      offset ,
  CommBuffer & buf )
{
  double item_buffer[2] ;

  for ( std::vector<const Entity *>::const_iterator
        i = ib ; i != ie ; ++i ) {
    const Entity & e = **i ;
    T * const data = reinterpret_cast<T*>( field_data(f_data,e) );

    item_buffer[0] = (double)( * field_data( f_index , e ) );
    item_buffer[1] = (double)( data[ offset ] );

    pack_data_2( item_buffer , buf );
  }
}

void pack_entity_data(
  const std::vector<const Entity *>::const_iterator ib ,
  const std::vector<const Entity *>::const_iterator ie ,
  const FileSchema::IndexField      & f_index ,
  const FieldBase                   & f_data ,
  const unsigned                      offset ,
  CommBuffer & buf )
{
  switch( f_data.numeric_type_ordinal() ) {
  case NumericEnum< signed char >::value :
    pack_entity_data<signed char>( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< unsigned char >::value :
    pack_entity_data<unsigned char>( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< signed short >::value :
    pack_entity_data< signed short >( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< unsigned short >::value :
    pack_entity_data< unsigned short >( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< signed int >::value :
    pack_entity_data< signed int >( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< unsigned int >::value :
    pack_entity_data< unsigned int >( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< signed long >::value :
    pack_entity_data< signed long >( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< unsigned long >::value :
    pack_entity_data< unsigned long >( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< float >::value :
    pack_entity_data< float >( ib,ie, f_index, f_data, offset, buf );
    break ;
  case NumericEnum< double >::value :
    pack_entity_data< double >( ib,ie, f_index, f_data, offset, buf );
    break ;
  default : break ;
  }
}

//----------------------------------------------------------------------

struct WriteNodeValues {

  FileOutput & output ;
  FieldIO * field ;
  unsigned     maximum ;
  bool         writer ;
  int          exo_error ;
  const char * exo_func ;
  std::vector<double> recv_buffer ;

  unsigned max() const { return maximum ; }

  void operator()( const int , const int ,
                   const std::vector<const Entity *>::const_iterator ,
                   const std::vector<const Entity *>::const_iterator );

  WriteNodeValues( FileOutput & , unsigned );

  void set_field( FieldIO & f ) { field = & f ; }
};

WriteNodeValues::WriteNodeValues(
  FileOutput & arg_output , unsigned max_buffer )
  : output( arg_output ),
    field( NULL ),
    maximum( 0 ),
    writer( false ),
    exo_error( 0 ),
    exo_func( NULL ),
    recv_buffer()
{
  const BulkData       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = SF.m_io_rank ;

  writer = p_write == p_rank ;

  // Two buffers for field data and one for index

  maximum = max_buffer / ( 3 * sizeof(double) );
}

void WriteNodeValues::operator()(
  const int index_begin ,
  const int index_end ,
  const std::vector<const Entity *>::const_iterator send_begin ,
  const std::vector<const Entity *>::const_iterator send_end )
{
  const BulkData       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned send_size = 2 * sizeof(double) * ( send_end - send_begin );

  double item_buffer[2] ;

  CommGather gather( p_comm , p_write , send_size );

  if ( field->m_field ) {
    pack_entity_data( send_begin , send_end , SF.m_field_index ,
                      * field->m_field , field->m_offset ,
                       gather.send_buffer() );
  }
  else {
    pack_entity_proc( send_begin, send_end,
                      SF.m_field_index , gather.send_buffer() );
  }

  gather.communicate();

  if ( writer ) {

    const int exo_id   = output.exo_id() ;
    const int exo_step = output.exo_step();
    const int number   = index_end - index_begin ;
    const unsigned n = number ;

    if ( recv_buffer.size() < n ) { recv_buffer.resize( n ); }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = gather.recv_buffer(p);
      while ( buf.remaining() ) {
        buf.unpack<double>( item_buffer , 2 );

        const unsigned index  = (unsigned) item_buffer[0] ;
        const unsigned offset = index - index_begin ;

        recv_buffer[ offset ] = item_buffer[1] ;
      }
    }

    exo_func = "ne_put_nodal_var_slab" ;
    exo_error = ne_put_nodal_var_slab( exo_id , exo_step, field->m_var_index ,
                                       index_begin , number ,
                                       & recv_buffer[0] );
  }
}

//----------------------------------------------------------------------

struct WriteElemValues {

  FileOutput & output ;
  FieldIO * field ;
  int          index_part ;
  unsigned     maximum ;
  bool         writer ;
  int          exo_ident ;
  int          exo_error ;
  const char * exo_func ;
  std::vector<double> recv_buffer ;

  unsigned max() const { return maximum ; }

  void operator()( const int , const int ,
                   const std::vector<const Entity *>::const_iterator ,
                   const std::vector<const Entity *>::const_iterator );

  WriteElemValues( FileOutput & , unsigned );

  void set_field( FieldIO & f , int index , int ident )
    { field = & f ; index_part = index ; exo_ident = ident ; }
};

WriteElemValues::WriteElemValues(
  FileOutput & arg_output , unsigned max_buffer )
  : output( arg_output ),
    field( NULL ),
    index_part( 0 ),
    maximum( 0 ),
    writer( false ),
    exo_ident( 0 ),
    exo_error( 0 ),
    exo_func( NULL ),
    recv_buffer()
{
  const BulkData       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = SF.m_io_rank ;

  writer = p_write == p_rank ;

  // Two buffers for field data and one for index

  maximum = max_buffer / ( 3 * sizeof(double) );
}

void WriteElemValues::operator()(
  const int index_begin ,
  const int index_end ,
  const std::vector<const Entity *>::const_iterator send_begin ,
  const std::vector<const Entity *>::const_iterator send_end )
{
  const BulkData       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned send_size = 2 * sizeof(double) * ( send_end - send_begin );

  double item_buffer[2] ;

  CommGather gather( p_comm , p_write , send_size );

  if ( field->m_field ) {
    pack_entity_data( send_begin, send_end, SF.m_field_index ,
                      * field->m_field , field->m_offset ,
                       gather.send_buffer() );
  }
  else {
    pack_entity_proc( send_begin, send_end,
                      SF.m_field_index , gather.send_buffer() );
  }

  gather.communicate();

  if ( writer ) {

    const int exo_id   = output.exo_id() ;
    const int exo_step = output.exo_step();
    const int number   = index_end - index_begin ;
    const unsigned n = number ;
    const int index_begin_part = 1 + index_begin - index_part ;

    if ( recv_buffer.size() < n ) { recv_buffer.resize( n ); }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = gather.recv_buffer(p);
      while ( buf.remaining() ) {
        buf.unpack<double>( item_buffer , 2 );

        const unsigned index  = (unsigned) item_buffer[0] ;
        const unsigned offset = index - index_begin ;

        recv_buffer[ offset ] = item_buffer[1] ;
      }
    }

    exo_func = "ne_put_elem_var_slab" ;
    exo_error = ne_put_elem_var_slab( exo_id , exo_step,
                                      field->m_var_index ,
                                      exo_ident ,
                                      index_begin_part , number ,
                                      & recv_buffer[0] );
  }
}

}

//----------------------------------------------------------------------

void FileOutput::write( double time_value )
{
  static const char method[] = "phdmesh::exodus::FieldIO::wrote" ;

  const BulkData       & M  = m_mesh ;
  const MetaData     & SM = M.mesh_meta_data();
  const FileSchema & FS = m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_rank = M.parallel_rank() ;
  const unsigned  p_write = FS.m_io_rank ;

  const std::vector<const FilePart *> & node_parts = FS.parts( Node );
  const std::vector<const FilePart *> & edge_parts = FS.parts( Edge );
  const std::vector<const FilePart *> & face_parts = FS.parts( Face );
  const std::vector<const FilePart *> & elem_parts = FS.parts( Element );

  const unsigned num_node_parts = node_parts.size();
  const unsigned num_edge_parts = edge_parts.size();
  const unsigned num_face_parts = face_parts.size();
  const unsigned num_elem_parts = elem_parts.size();

  const bool writer = p_write == p_rank ;

  Part & owns_part = SM.locally_owned_part();

  int exo_error = 0 ;
  const char * exo_func = NULL ;

  ++m_counter ;

  if ( writer ) {
    exo_func = "ex_put_time" ;
    exo_error = ex_put_time( m_exo_id , m_counter , & time_value );
  }

  // exo_error is NOT parallel consistent

  {
    const int num_nodes_global = m_global_counts[ Node ] ;

    WriteNodeValues work( *this , m_max_buffer );

    for ( std::vector<FieldIO>::iterator
          i =  m_field_node_universal.begin() ;
          i != m_field_node_universal.end() ; ++i ) {

      work.set_field( *i );

      const int index = 1 ;
      const int index_end = index + num_nodes_global ;

      std::vector< const Entity * > tmp ;
      get_entities( tmp , M , Node , owns_part );

      chunk( tmp , FS.m_field_index , index , index_end , work );
    }

    exo_func  = work.exo_func ;
    exo_error = work.exo_error ;

    broadcast( p_comm , p_write , & exo_error , 1 );
  }

  if ( ! exo_error ) {
    // Iteration: Variables, element blocks, element chunking

    const int * const global_num_elem_in_block =
      & m_global_counts[ 4 + num_node_parts + num_edge_parts + num_face_parts ];

    WriteElemValues work( *this , m_max_buffer );

    for ( std::vector<FieldIO>::iterator
          k =  m_field_elem.begin() ;
          k != m_field_elem.end() ; ++k ) {

      const FilePart * const part = k->m_part ;

      int index = 1 ;

      for ( unsigned j = 0 ; j < num_elem_parts ; ++j ) {

        const int index_end = index + global_num_elem_in_block[j] ;

        if ( part == NULL || part == elem_parts[j] ) {

          work.set_field( *k , index , elem_parts[j]->m_identifier );

          std::vector<const Entity *> tmp ;
          get_entities( tmp , M , Element, owns_part , elem_parts[j]->m_part );

          chunk( tmp , FS.m_field_index , index , index_end , work );

          if ( ! exo_error ) {
            exo_error = work.exo_error ;
            exo_func  = work.exo_func ;
          }
        }
        index = index_end ;
      }
    }

    exo_func  = work.exo_func ;
    exo_error = work.exo_error ;
    broadcast( p_comm , p_write , & exo_error , 1 );
  }

  if ( exo_error ) {
    std::ostringstream msg ;
    if ( p_rank == p_write ) {
      msg << method << " FAILED calling " << exo_func << " with " << exo_error ;
    }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Open and read file for meta-data

namespace {

std::pair< const CellTopology * , int >
map_from_exodus( const char * const element_type ,
                 const int node_count ,
                 const int space_dim )
{
  std::pair< const CellTopology * , int > result ;
  result.first = NULL ;
  result.second = 0 ;

  if ( 0 == strncasecmp( "CIRCLE" , element_type , 3 ) ) {
    result.second = 1 ;
  }
  else if ( 0 == strncasecmp( "SPHERE" , element_type , 3 ) ) {
    result.second = 1 ;
  }
  else if ( 0 == strncasecmp( "TRUSS" , element_type , 3 ) ) {
    if ( node_count == 2 ) {
      result.first  = cell_topology< Line<2> >();
      result.second = 1 ;
    }
    else if ( node_count == 3 ) {
      result.first  = cell_topology< Line<3> >();
      result.second = 1 ;
    }
  }
  else if ( 0 == strncasecmp( "BEAM" , element_type , 3 ) ) {
    if ( node_count == 2 ) {
      result.first  = cell_topology< Line<2> >();
      result.second = space_dim == 2 ? 3 : (
                      space_dim == 3 ? 7 : 0 );
    }
    else if ( node_count == 3 ) {
      result.first  = cell_topology< Line<3> >();
      result.second = space_dim == 2 ? 3 : (
                      space_dim == 3 ? 7 : 0 );
    }
  }
  else if ( 0 == strncasecmp( "SHELL" , element_type , 3 ) ) {
     if ( node_count == 2 ) {
       result.first  = cell_topology< ShellLine<2> >();
       result.second = 1 ;
     }
     else if ( node_count == 3 ) {
       result.first  = cell_topology< ShellLine<3> >();
       result.second = 1 ;
     }
     else if ( node_count == 4 ) {
       result.first  = cell_topology< ShellQuadrilateral<4> >();
       result.second = 1 ;
     }
     else if ( node_count == 8 ) {
       result.first  = cell_topology< ShellQuadrilateral<8> >();
       result.second = 1 ;
     }
     else if ( node_count == 9 ) {
       result.first  = cell_topology< ShellQuadrilateral<9> >();
       result.second = 1 ;
     }
  }
  else if ( 0 == strncasecmp( "QUAD" , element_type , 3 ) ) {
    if ( space_dim == 2 ) {
      if ( node_count == 4 ) {
        result.first = cell_topology< Quadrilateral<4> >();
      }
      else if ( node_count == 8 ) {
        result.first = cell_topology< Quadrilateral<8> >();
      }
      else if ( node_count == 9 ) {
        result.first = cell_topology< Quadrilateral<9> >();
      }
    }
    else if ( space_dim == 3 ) {
      if ( node_count == 4 ) {
        result.first  = cell_topology< ShellQuadrilateral<4> >();
        result.second = 1 ;
      }
      else if ( node_count == 8 ) {
        result.first  = cell_topology< ShellQuadrilateral<8> >();
        result.second = 1 ;
      }
      else if ( node_count == 9 ) {
        result.first  = cell_topology< ShellQuadrilateral<9> >();
        result.second = 1 ;
      }
    }
  }
  else if ( 0 == strncasecmp( "TRIANGLE" , element_type , 3 ) ) {
    if ( space_dim == 2 ) {
      if ( node_count == 3 ) {
        result.first = cell_topology< Triangle<3> >();
      }
      else if ( node_count == 6 ) {
        result.first = cell_topology< Triangle<6> >();
      }
    }
    else if ( space_dim == 3 ) {
      if ( node_count == 3 ) {
        result.first = cell_topology< ShellTriangle<3> >();
      }
      else if ( node_count == 6 ) {
        result.first = cell_topology< ShellTriangle<6> >();
      }
    }
  }
  else if ( 0 == strncasecmp( "PYRAMID" , element_type , 3 ) ) {
    if ( node_count == 5 ) {
      result.first = cell_topology< Pyramid<5> >();
    }
    else if ( node_count == 13 ) {
      result.first = cell_topology< Pyramid<13> >();
    }
    else if ( node_count == 14 ) {
      result.first = cell_topology< Pyramid<14> >();
    }
  }
  else if ( 0 == strncasecmp( "TETRA" , element_type , 3 ) ) {
    if ( node_count == 4 ) {
      result.first = cell_topology< Tetrahedron<4> >();
    }
    else if ( node_count == 10 ) {
      result.first = cell_topology< Tetrahedron<10> >();
    }
  }
  else if ( 0 == strncasecmp( "WEDGE" , element_type , 3 ) ) {
    if ( node_count == 6 ) {
      result.first = cell_topology< Wedge<6> >();
    }
    else if ( node_count == 15 ) {
      result.first = cell_topology< Wedge<15> >();
    }
    else if ( node_count == 18 ) {
      result.first = cell_topology< Wedge<18> >();
    }
  }
  else if ( 0 == strncasecmp( "HEX" , element_type , 3 ) ) {
    if ( node_count == 8 ) {
      result.first = cell_topology< Hexahedron<8> >();
    }
    else if ( node_count == 20 ) {
      result.first = cell_topology< Hexahedron<20> >();
    }
    else if ( node_count == 27 ) {
      result.first = cell_topology< Hexahedron<27> >();
    }
  }

  return result ;
}
 

enum { MAX_NUM_ATTR = 8 };

struct exo_elem_block_data {
  char name[ Maximum_Name_Length ];
  char type[ Maximum_Name_Length ];
  char attr[ MAX_NUM_ATTR ][ Maximum_Name_Length ];
  int  block_id ;
  int  num_nodes ;
  int  num_attr ;
  int  error ;
};

}

FileSchema::FileSchema(
  MetaData & arg_mesh_meta_data ,
  const FieldBase & arg_node_coordinates ,
  const FileSchema::AttributeField  & arg_elem_attributes ,
  const std::string     & arg_file_path ,
  ParallelMachine         arg_comm ,
  const unsigned          arg_reader_rank )
  : m_schema( arg_mesh_meta_data ),
    m_io_rank( arg_reader_rank ),
    m_dimension( arg_node_coordinates.max_size( Node ) ),
    m_field_node_coord( arg_node_coordinates ),
    m_field_elem_attr(  arg_elem_attributes ),
    m_field_index( exo_index( arg_mesh_meta_data ) )
{
  static const char method[] = "phdmesh::exodus::FileSchema::FileSchema" ;

  ParallelMachine p_comm = arg_comm ;
  const unsigned  p_rank = parallel_machine_rank( arg_comm );
  const unsigned  p_read = arg_reader_rank ;

  //--------------------------------------------------------------------

  verify_node_coordinate_field( arg_node_coordinates );

  //--------------------------------------------------------------------

  int exo_data[ 32 ];

  int exo_id = 0 ;
  int exo_error = 0 ;
  const char * exo_func = NULL ;

  const bool reader = p_rank == p_read ;

  if ( reader ) { // Open file:
    int comp_ws = sizeof(double) ;
    int io_ws   = 0 ;
    float version = 0 ;

    exo_func = "ex_open" ;
    exo_id = ex_open( arg_file_path.c_str() , EX_READ ,
                        & comp_ws , & io_ws , & version );

    if ( exo_id < 0 ) { exo_error = -1 ; }
  }

  // Get sizes:

  if ( reader && ! exo_error ) {
    char title[ MAX_LINE_LENGTH ];
    exo_func = "ex_get_init" ;
    exo_error = ex_get_init( exo_id , title ,
                             exo_data + 1 ,
                             exo_data + 2 ,
                             exo_data + 3 ,
                             exo_data + 4 ,
                             exo_data + 5 ,
                             exo_data + 6 );
  }

  exo_data[0] = exo_error ;

  broadcast( p_comm , p_read , exo_data , 7 );

  exo_error = exo_data[0] ;

  const int num_dim          = exo_data[1] ;
  // const int num_nodes_global = exo_data[2] ;
  // const int num_elems_global = exo_data[3] ;
  const int num_elem_blk     = exo_data[4] ;
  // const int num_node_sets    = exo_data[5] ;
  // const int num_side_sets    = exo_data[6] ;

  if ( num_dim != (int) m_dimension ) {
    std::ostringstream msg ;
    if ( reader ) {
      msg << method << " FAILED: incompatible spatial dimensions "
          << m_dimension << " != "
          << num_dim ;
    }
    throw std::runtime_error( msg.str() );
  }

  //--------------------------------------------------------------------
  // Get element blocks names and declare element parts

  if ( ! exo_error ) {

    std::vector<exo_elem_block_data> block_data( num_elem_blk );

    {
      char * begin = reinterpret_cast<char *>( & block_data[0] );
      char * end   = reinterpret_cast<char *>( & block_data[ num_elem_blk ] );
      const unsigned n = end - begin ;

      memset( begin , 0 , n );
    }

    if ( reader ) {
      const unsigned num_names =
        num_elem_blk < MAX_NUM_ATTR ? MAX_NUM_ATTR : num_elem_blk ;

      std::vector<char*> names( num_names );

      std::vector<int> ids( num_elem_blk );

      exo_func = "ex_get_elem_blk_ids" ;
      exo_error = ex_get_elem_blk_ids( exo_id , & ids[0] );

      if ( ! exo_error ) {

        for ( int i = 0 ; i < num_elem_blk ; ++i ) {
          block_data[i].block_id = ids[i] ;
          names[i] = block_data[i].name ;
        }

        exo_func = "ex_get_names" ;
        exo_error = ex_get_names( exo_id , EX_ELEM_BLOCK , & names[0] );

        if ( 0 < exo_error ) { // Names are not defined
          for ( int i = 0 ; i < num_elem_blk ; ++i ) {
            sprintf( block_data[i].name , "block_%d" , ids[i] );
          }
          exo_error = 0 ;
        }
      }

      for ( int i = 0 ; ! exo_error && i < num_elem_blk ; ++i ) {
        int num_elem_this_blk ; // discard for now

        exo_func = "ex_get_elem_block" ;
        exo_error = ex_get_elem_block( exo_id , ids[i] ,
                                       block_data[i].type ,
                                       & num_elem_this_blk ,
                                       & block_data[i].num_nodes ,
                                       & block_data[i].num_attr );

        if ( ! exo_error ) {
          for ( unsigned j = 0 ; j < MAX_NUM_ATTR ; ++j ) {
            names[i] = block_data[i].attr[j] ;
          }
          exo_func = "ex_get_elem_attr_names" ;
          exo_error = ex_get_elem_attr_names( exo_id , ids[i] , & names[0] );

          if ( 0 < exo_error ) {
            for ( unsigned j = 0 ; j < MAX_NUM_ATTR ; ++j ) {
              block_data[i].attr[j][0] = 0 ;
            }
            exo_error = 0 ;
          }
        }
      }

      block_data[0].error = exo_error ;
    }

    {
      char * begin = reinterpret_cast<char *>( & block_data[0] );
      char * end   = reinterpret_cast<char *>( & block_data[ num_elem_blk ] );
      const unsigned n = end - begin ;

      broadcast( p_comm , p_read , begin , n );

      exo_error = block_data[0].error ;
    }

    if ( ! exo_error ) {
      for ( int i = 0 ; i < num_elem_blk ; ++i ) {

        std::pair< const CellTopology * , int > elem_info =
          map_from_exodus( block_data[i].type ,
                           block_data[i].num_nodes ,
                           num_dim );

        Part & part = m_schema.declare_part( std::string(block_data[i].name) );

        const CellTopology * top = elem_info.first ;

        if ( NULL == top || elem_info.second != block_data[i].num_attr ) {
          std::ostringstream msg ;
          msg << "phdmesh::exodus::FileSchema FAILED, Read unknown element type from ExodusII : type = " ;
          msg << block_data[i].type ;
          msg << " , nodes = " ;
          msg << block_data[i].num_nodes ;
          msg << " , attr = " ;
          msg << block_data[i].num_attr ;
          throw std::runtime_error( msg.str() );
        }

        set_cell_topology( part , top );

        const FilePart * fp =
          internal_declare_part( m_schema ,
                                 part ,
                                 block_data[i].block_id ,
                                 Element , top ,
                                 m_dimension ,
                                 block_data[i].num_attr ,
                                 m_field_elem_attr );

        m_parts[ Element ].push_back( fp );
      }
    }
  }

  //--------------------------------------------------------------------

  if ( exo_error ) {
    std::ostringstream msg ;
    if ( reader ) {
      msg << method << " FAILED calling " << exo_func << " with " << exo_error ;
    }
    throw std::runtime_error( msg.str() );
  }

  if ( reader ) { ex_close( exo_id ); }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

// map index to processor:  p = ( p_size * i ) / n_total

unsigned index_processor_map( const unsigned p_size ,
                              const unsigned i ,
                              const unsigned n_total )
{
  return ( p_size * i ) / n_total ;
}

unsigned index_processor_count( const unsigned p ,
                                const unsigned p_size ,
                                const unsigned i ,
                                const unsigned n ,
                                const unsigned n_total )
{
  // The beginning and ending indices for the given processor
  // within the span of indices [i,i+n)
  unsigned i_beg = ( n_total * p ) / p_size ;
  unsigned i_end = ( n_total * ( p + 1 ) ) / p_size ;
  if ( i_beg < i ) { i_beg = i ; }
  if ( i + n < i_end ) { i_end = i + n ; }
  return i_beg < i_end ? ( i_end - i_beg ) : 0 ;
}

struct NodeData {
  double coord[3] ;
  int    ident ;
  int    index ;

  static unsigned size_of()
    {
      NodeData * tmp = NULL ;
      return reinterpret_cast<char*>(tmp+1) - reinterpret_cast<char*>(tmp);
    }
};

struct less_NodeData {
  bool operator()( const NodeData & lhs , const int rhs ) const
    { return lhs.index < rhs ; }
};

}

FileInput::~FileInput()
{
  const BulkData       & M  = m_mesh ;
  const FileSchema & FS = m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_read = FS.m_io_rank ;

  if ( p_rank == p_read ) { ex_close( m_exo_id ); }
}

FileInput::FileInput(
  const FileSchema  & arg_schema ,
        BulkData        & arg_mesh ,
  const std::string & arg_file_path ,
  const std::vector< const FieldBase * > & arg_fields )
  : m_schema( arg_schema ),
    m_mesh( arg_mesh ),
    m_exo_id( 0 ),
    m_counter( 0 ),
    m_max_buffer( 0x0200000 )
{
  static const char method[] = "phdmesh::exodus::FileInput::FileInput" ;

        BulkData       & M  = arg_mesh ;
  const MetaData     & SM = M.mesh_meta_data();
  const FileSchema & FS = arg_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_read = FS.m_io_rank ;

  const bool reader = p_read == p_rank ;

  const Part & universal_part = SM.universal_part();
        Part & owns_part      = SM.locally_owned_part();

  // const std::vector<const FilePart *> & node_parts = FS.parts( Node );
  // const std::vector<const FilePart *> & edge_parts = FS.parts( Edge );
  // const std::vector<const FilePart *> & face_parts = FS.parts( Face );
  const std::vector<const FilePart *> & elem_parts = FS.parts( Element );

  //--------------------------------------------------------------------

  std::vector< std::string > name_node_var ;
  std::vector< std::string > name_elem_var ;
  std::vector<int> exist_elem_var ; 

  for ( std::vector< const FieldBase * >::const_iterator
        i = arg_fields.begin() ; i != arg_fields.end() ; ++i ) {

    const FieldBase  & f = **i ;
    const FieldBase::Restriction & d = f.restriction( Node , universal_part );
    const unsigned f_num_dim = f.rank();

    if ( d.stride[0] ) {
      FieldIO tmp ;

      tmp.m_field     = & f ;
      tmp.m_part      = NULL ;
      tmp.m_offset    = 0 ;
      tmp.m_var_index = 0 ;

      const unsigned n = array_stride_size( f_num_dim, d.stride );

      for ( unsigned k = 0 ; k < n ; ++k ) {
        name_node_var.push_back( variable_name(Node,universal_part, f, k) );
        tmp.m_offset    = k ;
        tmp.m_var_index = name_node_var.size();
        m_field_node_universal.push_back( tmp );
      }
    }
  }

  for ( std::vector< const FieldBase * >::const_iterator
        i = arg_fields.begin() ; i != arg_fields.end() ; ++i ) {

    const FieldBase & f = **i ;

    bool is_element_var = false ;

    for ( std::vector<FieldBase::Restriction>::const_iterator
          j =  f.restrictions().begin() ;
          j != f.restrictions().end() && ! is_element_var ; ++j ) {

      is_element_var = entity_type( j->key ) == Element ;
    }
    if ( is_element_var ) {
      variable_add( Element , f , elem_parts , name_elem_var , m_field_elem );
    }
  }

  //--------------------------------------------------------------------
  // Sizes and counts

  int exo_data[ 32 ];

  int exo_error = 0 ;
  const char * exo_func = NULL ;

  if ( reader ) { // Open file:
    int comp_ws = sizeof(double) ;
    int io_ws   = 0 ;
    float version = 0 ;

    exo_func = "ex_open" ;
    m_exo_id = ex_open( arg_file_path.c_str() , EX_READ ,
                        & comp_ws , & io_ws , & version );

    if ( m_exo_id < 0 ) { exo_error = -1 ; }
  }

  // Get sizes:

  if ( reader && ! exo_error ) {
    char title[ MAX_LINE_LENGTH ];
    exo_func = "ex_get_init" ;
    exo_error = ex_get_init( m_exo_id , title ,
                             exo_data + 1 ,
                             exo_data + 2 ,
                             exo_data + 3 ,
                             exo_data + 4 ,
                             exo_data + 5 ,
                             exo_data + 6 );

    if ( 0 == exo_error ) {
      int tmp ;

      exo_func = "ne_get_n_node_num_map" ;
      exo_error = ne_get_n_node_num_map( m_exo_id, 1, 1, & tmp );

      exo_data[7] = 0 == exo_error ;

      if ( 0 < exo_error ) { exo_error = 0 ; }
    }

    if ( 0 == exo_error ) {
      int tmp ;

      exo_func = "ne_get_n_elem_num_map" ;
      exo_error = ne_get_n_elem_num_map( m_exo_id, 1, 1, & tmp );

      exo_data[8] = 0 == exo_error ;

      if ( 0 < exo_error ) { exo_error = 0 ; }
    }
  }

  exo_data[0] = exo_error ;

  broadcast( p_comm , p_read , exo_data , 9 );

  exo_error = exo_data[0] ;

  const int num_dim          = exo_data[1] ;
  const int num_nodes_global = exo_data[2] ;
  const int num_elems_global = exo_data[3] ;
  const int num_elem_blk     = exo_data[4] ;
  // const int num_node_sets    = exo_data[5] ;
  // const int num_side_sets    = exo_data[6] ;
  const bool has_node_identifiers = exo_data[7] ;
  const bool has_elem_identifiers = exo_data[8] ;

  if ( num_dim != (int) FS.m_dimension ||
       num_elem_blk != (int) elem_parts.size() ) {
    std::ostringstream msg ;
    if ( reader ) {
      msg << method << " FAILED: incompatible spatial dimension "
          << FS.m_dimension << " != "
          << num_dim
          << " OR element block count "
          << elem_parts.size() << " != "
          << num_elem_blk ;
    }
    throw std::runtime_error( msg.str() );
  }

  //--------------------------------------------------------------------
  // Read and broadcast element block sizes

  std::vector<int> elem_blk_counts( elem_parts.size() + 1 );

  if ( ! exo_error ) {

    if ( reader ) {

      for ( int i_blk = 0 ; ! exo_error &&
                            i_blk < num_elem_blk ; ++i_blk ) {
        const FilePart & part = * elem_parts[i_blk] ;
        int num_elem_this_blk ;
        int num_node_per_elem ;
        int num_attr_per_elem ;
        char elem_type[ Maximum_Name_Length ];

        exo_func = "ex_get_elem_block" ;
        exo_error = ex_get_elem_block( m_exo_id , part.m_identifier ,
                                       elem_type ,
                                       & num_elem_this_blk ,
                                       & num_node_per_elem ,
                                       & num_attr_per_elem );

        elem_blk_counts[i_blk] = num_elem_this_blk ;
      }
    }

    elem_blk_counts[ num_elem_blk ] = exo_error ;

    broadcast( p_comm, p_read, & elem_blk_counts[0], elem_blk_counts.size() );

    exo_error = elem_blk_counts[ num_elem_blk ];
  }

  //--------------------------------------------------------------------
  // Read and scatter element identifiers and relationivity.
  // Element data is partitioned by element block.
  // Determine size of element data and upper bound size of
  // the needed node map.

  // { elem-index , elem-identifier , { node-index } }

  unsigned size_elem_data_local = 0 ;
  unsigned size_needed_local_nodes = 0 ;
  
  if ( ! exo_error ) { // Sizing pass

    unsigned i = 0 ;

    for ( unsigned i_blk = 0 ; i_blk < elem_parts.size() ; ++i_blk ) {
      const FilePart & part = * elem_parts[i_blk] ;

      const unsigned num_elem_this_blk = elem_blk_counts[i_blk] ;
      const unsigned num_value = 2 + part.m_topology->node_count ;
      const unsigned num_items =
        index_processor_count( p_rank, p_size,
                               i, num_elem_this_blk, num_elems_global );
      i += num_elem_this_blk ;

      size_elem_data_local += num_value * num_items ;
      size_needed_local_nodes += part.m_topology->node_count * num_items ;
    }
  }

  std::vector<int> elem_data_local( size_elem_data_local );
  std::vector<int> needed_local_nodes( size_needed_local_nodes );

  size_elem_data_local = 0 ;
  size_needed_local_nodes = 0 ;

  if ( ! exo_error ) {

    std::vector<int> data ;
    std::vector<int> relation ;
    std::vector<int> identifiers ;
    std::vector<unsigned> send_size ;

    if ( reader ) { send_size.resize( p_size ); }

    // Load element identifiers and
    // relationivity node indices (not identifiers)

    int i = 0 ;

    for ( unsigned i_blk = 0 ; i_blk < elem_parts.size() ; ++i_blk ) {
      const FilePart & part = * elem_parts[i_blk] ;

      const unsigned num_elem_this_blk = elem_blk_counts[i_blk] ;
      const unsigned num_value = 2 + part.m_topology->node_count ;
      const unsigned index_part = 1 + i ;

      unsigned num_item_per_chunk =
        m_max_buffer / ( 2 * num_value * sizeof(int) );

      if ( num_elem_this_blk < num_item_per_chunk ) {
        num_item_per_chunk = num_elem_this_blk ;
      }

      if ( data.size() < num_item_per_chunk * num_value ) {
        data.resize( num_item_per_chunk * num_value );
      }

      if ( reader ) {
        if ( relation.size() < num_item_per_chunk * part.m_topology->node_count ) {
          relation.resize( num_item_per_chunk * part.m_topology->node_count );
        }
        if ( has_elem_identifiers && identifiers.size() < num_item_per_chunk ) {
          identifiers.resize( num_item_per_chunk );
        }
      }

      const int i_end = i + num_elem_this_blk ;

      while ( i < i_end ) {

        const int index_beg = 1 + i ;

        const int number = num_elem_this_blk < num_item_per_chunk ?
                           num_elem_this_blk : num_item_per_chunk ;

        const unsigned recv_count = num_value *
          index_processor_count( p_rank, p_size, i, number, num_elems_global );

        const unsigned recv_size = sizeof(int) * recv_count ;

        if ( reader ) {

          const unsigned i_last = i + number - 1 ;
          const unsigned p_beg = index_processor_map(p_size,i,num_elems_global);
          const unsigned p_end =
            1 + index_processor_map(p_size,i_last,num_elems_global);

          for ( unsigned p = 0 ;     p < p_beg ;  ++p ) { send_size[p] = 0 ; }
          for ( unsigned p = p_end ; p < p_size ; ++p ) { send_size[p] = 0 ; }

          for ( unsigned p = p_beg ; p < p_end ; ++p ) {
            send_size[p] = sizeof(int) * num_value *
              index_processor_count( p, p_size, i, number, num_elems_global );
          }

          // Read a 'chunk' worth of element identifiers and relationivity
          // Pack into the 'data' array for scattering.

          if ( 0 == exo_error ) {
            const int index_beg_part = 1 + index_beg - index_part ;
            exo_func = "ne_get_n_elem_conn" ;
            exo_error = ne_get_n_elem_conn( m_exo_id , part.m_identifier ,
                                            index_beg_part , number ,
                                            & relation[0] );
          }

          if ( 0 == exo_error && has_elem_identifiers ) {
            exo_func = "ne_get_n_elem_num_map" ; 
            exo_error = ne_get_n_elem_num_map( m_exo_id ,
                                               index_beg , number ,
                                               & identifiers[0] );
          }

          // Copy into scatter data array

          for ( int k = 0 ; k < number ; ++k ) {
            int * const elem_data = & data[ k * num_value ];
            const int * const conn_data = & relation[ k * part.m_topology->node_count ];
            elem_data[0] = index_beg + k ;      // Element index
            elem_data[1] = has_elem_identifiers
                           ? identifiers[k] : index_beg + k ; // Identifier
            for ( unsigned j = 0 ; j < part.m_topology->node_count ; ++j ) {
              elem_data[2+j] = conn_data[ j ];  // Element-node index
            }
          }
        }

        const unsigned * const p_send_size = reader ? & send_size[0] : NULL ;

        scatter( p_comm , p_read , p_send_size , recv_size , & data[0] );

        for ( unsigned k = 0 ; k < recv_count ; ) {
          elem_data_local[ size_elem_data_local++ ] = data[k++] ; // Index
          elem_data_local[ size_elem_data_local++ ] = data[k++] ; // Id
          for ( unsigned j = 0 ; j < part.m_topology->node_count ; ++j ) {
            const int node_index = data[k++] ;
            elem_data_local[ size_elem_data_local++ ] = node_index ;
            needed_local_nodes[ size_needed_local_nodes++ ] = node_index ;
          }
        }

        i += number ;
      }
    }

    { // Sort and unique the index node map, used to request node data.
      std::sort( needed_local_nodes.begin() , needed_local_nodes.end() );
      std::vector<int>::iterator iter =
        std::unique( needed_local_nodes.begin() , needed_local_nodes.end() );
      needed_local_nodes.erase( iter , needed_local_nodes.end() );
    }

    broadcast( p_comm, p_read, & exo_error, 1 );
  }

  //--------------------------------------------------------------------
  // Read and communicate node identifiers and coordinates as needed.

  std::vector<NodeData> node_data_local( needed_local_nodes.size() );

  if ( ! exo_error ) {

    const int num_items_per_chunk =
      m_max_buffer / ( 2 * NodeData::size_of() );

    std::vector<unsigned> send_size ;
    std::vector<int>      ident ;
    std::vector<double>   coord ;
    std::vector<NodeData> data ;

    if ( reader ) {
      send_size.resize( p_size );
      ident.resize( num_items_per_chunk );
      coord.resize( 3 * num_items_per_chunk );
    }

    unsigned count_node_data_local = 0 ;

    std::vector<int>::iterator iter_needed_node_beg =
      needed_local_nodes.begin();

    for ( int i = 0 ; i < num_nodes_global ; ) {
      const int index_beg = 1 + i ;
      int number = num_nodes_global - i ;
      if ( num_items_per_chunk < number ) { number = num_items_per_chunk ; }
      const int index_end = index_beg + number ;

      // Request node data in the span [i,i+number) from the reader.

     iter_needed_node_beg =
        std::lower_bound( iter_needed_node_beg ,
                          needed_local_nodes.end() , index_beg );

      std::vector<int>::iterator iter_needed_node_end =
        std::lower_bound( iter_needed_node_beg ,
                          needed_local_nodes.end() , index_end );

      const unsigned num_needed_node = std::distance( iter_needed_node_beg ,
                                                      iter_needed_node_end );

      const unsigned recv_size = num_needed_node * NodeData::size_of();

      CommGather node_request( p_comm, p_read, sizeof(int) * num_needed_node );

      node_request.send_buffer().pack<int>( & *iter_needed_node_beg ,
                                            num_needed_node );

      iter_needed_node_beg = iter_needed_node_end ;

      node_request.communicate();

      if ( reader ) {

        int    * const m = & ident[ 0 ] ;
        double * const x = & coord[ 0 ] ;
        double * const y = & coord[ number ] ;
        double * const z = & coord[ number * 2 ] ;

        exo_func = "ne_get_n_coord" ;
        exo_error = ne_get_n_coord( m_exo_id, index_beg, number, x, y, z );

        if ( ! exo_error && has_node_identifiers ) {
          exo_func = "ne_get_n_node_num_map" ;
          exo_error = ne_get_n_node_num_map( m_exo_id, index_beg, number, m );
        }

        unsigned total_num_send = 0 ;
        for ( unsigned p = 0 ; p < p_size ; ++p ) {
          CommBuffer & buf_request = node_request.recv_buffer(p);
          const unsigned num_send = buf_request.remaining() / sizeof(int);
          send_size[p] = num_send * NodeData::size_of();
          total_num_send += num_send ;
        }

        if ( data.size() < total_num_send ) { data.resize( total_num_send ); }

        total_num_send = 0 ;

        for ( unsigned p = 0 ; p < p_size ; ++p ) {
          CommBuffer & buf_request = node_request.recv_buffer(p);
          while ( buf_request.remaining() ) {
            int index_node_request ;
            buf_request.unpack<int>( index_node_request );
            const unsigned offset = index_node_request - index_beg ;

            data[ total_num_send ].coord[0] = x[ offset ];
            data[ total_num_send ].coord[1] = y[ offset ];
            data[ total_num_send ].coord[2] = z[ offset ];
            data[ total_num_send ].ident =
              has_node_identifiers ? m[offset] : index_node_request ;
            data[ total_num_send ].index = index_node_request ;
            ++total_num_send ;
          }
        }
      }
      else if ( data.size() < num_needed_node ) {
        data.resize( num_needed_node );
      }

      const unsigned * const p_send_size = reader ? & send_size[0] : NULL ;

      scatter( p_comm , p_read , p_send_size , recv_size , & data[0] );

      for ( unsigned k = 0 ; 0 == exo_error && k < num_needed_node ; ++k ) {
        node_data_local[ count_node_data_local ].coord[0] = data[k].coord[0] ;
        node_data_local[ count_node_data_local ].coord[1] = data[k].coord[1] ;
        node_data_local[ count_node_data_local ].coord[2] = data[k].coord[2] ;
        node_data_local[ count_node_data_local ].ident = data[k].ident ;
        node_data_local[ count_node_data_local ].index = data[k].index ;
        ++count_node_data_local ;
      }

      i += number ;
    }

    broadcast( p_comm , p_read , & exo_error , 1 );
  }

  if ( ! exo_error ) {

    PartSet entity_parts(2);
    entity_parts[0] = & owns_part ;

    // Now have all needed data to create nodes and elements
    // std::vector<int> elem_data_local ;
    // std::vector<NodeData> node_data_local ;

    size_elem_data_local = 0 ;
    unsigned i = 0 ;

    for ( unsigned i_blk = 0 ; i_blk < elem_parts.size() ; ++i_blk ) {
      const FilePart & part = * elem_parts[i_blk] ;

      entity_parts[1] = & part.m_part ;

      const unsigned num_elem_this_blk = elem_blk_counts[i_blk] ;
      const unsigned num_items =
        index_processor_count( p_rank, p_size,
                               i, num_elem_this_blk, num_elems_global );
      i += num_elem_this_blk ;

      for ( unsigned k = 0 ; k < num_items ; ++k ) {

        const int elem_index = elem_data_local[ size_elem_data_local++ ];
        const entity_id_type elem_ident =
            elem_data_local[ size_elem_data_local++ ] ;
        const entity_key_type elem_key = entity_key( Element , elem_ident );

        Entity & elem = M.declare_entity( elem_key , entity_parts );

        field_data( FS.m_field_index , elem )[0] = elem_index ;

        for ( unsigned j = 0 ; j < part.m_topology->node_count ; ++j ) {
          const int node_index = elem_data_local[ size_elem_data_local++ ];

          // Find the node data

          std::vector<NodeData>::iterator iter_node_data =
            std::lower_bound( node_data_local.begin() ,
                              node_data_local.end() ,
                              node_index ,
                              less_NodeData() );

          if ( iter_node_data == node_data_local.end() ||
               iter_node_data->index != node_index ) {
            std::ostringstream msg ;
            msg << method ;
            msg << " : FAILED to find node_index = " ;
            msg << node_index ;
            throw std::logic_error( msg.str() );
          }

          const entity_id_type node_ident = iter_node_data->ident ;
          const entity_key_type node_key = entity_key( Node , node_ident );

          Entity & node = M.declare_entity( node_key, entity_parts );

          M.declare_relation( elem , node , j );

          field_data( FS.m_field_index , node )[0] = node_index ;

          double * const node_coord =
            (double *) field_data(FS.m_field_node_coord, node);

          node_coord[0] = iter_node_data->coord[0] ;
          node_coord[1] = iter_node_data->coord[1] ;
          node_coord[2] = iter_node_data->coord[2] ;
        }
      }
    }
  }

  // Nodes and elements are created, discover sharing and generate aura

  comm_mesh_discover_sharing( M );
  comm_mesh_regenerate_aura( M );

  if ( exo_error ) {
    std::ostringstream msg ;
    if ( reader ) {
      msg << method << " FAILED calling " << exo_func << " with " << exo_error ;
    }
    throw std::runtime_error( msg.str() );
  }
}

}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#else

namespace phdmesh {
namespace exodus {

FileSchema::FileSchema(
  MetaData & arg_mesh_meta_data ,
  const FieldBase & arg_node_coordinates ,
  const FileSchema::AttributeField  & arg_elem_attributes ,
  const std::string     & ,
  ParallelMachine         ,
  const unsigned          arg_reader_rank )
  : m_schema( arg_mesh_meta_data ),
    m_io_rank( arg_reader_rank ),
    m_dimension( arg_node_coordinates.max_size( Node ) ),
    m_field_node_coord( arg_node_coordinates ),
    m_field_elem_attr(  arg_elem_attributes ),
    m_field_index( exo_index( arg_mesh_meta_data ) )
{
  verify_node_coordinate_field( arg_node_coordinates );
}

FileOutput::~FileOutput() {}

FileOutput::FileOutput(
  const FileSchema  & arg_schema ,
  const BulkData        & arg_mesh ,
  const std::string & ,
  const std::string & ,
  const bool ,
  const std::vector<const FieldBase * > & ,
  const int * const )
  : m_schema( arg_schema ),
    m_mesh( arg_mesh ),
    m_exo_id( 0 ),
    m_counter( 0 ),
    m_max_buffer( 0x0200000 )
  {}

void FileOutput::write( double ) {}

FileInput::~FileInput() { }

FileInput::FileInput(
  const FileSchema  & arg_schema ,
        BulkData        & arg_mesh ,
  const std::string & ,
  const std::vector< const FieldBase * > & )
  : m_schema( arg_schema ),
    m_mesh( arg_mesh ),
    m_exo_id( 0 ),
    m_counter( 0 ),
    m_max_buffer( 0x0200000 )
{}

}
}

#endif


