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

#include <stdexcept>
#include <sstream>
#include <algorithm>

#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>

#include <mesh/MetaData.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/FieldParallel.hpp>


namespace phdmesh {

// Heterogeneity?

bool communicate_field_data(
  const BulkData & mesh ,
  const std::vector<EntityProc> & domain ,
  const std::vector<EntityProc> & range ,
  const std::vector<const FieldBase *> & fields ,
  bool local_flag )
{
  if ( fields.empty() ) { return local_flag ; }

  const unsigned parallel_size = mesh.parallel_size();
  const unsigned parallel_rank = mesh.parallel_rank();
  const bool     asymmetric    = & domain != & range ;

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  std::vector<EntityProc>::const_iterator i ;

  for ( i = domain.begin() ; i != domain.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || parallel_rank == e.owner_rank() ) {
      unsigned e_size = 0 ;
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        e_size += field_data_size( f , e );
      } 
      send_size[ p ] += e_size ;
    }
  }

  for ( i = range.begin() ; i != range.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || p == e.owner_rank() ) {
      unsigned e_size = 0 ;
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        e_size += field_data_size( f , e );
      } 
      recv_size[ p ] += e_size ;
    }
  }

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const s_size = & send_size[0] ;
    const unsigned * const r_size = & recv_size[0] ;
    sparse.allocate_buffers( mesh.parallel(), parallel_size / 4 , s_size, r_size);
  }

  // Pack for send:

  for ( i = domain.begin() ; i != domain.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || parallel_rank == e.owner_rank() ) {
      CommBuffer & b = sparse.send_buffer( p );
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        const unsigned size = field_data_size( f , e );
        if ( size ) {
          unsigned char * ptr = (unsigned char *) field_data( f , e );
          b.pack<unsigned char>( ptr , size );
        }
      }
    }
  }

  // Communicate:

  sparse.communicate();

  // Unpack for recv:

  for ( i = range.begin() ; i != range.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || p == e.owner_rank() ) {
      CommBuffer & b = sparse.recv_buffer( p );
      for ( fi = fb ; fi != fe ; ++fi ) {
        const FieldBase & f = **fi ;
        const unsigned size = field_data_size( f , e );
        if ( size ) {
          unsigned char * ptr = (unsigned char *) field_data( f , e );
          b.unpack<unsigned char>( ptr , size );
        }
      }
    }
  }

  return local_flag ;
}

//----------------------------------------------------------------------

void communicate_field_data(
  ParallelMachine machine ,
  const std::vector<EntityProc> & shared ,
  const unsigned field_count ,
  const FieldBase * fields[] ,
  CommAll & sparse )
{
  const unsigned parallel_size = parallel_machine_size( machine );

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> msg_size( parallel_size , zero );

  size_t j ;
  std::vector<EntityProc>::const_iterator i ;

  for ( j = 0 ; j < field_count ; ++j ) {
    const FieldBase & f = * fields[j] ;
    for ( i = shared.begin() ; i != shared.end() ; ++i ) {
      Entity       & e = * i->first ;
      const unsigned p = i->second ;
      msg_size[ p ] += field_data_size( f , e );
    }
  }

  // Allocate send and receive buffers:

  {
    const unsigned * const s_size = & msg_size[0] ;
    sparse.allocate_buffers( machine, parallel_size / 4 , s_size, s_size);
  }

  // Pack for send:
 
  for ( j = 0 ; j < field_count ; ++j ) {
    const FieldBase & f = * fields[j] ;
    for ( i = shared.begin() ; i != shared.end() ; ++i ) {
      Entity       & e = * i->first ;
      const unsigned p = i->second ;

      CommBuffer & b = sparse.send_buffer( p );

      const unsigned size = field_data_size( f , e );

      if ( size ) {
        unsigned char * ptr = (unsigned char *) field_data( f , e );
        b.pack<unsigned char>( ptr , size );
      }
    }
  }

  // Communicate:

  sparse.communicate();
}

//----------------------------------------------------------------------

}

