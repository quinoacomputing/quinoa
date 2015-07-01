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

#include <algorithm>
#include <util/PairIter.hpp>
#include <util/ParallelComm.hpp>
#include <util/ParallelIndex.hpp>

namespace phdmesh {

namespace {

inline
unsigned proc_map( const unsigned proc_size ,
                   const ParallelIndex::key_type key )
{
  return ( key >> 8 ) % proc_size ;
}

struct EqualKeyProc {
  bool operator()( const ParallelIndex::KeyProc & lhs ,
                   const ParallelIndex::KeyProc & rhs ) const
    { return lhs == rhs ; }
};

void sort_unique( std::vector<ParallelIndex::KeyProc> & key_proc )
{
  std::vector<ParallelIndex::KeyProc>::iterator
   i = key_proc.begin() , j = key_proc.end();

  std::sort( i , j , ParallelIndex::LessKeyProc() );
  i = std::unique( i , j , EqualKeyProc() );
  key_proc.erase( i , j );
}

typedef
  PairIter< std::vector< ParallelIndex::KeyProc >::const_iterator >
    PairIterKeyProc ;

PairIterKeyProc
span( const std::vector< ParallelIndex::KeyProc > & key_proc ,
      const ParallelIndex::key_type key )
{
  std::vector< ParallelIndex::KeyProc >::const_iterator i = key_proc.begin();
  std::vector< ParallelIndex::KeyProc >::const_iterator j = key_proc.end();

  i = std::lower_bound( i , j , key , ParallelIndex::LessKeyProc() );

  for ( j = i ; j != key_proc.end() && j->first == key ; ++j );

  return PairIterKeyProc( i , j );
}

//----------------------------------------------------------------------

void pack_map( CommAll & all ,
               const std::vector<ParallelIndex::key_type> & local )
{
  const unsigned p_size = all.parallel_size();

  std::vector<ParallelIndex::key_type>::const_iterator i ;

  for ( i = local.begin() ; i != local.end() ; ++i ) {
    const ParallelIndex::key_type value = *i ;
    const unsigned      proc = proc_map( p_size , value );
    all.send_buffer( proc ).pack<ParallelIndex::key_type>( value );
  }
}

void unpack_map( CommAll & all ,
                 std::vector< ParallelIndex::KeyProc > & key_proc )
{
  const unsigned p_size = all.parallel_size();

  unsigned count = 0 ;
  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    count += all.recv_buffer( p ).capacity() / sizeof(ParallelIndex::key_type);
  }

  key_proc.clear();
  key_proc.reserve( count );

  ParallelIndex::KeyProc value ;

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    value.second = p ;
    while ( buf.remaining() ) {
      buf.unpack<ParallelIndex::key_type>( value.first );
      key_proc.push_back( value );
    }
  }
}

//----------------------------------------------------------------------

void pack_query( CommAll & all ,
                 const std::vector< ParallelIndex::KeyProc > & key_proc )
{
  ParallelIndex::key_type value[2] ;

  std::vector< ParallelIndex::KeyProc >::const_iterator i ;

  for ( i = key_proc.begin() ; i != key_proc.end() ; ) {
    value[0] = i->first ;

    const std::vector< ParallelIndex::KeyProc >::const_iterator i_beg = i ;

    for ( ; i != key_proc.end() && value[0] == i->first ; ++i );

    const std::vector< ParallelIndex::KeyProc >::const_iterator i_end = i ;
    
    for ( i = i_beg ; i != i_end ; ++i ) {
      CommBuffer & buf = all.send_buffer( i->second );

      std::vector< ParallelIndex::KeyProc >::const_iterator j ;

      for ( j = i_beg ; j != i_end ; ++j ) {
        if ( j != i ) {
          value[1] = j->second ;
          buf.pack<ParallelIndex::key_type>( value , 2 );
        }
      }
    }
  }
}

void unpack_query( CommAll & all ,
                   std::vector< ParallelIndex::KeyProc > & key_proc )
{
  const unsigned p_size = all.parallel_size();

  ParallelIndex::KeyProc entry ;
  ParallelIndex::key_type value[2] ;

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    while ( buf.remaining() ) {
      buf.unpack<ParallelIndex::key_type>( value , 2 );
      entry.first = value[0] ;
      entry.second = value[1] ;
      key_proc.push_back( entry );
    }
  }
}

//----------------------------------------------------------------------

void pack_query( CommAll & all ,
                 const std::vector< ParallelIndex::KeyProc > & key_proc_map ,
                 const std::vector< ParallelIndex::KeyProc > & query )
{
  ParallelIndex::key_type value[2] ;

  std::vector< ParallelIndex::KeyProc >::const_iterator i ;

  for ( i = query.begin() ; i != query.end() ; ) {
    value[0] = i->first ;

    const PairIterKeyProc s = span( key_proc_map , value[0] );

    for ( ; i != query.end() && value[0] == i->first ; ++i ) {
      CommBuffer & buf = all.send_buffer( i->second );

      std::vector< ParallelIndex::KeyProc >::const_iterator j ;

      for ( j = s.begin() ; j != s.end() ; ++j ) {
        value[1] = j->second ;
        buf.pack<ParallelIndex::key_type>( value , 2 );
      }
    }
  }
}

//----------------------------------------------------------------------

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

ParallelIndex::~ParallelIndex() {}

ParallelIndex::ParallelIndex(
  ParallelMachine arg_comm ,
  const std::vector<ParallelIndex::key_type> & arg_local )
  : m_comm( arg_comm ), m_key_proc()
{
  const unsigned p_size = parallel_machine_size( arg_comm );

  CommAll all( arg_comm );

  pack_map( all , arg_local );

  all.allocate_buffers( p_size / 4 , false );
  
  pack_map( all , arg_local );

  all.communicate();

  unpack_map( all , m_key_proc );

  sort_unique( m_key_proc );
}

void ParallelIndex::query(
  std::vector<ParallelIndex::KeyProc> & arg_global ) const
{
  const unsigned p_size = parallel_machine_size( m_comm );

  CommAll all( m_comm );

  pack_query( all , m_key_proc );

  all.allocate_buffers( p_size / 4 , false );

  pack_query( all , m_key_proc );

  all.communicate();

  unpack_query( all , arg_global );

  sort_unique( arg_global );
}

void ParallelIndex::query(
  const std::vector<ParallelIndex::key_type> & arg_local ,
  std::vector<ParallelIndex::KeyProc> & arg_global ) const
{
  const unsigned p_size = parallel_machine_size( m_comm );

  std::vector<ParallelIndex::KeyProc> tmp ;

  {
    CommAll all( m_comm );

    pack_map( all , arg_local );

    all.allocate_buffers( p_size / 4 , false );

    pack_map( all , arg_local );

    all.communicate();

    unpack_map( all , tmp ); // { ( key , querying_processor ) }

    sort_unique( tmp );
  }

  {
    CommAll all( m_comm );

    pack_query( all , m_key_proc , tmp );

    all.allocate_buffers( p_size / 4 , false );

    pack_query( all , m_key_proc , tmp );

    all.communicate();

    unpack_query( all , arg_global );

    sort_unique( arg_global );
  }
}

}

