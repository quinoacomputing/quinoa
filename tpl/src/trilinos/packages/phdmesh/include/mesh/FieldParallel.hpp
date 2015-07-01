/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2008) Sandia Corporation                     */
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

#ifndef phdmesh_FieldParallel_hpp
#define phdmesh_FieldParallel_hpp

//----------------------------------------------------------------------

#include <util/SimpleArrayOps.hpp>
#include <util/Parallel.hpp>
#include <util/ParallelComm.hpp>

#include <mesh/Types.hpp>
#include <mesh/Field.hpp>
#include <mesh/Entity.hpp>
#include <mesh/BulkData.hpp>

namespace phdmesh {

/** Communicate field data from domain to range.
 *  The fields array must be identical on all processors.
 *  All fields and mesh entities must belong to the same mesh.
 *  If symmetric ( & domain == & range) then from owned to not owned.
 */
bool communicate_field_data(
  const BulkData & mesh ,
  const std::vector<EntityProc> & domain ,
  const std::vector<EntityProc> & range ,
  const std::vector< const FieldBase *> & fields ,
  bool local_flag = false );

/** Communicate field data given symmetric communication plan.  */
void communicate_field_data(
  ParallelMachine ,
  const std::vector<EntityProc> & ,
  const unsigned field_count ,
  const FieldBase * fields[] ,
  CommAll & sparse );

//----------------------------------------------------------------------

namespace {

//----------------------------------------------------------------------
/** Parallel reduction of shared entities' field data.
 *  In the anonymous namespace to avoid redundant link symbols.
 *
 *  example usage:
 *    parallel_reduce( mesh , sum( field ) );
 *
 *  where the operations are: sum, max, min
 */
template< class OpField >
void parallel_reduce( const BulkData & mesh ,
                      const OpField      & op )
{
  const FieldBase * fields[1] = { & op.field };

  const std::vector<EntityProc> & shared = mesh.shared_entities();

  CommAll sparse ;

  communicate_field_data( mesh.parallel(), shared, 1, fields, sparse );

  op( shared , sparse );
}

/** Parallel reduction of shared entities' field data.
 *  In the anonymous namespace to avoid redundant link symbols.
 *
 *  example usage:
 *    parallel_reduce( mesh , sum( fieldA ) , max( fieldB ) );
 */
template< class OpField1 , class OpField2 >
void parallel_reduce( const BulkData & mesh ,
                      const OpField1     & op1 ,
                      const OpField2     & op2 )
{
  const FieldBase * fields[2] = { & op1.field , & op2.field };

  const std::vector<EntityProc> & shared = mesh.shared_entities();

  CommAll sparse ;

  communicate_field_data( mesh.parallel(), shared, 2, fields, sparse );

  op1( shared , sparse );
  op2( shared , sparse );
}

//----------------------------------------------------------------------

template< class ReduceOp ,
          class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
struct ParallelReduceField {
  typedef Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> field_type ;

  const field_type & field ;

  ParallelReduceField( const field_type & f ) : field(f) {}
  ParallelReduceField( const ParallelReduceField & p ) : field(p.field) {}

  void operator()( const std::vector<EntityProc> & shared ,
                   CommAll & sparse ) const ;

private:
  ParallelReduceField & operator = ( const ParallelReduceField & );
};

template< class ReduceOp ,
          class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
void ParallelReduceField< ReduceOp , Type ,  Tag1,  Tag2,  Tag3 ,
                                     Tag4 ,  Tag5,  Tag6,  Tag7 >::
  operator()( const std::vector<EntityProc> & shared , CommAll & sparse ) const
{
  typedef EntityArray< field_type > array_type ;

  for ( std::vector<EntityProc>::const_iterator
        i = shared.begin(); i != shared.end() ; ++i ) {

    array_type array( field , * i->first );
    Type * ptr           = array.contiguous_data();
    Type * const ptr_end = ptr + array.size();

    CommBuffer & b = sparse.recv_buffer( i->second );

    for ( ; ptr < ptr_end ; ++ptr ) {
      Type tmp ;
      b.template unpack<unsigned char>( (unsigned char *)(&tmp), sizeof(Type) );
      ReduceOp::op( ptr , & tmp );
    }
  }
}

}

//----------------------------------------------------------------------

template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<impl::Sum<Type,1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
sum( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f )
{
  return ParallelReduceField<impl::Sum<Type,1>,
                             Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f );
}

template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<impl::Max<Type,1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
max( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f )
{
  return ParallelReduceField<impl::Max<Type,1>,
                             Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f );
}

template< class Type , class Tag1, class Tag2, class Tag3 ,
          class Tag4 , class Tag5, class Tag6, class Tag7 >
ParallelReduceField<impl::Min<Type,1>,Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
inline
min( const Field<Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & f )
{
  return ParallelReduceField<impl::Min<Type,1>,
                             Type,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>( f );
}

}

#endif

