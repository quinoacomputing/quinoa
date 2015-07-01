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

#include <stdlib.h>
#include <memory.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <mesh/Kernel.hpp>
#include <mesh/Entity.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/FieldData.hpp>

namespace phdmesh {

namespace {

void memory_copy( unsigned char * dst , unsigned char * src , unsigned n )
{ memcpy( dst , src , n ); }

void memory_zero( unsigned char * dst , unsigned n )
{ memset( dst , 0 , n ); }

}

//----------------------------------------------------------------------
// KernelKey key = ( part-count , { part-ordinals } , counter )
//  key[ key[0] ] == counter

namespace {

inline
unsigned kernel_counter( const unsigned * const key )
{ return key[ *key ]; }

// The part count and parts are equal
bool kernel_part_equal( const unsigned * lhs , const unsigned * rhs )
{
  bool result = true ;
  {
    const unsigned * const end_lhs = lhs + *lhs ;
    while ( result && end_lhs != lhs ) {
      result = *lhs == *rhs ;
      ++lhs ; ++rhs ;
    }
  }
  return result ;
}

}

//----------------------------------------------------------------------

bool Kernel::member( const Part & part ) const
{
  const unsigned * const key_base = key();
  const unsigned * const i_beg = key_base + 1 ;
  const unsigned * const i_end = key_base + key_base[0] ;

  const unsigned ord = part.mesh_meta_data_ordinal();
  const unsigned * const i = std::lower_bound( i_beg , i_end , ord );

  return i_end != i && ord == *i ;
}

bool Kernel::member_all( const std::vector<Part*> & parts ) const
{
  const unsigned * const key_base = key();
  const unsigned * const i_beg = key_base + 1 ;
  const unsigned * const i_end = key_base + key_base[0] ;

  const std::vector<Part*>::const_iterator ip_end = parts.end();
        std::vector<Part*>::const_iterator ip     = parts.begin() ; 

  bool result_all = true ;

  for ( ; result_all && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_all = i_end != i && ord == *i ;
  }
  return result_all ;
}

bool Kernel::member_any( const std::vector<Part*> & parts ) const
{
  const unsigned * const key_base = key();
  const unsigned * const i_beg = key_base + 1 ;
  const unsigned * const i_end = key_base + key_base[0] ;

  const std::vector<Part*>::const_iterator ip_end = parts.end();
        std::vector<Part*>::const_iterator ip     = parts.begin() ; 

  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_none = i_end == i || ord != *i ;
  }
  return ! result_none ;
}

//----------------------------------------------------------------------
// The part count and part ordinals are less
bool KernelLess::operator()( const unsigned * lhs ,
                             const unsigned * rhs ) const
{
  const unsigned * const last_lhs = lhs + ( *lhs < *rhs ? *lhs : *rhs );
  while ( last_lhs != lhs && *lhs == *rhs ) { ++lhs ; ++rhs ; }
  return *lhs < *rhs ;
}

bool Kernel::has_superset( const Part & p ) const
{
  const unsigned ordinal = p.mesh_meta_data_ordinal();

  std::pair<const unsigned *, const unsigned *> 
    part_ord = superset_part_ordinals();

  part_ord.first =
    std::lower_bound( part_ord.first , part_ord.second , ordinal );

  return part_ord.first < part_ord.second && ordinal == *part_ord.first ;
}

bool Kernel::has_superset( const PartSet & ps ) const
{
  std::pair<const unsigned *, const unsigned *> 
    part_ord = superset_part_ordinals();

  bool result = ! ps.empty();

  for ( PartSet::const_iterator
        i = ps.begin() ; result && i != ps.end() ; ++i ) {

    const unsigned ordinal = (*i)->mesh_meta_data_ordinal();

    part_ord.first =
      std::lower_bound( part_ord.first , part_ord.second , ordinal );

    result = part_ord.first < part_ord.second && ordinal == *part_ord.first ;
  }
  return result ;
}

void Kernel::supersets( PartSet & ps ) const
{
  const MetaData & mesh_meta_data = m_mesh.mesh_meta_data();

  std::pair<const unsigned *, const unsigned *> 
    part_ord = superset_part_ordinals();

  ps.resize( part_ord.second - part_ord.first );

  for ( unsigned i = 0 ;
        part_ord.first < part_ord.second ; ++(part_ord.first) , ++i ) {
    ps[i] = & mesh_meta_data.get_part( * part_ord.first );
  }
}

//----------------------------------------------------------------------

bool field_data_valid( const FieldBase & f ,
                       const Kernel & k ,
                       unsigned ord ,
                       const char * required_by )
{
  const MetaData * const k_mesh_meta_data = & k.mesh().mesh_meta_data();
  const MetaData * const f_mesh_meta_data = & f.mesh_meta_data();
  const bool ok_mesh_meta_data  = k_mesh_meta_data == f_mesh_meta_data ;
  const bool ok_ord     = ord < k.size() ;
  const bool exists     = ok_mesh_meta_data && ok_ord && NULL != field_data( f , k );

  if ( required_by && ! exists ) {
    std::ostringstream msg ;
    msg << "phdmesh::field_data_valid( " ;
    msg << f ;
    msg << " , " ;
    msg << k ;
    msg << " , " ;
    msg << ord ;
    msg << " , " ;
    msg << required_by ;
    msg << " ) FAILED with " ;
    if ( ! ok_mesh_meta_data ) {
      msg << " different MetaData" ;
    }
    else if ( ! ok_ord ) {
      msg << " Ordinal " ;
      msg << ord ;
      msg << " >= " ;
      msg << " size " ;
      msg << k.size();
    }
    else {
      msg << " no data" ;
    }
    throw std::runtime_error( msg.str() );
  }

  return exists ;
}

//----------------------------------------------------------------------

Kernel::Kernel( BulkData & arg_mesh ,
                EntityType arg_type ,
                const unsigned * arg_key )
: SetvMember<const unsigned * const>( arg_key ),
  m_mesh( arg_mesh ) ,
  m_entity_type( arg_type ),
  m_size( 0 ),
  m_capacity( 0 ),
  m_alloc_size( 0 ),
  m_kernel(),
  m_field_map( NULL ),
  m_entities( NULL )
{}


//----------------------------------------------------------------------

void Kernel::zero_fields( Kernel & k_dst , unsigned i_dst )
{
  const std::vector<FieldBase*> & field_set =
    k_dst.mesh().mesh_meta_data().get_fields();

  unsigned char * const p = reinterpret_cast<unsigned char*>(k_dst.m_entities);
  const DataMap *       i = k_dst.m_field_map ;
  const DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i ) {
    if ( i->m_size ) {
      memory_zero( p + i->m_base + i->m_size * i_dst , i->m_size );
    }
  }
}

void Kernel::copy_fields( Kernel & k_dst , unsigned i_dst ,
                          Kernel & k_src , unsigned i_src )
{
  static const char method[] = "phdmesh::Kernel::copy_fields" ;

  const std::vector<FieldBase*> & field_set =
    k_dst.mesh().mesh_meta_data().get_fields();

  unsigned char * const s = reinterpret_cast<unsigned char*>(k_src.m_entities);
  unsigned char * const d = reinterpret_cast<unsigned char*>(k_dst.m_entities);
  DataMap *       j = k_src.m_field_map ;
  DataMap *       i = k_dst.m_field_map ;
  DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i , ++j ) {

    if ( i->m_size ) {
      if ( j->m_size ) {
        if ( i->m_size == j->m_size ) {
          memory_copy( d + i->m_base + i->m_size * i_dst ,
                       s + j->m_base + j->m_size * i_src , i->m_size );
        }
        else {
          std::ostringstream msg ;
          msg << method ;
          msg << " FAILED WITH INCOMPATIBLE FIELD SIZES" ;
          throw std::runtime_error( msg.str() );
        }
      }
      else {
        memory_zero( d + i->m_base + i->m_size * i_dst , i->m_size );
      }
    }
  }
}

//----------------------------------------------------------------------

namespace {

inline unsigned align( unsigned nb )
{
  enum { BYTE_ALIGN = 16 };
  const unsigned gap = nb % BYTE_ALIGN ;
  if ( gap ) { nb += BYTE_ALIGN - gap ; }
  return nb ;
}

struct DimLess {
  bool operator()( const FieldBase::Restriction & lhs ,
                   const unsigned & rhs ) const
    { return lhs < rhs ; }
};

const FieldBase::Restriction & dimension( const FieldBase & field ,
                                  EntityType etype ,
                                  const unsigned num_part_ord ,
                                  const unsigned part_ord[] ,
                                  const char * const method )
{
  static const FieldBase::Restriction empty ;

  const FieldBase::Restriction * dim = & empty ;

  const std::vector<FieldBase::Restriction> & dim_map = field.restrictions();
  const std::vector<FieldBase::Restriction>::const_iterator iend = dim_map.end();
        std::vector<FieldBase::Restriction>::const_iterator ibeg = dim_map.begin();

  for ( unsigned i = 0 ; i < num_part_ord && iend != ibeg ; ++i ) {

    const unsigned key = FieldBase::Restriction::key_value( etype , part_ord[i] );

    ibeg = std::lower_bound( ibeg , iend , key , DimLess() );

    if ( iend != ibeg && ibeg->key == key ) {
      if ( dim == & empty ) { dim = & *ibeg ; }

      if ( NotEqual< MaximumFieldDimension >( ibeg->stride , dim->stride ) ) {

        Part & p_old = field.mesh_meta_data().get_part( ibeg->ordinal() );
        Part & p_new = field.mesh_meta_data().get_part( dim->ordinal() );

        std::ostringstream msg ;
        msg << method ;
        msg << " FAILED WITH INCOMPATIBLE DIMENSIONS FOR " ;
        msg << field ;
        msg << " Part[" << p_old.name() ;
        msg << "] and Part[" << p_new.name() ;
        msg << "]" ;
     
        throw std::runtime_error( msg.str() );
      }
    }
  }

  return *dim ;
}

}

//----------------------------------------------------------------------

void Kernel::update_state()
{
  if ( 0 == kernel_counter( key() ) ) {

    const MetaData & S = m_mesh.mesh_meta_data();
    const std::vector<FieldBase*> & field_set = S.get_fields();

    for ( unsigned i = 0 ; i < field_set.size() ; ) {

      DataMap * const tmp = m_field_map + i ;
      const FieldBase & field = * field_set[i] ;
      const unsigned num_state = field.number_of_states();
      i += num_state ;

      if ( 1 < num_state && tmp->m_size ) {
        unsigned offset[ MaximumFieldStates ] ;

        for ( unsigned j = 0 ; j < num_state ; ++j ) {
          offset[j] = tmp[j].m_base ;
        }

        for ( unsigned j = 0 ; j < num_state ; ++j ) {
          const unsigned j_new = ( j + num_state - 1 ) % num_state ;
          tmp[j_new].m_base = offset[j] ;
        }
      }
    }
  }
}

//----------------------------------------------------------------------

namespace {

void * local_malloc( size_t n )
{
  void * const ptr = malloc( n );

  if ( NULL == ptr ) {
    std::ostringstream msg ;
    msg << "phdmesh::Kernel FAILED malloc( " << n << " )" ;
    throw std::runtime_error( msg.str() );
  }

  return ptr ;
}

}

Kernel::~Kernel()
{}

void BulkData::destroy_kernel( KernelSet::iterator ik )
{
  Kernel * const k = & * ik ;

  KernelSet & kernel_set = m_kernels[ k->entity_type() ];

  kernel_set.remove( *ik ); // 'ik' is now invalidated

  if ( 0 == kernel_counter( k->key() ) ) {
    free( (void*) k->m_field_map );
  }
  k->m_field_map = NULL ;

  k->~Kernel();

  free( (void*) k );
}

//----------------------------------------------------------------------
// The input part ordinals are complete and contain all supersets.

KernelSet::iterator
BulkData::declare_kernel( const EntityType arg_entity_type ,
                              const unsigned part_size ,
                              const unsigned part_ord[] )
{
  enum { KEY_TMP_BUFFER_SIZE = 32 };

  static const char method[] = "phdmesh::Kernel" ;

  const unsigned max = ~((unsigned) 0);

  KernelSet & kernel_set = m_kernels[ arg_entity_type ];

  const std::vector< FieldBase * > & field_set = m_mesh_meta_data.get_fields();

  const unsigned num_fields = field_set.size();

  //----------------------------------
  // For performance try not to allocate a temporary.

  unsigned key_tmp_buffer[ KEY_TMP_BUFFER_SIZE ];

  std::vector<unsigned> key_tmp_vector ;

  const unsigned key_size = 2 + part_size ;

  unsigned * const key =
    ( key_size <= KEY_TMP_BUFFER_SIZE )
    ? key_tmp_buffer
    : ( key_tmp_vector.resize( key_size ) , & key_tmp_vector[0] );

  //----------------------------------
  // Key layout:
  // { part_size + 1 , { part_ordinals } , family_count }
  // Thus family_count = key[ key[0] ]
  //
  // for upper bound search use the maximum key.

  key[ key[0] = part_size + 1 ] = max ;

  {
    unsigned * const k = key + 1 ;
    for ( unsigned i = 0 ; i < part_size ; ++i ) { k[i] = part_ord[i] ; }
  }

  //----------------------------------
  // Look for the last kernel in this family:

  const unsigned * const key_const = key ;

  KernelSet::iterator ik = kernel_set.upper_bound( key_const );

  Kernel * const last_kernel =
    ( ik != kernel_set.begin() ) && kernel_part_equal( (--ik)->key() , key )
    ? ( & * ik ) : NULL ;

  Kernel          * kernel    = NULL ;
  Kernel::DataMap * field_map = NULL ;

  if ( last_kernel == NULL ) { // First kernel in this family
    key[ key[0] ] = 0 ; // Set the key's family count to zero
  }
  else { // Last kernel present, can it hold one more entity?

    field_map = last_kernel->m_field_map ;

    const unsigned last_count = last_kernel->key()[ key[0] ];

    const unsigned cap = last_kernel->capacity();

    if ( last_kernel->size() < cap ) {
      kernel = last_kernel ;
    }
    else if ( last_count < max ) {
      key[ key[0] ] = 1 + last_count ; // Increment the key's family count.
    }
    else {
      // ERROR insane number of kernels!
      std::string msg ;
      msg.append( method );
      msg.append( " FAILED due to impossibly large number of kernels" );
      throw std::logic_error( msg );
    }
  }

  //----------------------------------
  // Family's field map does not exist, create it:

  if ( NULL == field_map ) {

    field_map = reinterpret_cast<Kernel::DataMap*>(
                local_malloc( sizeof(Kernel::DataMap) * ( num_fields + 1 )));

    unsigned value_offset = 0 ;

    value_offset += align( sizeof(Entity*) * m_kernel_capacity );

    for ( unsigned i = 0 ; i < num_fields ; ++i ) {
      const FieldBase  & field = * field_set[i] ;

      unsigned value_size = 0 ;

      const FieldBase::Restriction & dim =
        dimension( field, arg_entity_type, key[0] - 1, key + 1, method);

      if ( dim.stride[0] ) { // Exists

        const unsigned scalar_size =
          NumericEnum<void>::size( field.numeric_type_ordinal() );

        const unsigned field_num_dim = field.rank();

        value_size = scalar_size *
          ( field_num_dim ? dim.stride[ field_num_dim - 1 ] : 1 );
      }

      field_map[i].m_base = value_offset ;
      field_map[i].m_size = value_size ;
      field_map[i].m_stride = dim.stride ;

      value_offset += align( value_size * m_kernel_capacity );
    }
    field_map[ num_fields ].m_base  = value_offset ;
    field_map[ num_fields ].m_size = 0 ;
    field_map[ num_fields ].m_stride = NULL ;
  }

  //----------------------------------

  if ( NULL == kernel ) {

    // Required kernel does not exist, must allocate and insert
    //
    // Allocation size:
    //   sizeof(Kernel) +
    //   key_size * sizeof(unsigned) +
    //   sizeof(Entity*) * capacity() +
    //   sum[number_of_fields]( fieldsize * capacity )

    const unsigned size = align( sizeof(Kernel) ) +
                          align( sizeof(unsigned) * key_size ) +
                          field_map[ num_fields ].m_base ;

    // All fields checked and sized, Ready to allocate

    unsigned char * ptr = (unsigned char *) local_malloc( size );

    kernel = (Kernel *) ptr ; ptr += align( sizeof(Kernel) );

    {
      unsigned * new_key = (unsigned *) ptr ;

      ptr += align( sizeof(unsigned) * key_size );

      for ( unsigned i = 0 ; i < key_size ; ++i ) { new_key[i] = key[i] ; }

      new(kernel) Kernel( *this , arg_entity_type , new_key );
    }

    kernel->m_size       = 0 ;
    kernel->m_capacity   = m_kernel_capacity ;
    kernel->m_alloc_size = size ;
    kernel->m_field_map  = field_map ;
    kernel->m_entities   = (Entity **) ptr ;

    std::pair<KernelSet::iterator,bool> result = kernel_set.insert( kernel );

    if ( ! result.second ) {
      std::string msg ;
      msg.append( method );
      msg.append( " FAILED INSERTION" );
      throw std::logic_error( msg );
    }

    ik = result.first ;

    KernelSet::iterator
      first_kernel = last_kernel ? last_kernel->m_kernel : ik ;

    first_kernel->m_kernel = ik ; // New last kernel

    ik->m_kernel = first_kernel ;
  }

  //----------------------------------

  return ik ;
}

//----------------------------------------------------------------------

void BulkData::remove_entity( KernelSet::iterator ik , unsigned i )
{
  const KernelSet::iterator first =
    kernel_counter( ik->key() ) ? ik->m_kernel : ik ;

  const KernelSet::iterator last = first->m_kernel ;

  // Only move if not the last entity being removed

  if ( last != ik || ik->m_size != i + 1 ) {

    // Not the same kernel or not the last entity

    // Copy last entity in last to ik slot i

    Entity * const entity = last->m_entities[ last->m_size - 1 ];

    Kernel::copy_fields( *ik , i , *last , last->m_size - 1 );

    ik->m_entities[i]    = entity ;
    entity->m_kernel     = ik ;
    entity->m_kernel_ord = i ;

    // Entity field data has relocated

    internal_propagate_relocation( *entity );
  }

  --( last->m_size );

  if ( last->m_size != 0 ) {
    last->m_entities[ last->m_size ] = NULL ;
  }
  else {

    if ( first != last ) { --( first->m_kernel ); }

    destroy_kernel( last );
  }
}

//----------------------------------------------------------------------

namespace {

struct LessEntityPointer {
  bool operator()( const Entity * const lhs , const Entity * const rhs ) const
    { return lhs->key() < rhs->key() ; }
};

}

bool BulkData::internal_sort_kernel_entities()
{
  bool change = false ;

  for ( unsigned type = 0 ; type < EntityTypeEnd ; ++type ) {
    KernelSet & ks = m_kernels[type] ;

    // bk == first kernel in the family
    KernelSet::iterator bk = ks.begin();

    while ( bk != ks.end() ) {

      KernelSet::iterator ik_vacant = bk->m_kernel ; // Last kernel
      unsigned            ie_vacant = ik_vacant->size();

      if ( ik_vacant->capacity() <= ie_vacant ) {
        // Have to create a kernel just for the scratch space...
        const unsigned * const key = bk->key();
        ik_vacant = declare_kernel( (EntityType) type, key[0]-1 , key+1 );
        ie_vacant = 0 ;
      }

      ik_vacant->m_entities[ ie_vacant ] = NULL ;

      // ek == end kernel for the family
      KernelSet::iterator ek = bk->m_kernel ; ++ek ;

      unsigned count = 0 ;
      for ( KernelSet::iterator ik = bk ; ik != ek ; ++ik ) {
        count += ik->size();
      }

      std::vector<Entity*> entities( count );

      std::vector<Entity*>::iterator j = entities.begin();

      for ( KernelSet::iterator ik = bk ; ik != ek ; ++ik ) {
        const unsigned n = ik->size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          *j = ik->m_entities[i] ;
        }
      }

      std::sort( entities.begin() , entities.end() , LessEntityPointer() );

      j = entities.begin();

      bool change_this_family = false ;

      for ( KernelSet::iterator ik = bk ; ik != ek ; ++ik ) {
        const unsigned n = ik->size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          Entity * const current = ik->m_entities[i] ;

          if ( current != *j ) {

            if ( current ) {
              // Move current entity to the vacant spot
              Kernel::copy_fields( *ik_vacant , ie_vacant , *ik, i );
              current->m_kernel     = ik_vacant ;
              current->m_kernel_ord = ie_vacant ;
              ik_vacant->m_entities[ ie_vacant ] = current ;
            }

            // Set the vacant spot to where the required entity is now.
            ik_vacant = (*j)->m_kernel ;
            ie_vacant = (*j)->m_kernel_ord ;
            ik_vacant->m_entities[ ie_vacant ] = NULL ;

            // Move required entity to the required spot
            Kernel::copy_fields( *ik, i, *ik_vacant , ie_vacant );
            (*j)->m_kernel     = ik ;
            (*j)->m_kernel_ord = i ;
            ik->m_entities[i] = *j ;

            change_this_family = true ;
          }

          // Once a change has occured then need to propagate the
          // relocation for the remainder of the family.
          // This allows the propagation to be performed once per
          // entity as opposed to both times the entity is moved.

          if ( change_this_family ) { internal_propagate_relocation( **j ); }
        }
      }

      if ( change_this_family ) { change = true ; }

      bk = ek ;
    }
  }

  return change ;
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const Kernel & k )
{
  const char * const entity_name = entity_type_name( k.entity_type() );

  PartSet parts ; k.supersets( parts );

  s << "Kernel( " << entity_name << " : " ;
  for ( PartSet::iterator i = parts.begin() ; i != parts.end() ; ++i ) {
    s << (*i)->name() << " " ;
  }
  s << ")" ;

  return s ;
}


std::ostream &
Kernel::print( std::ostream & os , const std::string & lead ) const
{
  const MetaData & mesh_meta_data = m_mesh.mesh_meta_data();
  const char * const entity_name = entity_type_name( m_entity_type );

  const unsigned * k = key();
  const unsigned n = *k ; ++k ;

  os << lead
     << "Kernel(" << entity_name << " : " ;
  for ( unsigned i = 1 ; i < n ; ++i , ++k ) {
    const std::string & name = mesh_meta_data.get_part( *k ).name(); os << " " << name ;
  }

  os << " " << *k << " ){ "
     << m_size << " of " << m_capacity << " }"
     << std::endl << lead
     << "  members {" ;

  for ( unsigned j = 0 ; j < m_size ; ++j ) {
    const long id = m_entities[j]->identifier(); os << " " << id ;
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

} // namespace phdmesh

