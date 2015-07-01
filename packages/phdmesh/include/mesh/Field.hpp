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

#ifndef phdmesh_Field_hpp
#define phdmesh_Field_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <string>
#include <vector>

#include <util/NumericEnum.hpp>
#include <util/SimpleArrayOps.hpp>
#include <util/CSet.hpp>

#include <mesh/Types.hpp>
#include <mesh/FieldTraits.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

/** Print the text name for a field, depending on the number of states.  */
std::ostream & operator << ( std::ostream & , const FieldBase & );

/** Print field and field dimension map entries on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const , const FieldBase & );

//----------------------------------------------------------------------
/** Field value states.
 *    Default state is 'None' or 'New' or 'N+1'
 */
enum FieldState {
  /* Exactly one state  */    StateNone = 0 ,
  /* Exactly two states */    StateNew  = 0 ,
                              StateNP1  = 0  /* N+1 */ ,
                              StateOld  = 1 ,
  /* Three or more states */  StateN    = 1 ,
                              StateNM1  = 2  /* N-1 */ ,
                              StateNM2  = 3  /* N-2 */ ,
                              StateNM3  = 4  /* N-3 */ ,
                              StateNM4  = 5  /* N-4 */ };

enum { MaximumFieldStates = 6 };

const char * field_state_name( FieldState );

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<>
class Field< void , void , void , void , void , void , void , void > {
public:

  MetaData & mesh_meta_data() const { return m_mesh_meta_data ; }

  unsigned mesh_meta_data_ordinal() const { return m_mesh_meta_data_ordinal ; }

  const std::string & name() const { return m_name ; }

  template<typename NumType> bool type_is() const
    { return m_scalar_type == NumericEnum<NumType>::value ; }

  unsigned numeric_type_ordinal() const { return m_scalar_type ; }

  unsigned number_of_states() const { return m_num_states ; }

  FieldState state() const { return m_this_state ; }

  unsigned rank() const { return m_rank ; }

  const ArrayDimTag * const * dimension_tags() const
    { return m_dim_tags ; }

  unsigned max_size( EntityType ) const ;

  //----------------------------------------

  template<class A>
  const A * attribute() const { return m_cset.template get<A>(); }

  //----------------------------------------
  /** An internal data structure that should never need to be
   *  used by a user of the phdBulkData package.
   */
  struct Restriction {
    unsigned key ; /* ( Entity type , part ordinal ) */
    unsigned stride[ MaximumFieldDimension ];

    Restriction();
    Restriction( const Restriction & rhs );
    Restriction & operator = ( const Restriction & rhs );

    Restriction( EntityType t , unsigned );
    EntityType type() const ;
    unsigned   ordinal() const ;

    static unsigned key_value( EntityType , unsigned );

    bool operator < ( const Restriction & rhs ) const { return key < rhs.key ; }
    bool operator < ( const unsigned rhs_key ) const { return key < rhs_key ; }

  private:
    enum { ord_digits = std::numeric_limits<unsigned>::digits -
                        entity_key_type_digits };
  };

  /** Volatile until the mesh_meta_data is committed */
  const std::vector<Restriction> & restrictions() const ;

  /** Volatile until the mesh_meta_data is committed */
  const Restriction & restriction( EntityType , const Part & ) const ;

  //----------------------------------------

private:

  friend class MetaData ;

  template< typename Scalar , class Tag1 , class Tag2 ,
                              class Tag3 , class Tag4 ,
                              class Tag5 , class Tag6 ,
                              class Tag7 >
    friend class phdmesh::Field ;

  ~Field();

  Field();
  Field( const Field & );
  Field & operator = ( const Field & );

  Field( MetaData & ,
         const std::string & ,
         unsigned scalar_type ,
         unsigned rank ,
         const ArrayDimTag * const * ,
         unsigned number_of_states ,
         FieldState );

  std::vector<Restriction> & restrictions();

  CSet        m_cset ;
  std::string m_name ;
  MetaData & m_mesh_meta_data ;   // in which this field resides
  unsigned    m_mesh_meta_data_ordinal ; // Ordinal in the field set
  unsigned    m_scalar_type ;    // Ordinal in FieldTypes
  unsigned    m_rank ;           // Number of dimensions
  unsigned    m_num_states ;     // Number of states of this field
  FieldState  m_this_state ;     // Field state of this field

  std::vector<Restriction> m_dim_map ; // Only valid on StateNone
  Field             * m_field_states[ MaximumFieldStates ];
  const ArrayDimTag * m_dim_tags[ MaximumFieldDimension ];
};

//----------------------------------------------------------------------

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
                            class Tag5 , class Tag6 , class Tag7 >
class Field : public FieldBase {
public:

  /** Query this field for a given field state. */
  const Field & operator[]( FieldState state ) const
    { return static_cast<Field &>( * FieldBase::m_field_states[state] ); }

  Field & operator[]( FieldState state )
    { return static_cast<Field &>( * FieldBase::m_field_states[state] ); }

private:

  ~Field();
  Field();
  Field( const Field & );
  Field & operator = ( const Field & );
};

}

//----------------------------------------------------------------------
/** \cond */
// Internal implementation details to follow.
//----------------------------------------------------------------------

namespace phdmesh {

template< typename Scalar >
struct FieldTraits< Field<Scalar,void,void,void,void,void,void,void> >
{
  typedef Scalar data_type ;
  typedef void   tag1 ;
  typedef void   tag2 ;
  typedef void   tag3 ;
  typedef void   tag4 ;
  typedef void   tag5 ;
  typedef void   tag6 ;
  typedef void   tag7 ;

  enum { Numeric = NumericEnum< data_type >::value };
  enum { Rank = 0 };

  unsigned * const stride ;

  FieldTraits() : stride(NULL) {}

  static const ArrayDimTag * const * tags()
    { return array_dim_tags<void,void,void,void,void,void,void,void>(); }
};


template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct FieldTraits< Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
{
private:
  typedef
    Array<void,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> atype ;
public:
  typedef Scalar data_type ;
  typedef Tag1   tag1 ;
  typedef Tag2   tag2 ;
  typedef Tag3   tag3 ;
  typedef Tag4   tag4 ;
  typedef Tag5   tag5 ;
  typedef Tag6   tag6 ;
  typedef Tag7   tag7 ;

  enum { Numeric = NumericEnum< data_type >::value };
  enum { Rank = atype::Rank };

  unsigned stride[ Rank ];

  FieldTraits()
    { atype::assign( stride ); }

  FieldTraits( const unsigned & n7 )
    { atype::assign( stride , n7 ); }

  FieldTraits( const unsigned & n6 , const unsigned & n7 )
    { atype::assign( stride , n6 , n7 ); }

  FieldTraits( const unsigned & n5 , const unsigned & n6 , const unsigned & n7 )
    { atype::assign( stride , n5 , n6 , n7 ); }

  FieldTraits( const unsigned & n4 , const unsigned & n5 ,
               const unsigned & n6 , const unsigned & n7 )
    { atype::assign( stride , n4 , n5 , n6 , n7 ); }

  FieldTraits( const unsigned & n3 , const unsigned & n4 , const unsigned & n5 ,
               const unsigned & n6 , const unsigned & n7 )
    { atype::assign( stride , n3 , n4 , n5 , n6 , n7 ); }

  FieldTraits( const unsigned & n2 , const unsigned & n3 ,
               const unsigned & n4 , const unsigned & n5 ,
               const unsigned & n6 , const unsigned & n7 )
    { atype::assign( stride , n2 , n3 , n4 , n5 , n6 , n7 ); }

  FieldTraits( const unsigned & n1 , const unsigned & n2 , const unsigned & n3 ,
               const unsigned & n4 , const unsigned & n5 ,
               const unsigned & n6 , const unsigned & n7 )
    { atype::assign( stride , n1 , n2 , n3 , n4 , n5 , n6 , n7 ); }

  static const ArrayDimTag * const * tags()
    { return array_dim_tags<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>(); }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

inline
FieldBase::Restriction::Restriction()
  : key(0) { Copy<MaximumFieldDimension>( stride , 0u ); }

inline
FieldBase::Restriction::Restriction( const FieldBase::Restriction & rhs )
  : key( rhs.key ) { Copy< MaximumFieldDimension >( stride , rhs.stride ); }

inline
FieldBase::Restriction &
FieldBase::Restriction::operator = ( const FieldBase::Restriction & rhs )
  {
    key = rhs.key ;
    Copy< MaximumFieldDimension >( stride , rhs.stride );
    return *this ;
  }

inline
unsigned FieldBase::Restriction::key_value( EntityType t , unsigned ord )
{ return ( ((unsigned)t) << ord_digits ) | ord ; }

inline
FieldBase::Restriction::Restriction( EntityType t , unsigned ord )
  : key( key_value(t,ord) )
    { Copy< MaximumFieldDimension >( stride , 0u ); }

inline
EntityType
FieldBase::Restriction::type() const
{ return EntityType( key >> ord_digits ); }

inline
unsigned FieldBase::Restriction::ordinal() const
{
  enum { mask = ~( ~(0u) << ord_digits ) };
  return key & mask ;
}

//----------------------------------------------------------------------
/** A field relation for fields that are a pointers to other fields.
 *
 *  If Entity 'e1' has a relation to Entity 'e2' that is in the
 *     domain of the relation stencil 'm_function' AND
 *     field 'm_root'   has a pointer scalar type 'T *' AND
 *     field 'm_target' has a scalar type 'T' AND
 *     field 'm_root'   has field data for Entity 'e1' AND
 *     field 'm_target' has field data for Entity 'e2' AND
 *     the 'e1' to 'e2' relation identifier is within the
 *     the 'm_root' field data size
 *   then
 *      field_data(*m_root,e1)[index] == field_data(*m_target,e2)
 *      where index = (*m_function)( e1.entity_type() ,
 *                                   e2.entity_type() ,
 *                                   relation_identifier_e1_to_e2 ,
 *                                   relation_kind_e1_to_e2 )
 *
 *  This data structure is used internally and should never need to be
 *  used by a user of the phdBulkData package.
 */
struct FieldRelation {
  FieldBase          * m_root ;
  FieldBase          * m_target ;
  relation_stencil_ptr m_function ;

  FieldRelation() : m_root( NULL ), m_target( NULL ), m_function( NULL ) {}

  FieldRelation( const FieldRelation & rhs )
    : m_root( rhs.m_root ),
      m_target( rhs.m_target ),
      m_function( rhs.m_function ) {}

  FieldRelation & operator = ( const FieldRelation & rhs )
    {
      m_root = rhs.m_root ;
      m_target = rhs.m_target ;
      m_function = rhs.m_function ;
      return *this ;
    }
};

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

