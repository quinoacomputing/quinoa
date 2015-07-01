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

#ifndef phdmesh_FieldData_hpp
#define phdmesh_FieldData_hpp

//----------------------------------------------------------------------

#include <mesh/Field.hpp>
#include <mesh/Kernel.hpp>
#include <mesh/Entity.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** Check validity of field and kernel: compatibility and existence.
 *  For performance none of the remaining field_data functions have
 *  internal validity checks.
 */
bool field_data_valid( const FieldBase & f ,
                       const Kernel & k ,
                       unsigned ord = 0,
                       const char * required_by = NULL );

inline
bool field_data_valid( const FieldBase & f ,
                       const Entity & e ,
                       const char * required_by = NULL )
{ return field_data_valid( f, e.kernel(), e.kernel_ordinal(), required_by ); }

//----------------------------------------------------------------------

inline
unsigned field_data_size( const FieldBase & f , const Kernel & k )
{
  const Kernel::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];
  return pd.m_size ;
}

inline
unsigned field_data_size( const FieldBase & f , const Entity & e )
{ return field_data_size( f , e.kernel() ); }

//----------------------------------------------------------------------

template< class field_type >
inline
typename FieldTraits< field_type >::data_type *
field_data( const field_type & f , const Kernel & k )
{
  typedef unsigned char * byte_p ;
  typedef typename FieldTraits< field_type >::data_type * data_p ;

  data_p ptr = NULL ;

  {
    const Kernel::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];

    if ( pd.m_size ) {
      ptr = (data_p)( ((byte_p)(k.m_entities)) + pd.m_base );
    }
  }
  return ptr ;
}

template< class field_type >
inline
typename FieldTraits< field_type >::data_type *
field_data( const field_type & f , const Entity & e )
{
  typedef unsigned char * byte_p ;
  typedef typename FieldTraits< field_type >::data_type * data_p ;

  data_p ptr = NULL ;

  {
    const Kernel & k = e.kernel();
    const unsigned i = e.kernel_ordinal();

    const Kernel::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];

    if ( pd.m_size ) {
      ptr = (data_p)( ((byte_p)(k.m_entities)) + pd.m_base + pd.m_size * i );
    }
  }
  return ptr ;
}

//----------------------------------------------------------------------

template< typename Scalar >
struct EntityArray< Field<Scalar,void,void,void,void,void,void,void> >
  : public Array<Scalar,RankZero,void,void,void,void,void,void,void>
{
  typedef Field<Scalar,void,void,void,void,void,void,void> field_type ;
  typedef Array<Scalar,RankZero,void,void,void,void,void,void,void> array_type ;

  EntityArray( const field_type & f , const Entity & e )
    {
      typedef unsigned char * byte_p ;
      typedef Scalar * data_p ;
      enum { R = array_type::Rank };

      const Kernel & k = e.kernel();

      const Kernel::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];

      if ( pd.m_size ) {
        array_type::m_ptr = (data_p)(
                (byte_p)(k.m_entities) + pd.m_base +
                                         pd.m_size * e.kernel_ordinal() );
      }
    }
private:
  EntityArray();
  EntityArray( const EntityArray & );
  EntityArray & operator = ( const EntityArray & );
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct EntityArray< Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
  : public Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
{
  typedef Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> field_type ;
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> array_type ;

  EntityArray( const field_type & f , const Entity & e )
    {
      typedef unsigned char * byte_p ;
      typedef Scalar * data_p ;
      enum { R = array_type::Rank };

      const Kernel & k = e.kernel();

      const Kernel::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];

      if ( pd.m_size ) {
        array_type::m_ptr = (data_p)(
                (byte_p)(k.m_entities) + pd.m_base +
                                         pd.m_size * e.kernel_ordinal() );
        Copy< R >( array_type::m_stride , pd.m_stride );
      }
    }
private:
  EntityArray();
  EntityArray( const EntityArray & );
  EntityArray & operator = ( const EntityArray & );
};

template< typename Scalar >
struct KernelArray< Field<Scalar,void,void,void,void,void,void,void> >
  : public
      Array<Scalar,FortranOrder,EntityDimension,void,void,void,void,void,void> 
{
  typedef Field<Scalar,void,void,void,void,void,void,void> field_type ;
  typedef
    Array<Scalar,FortranOrder,EntityDimension,void,void,void,void,void,void> 
      array_type ;

  KernelArray( const field_type & f , const Kernel & k )
    {
      typedef unsigned char * byte_p ;
      typedef Scalar * data_p ;
      enum { R = array_type::Rank };

      const Kernel::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];

      if ( pd.m_size ) {
        array_type::m_ptr = (data_p)( (byte_p)(k.m_entities) + pd.m_base );
        array_type::m_stride[0] = k.size();
      }
    }
private:
  KernelArray();
  KernelArray( const KernelArray & );
  KernelArray & operator = ( const KernelArray & );
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct KernelArray< Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
  : public ArrayAppend<
      Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> ,
      EntityDimension >::type
{
  typedef Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> field_type ;
  typedef typename ArrayAppend<
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> ,
    EntityDimension >::type array_type ;

  KernelArray( const field_type & f , const Kernel & k )
    {
      typedef unsigned char * byte_p ;
      typedef Scalar * data_p ;
      enum { R = array_type::Rank };

      const Kernel::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];

      if ( pd.m_size ) {
        array_type::m_ptr = (data_p)( (byte_p)(k.m_entities) + pd.m_base );

        Copy< R - 1 >( array_type::m_stride , pd.m_stride );
        array_type::m_stride[R-1] = array_type::m_stride[R-2] * k.size();
      }
    }
private:
  KernelArray();
  KernelArray( const KernelArray & );
  KernelArray & operator = ( const KernelArray & );
};

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

