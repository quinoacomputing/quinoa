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
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   October 2007
 */

#ifndef util_NamedValue_hpp
#define util_NamedValue_hpp

#include <typeinfo>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cstring>
#include <util/TypeName.hpp>

using std::memcpy;

namespace phdmesh {

template< typename T = void > class NamedValue ;

class NamedValueSet ;

std::istream & operator >> ( std::istream & s , NamedValueSet & v );
std::ostream & operator << ( std::ostream & s , const NamedValueSet & v );

}

//----------------------------------------------------------------------

namespace phdmesh {

class NamedValueSet {
public:
  const std::vector< NamedValue<void> * > get() const { return m_members ; }

  NamedValue<void> *
    find( const std::string & s , const char sep = '.' ) const ;

  NamedValue<void> * insert( NamedValue<void> * );

  void remove( NamedValue<void> * );
  void clear();

  ~NamedValueSet();
  NamedValueSet();
private:
  NamedValueSet( const NamedValueSet & );
  NamedValueSet & operator = ( const NamedValueSet & );
  std::vector< NamedValue<void> * > m_members ;
};

//----------------------------------------------------------------------
/**
 *  \brief  Base class for references provides access to anonymous type.
 */
template<>
class NamedValue<void> {
public:
  const std::string      name ;
  const std::type_info & type ;

  //----------------------------------

  virtual     unsigned get_max() const = 0 ;
  virtual     unsigned put_max() const = 0 ;

  virtual       void * put_void( unsigned ) = 0 ;
  virtual const void * get_void( unsigned ) const = 0 ;

  virtual void     tell(  std::ostream & ) const = 0 ;
  virtual void     write( std::ostream & ) const = 0 ;
  virtual unsigned read(  std::istream & ) = 0 ;

  virtual ~NamedValue();

  //----------------------------------

  template<typename T>
  const T & get( unsigned i = 0 ) const
    {
      const std::type_info & t = typeid(T);
      const void * const p = get_void(i);
      if ( t != type || NULL == p ) { get_throw(t,i); }
      return * ( (const T *) p );
    }

  template<typename T>
  T & put( unsigned i = 0 )
    {
      const std::type_info & t = typeid(T);
      void * const p = put_void(i);
      if ( t != type || NULL == p ) { put_throw(t,i); }
      return * ( (T *) p );
    }

  //----------------------------------

  /** Pack referenced value into a buffer, return the number of bytes */
  virtual unsigned pack( void * ) const = 0 ;

  /** Unpack referenced value from a buffer, return the number of bytes */
  virtual unsigned unpack( void * ) = 0 ;

  //----------------------------------

  const std::vector< NamedValueSet * > references() const ;

protected:

  NamedValue( const std::string & n , const std::type_info & t )
    : name(n), type(t) {}

private:

  friend class NamedValueSet ;

  std::vector< NamedValueSet * > m_holders ;

  void get_throw( const std::type_info & , unsigned ) const ;
  void put_throw( const std::type_info & , unsigned ) const ;

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------

namespace {

template< typename T >
unsigned read_array( std::istream & s , T * const v , const unsigned n )
{
  unsigned i = 0 ;
  while ( i < n && ( s >> v[i] ) ) { ++i ; }
  s.clear( s.rdstate() & ~std::ios::failbit );
  return i ;
}

template< typename T >
unsigned read_vector( std::istream & s , std::vector<T> & v )
{
  const unsigned n = v.size();
  unsigned i = 0 ;
  for ( T tmp ; s >> tmp ; ++i ) {
    if ( i < n ) { v[i] = tmp ; } 
    else { v.push_back( tmp ); }
  }
  s.clear( s.rdstate() & ~std::ios::failbit );
  return i ;
}

template< typename T >
void write_array( std::ostream & s , const T * const v , const unsigned n )
{
  for ( unsigned i = 0 ; i < n ; ++i ) {
    if ( i ) { s << " " ; }
    s << v[i] ;
  }
}

template< typename T >
unsigned pack_array( void * b , const T * const p , unsigned n )
{
  n = ((const unsigned char *)(p+n)) - ((const unsigned char *)(p));
  if ( b ) { memcpy( b , p , n ); }
  return n ;
}

template< typename T >
unsigned unpack_array( const void * b , T * const p , unsigned n )
{
  if ( b ) {
    n = ((unsigned char *)(p+n)) - ((unsigned char *)(p));
    memcpy( p , b , n );
  }
  else {
    n = 0 ;
  }
  return n ;
}

template< typename T >
unsigned pack_vector( void * b , const std::vector<T> & v )
{
  const unsigned n = v.size();
  if ( b ) {
    memcpy( b , & n , sizeof(unsigned) );
    b = ((unsigned char *)b) + sizeof(unsigned);
  }
  return sizeof(unsigned) + pack_array<T>( b , & v[0] , n );
}

template< typename T >
unsigned unpack_vector( void * b , std::vector<T> & v )
{
  unsigned n = 0 ;
  if ( b ) {
    memcpy( & n , b , sizeof(unsigned) );
    b = ((unsigned char *)b) + sizeof(unsigned);
    if ( v.size() < n ) { v.resize(n); }
    n = sizeof(unsigned) + unpack_array<T>( b , & v[0] , n );
  }
  return n ;
}

unsigned pack_value( void * b , const std::string & s )
{
  const unsigned n = s.size() + 1 ;
  if ( b ) { memcpy( b , s.c_str() , n ); }
  return n ;
}

unsigned unpack_value( void * b , std::string & s )
{
  unsigned n = 0 ;
  if ( b ) { s.assign( (char *) b ); n = s.size() + 1 ; }
  return n ;
}

template< typename T >
unsigned pack_value( void * b , const T & v )
{
  const unsigned n = sizeof(T);
  if ( b ) { memcpy( b , & v , n ); }
  return n ;
}

template< typename T >
unsigned unpack_value( void * b , T & v )
{
  unsigned n = 0 ;
  if ( b ) { memcpy( & v , b , n = sizeof(T) ); }
  return n ;
}

}

//----------------------------------------------------------------------
/**
 * \brief  NamedValue to an ordinary value.
 */
template< typename T >
class NamedValue : public NamedValue<void> {
public:

  T value ;

  ~NamedValue() {}

  /** \brief  Constructor */
  NamedValue( const std::string & n ) : NamedValue<void>( n , typeid(T) ) {}
  NamedValue( const std::string & n , const T & v )
    : NamedValue<void>( n , typeid(T) ), value(v) {}

  /** \brief  Tell the type to a stream */
  void tell( std::ostream & s ) const
    { s << name << " = " << TypeName<T>::value(); }

  /** \brief  Write value to stream */
  void write( std::ostream & s ) const { s << value ; }

  /** \brief  Read value from stream */
  unsigned read( std::istream & s ) { return s >> value ? 1 : 0 ; }

  unsigned pack(   void * b ) const { return pack_value( b , value ); }
  unsigned unpack( void * b )       { return unpack_value( b , value ); }

private:

        unsigned get_max() const { return 1 ; }
        unsigned put_max() const { return 1 ; }
        void * put_void( unsigned i )       { return i ? NULL : & value ; }
  const void * get_void( unsigned i ) const { return i ? NULL : & value ; }

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------
/**
 * \brief  NamedValue to an ordinary reference.
 */
template< typename T >
class NamedValue< T & > : public NamedValue<void> {
public:

  T & ref ;

  ~NamedValue() {}

  /** \brief  Constructor */
  NamedValue( const std::string & n , T & v )
    : NamedValue<void>( n , typeid(T) ), ref(v) {}

  /** \brief  Tell the type to a stream */
  void tell( std::ostream & s ) const
    { s << name << " = " << TypeName<T &>::value(); }

  /** \brief  Write value to stream */
  void write( std::ostream & s ) const { s << ref ; }

  /** \brief  Read value from stream */
  unsigned read( std::istream & s ) { return s >> ref ? 1 : 0 ; }

  unsigned pack(   void * b ) const { return pack_value( b , ref ); }
  unsigned unpack( void * b )       { return unpack_value( b , ref ); }

private:

        unsigned get_max() const { return 1 ; }
        unsigned put_max() const { return 1 ; }
        void * put_void( unsigned i )       { return i ? NULL : & ref ; }
  const void * get_void( unsigned i ) const { return i ? NULL : & ref ; }

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------
/**
 * \brief  NamedValue to an ordinary constant value.
 */
template< typename T >
class NamedValue< const T & > : public NamedValue<void> {
public:
  const T & ref ;

  NamedValue( const std::string & n , const T & v )
    : NamedValue<void>( n , typeid(T) ), ref(v) {}

  ~NamedValue() {}

  void write( std::ostream & s ) const { s << ref ; }

  unsigned read( std::istream & ) { return 0 ; }

  void tell( std::ostream & s ) const
    { s << name << " = " << TypeName<const T &>::value(); }

  unsigned pack( void * b ) const { return pack_value( b , ref ); }

  unsigned unpack( void * b ) { return 0 ; }

private:

        unsigned get_max() const { return 1 ; }
        unsigned put_max() const { return 0 ; }
        void * put_void( unsigned i ) { return NULL ; }
  const void * get_void( unsigned i ) const { return i ? NULL : & ref ; }

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------

/**
 * \brief  NamedValue to a fixed size array.
 */
template< typename T , unsigned N >
class NamedValue< T[N] > : public NamedValue<void> {
public:
  T value[N] ;

  NamedValue( const std::string & n ) : NamedValue<void>( n , typeid(T) ) {}

  ~NamedValue() {}

  void write( std::ostream & s ) const { write_array<T>(s,value,N); }

  unsigned read( std::istream & s ) { return read_array<T>(s,value,N); }

  void tell( std::ostream & s ) const
    { s << name << " = " << TypeName<T[N]>::value(); }

  unsigned pack(   void * b ) const { return pack_array<T>( b , value, N); }
  unsigned unpack( void * b )       { return unpack_array<T>(b, value, N); }

private:

        unsigned get_max() const { return N ; }
        unsigned put_max() const { return N ; }
        void * put_void( unsigned i )       { return i ? NULL : value + i ; }
  const void * get_void( unsigned i ) const { return i ? NULL : value + i ; }

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------
/**
 * \brief  NamedValue to a fixed size array of ordinary values.
 */
template< typename T >
class NamedValue< T * > : public NamedValue<void> {
public:
  T * const      ref ;
  const unsigned size ;

  NamedValue( const std::string & n , T * v , unsigned s )
    : NamedValue<void>( n , typeid(T) ), ref(v), size(s) {}

  ~NamedValue() {}

  void write( std::ostream & s ) const { write_array<T>(s,ref,size); }

  unsigned read( std::istream & s )
    { return read_array<T>(s,ref,size); }

  void tell( std::ostream & s ) const
    { s << name << " = " << type_name_array<T>(size) << " *" ; }

  unsigned pack(   void * b ) const { return pack_array<T>(  b,ref,size); }
  unsigned unpack( void * b )       { return unpack_array<T>(b,ref,size); }

private:

        unsigned get_max() const { return size ; }
        unsigned put_max() const { return size ; }
        void * put_void(unsigned i)       {return i < size ? ref+i : NULL;}
  const void * get_void(unsigned i) const {return i < size ? ref+i : NULL;}

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};


/**
 * \brief  NamedValue to a const fixed size array of ordinary values .
 */
template< typename T >
class NamedValue< const T * > : public NamedValue<void> {
public:
  const T * const ref ;
  const unsigned  size ;

  NamedValue( const std::string & n , const T * v , unsigned s )
    : NamedValue<void>( n , typeid(T) ), ref(v), size(s) {}

  ~NamedValue() {}

  void write( std::ostream & s ) const { write_array<T>(s,ref,size); }
  unsigned read( std::istream & ){ return 0 ; }

  void tell( std::ostream & s ) const
    { s << name << " = " << type_name_array<const T>(size) << " *" ; }

  unsigned pack( void * b ) const { return pack_array<T>(b,ref,size); }
  unsigned unpack( void * b ) { return 0 ; }

private:

        unsigned get_max() const { return size ; }
        unsigned put_max() const { return 0 ; }
        void * put_void(unsigned ) {return NULL ; }
  const void * get_void(unsigned i) const {return i < size ? ref+i : NULL;}

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------
/**
 *  \brief  NamedValue to a std::vector
 */
template< typename T >
class NamedValue< std::vector<T> > : public NamedValue<void> {
public:
  std::vector<T> value ;

  NamedValue( const std::string & n ) : NamedValue<void>( n , typeid(T) ) {}

  NamedValue( const std::string & n , const std::vector<T> & v )
    : NamedValue<void>( n , typeid(T) ), value(v) {}

  ~NamedValue() {}

  void tell( std::ostream & s ) const
    { s << name << " = " << type_name_vector<T>( value.size() ); }

  void write( std::ostream & s ) const
    { write_array<T>( s , & value[0] , value.size() ); }

  unsigned read( std::istream & s )
    { return read_vector<T>( s , value ); }

  unsigned pack( void * b ) const { return pack_vector<T>( b , value ); }

  unsigned unpack( void * b ) { return unpack_vector<T>( b , value ); }

private:

  unsigned get_max() const { return value.size(); }
  unsigned put_max() const { return std::numeric_limits<unsigned>::max(); }

  void * put_void( unsigned i )
    {
      if ( value.size() <= i ) { value.resize(i+1); }
      return & value[i] ;
    }

  const void * get_void( unsigned i ) const
    { return i < value.size() ? & value[i] : NULL ; }

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

template< typename T >
class NamedValue< std::vector<T> & > : public NamedValue<void> {
public:
  std::vector<T> & ref ;

  NamedValue( const std::string & n , std::vector<T> & v )
    : NamedValue<void>( n , typeid(T) ), ref(v) {}

  ~NamedValue() {}

  void tell( std::ostream & s ) const
    { s << name << " = " << type_name_vector<T>( ref.size() ) << " &" ; }

  void write( std::ostream & s ) const
    { write_array<T>( s , & ref[0] , ref.size() ); }

  unsigned read( std::istream & s )
    { return read_vector<T>( s , ref ); }

  unsigned pack( void * b ) const { return pack_vector<T>( b , ref ); }

  unsigned unpack( void * b ) { return unpack_vector<T>( b , ref ); }

private:

  unsigned get_max() const { return ref.size(); }
  unsigned put_max() const { return std::numeric_limits<unsigned>::max(); }

  void * put_void( unsigned i )
    {
      if ( ref.size() <= i ) { ref.resize(i+1); }
      return & ref[i] ;
    }

  const void * get_void( unsigned i ) const
    { return i < ref.size() ? & ref[i] : NULL ; }

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------
/**
 *  \brief  NamedValue to a const std::vector
 */
template< typename T >
class NamedValue< const std::vector<T> & > : public NamedValue<void> {
public:
  const std::vector<T> & ref ;

  explicit NamedValue( const std::string & n , const std::vector<T> & arg )
    : NamedValue<void>( n , typeid(T) ), ref(arg) {}

  ~NamedValue() {}

  void tell( std::ostream & s ) const
    { s << name << " = " << type_name_vector<const T>( ref.size() ) << " &" ; }

  void write( std::ostream & s ) const
    { write_array<T>( s , & ref[0] , ref.size() ); }

  unsigned read( std::istream & ) { return 0 ; }

  unsigned pack( void * b ) const { return pack_vector<T>( b , ref ); }

  unsigned unpack( void * b ) { return 0 ; }

private:

  unsigned get_max() const { return ref.size(); }
  unsigned put_max() const { return 0 ; }

  void * put_void( unsigned ) { return NULL ; }

  const void * get_void( unsigned i ) const
    { return i < ref.size() ? & ref[i] : NULL ; }

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<>
class NamedValue< NamedValueSet > : public NamedValue<void> {
public:
  NamedValueSet value ;

  ~NamedValue();

  NamedValue( const std::string & n )
    : NamedValue<void>( n , typeid(NamedValueSet) ) {}

  void tell( std::ostream & ) const ;

  void     write( std::ostream & s ) const { s << value ; }
  unsigned read(  std::istream & s ) { return s >> value ? 1 : 0 ; }

  //----------------------------------

  /** Pack referenced value into a buffer, return the number of bytes */
  unsigned pack( void * ) const ;

  /** Unpack referenced value from a buffer, return the number of bytes */
  unsigned unpack( void * );

private:

        unsigned get_max() const { return 1 ; }
        unsigned put_max() const { return 1 ; }
        void * put_void( unsigned i )       { return i ? NULL : & value ; }
  const void * get_void( unsigned i ) const { return i ? NULL : & value ; }


  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};


} // namespace phdmesh

//----------------------------------------------------------------------

#endif

