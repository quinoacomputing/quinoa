//******************************************************************************
/*!
  \file      src/Base/Factory.h
  \author    J. Bakosi
  \date      Sat 08 Mar 2014 04:25:03 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Factory utils
  \details   Factory utils
*/
//******************************************************************************
#ifndef Factory_h
#define Factory_h

#include <list>
#include <functional>

#include <boost/functional/factory.hpp>

#include <Exception.h>

namespace tk {

template< class C, class Key, class Factory, typename... ConstructorArgs >
void add( Factory& factory, const Key& key, ConstructorArgs&&... args ) {
  factory[ key ] = std::bind( boost::factory< C* >(),
                              std::forward< ConstructorArgs >( args )... );
}

//! Register class into factory with given key
template< class C, class Key, class Factory, typename... ConstructorArgs >
void record( Factory& factory, const Key& key, ConstructorArgs&&... args ) {
  add< C >( factory, key, std::forward< ConstructorArgs >( args )... );
}

//! Instantiate object from factory
//! Concept: Factory must have a mapped_value which must have a result_type ptr,
//! e.g., std::map< Key, std::function< Obj*() > >
template< class Factory, class Key,
          class Obj = typename std::remove_pointer<
                        typename Factory::mapped_type::result_type >::type >
std::unique_ptr< Obj > instantiate( const Factory& factory, const Key& key ) {
  auto it = factory.find( key );
  Assert( it != end( factory ), tk::ExceptType::FATAL,
          "No such object registered in factory" );
  return std::unique_ptr< Obj >( it->second() );
}

} // tk::

#endif // Factory_h
