//******************************************************************************
/*!
  \file      src/Base/Factory.h
  \author    J. Bakosi
  \date      Thu 13 Feb 2014 09:18:43 PM CET
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Factory utils
  \details   Factory utils
*/
//******************************************************************************
#ifndef Factory_h
#define Factory_h

#include <list>

#include <boost/functional/factory.hpp>

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

} // tk::

#endif // Factory_h
