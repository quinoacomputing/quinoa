//******************************************************************************
/*!
  \file      src/Base/Factory.h
  \author    J. Bakosi
  \date      Thu 06 Feb 2014 08:55:59 PM MST
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

//! Register option into factory with enum key
template< class C, class F, class O, typename E, typename... Args >
void regist( F& f, std::list<E>& reg, const O& o, E e, Args&&... args ) {
  reg.push_back( o.template add<C>( f, e, std::forward<Args>(args)... ) );
}

template< class C, class Key, class Factory, typename... ConstructorArgs >
void add( Factory& factory, const Key& key, ConstructorArgs&&... args ) {
  factory[ key ] = std::bind( boost::factory< C* >(),
                              std::forward< ConstructorArgs >( args )... );
}

template< class C, class Key, class Factory, typename... ConstructorArgs >
void regSDE( Factory& factory, const Key& key, ConstructorArgs&&... args ) {
  add< C >( factory, key, std::forward< ConstructorArgs >( args )... );
}

} // tk::

#endif // Factory_h
