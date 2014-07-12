//******************************************************************************
/*!
  \file      src/Base/Factory.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:06:47 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Factory utils
  \details   Factory utils
*/
//******************************************************************************
#ifndef Factory_h
#define Factory_h

#include <list>
#include <functional>

#include <boost/functional/factory.hpp>
#include <boost/functional/value_factory.hpp>

#include <Exception.h>

namespace tk {

//! Register class into factory with given key
template< class C, class Key, class Factory, typename... ConstructorArgs >
void record( Factory& f, const Key& key, ConstructorArgs&&... args ) {
  f.emplace( key,
             std::bind( boost::factory< C* >(),
                        std::forward< ConstructorArgs >( args )... ) );
}

//! Instantiate object from factory
//! Concept: Factory must have a mapped_value which must have a result_type ptr,
//! e.g., std::map< Key, std::function< Obj*() > >
template< class Factory, class Key,
          class Obj = typename std::remove_pointer<
                        typename Factory::mapped_type::result_type >::type >
std::unique_ptr< Obj > instantiate( const Factory& f, const Key& key ) {
  const auto it = f.find( key );
  Assert( it != end( f ), "No such object registered in factory" );
  return std::unique_ptr< Obj >( it->second() );
}

//! Register model class of host into factory with given key
template< class Host, class ModelConstructor, class Factory, class Key,
          typename... ModelConstrArgs >
void recordModel( Factory& f, const Key& key, ModelConstrArgs&&... args ) {
  // Bind model constructor to its arguments
  std::function< ModelConstructor() > c =
    std::bind( boost::value_factory< ModelConstructor >(),
               std::forward< ModelConstrArgs >( args )... );
  // Bind host to std::function of model constructor and place in factory
  f.emplace( key, std::bind( boost::value_factory< Host >(), std::move(c) ) );
}

//! Register Charm++ model class of host into factory with given key. We bind a
//! host constructor to its arguments of which the first one is a std::function
//! holding a model constructor type (modeling, i.e., used polymorhically with
//! host), the second one is a key followed by an optional number of others
//! (possibly zero) with arbitrary types. Note that the model constructor is
//! nullptr and only used to forward its type to the call site inside
//! std::function. The host constructor function is then placed into the factory.
//! This is because Charm++ chares do not explicitly invoke constructors,
//! only call ckNew() on their proxy, which requires all constructor arguments
//! to be present and forwarded to the actual constructor that is only called at
//! a later point in time. This can then be used by those constructors of hosts
//! that invoke the model constructors' proxies' ckNew() and ignore the
//! std::function. See, e.g., rngtest::Battery().
template< class Host, class ModelConstructor, class Factory, class Key,
          typename... ModelConstrArgs >
void recordCharmModel( Factory& f, const Key& key, ModelConstrArgs&&... args ) {
  f.emplace( key, std::bind( boost::value_factory< Host >(),
                             std::function< ModelConstructor() >(),
                             key, std::forward< ModelConstrArgs >( args )...) );
}

} // tk::

#endif // Factory_h
