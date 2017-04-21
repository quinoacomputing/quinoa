// *****************************************************************************
/*!
  \file      src/Base/Factory.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Factory utilities
  \details   Factory utilities. The functions defined in this file help
    interfacing with object factories. For a short introduction on what
    factories are good for, see
    http://www.boost.org/doc/libs/release/libs/functional/factory.
*/
// *****************************************************************************
#ifndef Factory_h
#define Factory_h

#include <list>
#include <functional>

#include "NoWarning/bind.h"
#include "NoWarning/factory.h"
#include "NoWarning/value_factory.h"

#include "Exception.h"

namespace tk {

//! Register class into factory with given key. This is used to register a
//! derived-class object's constructor (deriving from some base class) to a
//! factory. The factory itself is a std::map< key, std::function< Child*() > >,
//! i.e., an associative container, associating some key to a std::function
//! object holding a pointer of Child's base class constructor. The constructor
//! and its bound arguments are stored via boost::factory, which, in this
//! use-case, yields the correct function object of type Base constructor pointer
//! and thus facilitates runtime polymorphism. This function works in conjunction
//! with boost::factory, i.e., uses reference semantics (works with storing
//! pointers of objects). For a simple example on how to use this function, see
//! UnitTest/tests/Base/Factory.h.
//! \param[in] f Factory to register to (std::map with value using reference
//!   semantics)
//! \param[in] key Key used to identify the entry in the factory
//! \param[in] args Variable number of arguments to pass to the constructor
//!   being registered. Note that the constructor arguments are only bound to
//!   the constructor and stored in the factory (an std::map with value using
//!   reference semantics). The object is not instantiated here, i.e., the
//!   constructor is not called here. The object can be instantiated by function
//!   instantiate. \see instantiate
//! \author J. Bakosi
template< class C, class Key, class Factory, typename... ConstructorArgs >
void record( Factory& f, const Key& key, ConstructorArgs&&... args ) {
  f.emplace( key,
             std::bind( boost::factory< C* >(),
                        std::forward< ConstructorArgs >( args )... ) );
}

//! Instantiate object from factory. Factory must have a mapped_value which must
//! have a result_type ptr, e.g., std::map< Key, std::function< Obj*() > >. This
//! wrapper function can be used to instantiate an derived-class object from a
//! factory, repeatedly filled with wrapper function 'record' above. The
//! factory, as described in the documentation of 'record', stores base class
//! pointers in an associative container, thereby facilitating runtime
//! polymorphism and a simple lookup-and-instantiate-style object creation. The
//! object instantiated is of type Child class. This function works in
//! conjunction with boost::factory, i.e., uses reference semantics (works with
//! storing pointers of objects). For a simple example on how to
//! use this function, see UnitTest/tests/Base/Factory.h.
//! \param[in] f Factory to instantiate object from (std::map with value using
//!   reference semantics)
//! \param[in] key Key used to identify the object to instantiate from factory
//! \return std::unique_ptr pointing to the object instantiated from factory
//! \see record
//! \author J. Bakosi
template< class Factory, class Key,
          class Obj = typename std::remove_pointer<
                        typename Factory::mapped_type::result_type >::type >
std::unique_ptr< Obj > instantiate( const Factory& f, const Key& key ) {
  const auto it = f.find( key );
  Assert( it != end( f ), "No such object registered in factory" );
  return std::unique_ptr< Obj >( it->second() );
}

//! Register "model" class of "host" into factory with given key. This wrapper
//! can be used to in a similar manner to 'record', but uses
//! boost::value_factory to bind the model object constructor to its arguments
//! and place it in the associative container storing host class objects. The
//! container is thus of type std::map< key, std::function< T() > >, i.e.,
//! associating a key to a function holding a constructor (and not its
//! pointer). Runtime polymorphism here is realized entirely within the "base"
//! class. See walker::DiffEq in DiffEq/DiffEq.h for an example and more
//! information on runtime polymorphism without client-side inheritance. As a
//! result, this wrapper works with factories that use value semantics, as opposed
//! to 'record' and instantiate which work with reference semantics factories.
//! In order to differentiate between runtime polymorphic classes using
//! reference semantics, consistent with classes realizing runtime polymorphism
//! without client-side inheritance, we call Host as the "Base" class and Model
//! as the "derived" (or child) class. This wrapper function works in
//! conjunction with boost::value_factory, i.e., uses value semantics (works
//! with storing objects instead of object pointers). For a simple example on
//! how to use this function, see UnitTest/tests/Base/Factory.h.
//! \param[in] f Factory to register to (std::map with value using value
//!   semantics)
//! \param[in] key Key used to identify the entry in the factory
//! \param[in] args Variable number of arguments to pass to the constructor
//!   being registered. Note that the constructor arguments are only bound to
//!   the constructor and stored in the factory (an std::map with value using
//!   value semantics). The object is not instantiated here, i.e., the
//!   constructor is not called here. The object can be instantiated by simply
//!   calling the function call operator () on the mapped value. For an example,
//!   RNGStack::selected() in RNG/RNGStack.C.
//! \author J. Bakosi
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

//! Register model class of host into factory with given key using late binding.
//! This variant of 'record' is very similar to 'recordModel', but registers a
//! model class constructor to a factory with late binding of the constructor
//! argument. Late binding allows specifying the constructor argument at the
//! time when the object is instantiated instead of at the time when it is
//! registered. This has all the benefits of using a factory and allows passing
//! information into the model object only when it is available. The late bind
//! is facilitated via boost::bind instead of std::bind using a placeholder,
//! _1, which stands for the first argument (bound later, i.e., not here). The
//! value of the model constructor argument is then not used here, only its
//! type, used to perform the late binding. The binding happens to both the
//! model constructor via std::function (passed to the host constructor) as well
//! as explicitly to the host constructor. Prescribing late binding to the model
//! constructor ensures that the compiler requires the argument to the model
//! constructor, i.e., ensures that the host constructor is required to pass the
//! argument to the model constructor. Prescribing late binding to the host
//! constructor puts in the actual request that an argument (with the correct
//! type) must be passed to the host constructor at instantiate time, which then
//! will forward it to the model constructor. See also, for example,
//! walker::DiffEq's corresponding constructor. An example of client-side code
//! is in walker::DiffEqStack::registerDiffEq for registration into factory, and
//! DiffEqStack::createDiffEq for instantiation late-passing the argument.
//! \param[in] f Factory to register to (std::map with value using value
//!   semantics)
//! \param[in] key Key used to identify the entry in the factory
//! \warning Only works with a single constructor argument
//! \author J. Bakosi
template< class Host, class ModelConstructor, class Factory, class Key,
          typename ModelConstrArg >
void recordModelLate( Factory& f, const Key& key, ModelConstrArg ) {
  // Prescribe late binding the model constructor to its single argument
  std::function< ModelConstructor(const ModelConstrArg&) > c =
    boost::bind( boost::value_factory< ModelConstructor >(), boost::arg<1>() );
  // Bind host to std::function of model constructor and place in factory and
  // also explicitly bind single model constructor argument to host constructor
  f.emplace( key,
    boost::bind( boost::value_factory< Host >(), std::move(c), boost::arg<1>() )
  );
}

//! Register Charm++ model class of host into factory with given key. We bind a
//! host constructor to its arguments of which the first one is a std::function
//! holding a model constructor type (modeling, i.e., used polymorhically with
//! host), followed by an optional number of others (possibly zero) with
//! arbitrary types. Note that the model constructor is a nullptr (default-
//! constructed) and only used to forward its type to the call site inside
//! std::function. The host constructor function is then placed into the
//! factory. This is because Charm++ chares do not explicitly invoke
//! constructors, only call ckNew() on their proxy, which requires all
//! constructor arguments to be present and forwarded to the actual constructor
//! that is only called at a later point in time. This can then be used by those
//! constructors of hosts that invoke the model constructors' proxies' ckNew()
//! and ignore the std::function. See, e.g., rngtest::Battery() and the
//! associated unit tests in UnitTest/tests/Base/Factory.h.
//! \param[in] f Factory to register to (std::map with value using value
//!   semantics)
//! \param[in] key Key used to identify the entry in the factory
//! \param[in] args Variable number of arguments to pass to the constructor
//!   being registered.
//! \author J. Bakosi
template< class Host, class ModelConstructor, class Factory, class Key,
          typename... ModelConstrArgs >
void recordCharmModel( Factory& f, const Key& key, ModelConstrArgs&&... args ) {
  f.emplace( key, std::bind( boost::value_factory< Host >(),
                             std::function< ModelConstructor() >(), // nullptr
                             std::forward< ModelConstrArgs >( args )...) );
}

} // tk::

#endif // Factory_h
