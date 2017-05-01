// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestFactory.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/Factory.h
  \details   Unit tests for Base/Factory.h
*/
// *****************************************************************************
#ifndef test_Factory_h
#define test_Factory_h

#include <unistd.h>

#include "NoWarning/tut.h"

#include "Macro.h"
#include "Make_unique.h"
#include "Factory.h"

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

} // unittest::

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct Factory_common {

  // For testing boost::factory (runtime polymorphism using reference semantics)
  struct Base {
    Base( std::string init ) : type( init ) {}
    std::string type;
  };
  struct Child : Base {
    Child() : Base( "def" ) {}
    Child( int ) : Base( "int" ) {}
  };
  using Factory = std::map< int, std::function< Base*() > >;

  //! For testing boost::value_factory (runtime polymorphism using value
  //! semantics). The struct below uses runtime polymorphism without client-side
  //! inheritance: inheritance is confined to the internals of the struct below,
  //! inivisble to client-code. The struct exclusively deals with ownership
  //! enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
  //! https://github.com/sean-parent/sean-parent.github.com/wiki/
  //! Papers-and-Presentations
  struct VBase {
    //! Constructor taking an object modeling Concept (see below). The object
    //! of class T was pre-constructed.
    template< typename T >
    explicit VBase( T x ) :
      self( tk::make_unique< Model<T> >( std::move(x) ) ),
      ctor( "val" ),
      assg() {}

    //! Constructor taking a function pointer to a constructor of an object
    //! modeling Concept (see below). Passing std::function allows late
    //! execution of the constructor, i.e., as late as inside this struct's
    //! constructor, and thus usage from a factory.
    template< typename T >
    explicit VBase( std::function<T()> x ) :
      self( tk::make_unique< Model<T> >( std::move(x()) ) ),
      ctor( "fun" ),
      assg() {}

    //! Public interface to querying the child constructor type invoked
    std::string Type() const { return self->Type(); }
    //! Copy assignment
    VBase& operator=( const VBase& x )
    { VBase tmp(x); *this = std::move(tmp); assg = "cpy"; return *this; }
    //! Copy constructor
    VBase( const VBase& x ) : self( x.self->copy() ), ctor( "cpy" ), assg() {}
    //! Move assignment: could be default
    VBase& operator=( VBase&& x ) noexcept {
      self = std::move( x.self );
      ctor = std::move( x.ctor );
      assg = "mov";
      return *this;
    }
    //! Move constructor: could be default, but in terms of move assignment
    VBase( VBase&& x ) noexcept { *this = std::move(x); ctor = "mov"; }

    //! Accessor to ctor type
    std::string ctortype() const { return ctor; }
    //! Accessor to assignment type
    std::string assgntype() const { return assg; }

    //! Concept is a pure virtual base struct specifying the requirements of
    //! polymorphic objects deriving from it
    struct Concept {
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual std::string Type() const = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      std::string Type() const override { return data.Type(); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< VBase pointer used polymorphically
    std::string ctor;
    std::string assg;
  };

  //! Child struct used polymorphically with VBase
  struct VChild {
    VChild() : type( "def" ) {}
    VChild( int ) : type( "int" ) {}
    std::string Type() const { return type; }
    std::string type;
  };

  using ValueFactory = std::map< int, std::function< VBase() > >;
}; // Factory_common

//! Test group shortcuts
using Factory_group = test_group< Factory_common, MAX_TESTS_IN_GROUP >;
using Factory_object = Factory_group::object;

//! Define test group
static Factory_group Factory( "Base/Factory" );

//! Test definitions for group

//! Test if tk::record correctly registers a default child constructor in actory
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 1 >() {
  set_test_name( "record default child ctor in factory" );

  Factory_common::Factory f;
  tk::record< Child >( f, 2 );  // key: 2
  ensure_equals( "recorded 1 default ctor in factory", f.size(), 1UL );

  const auto it = f.find( 2 );
  if ( it != f.end() )
    ensure_equals( "ctor() invoked from factory", it->second()->type, "def" );
  else
    fail( "cannot find key in factory" );
}

//! Test if tk::record correctly registers a child constructor in factory
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 2 >() {
  set_test_name( "record child ctor(int) in factory" );

  Factory_common::Factory f;
  tk::record< Child >( f, 2, 23 );  // key: 2, ctor arg: 23
  ensure_equals( "recorded 1 ctor in factory", f.size(), 1UL );

  const auto it = f.find( 2 );
  if ( it != f.end() )
    ensure_equals( "ctor(int) invoked from factory",
                   it->second()->type, "int" );
  else
    fail( "cannot find key in factory" );
}

//! Test if instantiate can correctly instantiates objects from factory
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 3 >() {
  set_test_name( "instantiate objects from factory" );

  Factory_common::Factory f;
  tk::record< Child >( f, 1 );      // key: 1, default ctor
  tk::record< Child >( f, 2, 23 );  // key: 2, ctor arg: 23
  ensure_equals( "recorded 2 ctors in factory", f.size(), 2UL );

  if ( f.find(1) != f.end() )
    ensure_equals( "default-ctor object instantiated from factory",
                   tk::instantiate(f,1)->type, "def" );
  else
    fail( "cannot find key in factory" );

  if ( f.find(2) != f.end() )
    ensure_equals( "ctor(int)-object instantiated from factory",
                   tk::instantiate(f,2)->type, "int" );
  else
    fail( "cannot find key in factory" );
}

//! Test if instantiate correctly throws an exception if key does not exist in
//! factory (only test in DEBUG mode, RELEASE would result in undefined
//! behavior)
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 4 >() {
  set_test_name( "instantiate non-existent object from factory" );

  // Only test in debug mode. The code below would result in undefined behavior
  // in non-debug mode, so we don't test for it here, but skip the test. See
  // also http://stackoverflow.com/a/3829185.
  #ifndef NDEBUG
  Factory_common::Factory f;
  tk::record< Child >( f, 1 );      // key: 1, default ctor

  try {
    tk::instantiate(f,2);
    fail( "tk::instantiate with non-existent key must throw if DEBUG mode" );
  }
  catch( tk::Exception& ) {
    // exception thrown, test ok
  }
  #else
  skip( "in RELEASE mode, would yield undefined behavior" );
  #endif
}

//! Test if tk::recordModel correctly registers a default child constructor in
//! value factory
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 5 >() {
  set_test_name( "record default child ctor in value factory" );

  ValueFactory f;
  tk::recordModel< VBase, VChild >( f, 2 );  // key: 2
  ensure_equals( "recorded 1 default ctor in value factory", f.size(), 1UL );

  const auto it = f.find( 2 );
  if ( it != f.end() )
    ensure_equals( "ctor() invoked from value factory",
                   (it->second()).Type(), "def" );
  else
    fail( "cannot find key in factory" );
}

//! Test if tk::record correctly registers a child constructor in value factory
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 6 >() {
  set_test_name( "record child ctor(int) in value factory" );

  ValueFactory f;
  tk::recordModel< VBase, VChild >( f, 2, 23 );  // key: 2, ctor arg: 23
  ensure_equals( "recorded 1 ctor in value factory", f.size(), 1UL );

  const auto it = f.find( 2 );
  if ( it != f.end() )
    ensure_equals( "ctor(int) invoked from value factory",
                   (it->second()).Type(), "int" );
  else
    fail( "cannot find key in factory" );
}

//! Test value constructor of VBase. This test does not necessarily test
//! functionality in Base/Factory.h, but functionality in the class, used for
//! concept-based polymorphism, used in conjunction with the functions in
//! Base/Factory. It also provides a simple example of how to use the
//! functionality tested.
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 7 >() {
  set_test_name( "val-ctor of concept-based polymorphic base" );

  auto c = VChild();
  VBase a( c );

  ensure_equals( "value constructor used to instantiate concept-based "
                 "polymorphic base via child", a.ctortype(), "val" );
}

//! Test std::function constructor of VBase. This test does not necessarily test
//! functionality in Base/Factory.h, but functionality in the class, used for
//! concept-based polymorphism, used in conjunction with the functions in
//! Base/Factory. It also provides a simple example of how to use the
//! functionality tested.
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 8 >() {
  set_test_name( "func-ctor of concept-based polymorphic base" );

  std::function< VChild() > f = boost::value_factory< VChild >();
  VBase a( f );
  ensure_equals( "std::function constructor used to instantiate concept-based "
                 "polymorphic base via child", a.ctortype(), "fun" );
}

//! Test copy cnstructor of VBase. This test does not necessarily test
//! functionality in Base/Factory.h, but functionality in the class, used for
//! concept-based polymorphism, used in conjunction with the functions in
//! Base/Factory. It also provides a simple example of how to use the
//! functionality tested.
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 9 >() {
  set_test_name( "copy-ctor of concept-based polymorphic base" );

  auto c = VChild();
  VBase a( c );
  std::vector< VBase > v;
  v.push_back( a );
  ensure_equals( "copy constructor used to push_back a concept-based "
                 "polymorphic base to a std::vector", v[0].ctortype(), "cpy" );
}

//! Test move cnstructor of VBase. This test does not necessarily test
//! functionality in Base/Factory.h, but functionality in the class, used for
//! concept-based polymorphism, used in conjunction with the functions in
//! Base/Factory. It also provides a simple example of how to use the
//! functionality tested.
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 10 >() {
  set_test_name( "move-ctor of concept-based polymorphic base" );

  auto c = VChild();
  std::vector< VBase > v;
  v.emplace_back( VBase( c ) );
  ensure_equals( "move constructor used to emplace_back a concept-based "
                 "polymorphic base to a std::vector", v[0].ctortype(), "mov" );
}

//! Test copy assignment of VBase. This test does not necessarily test
//! functionality in Base/Factory.h, but functionality in the class, used for
//! concept-based polymorphism, used in conjunction with the functions in
//! Base/Factory. It also provides a simple example of how to use the
//! functionality tested.
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 11 >() {
  set_test_name( "copy-assgn of concept-based polymorphic base" );

  auto c1 = VChild();
  auto c2 = VChild( 23 );
  VBase b1( c1 );
  VBase b2( c2 );
  b2 = b1;
  ensure_equals( "copy assignment used to copy a concept-based polymorphic "
                 "base", b2.assgntype(), "cpy" );
}

//! Test move assignment of VBase. This test does not necessarily test
//! functionality in Base/Factory.h, but functionality in the class, used for
//! concept-based polymorphism, used in conjunction with the functions in
//! Base/Factory. It also provides a simple example of how to use the
//! functionality tested.
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 12 >() {
  set_test_name( "move-assgn of concept-based polymorphic base" );

  auto c1 = VChild();
  auto c2 = VChild( 23 );
  VBase b1( c1 );
  VBase b2( c2 );
  b2 = std::move( b1 );
  ensure_equals( "move assignment used to move a concept-based polymorphic "
                 "base", b2.assgntype(), "mov" );
  ensure( "move assignment assigns std::unique_ptr in concept-based polymorphic"
          " base", b2.self.get() != nullptr );
  ensure( "move assignment leaves its std::unique_ptr as nullptr in "
          "concept-based polymorphic base", b1.self == nullptr );
}

//! For testing boost::value_factory (runtime polymorphism using value
//! semantics) with Charm++ children. Not in Factory_common so Charm++ can find
//! it. The struct below uses runtime polymorphism without client-side
//! inheritance: inheritance is confined to the internals of the struct below,
//! inivisble to client-code. The struct exclusively deals with ownership
//! enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//! https://github.com/sean-parent/sean-parent.github.com/wiki/
//! Papers-and-Presentations
struct VBase {
  //! Constructor taking a function pointer to a constructor of an object
  //! modeling Concept (see below). Passing std::function allows late
  //! execution of the constructor of T, i.e., at some future time, and thus
  //! usage from a factory. Note that the value of the first function
  //! argument, std::function<T()>, is not used here, but its constructor
  //! type, T, is used to enable the compiler to deduce the model constructor
  //! type, used to create its Charm proxy, defined by T::Proxy. The actual
  //! constructor of T is not called here but at some future time by the
  //! Charm++ runtime system, here only an asynchrounous ckNew() is called,
  //! i.e., a message (or request) for a future call to T's constructor. This
  //! overload is only enabled for Charm++ chare objects defining typedef
  //! 'Proxy', which must define the Charm++ proxy. All optional constructor
  //! arguments are forwarded to ckNew() and thus to T's constructor. If
  //! it was somehow possible to obtain all bound arguments' types and values
  //! from an already-bound std::function, we could use those instead of
  //! having to explicitly forward the model constructor arguments via this
  //! host constructor. See also tk::recordCharmModel().
  template< typename T, typename... ConstrArgs,
    typename std::enable_if< tk::HasTypedefProxy<T>::value, int >::type = 0 >
  explicit VBase( std::function<T()> c, ConstrArgs... args ) :
    self( tk::make_unique< Model< typename T::Proxy > >
         (std::move(T::Proxy::ckNew(std::forward<ConstrArgs>(args)...))) ) {
    Assert( c == nullptr, "std::function argument to VBase Charm "
                          "constructor must be nullptr" );
    #ifdef NDEBUG
    IGNORE(c);
    #endif
  }
  //! Copy assignment
  VBase& operator=( const VBase& x )
  { VBase tmp(x); *this = std::move(tmp); return *this; }
  //! Copy constructor
  VBase( const VBase& x ) : self( x.self->copy() ) {}
  //! Move assignment
  VBase& operator=( VBase&& ) noexcept = default;
  //! Move constructor
  VBase( VBase&& ) noexcept = default;
  //! Concept is a pure virtual base struct specifying the requirements of
  //! polymorphic objects deriving from it
  struct Concept {
    Concept() = default;
    Concept( const Concept& ) = default;
    virtual ~Concept() = default;
    virtual Concept* copy() const = 0;
  };
  //! Model models the Concept above by deriving from it and overriding the
  //! the virtual functions required by Concept
  template< typename T >
  struct Model : Concept {
    Model( T x ) : data( std::move(x) ) {}
    Concept* copy() const { return new Model( *this ); }
    T data;
  };
  std::unique_ptr< Concept > self;    //!< VBase pointer used polymorphically
};

//! Child struct used polymorphically with VBase (Charm++ chare)
//! (not in Factory_common, so Charm++ can find it)
class CharmChild : public CBase_CharmChild {
  public:
  using Proxy = CProxy_CharmChild;
  CharmChild() {
    // If we got here, the second part of this test succeeded. Construct and
    // send back a new test result, with tag "2", signaling the second part.
    tut::test_result tr( "Base/Factory", 7, 
                         "default Charm++ child ctor in val-factory 2",
                          tut::test_result::result_type::ok );
    unittest::g_suiteProxy.evaluate(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }
  CharmChild( tk::real ) {
    // If we got here, the second part of this test succeeded. Construct and
    // send back a new test result, with tag "2", signaling the second part.
    tut::test_result tr( "Base/Factory", 8, 
                         "Charm++ child ctor(real) in value factory 2",
                          tut::test_result::result_type::ok );
    unittest::g_suiteProxy.evaluate(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }
};

//! Value factory for Charm++ tests
//! (not in Factory_common, so Charm++ can find it)
using ValueFactory = std::map< int, std::function< VBase() > >;

//! \brief Test if tk::recordCharmModel correctly registers a default child
//!   constructor in value factory
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both triggers a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 13 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "default Charm++ child ctor in val-factory 1" );

  tut::ValueFactory f;
  tk::recordCharmModel< tut::VBase, tut::CharmChild >( f, 2 );  // key: 2
  ensure_equals( "recorded 1 Charm++ child def-ctor in value factory",
                 f.size(), 1UL );

  const auto it = f.find( 2 );
  if ( it != f.end() )
    it->second();      // fire up Charm++ chare
  else
    fail( "cannot find key in factory" );
}

//! \brief Test if tk::recordCharmModel correctly registers a child constructor
//!   in value factory
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both triggers a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
//! \author J. Bakosi
template<> template<>
void Factory_object::test< 14 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm++ child ctor(real) in value factory 1" );

  tut::ValueFactory f;
  // key: 2, ctor arg: 23.2
  tk::recordCharmModel< tut::VBase, tut::CharmChild >( f, 2, 23.2 );
  ensure_equals( "recorded 1 Charm++ child ctor(real) in value factory",
                 f.size(), 1UL );

  const auto it = f.find( 2 );
  if ( it != f.end() )
    it->second();      // fire up Charm++ chare
  else
    fail( "cannot find key in factory" );
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif // test_Factory_h
