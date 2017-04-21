// *****************************************************************************
/*!
  \file      src/RNGTest/StatTest.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator statistical test
  \details   This file defines a generic random number generator statistical
    test class. The class uses runtime polymorphism without client-side
    inheritance: inheritance is confined to the internals of the class,
    inivisble to client-code. The class exclusively deals with ownership
    enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
    https://github.com/sean-parent/
    sean-parent.github.com/wiki/Papers-and-Presentations.
*/
// *****************************************************************************
#ifndef StatTest_h
#define StatTest_h

#include <functional>

#include "NoWarning/charm++.h"

#include "Macro.h"
#include "Make_unique.h"
#include "CharmUtil.h"
#include "Options/RNG.h"

namespace rngtest {

//! \brief Random number generator statistical test
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   inivisble to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a Battery,
//!   see rngtest::TestU01.
//! \author J. Bakosi
class StatTest {

  public:
    //! \brief Constructor taking an object modeling Concept
    //! \details The object of class T comes pre-constructed.
    //! \param[in] x Instantiated object of type T given by the template
    //!   argument.
    template< typename T >
    explicit StatTest( T x ) :
      self( tk::make_unique< Model<T> >( std::move(x) ) ){}

    //! \brief Constructor taking a std::function holding a constructor bound to
    //!   its arguments of an object modeling Concept (see below)
    //! \details Passing std::function allows late execution of the constructor
    //!   of T, i.e., as late as inside this class' constructor, and thus usage
    //!   from a factory. Object of T is constructed here. This overload is
    //!   disabled for Charm++ chare objects defining typedef 'Proxy', see also
    //!   below.
    //! \param[in] x Function pointer to a constructor bound of an object
    //!   modeling Concept
    template< typename T,
      typename std::enable_if< !tk::HasTypedefProxy<T>::value, int >::type = 0 >
    explicit StatTest( std::function<T()> x ) :
      self( tk::make_unique< Model<T> >( std::move(x()) ) ) {}

    //! \brief Constructor taking a function pointer to a constructor of an
    //!    object modeling Concept
    //! \details Passing std::function allows late execution of the constructor
    //!   of T, i.e., at some future time, and thus usage from a factory. Note
    //!   that the value of the first function argument, std::function<T()>, is
    //!   not used here, but its constructor type, T, is used to enable the
    //!   compiler to deduce the model constructor type, used to create its
    //!   Charm proxy, defined by T::Proxy. The actual constructor of T is not
    //!   called here but at some future time by the Charm++ runtime system,
    //!   here only an asynchrounous ckNew() is called, i.e., a message (or
    //!   request) for a future call to T's constructor. This overload is only
    //!   enabled for Charm++ chare objects defining typedef 'Proxy', which must
    //!   define the Charm++ proxy. All optional constructor arguments are
    //!   forwarded to ckNew() and thus to T's constructor. If it was somehow
    //!   possible to obtain all bound arguments' types and values from an
    //!   already-bound std::function, we could use those instead of having to
    //!   explicitly forward the model constructor arguments via this host
    //!   constructor.
    //! \param[in] c Function pointer to a constructor of an object modeling
    //!    Concept
    //! \param[in] args Constructor arguments
    //! \see See also tk::recordCharmModel().
    template< typename T, typename... ConstrArgs,
      typename std::enable_if< tk::HasTypedefProxy<T>::value, int >::type = 0 >
    explicit StatTest( std::function<T()> c, ConstrArgs... args ) :
      self( tk::make_unique< Model< typename T::Proxy > >
            (std::move(T::Proxy::ckNew(std::forward<ConstrArgs>(args)...))) ) {
      Assert( c == nullptr, "std::function argument to StatTest Charm "
                            "constructor must be nullptr" );
      #ifdef NDEBUG
      IGNORE(c);
      #endif
    }

    //! Public interface to contribute number of results/test, i.e., p-values
    void npval() const { self->npval(); }

    //! Public interface to contribute test name(s)
    void names() const { self->names(); }

    //! Public interface to running a test
    void run() const { self->run(); }

    //! Public interface to contributing a test's run time measured in seconds
    void time() const { self->time(); }

    //! Copy assignment
    StatTest& operator=( const StatTest& x )
    { StatTest tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    StatTest( const StatTest& x ) : self( x.self->copy() ) {}
    //! Move assignment
    StatTest& operator=( StatTest&& ) noexcept = default;
    //! Move constructor
    StatTest( StatTest&& ) noexcept = default;

  private:
    //! Concept is a pure virtual base class specifying the requirements of
    //! polymorphic objects deriving from it
    struct Concept {
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void npval() = 0;
      virtual void names() = 0;
      virtual void run() = 0;
      virtual void time() = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      void npval() override { data.npval(); }
      void names() override { data.names(); }
      void run() override { data.run(); }
      void time() override { data.time(); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // rngtest::

#endif // StatTest_h
