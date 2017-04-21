// *****************************************************************************
/*!
  \file      src/RNGTest/Battery.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator test harness
  \details   This file defines a generic random number generator test harness
    class. The class uses runtime polymorphism without client-side inheritance:
    inheritance is confined to the internals of the class, inivisble to
    client-code. The class exclusively deals with ownership enabling client-side
    value semantics. Credit goes to Sean Parent at Adobe:
    https://github.com/sean-parent/
    sean-parent.github.com/wiki/Papers-and-Presentations.
*/
// *****************************************************************************
#ifndef Battery_h
#define Battery_h

#include <functional>

#include "NoWarning/charm++.h"

#include "Macro.h"
#include "Has.h"
#include "Make_unique.h"
#include "CharmUtil.h"

namespace rngtest {

//! \brief Battery
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   inivisble to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a Battery,
//!   see rngtest::TestU01Suite.
//! \author J. Bakosi
class Battery {

  public:
    //! \brief Constructor taking an object modeling Concept
    //! \details The object of class T comes pre-constructed.
    //! \param[in] x Instantiated object of type T given by the template
    //!   argument.
    template< typename T > explicit Battery( T x ) :
      self( tk::make_unique< Model<T> >( std::move(x) ) ) {}

    //! \brief Constructor taking a std::function holding a constructor bound to
    //!   its arguments of an object modeling Concept.
    //! \details Passing std::function allows late execution of the constructor
    //!   of T (given by the template argument), i.e., as late as inside this
    //!   class' constructor, and thus usage from a factory. Object of T is
    //!   constructed here. This overload is disabled for Charm++ chare objects
    //!   defining typedef 'Proxy', see also below.
    template< typename T,
      typename std::enable_if< !tk::HasTypedefProxy<T>::value, int >::type = 0 >
    explicit Battery( std::function<T()> x ) :
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
    //!    Concept.
    //! \param[in] args Constructor arguments
    //! \see See also tk::recordCharmModel().
    template< typename T, typename... ConstrArgs,
      typename std::enable_if< tk::HasTypedefProxy<T>::value, int >::type = 0 >
    explicit Battery( std::function<T()> c, ConstrArgs... args ) :
      self( tk::make_unique< Model< typename T::Proxy > >
            (std::move(T::Proxy::ckNew(std::forward<ConstrArgs>(args)...))) ) {
      Assert( c == nullptr, "std::function argument to Battery Charm++ "
                            "constructor must be nullptr" );
      #ifdef NDEBUG
      IGNORE(c);
      #endif
    }

    //! Public interface to evaluating a statistical test
    void evaluate( std::vector< std::vector< std::string > > status ) const
    { self->evaluate( status ); }

    //! Public interface to collecting the number of statistics from a test
    void npval( std::size_t n ) const { self->npval( n ); }

    //! Public interface to collecting test name(s) from a test
    void names( std::vector< std::string > n ) const { self->names( n ); }

    //! Copy assignment
    Battery& operator=( const Battery& x )
    { Battery tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    Battery( const Battery& x ) : self( x.self->copy() ) {}
    //! Move assignment
    Battery& operator=( Battery&& ) noexcept = default;
    //! Move constructor
    Battery( Battery&& ) noexcept = default;

  private:
    //! Concept is a pure virtual base class specifying the requirements of
    //! polymorphic objects deriving from it
    struct Concept {
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void evaluate( std::vector< std::vector< std::string > >
                               status ) = 0;
      virtual void npval( std::size_t n ) = 0;
      virtual void names( std::vector< std::string > n ) = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      void evaluate( std::vector< std::vector< std::string > > status )
        override { data.evaluate( status ); }
      void npval( std::size_t n ) override { data.npval( n ); }
      void names( std::vector< std::string > n ) override { data.names( n ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // rngtest::

#endif // Battery_h
