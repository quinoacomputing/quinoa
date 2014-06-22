//******************************************************************************
/*!
  \file      src/RNGTest/StatTest.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 05:24:24 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical test base
  \details   Statistical test base
*/
//******************************************************************************
#ifndef StatTest_h
#define StatTest_h

#include <functional>

#include <charm++.h>

#include <make_unique.h>
#include <CharmUtil.h>
#include <Options/RNG.h>

namespace rngtest {

//! Statistical test. The class below uses runtime polymorphism without
//! client-side inheritance: inheritance is confined to the internals of the
//! class below, inivisble to client-code. The class exclusively deals with
//! ownership enabling client-side value semantics. Credit goes to Sean Parent
//! at Adobe: https://github.com/sean-parent/sean-parent.github.com/wiki/
//! Papers-and-Presentations
class StatTest {

  public:
    //! Constructor taking an object modeling Concept (see below). The object
    //! of class T was pre-constructed.
    template< typename T >
    explicit StatTest( T x ) :
      self( tk::make_unique< Model<T> >( std::move(x) ) ){}

    //! Constructor taking a std::function holding a constructor bound to its
    //! arguments of an object modeling Concept (see below). Passing
    //! std::function allows late execution of the constructor of T, i.e., as
    //! late as inside this class' constructor, and thus usage from a factory.
    //! Object of T is constructed here. This overload is disabled for Charm++
    //! chare objects defining typedef 'Proxy', see also below.
    template< typename T,
      typename std::enable_if< !tk::HasProxy<T>::value, int >::type = 0 >
    explicit StatTest( std::function<T()> x ) :
      self( tk::make_unique< Model<T> >( std::move(x()) ) ) {}

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
      typename std::enable_if< tk::HasProxy<T>::value, int >::type = 0 >
    explicit StatTest( std::function<T()> c, ConstrArgs... args ) :
      self( tk::make_unique< Model< typename T::Proxy > >
            (std::move(T::Proxy::ckNew(std::forward<ConstrArgs>(args)...))) ) {
      Assert( c == nullptr, tk::ExceptType::FATAL,
             "std::function arg to StatTest Charm constructor must be nullptr" );
    }

    //! Public interface to running the test
    void run( std::size_t id ) const { self->run(id); }

    //! Public interface to contribute number of results/test, i.e., p-values
    void npval( CkFuture f ) const { self->npval(f); }

    //! Public interface to test name
    void name( CkFuture f ) const { self->name(f); }

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
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void run( std::size_t id ) = 0;
      virtual void npval( CkFuture f ) = 0;
      virtual void name( CkFuture f ) = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const { return new Model( *this ); }
      void run( std::size_t id ) override { data.run(id); }
      void npval( CkFuture f ) override { data.npval(f); }
      void name( CkFuture f ) override { data.name(f); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // rngtest::

#endif // StatTest_h
