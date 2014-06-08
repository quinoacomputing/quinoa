//******************************************************************************
/*!
  \file      src/RNGTest/Battery.h
  \author    J. Bakosi
  \date      Fri 06 Jun 2014 11:48:21 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Battery
  \details   Battery
*/
//******************************************************************************
#ifndef Battery_h
#define Battery_h

#include <make_unique.h>
#include <RNGTestPrint.h>

namespace rngtest {

//! Battery. The class below uses runtime polymorphism without client-side
//! inheritance: inheritance is confined to the internals of the class below,
//! inivisble to client-code. The class exclusively deals with ownership
//! enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//! https://github.com/sean-parent/sean-parent.github.com/wiki/
//! Papers-and-Presentations
class Battery {

  public:
    //! Constructor taking and object modeling Concept (see below)
    template< typename T >
    explicit Battery( T x ) :
      self( tk::make_unique< Model<T> >( std::move(x) ) ) {}

    //! Constructor taking a function pointer to a constructor of an object
    //! modeling Concept (see below). Passing std::function allows late
    //! execution of the constructor, i.e., as late as inside this class'
    //! constructor, and thus usage from a factory.
    template< typename T >
    explicit Battery( std::function<T()> x ) :
      self( tk::make_unique< Model<T> >( std::move(x()) ) ) {}

    //! Public interface to running a battery of RNG tests
    void run() const { self->run(); }

    //! Public interface to printing list of registered statistical tests
    void print( const RNGTestPrint& p ) const { self->print(p); }

    //! Public interface to evaluating a statistical test
    void evaluate( std::size_t id ) const { self->evaluate(id); }

    //! Public interface to returning number of statistical tests in battery
    std::size_t ntest() const { return self->ntest(); }

    //! Public interface to returning number of statistics produced by battery
    std::size_t nstat() const { return self->nstat(); }

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
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void run() const = 0;
      virtual void print( const RNGTestPrint& p ) const = 0;
      virtual void evaluate( std::size_t id ) const = 0;
      virtual std::size_t ntest() const = 0;
      virtual std::size_t nstat() const = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const { return new Model( *this ); }
      void run() const override { data.run(); }
      void print( const RNGTestPrint& p ) const override { data.print(p); }
      void evaluate( std::size_t id ) const override { data.evaluate(id); }
      std::size_t ntest() const override { return data.ntest(); }
      std::size_t nstat() const override { return data.nstat(); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // rngtest::

#endif // Battery_h
