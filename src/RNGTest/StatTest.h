//******************************************************************************
/*!
  \file      src/RNGTest/StatTest.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 06:33:51 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical test base
  \details   Statistical test base
*/
//******************************************************************************
#ifndef StatTest_h
#define StatTest_h

#include <functional>

#include <make_unique.h>
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
    //! Constructor taking an object modeling Concept (see below)
    template< typename T > explicit StatTest( T x ) :
      self( tk::make_unique< Model<T> >( std::move(x) ) ){}

    //! Constructor taking a function pointer to a constructor of an object
    //! modeling Concept (see below). Passing std::function allows late
    //! execution of the constructor, i.e., as late as inside this class'
    //! constructor, and thus usage from a factory.
    template< typename T >
    explicit StatTest( std::function<T()> x ) :
      self( tk::make_unique< Model<T> >( std::move(x()) ) ) {}

    //! Public interface to running the test
    void run( std::size_t id ) const { self->run(id); }

    //! Public interface to accessing test name
    const std::string& name( std::size_t i ) const { return self->name(i); }

    //! Public interface to accessing the number of results/test
    std::size_t nstat() const { return self->nstat(); }

    //! Public interface to accessing the RNG enum
    const tk::ctr::RNGType& rng() const { return self->rng(); }

    //! Public interface to accessing the RNG id
    std::size_t id() const { return self->id(); }

    //! Public interface to querying whether the test is failed
    bool fail( std::size_t p ) const { return self->fail(p); }

    //! Public interface to accessing the number of failed tests
    std::size_t nfail() const { return self->nfail(); }

    //! Public interface to accessing p-value as double
    double pval( std::size_t p ) const { return self->pval(p); }

    //! Public interface to accessing p-value as std::string
    std::string pvalstr( std::size_t p ) const { return self->pvalstr(p); }

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
      virtual const std::string& name( std::size_t i ) const =  0;
      virtual std::size_t nstat() const = 0;
      virtual const tk::ctr::RNGType& rng() const = 0;
      virtual std::size_t id() const = 0;
      virtual bool fail( std::size_t p ) const = 0;
      virtual std::size_t nfail() const = 0;
      virtual double pval( std::size_t p ) const = 0;
      virtual std::string pvalstr( std::size_t p ) const = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const { return new Model( *this ); }
      void run( std::size_t id ) override { data.run(id); }
      const std::string& name( std::size_t i ) const override {
        return data.name(i); }
      std::size_t nstat() const override { return data.nstat(); }
      const tk::ctr::RNGType& rng() const override { return data.rng(); }
      std::size_t id() const override { return data.id(); }
      bool fail( std::size_t p ) const override { return data.fail(p); }
      std::size_t nfail() const override { return data.nfail(); }
      double pval( std::size_t p ) const override { return data.pval(p); }
      std::string pvalstr( std::size_t p ) const override {
        return data.pvalstr(p); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // rngtest::

#endif // StatTest_h
