//******************************************************************************
/*!
  \file      src/DiffEq/DiffEq.h
  \author    J. Bakosi
  \date      Mon 18 Aug 2014 03:34:14 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Differential equation
  \details   Differential equation
*/
//******************************************************************************
#ifndef DiffEq_h
#define DiffEq_h

#include <string>
#include <functional>

#include <Types.h>
#include <make_unique.h>
#include <ParticleProperties.h>

namespace quinoa {

//! Differential equation. The class below uses runtime polymorphism without
//! client-side inheritance: inheritance is confined to the internals of the
//! class below, inivisble to client-code. The class exclusively deals with
//! ownership enabling client-side value semantics. Credit goes to Sean Parent
//! at Adobe: https://github.com/sean-parent/sean-parent.github.com/wiki/
//! Papers-and-Presentations
class DiffEq {

  public:
    //! Constructor taking an object modeling Concept (see below). The object
    //! of class T comes pre-constructed.
    template< typename T > explicit DiffEq( T x ) :
      self( tk::make_unique< Model<T> >( std::move(x) ) ) {}

    //! Constructor taking a function pointer to a constructor of an object
    //! modeling Concept (see below). Passing std::function allows late
    //! execution of the constructor, i.e., as late as inside this class'
    //! constructor, and thus usage from a factory. Note that there are at least
    //! two different ways of using this constructor: (1) Bind T's constructor
    //! arguments and place it in std::function<T()> and passing no arguments as
    //! args.... This case then instantiates the model via its constructor and
    //! stores in here. (2) Bind a single placeholder argument to T's
    //! constructor and pass it in as host's args..., which then forwards it to
    //! model's constructor. This allows late binding, i.e., binding the
    //! argument only here. See also the wrapper recordModel, which does (1),
    //! and recordModelLate, which does (2) defined in Base/Factory.h.
    template< typename T, typename...Args >
    explicit DiffEq( std::function<T(Args...)> x, Args... args ) :
      self( tk::make_unique< Model<T> >( std::move(x(args...)) ) ) {}

    //! Public interface to setting the initial conditions
    void initialize( ParProps& particles ) const
    { self->initialize( particles ); }

    //! Public interface to advancing particles in time
    void advance( ParProps& particles, int stream, tk::real dt ) const
    { self->advance( particles, stream, dt ); }

    //! Copy assignment
    DiffEq& operator=( const DiffEq& x )
    { DiffEq tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    DiffEq( const DiffEq& x ) : self( x.self->copy() ) {}
    //! Move assignment
    DiffEq& operator=( DiffEq&& ) noexcept = default;
    //! Move constructor
    DiffEq( DiffEq&& ) noexcept = default;

  private:
    //! Concept is a pure virtual base class specifying the requirements of
    //! polymorphic objects deriving from it
    struct Concept {
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void initialize( ParProps& ) const = 0;
      virtual void advance( ParProps&, int, tk::real ) const = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const { return new Model( *this ); }
      void initialize( ParProps& particles ) const override
      { data.initialize( particles ); }
      void advance( ParProps& particles, int stream, tk::real dt ) const
      override { data.advance( particles, stream, dt ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // quinoa::

#endif // DiffEq_h
