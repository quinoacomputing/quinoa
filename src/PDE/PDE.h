//******************************************************************************
/*!
  \file      src/PDE/PDE.h
  \author    J. Bakosi
  \date      Wed 03 Feb 2016 03:05:20 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Partial differential equation
  \details   This file defines a generic partial differential equation class.
    The class uses runtime polymorphism without client-side inheritance:
    inheritance is confined to the internals of the class, inivisble to
    client-code. The class exclusively deals with ownership enabling client-side
    value semantics. Credit goes to Sean Parent at Adobe:
    https://github.com/sean-parent/sean-parent.github.com/wiki/
    Papers-and-Presentations.
*/
//******************************************************************************
#ifndef PDE_h
#define PDE_h

#include <string>
#include <functional>

#include "Types.h"
#include "Make_unique.h"
#include "MeshNodes.h"

namespace inciter {

//! \brief Partial differential equation
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   inivisble to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a PDE,
//!   see inciter::Euler.
//! \author J. Bakosi
class PDE {

  public:
    //! \brief Constructor taking an object modeling Concept.
    //! \details The object of class T comes pre-constructed.
    //! \param[in] x Instantiated object of type T given by the template
    //!   argument.
    template< typename T > explicit PDE( T x ) :
      self( tk::make_unique< Model<T> >( std::move(x) ) ) {}

    //! \brief Constructor taking a function pointer to a constructor of an
    //!   object modeling Concept.
    //! \details Passing std::function allows late execution of the constructor,
    //!   i.e., as late as inside this class' constructor, and thus usage from
    //!   a factory. Note that there are at least two different ways of using
    //!   this constructor:
    //!   - Bind T's constructor arguments and place it in std::function<T()>
    //!   and passing no arguments as args.... This case then instantiates the
    //!   model via its constructor and stores it in here.
    //!   - Bind a single placeholder argument to T's constructor and pass it in
    //!   as host's args..., which then forwards it to model's constructor. This
    //!   allows late binding, i.e., binding the argument only here.
    //! \see See also the wrapper tk::recordModel() which does the former and
    //!   tk::recordModelLate() which does the latter, both defined in
    //!   src/Base/Factory.h.
    //! \param[in] x Function pointer to a constructor of an object modeling
    //!    Concept.
    //! \param[in] args Zero or more constructor arguments
    template< typename T, typename...Args >
    explicit PDE( std::function<T(Args...)> x, Args... args ) :
      self( tk::make_unique< Model<T> >( std::move(x(args...)) ) ) {}

    //! Public interface to setting the initial conditions for the diff eq
    void initialize( tk::MeshNodes& unk ) const
    { self->initialize( unk ); }

    //! Public interface to advancing the PDE in time
    void advance( tk::MeshNodes& unk, tk::real dt, tk::real t ) const
    { self->advance( unk, dt, t ); }

    //! Copy assignment
    PDE& operator=( const PDE& x )
    { PDE tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    PDE( const PDE& x ) : self( x.self->copy() ) {}
    //! Move assignment
    PDE& operator=( PDE&& ) noexcept = default;
    //! Move constructor
    PDE( PDE&& ) noexcept = default;

  private:
    //! \brief Concept is a pure virtual base class specifying the requirements
    //!   of polymorphic objects deriving from it
    struct Concept {
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void initialize( tk::MeshNodes& ) = 0;
      virtual void advance( tk::MeshNodes&, tk::real, tk::real ) = 0;
    };

    //! \brief Model models the Concept above by deriving from it and overriding
    //!   the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      void initialize( tk::MeshNodes& unk ) override { data.initialize( unk ); }
      void advance( tk::MeshNodes& unk, tk::real dt, tk::real t ) override
      { data.advance( unk, dt, t ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // inciter::

#endif // PDE_h
