// *****************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Driver base class declaration
  \details   Driver base class declaration. This class used as a base for
     various drivers. Note that "base" is not in a classical OOP sense. See also
     in-line code documentation.
*/
// *****************************************************************************
#ifndef Driver_h
#define Driver_h

#include "Make_unique.h"

namespace tk {

//! Driver. The class below uses runtime polymorphism without client-side
//! inheritance: inheritance is confined to the internals of the class below,
//! inivisble to client-code. The class exclusively deals with ownership
//! enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//! https://github.com/sean-parent/sean-parent.github.com/wiki/
//! Papers-and-Presentations
class Driver {

  public:
    //! Constructor taking an object modeling Concept (see below)
    template< typename T >
    explicit Driver( T x ) : self( make_unique< Model<T> >( std::move(x) ) ) {}

    //! Constructor taking a function pointer to a constructor of an object
    //! modeling Concept (see below). Passing std::function allows late
    //! execution of the constructor, i.e., as late as inside this class'
    //! constructor, and thus usage from a factory.
    template< typename T >
    explicit Driver( std::function<T()> x ) :
      self( make_unique< Model<T> >( std::move(x()) ) ) {}

    //! Public interface to execute
    void execute() const { self->execute(); }

    //! Copy assignment
    Driver& operator=( const Driver& x )
    { Driver tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    Driver( const Driver& x ) : self( x.self->copy() ) {}
    //! Move assignment
    Driver& operator=( Driver&& ) noexcept = default;
    //! Move constructor
    Driver( Driver&& ) noexcept = default;

  private:
    //! Concept is a pure virtual base class specifying the requirements of
    //! polymorphic objects deriving from it
    struct Concept {
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void execute() const = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const { return new Model( *this ); }
      void execute() const override { data.execute(); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // tk::

#endif // Driver_h
