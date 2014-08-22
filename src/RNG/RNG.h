//******************************************************************************
/*!
  \file      src/RNG/RNG.h
  \author    J. Bakosi
  \date      Tue 12 Aug 2014 05:04:15 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Random number generator base
  \details   Random number generator base
*/
//******************************************************************************
#ifndef RNG_h
#define RNG_h

#include <functional>

#include <make_unique.h>

namespace tk {

//! Random number generator. The class below uses runtime polymorphism without
//! client-side inheritance: inheritance is confined to the internals of the
//! class below, inivisble to client-code. The class exclusively deals with
//! ownership enabling client-side value semantics. Credit goes to Sean Parent
//! at Adobe: https://github.com/sean-parent/sean-parent.github.com/wiki/
//! Papers-and-Presentations
class RNG {

  public:
    //! Constructor taking an object modeling Concept (see below)
    template< typename T >
    explicit RNG( T x ) : self( make_unique< Model<T> >( std::move(x) ) ) {}

    //! Constructor taking a function pointer to a constructor of an object
    //! modeling Concept (see below). Passing std::function allows late
    //! execution of the constructor, i.e., as late as inside this class'
    //! constructor, and thus usage from a factory.
    template< typename T >
    explicit RNG( std::function<T()> x ) :
      self( make_unique< Model<T> >( std::move(x()) ) ) {}

    //! Public interface to uniform RNG
    void uniform( int stream, int num, double* r ) const
    { self->uniform( stream, num, r ); }

    //! Public interface to Gaussian RNG
    void gaussian( int stream, int num, double* r ) const
    { self->gaussian( stream, num, r ); }

    //! Copy assignment
    RNG& operator=( const RNG& x )
    { RNG tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    RNG( const RNG& x ) : self( x.self->copy() ) {}
    //! Move assignment
    RNG& operator=( RNG&& ) noexcept = default;
    //! Move constructor
    RNG( RNG&& ) noexcept = default;

  private:
    //! Concept is a pure virtual base class specifying the requirements of
    //! polymorphic objects deriving from it
    struct Concept {
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void uniform( int, int, double* ) const = 0;
      virtual void gaussian( int, int, double* ) const = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const { return new Model( *this ); }
      void uniform( int stream, int num, double* r ) const override
      { data.uniform( stream, num, r ); }
      void gaussian( int stream, int num, double* r ) const override
      { data.gaussian( stream, num, r ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // namespace tk

#endif // RNG_h
