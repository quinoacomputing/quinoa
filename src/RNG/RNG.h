// *****************************************************************************
/*!
  \file      src/RNG/RNG.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator
  \details   This file defines a generic random number generator class. The
    class uses runtime polymorphism without client-side inheritance: inheritance
    is confined to the internals of the class, inivisble to client-code. The
    class exclusively deals with ownership enabling client-side value semantics.
    Credit goes to Sean Parent at Adobe: https://github.com/sean-parent/
    sean-parent.github.com/wiki/Papers-and-Presentations.
*/
// *****************************************************************************
#ifndef RNG_h
#define RNG_h

#include <functional>

#include "Make_unique.h"
#include "Keywords.h"

namespace tk {

//! \brief Random number generator
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   inivisble to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a RNG, see
//!   see tk::MKLRNG or tk::RNGSSE.
//! \author J. Bakosi
class RNG {

    using ncomp_t = kw::ncomp::info::expect::type;    

  public:
    //! \brief Constructor taking an object modeling Concept
    //! \details The object of class T comes pre-constructed.
    //! \param[in] x Instantiated object of type T given by the template
    //!   argument.
    template< typename T >
    explicit RNG( T x ) : self( make_unique< Model<T> >( std::move(x) ) ) {}

    //! \brief Constructor taking a function pointer to a constructor of an
    //!   object modeling Concept
    //! \details Passing std::function allows late execution of the constructor,
    //!   i.e., as late as inside this class' constructor, and thus usage from a
    //!   factory.
    //! \param[in] x Function pointer to a constructor of an object modeling
    //!   Concept
    template< typename T >
    explicit RNG( std::function<T()> x ) :
      self( make_unique< Model<T> >( std::move(x()) ) ) {}

    //! Public interface to uniform RNG
    void uniform( int stream, ncomp_t num, double* r ) const
    { self->uniform( stream, num, r ); }

    //! Public interface to Gaussian RNG
    void gaussian( int stream, ncomp_t num, double* r ) const
    { self->gaussian( stream, num, r ); }

    //! Public interface to beta RNG
    void beta( int stream, ncomp_t num, double p, double q, double a, double b,
               double* r ) const
    { self->beta( stream, num, p, q, a, b, r ); }

    //! Public interface to number of threads accessor
    std::size_t nthreads() const noexcept { return self->nthreads(); }

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
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void uniform( int, ncomp_t, double* ) const = 0;
      virtual void gaussian( int, ncomp_t, double* ) const = 0;
      virtual void beta( int, ncomp_t, double, double, double, double, double* )
      const = 0;
      virtual std::size_t nthreads() const noexcept = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      void uniform( int stream, ncomp_t num, double* r ) const override
      { data.uniform( stream, num, r ); }
      void gaussian( int stream, ncomp_t num, double* r ) const override
      { data.gaussian( stream, num, r ); }
      void beta( int stream, ncomp_t num, double p, double q, double a,
                 double b, double* r ) const override
      { data.beta( stream, num, p, q, a, b, r ); }
      std::size_t nthreads() const noexcept override { return data.nthreads(); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // namespace tk

#endif // RNG_h
