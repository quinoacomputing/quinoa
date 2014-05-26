//******************************************************************************
/*!
  \file      src/RNG/RNG.h
  \author    J. Bakosi
  \date      Mon 26 May 2014 04:35:44 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator base
  \details   Random number generator base
*/
//******************************************************************************
#ifndef RNG_h
#define RNG_h

#include <make_unique.h>
#include <Types.h>

namespace tk {

//! Random number generator. The class below enables runtime polymorphism
//! without client-side inheritance: inheritance is confined to the internals of
//! the class below, inivisble to client-code. Credit goes to Sean Parent at
//! Adobe: https://github.com/sean-parent/sean-parent.github.com/wiki/
//! Papers-and-Presentations
class RNG {

  public:
    //! Constructor taking an object modeling Concept (see below)
    template< typename T >
    explicit RNG( T x ) : self( tk::make_unique< Model<T> >( std::move(x) ) ) {}

    //! Constructor taking a function pointer to a constructor of an object
    //! modeling Concept (see below). Passing std::function allows late
    //! execution and thus usage from a factory.
    template< typename T >
    explicit RNG( std::function<T*()> x ) :
      self( tk::make_unique< Model<T*> >( std::move(x()) ) ) {}

    //! Public interface to uniform RNG
    void uniform( int tid, int num, double* r ) const {
      self->uniform( tid, num, r );
    }

    //! Public interface to Gaussian RNG
    void gaussian( int tid, int num, double* r ) const {
      self->gaussian( tid, num, r );
    }

  private:
    //! Don't permit copy constructor
    RNG(const RNG&) = delete;
    //! Don't permit copy assigment
    RNG& operator=(const RNG&) = delete;
    //! Don't permit move constructor
    RNG(RNG&&) = delete;
    //! Don't permit move assigment
    RNG& operator=(RNG&&) = delete;

    //! Concept is a pure virtual base class specifying the requirements of
    //! polymorphic objects deriving from it
    struct Concept {
      virtual ~Concept() = default;
      virtual void uniform( int, int, double* ) const = 0;
      virtual void gaussian( int, int, double* ) const = 0;
    };

    //! Model models the Concept above by deriving from it and overriding the
    //! the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      void uniform( int tid, int num, double* r ) const override {
        data->uniform( tid, num, r );
      }
      void gaussian( int tid, int num, double* r ) const override {
        data->gaussian( tid, num, r );
      }
      T data;
    };

    std::unique_ptr< Concept > self;
};

} // namespace tk

#endif // RNG_h
