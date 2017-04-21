// *****************************************************************************
/*!
  \file      src/PDE/PDE.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Partial differential equation
  \details   This file defines a generic partial differential equation class.
    The class uses runtime polymorphism without client-side inheritance:
    inheritance is confined to the internals of the class, inivisble to
    client-code. The class exclusively deals with ownership enabling client-side
    value semantics. Credit goes to Sean Parent at Adobe:
    https://github.com/sean-parent/sean-parent.github.com/wiki/
    Papers-and-Presentations.
*/
// *****************************************************************************
#ifndef PDE_h
#define PDE_h

#include <array>
#include <string>
#include <vector>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include "Types.h"
#include "Make_unique.h"
#include "Fields.h"

namespace inciter {

//! \brief Partial differential equation
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   inivisble to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a PDE,
//!   see inciter::CompFlow.
//! \author J. Bakosi
class PDE {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

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
    explicit PDE( std::function<T(Args...)> x, Args&&... args ) :
      self( tk::make_unique< Model<T> >(
              std::move( x( std::forward<Args>(args)... ) ) ) ) {}

    //! Public interface to setting the initial conditions for the diff eq
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::Fields& unk,
                     tk::real t,
                     const std::vector< std::size_t >& gid,
                     const std::unordered_map< std::size_t,
                            std::vector< std::pair< bool, tk::real > > >& bc )
    const { self->initialize( coord, unk, t, gid, bc ); }

    //! Public interface to computing the left-hand side matrix for the diff eq
    void lhs( const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& psup,
              tk::Fields& lhsd,
              tk::Fields& lhso ) const
    { self->lhs( coord, inpoel, psup, lhsd, lhso ); }

    //! Public interface to computing the right-hand side vector for the diff eq
    void rhs( tk::real mult,
              tk::real deltat,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const tk::Fields& U,
              tk::Fields& R ) const
    { self->rhs( mult, deltat, coord, inpoel, U, R ); }

    //! Public interface for computing the minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const tk::Fields& U ) const
    { return self->dt( coord, inpoel, U ); }

    //! Public interface for computing the transport velocity at nodes
    std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields& U,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::array< std::size_t, 4 >& N ) const
    { return self->velocity( U, coord, N ); }

    //! \brief Public interface for collecting all side set IDs the user has
    //!   configured for all components of a PDE system
    void side( std::unordered_set< int >& conf ) const { self->side( conf ); }

    //! \brief Public interface for querying Dirichlet boundary condition values
    //!  set by the user on a given side set for all components in a PDE system
    std::vector< std::pair< bool, tk::real > >
    dirbc( int sideset ) const { return self->dirbc( sideset ); }

    //! Public interface to returning field output labels
    std::vector< std::string > names() const { return self->names(); }

    //! Public interface to returning field output
    std::vector< std::vector< tk::real > > output(
      tk::real t,
      const std::array< std::vector< tk::real >, 3 >& coord,
      tk::Fields& U ) const
    { return self->output( t, coord, U ); }

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
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void initialize( const std::array< std::vector< tk::real >, 3 >&,
                               tk::Fields&,
                               tk::real,
                               const std::vector< std::size_t >&,
                               const std::unordered_map< std::size_t,
                                 std::vector< std::pair< bool, tk::real > > >& )
        const = 0;
      virtual void lhs( const std::array< std::vector< tk::real >, 3 >&,
                        const std::vector< std::size_t >&,
                        const std::pair< std::vector< std::size_t >,
                                         std::vector< std::size_t > >&,
                        tk::Fields&, tk::Fields& ) const = 0;
      virtual void rhs( tk::real,
                        tk::real,
                        const std::array< std::vector< tk::real >, 3 >&,
                        const std::vector< std::size_t >&,
                        const tk::Fields&,
                        tk::Fields& ) const = 0;
      virtual tk::real dt( const std::array< std::vector< tk::real >, 3 >&,
                           const std::vector< std::size_t >&,
                           const tk::Fields& ) const = 0;
      virtual std::array< std::array< tk::real, 4 >, 3 > velocity(
        const tk::Fields&,
        const std::array< std::vector< tk::real >, 3 >&,
        const std::array< std::size_t, 4 >&  ) const = 0;
      virtual void side( std::unordered_set< int >& conf ) const = 0;
      virtual std::vector< std::pair< bool, tk::real > > dirbc( int ) const = 0;
      virtual std::vector< std::string > names() const = 0;
      virtual std::vector< std::vector< tk::real > > output(
        tk::real,
        const std::array< std::vector< tk::real >, 3 >&,
        tk::Fields& ) const = 0;
    };

    //! \brief Model models the Concept above by deriving from it and overriding
    //!   the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                       tk::Fields& unk,
                       tk::real t,
                       const std::vector< std::size_t >& gid,
                       const std::unordered_map< std::size_t,
                             std::vector< std::pair< bool, tk::real > > >& bc )
      const override { data.initialize( coord, unk, t, gid, bc ); }
      void lhs( const std::array< std::vector< tk::real >, 3 >& coord,
                const std::vector< std::size_t >& inpoel,
                const std::pair< std::vector< std::size_t >,
                                 std::vector< std::size_t > >& psup,
                tk::Fields& lhsd, tk::Fields& lhso ) const override
      { data.lhs( coord, inpoel, psup, lhsd, lhso ); }
      void rhs( tk::real mult,
                tk::real deltat,
                const std::array< std::vector< tk::real >, 3 >& coord,
                const std::vector< std::size_t >& inpoel,
                const tk::Fields& U,
                tk::Fields& R ) const override
      { data.rhs( mult, deltat, coord, inpoel, U, R ); }
      tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                   const std::vector< std::size_t >& inpoel,
                   const tk::Fields& U ) const override
      { return data.dt( coord, inpoel, U ); }
      std::array< std::array< tk::real, 4 >, 3 > velocity(
        const tk::Fields& U,
        const std::array< std::vector< tk::real >, 3 >& coord,
        const std::array< std::size_t, 4 >& N  ) const override
      { return data.velocity( U, coord, N ); }
      void side( std::unordered_set< int >& conf ) const override
      { data.side( conf ); }
      std::vector< std::pair< bool, tk::real > > dirbc( int sideset ) const
      override { return data.dirbc( sideset ); }
      std::vector< std::string > names() const override { return data.names(); }
      std::vector< std::vector< tk::real > > output(
        tk::real t,
        const std::array< std::vector< tk::real >, 3 >& coord,
        tk::Fields& U ) const override { return data.output( t, coord, U ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // inciter::

#endif // PDE_h
