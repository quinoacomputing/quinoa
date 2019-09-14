// *****************************************************************************
/*!
  \file      src/PDE/DGPDE.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Partial differential equation base for discontinuous Galerkin PDEs
  \details   This file defines a generic partial differential equation (PDE)
    class for PDEs that use discontinuous Galerkin spatial discretization.
    The class uses runtime polymorphism without client-side inheritance:
    inheritance is confined to the internals of the class, invisible to
    client-code. The class exclusively deals with ownership enabling client-side
    value semantics. Credit goes to Sean Parent at Adobe:
    https://github.com/sean-parent/sean-parent.github.com/wiki/
    Papers-and-Presentations.
*/
// *****************************************************************************
#ifndef DGPDE_h
#define DGPDE_h

#include <array>
#include <string>
#include <vector>
#include <memory>
#include <unordered_set>
#include <functional>

#include "Types.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"

namespace inciter {

//! \brief Partial differential equation base for discontinuous Galerkin PDEs
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   invisible to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a DGPDE,
//!   see inciter::CompFlow.
class DGPDE {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Default constructor taking no arguments for Charm++
    explicit DGPDE() = default;

    //! \brief Constructor taking an object modeling Concept.
    //! \details The object of class T comes pre-constructed.
    //! \param[in] x Instantiated object of type T given by the template
    //!   argument.
    template< typename T > explicit DGPDE( T x ) :
      self( std::make_unique< Model<T> >( std::move(x) ) ) {}

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
    explicit DGPDE( std::function<T(Args...)> x, Args&&... args ) :
      self( std::make_unique< Model<T> >(
              std::move( x( std::forward<Args>(args)... ) ) ) ) {}

    //! Public interface to find number of primitive quantities for the diff eq
    std::size_t nprim() const
    { return self->nprim(); }

    //! Public interface to setting the initial conditions for the diff eq
    void initialize( const tk::Fields& L,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord,
                     tk::Fields& unk,
                     tk::real t,
                     const std::size_t nielem ) const
    { self->initialize( L, inpoel, coord, unk, t, nielem ); }

    //! Public interface to computing the left-hand side matrix for the diff eq
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const
    { self->lhs( geoElem, l ); }

    //! Public interface to updating the primitives for the diff eq
    void updatePrimitives( const tk::Fields& unk,
                           tk::Fields& prim,
                           std::size_t nielem ) const
    { self->updatePrimitives( unk, prim, nielem ); }

    //! Public interface to reconstructing the second-order solution
    void reconstruct( tk::real t,
                      const tk::Fields& geoFace,
                      const tk::Fields& geoElem,
                      const inciter::FaceData& fd,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      tk::Fields& U,
                      tk::Fields& P ) const
    {
      self->reconstruct( t, geoFace, geoElem, fd, inpoel, coord, U, P );
    }

    //! Public interface to limiting the second-order solution
    void limit( tk::real t,
                const tk::Fields& geoFace,
                const tk::Fields& geoElem,
                const inciter::FaceData& fd,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& ndofel,
                tk::Fields& U,
                tk::Fields& P ) const
    {
      self->limit( t, geoFace, geoElem, fd, inpoel, coord, ndofel, U, P );
    }

    //! Public interface to computing the P1 right-hand side vector
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              const tk::Fields& P,
              const std::vector< std::size_t >& ndofel,
              tk::Fields& R ) const
    {
      self->rhs( t, geoFace, geoElem, fd, inpoel, coord, U, P, ndofel, R );
    }

    //! Public interface for computing the minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const inciter::FaceData& fd,
                 const tk::Fields& geoFace,
                 const tk::Fields& geoElem,
                 const std::vector< std::size_t >& ndofel,
                 const tk::Fields& U ) const
    { return self->dt( coord, inpoel, fd, geoFace, geoElem, ndofel, U ); }

    //! \brief Public interface for collecting all side set IDs the user has
    //!   configured for all components of a PDE system
    void side( std::unordered_set< int >& conf ) const { self->side( conf ); }

    //! Public interface to returning field output labels
    std::vector< std::string > fieldNames() const { return self->fieldNames(); }

    //! Public interface to returning variable names
    std::vector< std::string > names() const { return self->names(); }

    //! Public interface to returning field output
    std::vector< std::vector< tk::real > > fieldOutput(
      tk::real t,
      const tk::Fields& geoElem,
      tk::Fields& U ) const
    { return self->fieldOutput( t, geoElem, U ); }

    //! Public interface to returning nodal field output
    std::vector< std::vector< tk::real > > avgElemToNode(
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const tk::Fields& geoElem,
      const tk::Fields& U ) const
    { return self->avgElemToNode( inpoel, coord, geoElem, U ); }

    //! Public interface to returning analytic solution
    std::vector< tk::real >
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return self->analyticSolution( xi, yi, zi, t ); }

    //! Copy assignment
    DGPDE& operator=( const DGPDE& x )
    { DGPDE tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    DGPDE( const DGPDE& x ) : self( x.self->copy() ) {}
    //! Move assignment
    DGPDE& operator=( DGPDE&& ) noexcept = default;
    //! Move constructor
    DGPDE( DGPDE&& ) noexcept = default;

  private:
    //! \brief Concept is a pure virtual base class specifying the requirements
    //!   of polymorphic objects deriving from it
    struct Concept {
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual std::size_t nprim() const = 0;
      virtual void initialize( const tk::Fields&,
                               const std::vector< std::size_t >&,
                               const tk::UnsMesh::Coords&,
                               tk::Fields&,
                               tk::real,
                               const std::size_t nielem ) const = 0;
      virtual void lhs( const tk::Fields&, tk::Fields& ) const = 0;
      virtual void updatePrimitives( const tk::Fields&,
                                     tk::Fields&,
                                     std::size_t ) const = 0;
      virtual void reconstruct( tk::real,
                                const tk::Fields&,
                                const tk::Fields&,
                                const inciter::FaceData&,
                                const std::vector< std::size_t >&,
                                const tk::UnsMesh::Coords&,
                                tk::Fields&,
                                tk::Fields& ) const = 0;
      virtual void limit( tk::real,
                          const tk::Fields&,
                          const tk::Fields&,
                          const inciter::FaceData&,
                          const std::vector< std::size_t >&,
                          const tk::UnsMesh::Coords&,
                          const std::vector< std::size_t >&,
                          tk::Fields&,
                          tk::Fields& ) const = 0;
      virtual void rhs( tk::real,
                        const tk::Fields&,
                        const tk::Fields&,
                        const inciter::FaceData&,
                        const std::vector< std::size_t >&,
                        const tk::UnsMesh::Coords&,
                        const tk::Fields&,
                        const tk::Fields&,
                        const std::vector< std::size_t >&,
                        tk::Fields& ) const = 0;
      virtual tk::real dt( const std::array< std::vector< tk::real >, 3 >&,
                           const std::vector< std::size_t >&,
                           const inciter::FaceData&,
                           const tk::Fields&,
                           const tk::Fields&,
                           const std::vector< std::size_t >&,
                           const tk::Fields& ) const = 0;
      virtual void side( std::unordered_set< int >& conf ) const = 0;
      virtual std::vector< std::string > fieldNames() const = 0;
      virtual std::vector< std::string > names() const = 0;
      virtual std::vector< std::vector< tk::real > > fieldOutput(
        tk::real,
        const tk::Fields&,
        tk::Fields& ) const = 0;
      virtual std::vector< std::vector< tk::real > > avgElemToNode(
        const std::vector< std::size_t >&,
        const tk::UnsMesh::Coords&,
        const tk::Fields&,
        const tk::Fields& ) const = 0;
      virtual std::vector< tk::real > analyticSolution(
        tk::real xi, tk::real yi, tk::real zi, tk::real t ) const = 0;
    };

    //! \brief Model models the Concept above by deriving from it and overriding
    //!   the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      explicit Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      std::size_t nprim() const override
      { return data.nprim(); }
      void initialize( const tk::Fields& L,
                       const std::vector< std::size_t >& inpoel,
                       const tk::UnsMesh::Coords& coord,
                       tk::Fields& unk,
                       tk::real t,
                       const std::size_t nielem )
      const override { data.initialize( L, inpoel, coord, unk, t, nielem ); }
      void lhs( const tk::Fields& geoElem, tk::Fields& l ) const override
      { data.lhs( geoElem, l ); }
      void updatePrimitives( const tk::Fields& unk,
                             tk::Fields& prim,
                             std::size_t nielem )
      const override { data.updatePrimitives( unk, prim, nielem ); }
      void reconstruct( tk::real t,
                        const tk::Fields& geoFace,
                        const tk::Fields& geoElem,
                        const inciter::FaceData& fd,
                        const std::vector< std::size_t >& inpoel,
                        const tk::UnsMesh::Coords& coord,
                        tk::Fields& U,
                        tk::Fields& P ) const override
      {
        data.reconstruct( t, geoFace, geoElem, fd, inpoel, coord, U, P );
      }
      void limit( tk::real t,
                  const tk::Fields& geoFace,
                  const tk::Fields& geoElem,
                  const inciter::FaceData& fd,
                  const std::vector< std::size_t >& inpoel,
                  const tk::UnsMesh::Coords& coord,
                  const std::vector< std::size_t >& ndofel,
                  tk::Fields& U,
                  tk::Fields& P ) const override
      {
        data.limit( t, geoFace, geoElem, fd, inpoel, coord, ndofel, U, P );
      }
      void rhs( tk::real t,
                const tk::Fields& geoFace,
                const tk::Fields& geoElem,
                const inciter::FaceData& fd,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const tk::Fields& U,
                const tk::Fields& P,
                const std::vector< std::size_t >& ndofel,
                tk::Fields& R ) const override
      {
        data.rhs( t, geoFace, geoElem, fd, inpoel, coord, U, P, ndofel, R );
      }
      tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                   const std::vector< std::size_t >& inpoel,
                   const inciter::FaceData& fd,
                   const tk::Fields& geoFace,
                   const tk::Fields& geoElem,
                   const std::vector< std::size_t >& ndofel,
                   const tk::Fields& U ) const override
      { return data.dt( coord, inpoel, fd, geoFace, geoElem, ndofel, U ); }
      void side( std::unordered_set< int >& conf ) const override
      { data.side( conf ); }
      std::vector< std::string > fieldNames() const override
      { return data.fieldNames(); }
      std::vector< std::string > names() const override
      { return data.names(); }
      std::vector< std::vector< tk::real > > fieldOutput(
        tk::real t,
        const tk::Fields& geoElem,
        tk::Fields& U ) const override
      { return data.fieldOutput( t, geoElem, U ); }
      std::vector< std::vector< tk::real > > avgElemToNode(
        const std::vector< std::size_t >& inpoel,
        const tk::UnsMesh::Coords& coord,
        const tk::Fields& geoElem,
        const tk::Fields& U ) const override
      { return data.avgElemToNode( inpoel, coord, geoElem, U ); }
      std::vector< tk::real >
      analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t )
       const override { return data.analyticSolution( xi, yi, zi, t ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // inciter::

#endif // DGPDE_h
