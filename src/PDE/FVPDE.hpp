// *****************************************************************************
/*!
  \file      src/PDE/FVPDE.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Partial differential equation base for finite volume PDEs
  \details   This file defines a generic partial differential equation (PDE)
    class for PDEs that use finite volume spatial discretization.
    The class uses runtime polymorphism without client-side inheritance:
    inheritance is confined to the internals of the class, invisible to
    client-code. The class exclusively deals with ownership enabling client-side
    value semantics. Credit goes to Sean Parent at Adobe:
    https://github.com/sean-parent/sean-parent.github.com/wiki/
    Papers-and-Presentations.
*/
// *****************************************************************************
#ifndef FVPDE_h
#define FVPDE_h

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
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "FunctionPrototypes.hpp"
#include "History.hpp"

namespace inciter {

extern ctr::New2InputDeck g_inputdeck;

using ncomp_t = kw::ncomp::info::expect::type;

//! \brief Partial differential equation base for discontinuous Galerkin PDEs
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   invisible to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a FVPDE,
//!   see inciter::CompFlow.
class FVPDE {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Default constructor taking no arguments for Charm++
    explicit FVPDE() = default;

    //! \brief Constructor taking an object modeling Concept.
    //! \details The object of class T comes pre-constructed.
    //! \param[in] x Instantiated object of type T given by the template
    //!   argument.
    template< typename T > explicit FVPDE( T x ) :
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
    explicit FVPDE( std::function<T(Args...)> x, Args&&... args ) :
      self( std::make_unique< Model<T> >(
              std::move( x( std::forward<Args>(args)... ) ) ) ) {}

    //! Public interface to find number of primitive quantities for the diff eq
    std::size_t nprim() const
    { return self->nprim(); }

    //! Public interface to find number of materials for the diff eq
    std::size_t nmat() const
    { return self->nmat(); }

    //! Public interface to determine elements that lie inside the IC box
    void IcBoxElems( const tk::Fields& geoElem,
      std::size_t nielem,
      std::vector< std::unordered_set< std::size_t > >& inbox ) const
    { self->IcBoxElems( geoElem, nielem, inbox ); }

    //! Public interface to setting the initial conditions for the diff eq
    void initialize(
      const tk::Fields& L,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const std::vector< std::unordered_set< std::size_t > >& inbox,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&
        elemblkid,
      tk::Fields& unk,
      tk::real t,
      const std::size_t nielem ) const
    { self->initialize( L, inpoel, coord, inbox, elemblkid, unk, t, nielem ); }

    //! Public interface to computing the left-hand side matrix for the diff eq
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const
    { self->lhs( geoElem, l ); }

    //! Public interface to updating the primitives for the diff eq
    void updatePrimitives( const tk::Fields& unk,
                           tk::Fields& prim,
                           std::size_t nielem ) const
    { self->updatePrimitives( unk, prim, nielem ); }

    //! Public interface to cleaning up trace materials for the diff eq
    void cleanTraceMaterial( tk::real t,
                             const tk::Fields& geoElem,
                             tk::Fields& unk,
                             tk::Fields& prim,
                             std::size_t nielem ) const
    { self->cleanTraceMaterial( t, geoElem, unk, prim, nielem ); }

    //! Public interface to reconstructing the second-order solution
    void reconstruct( const tk::Fields& geoElem,
                      const inciter::FaceData& fd,
                      const std::map< std::size_t, std::vector< std::size_t > >&
                        esup,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      tk::Fields& U,
                      tk::Fields& P ) const
    {
      self->reconstruct( geoElem, fd, esup, inpoel, coord, U, P );
    }

    //! Public interface to limiting the second-order solution
    void limit( const tk::Fields& geoFace,
                const inciter::FaceData& fd,
                const std::map< std::size_t, std::vector< std::size_t > >& esup,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< int >& srcFlag,
                tk::Fields& U,
                tk::Fields& P ) const
    {
      self->limit( geoFace, fd, esup, inpoel, coord, srcFlag, U, P );
    }

    //! Public interface to computing the P1 right-hand side vector
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const tk::UnsMesh::Coords& coord,
              const std::unordered_map< std::size_t, std::set< std::size_t > >&
                elemblkid,
              const tk::Fields& U,
              const tk::Fields& P,
              tk::Fields& R,
              std::vector< int >& srcFlag ) const
    {
      self->rhs( t, geoFace, geoElem, fd, inpoel, coord, elemblkid, U, P, R,
        srcFlag );
    }

    //! Public interface for computing the minimum time step size
    tk::real dt( const inciter::FaceData& fd,
                 const tk::Fields& geoFace,
                 const tk::Fields& geoElem,
                 const tk::Fields& U,
                 const tk::Fields& P,
                 const std::size_t nielem,
                 const std::vector< int >& srcFlag,
                 std::vector< tk::real >& local_dte ) const
    { return self->dt( fd, geoFace, geoElem, U, P, nielem, srcFlag, local_dte );
    }

    //! Public interface to returning analytic field output labels
    std::vector< std::string > analyticFieldNames() const
    { return self->analyticFieldNames(); }

    //! Public interface to returning surface output labels
    std::vector< std::string > surfNames() const { return self->surfNames(); }

    //! Public interface to returning time history field output labels
    std::vector< std::string > histNames() const { return self->histNames(); }

    //! Public interface to returning variable names
    std::vector< std::string > names() const { return self->names(); }

    //! Public interface to returning surface field output
    std::vector< std::vector< tk::real > >
    surfOutput( const inciter::FaceData& fd,
                const tk::Fields& U,
                const tk::Fields& P ) const
    { return self->surfOutput( fd, U, P ); }

    //! Public interface to return point history output
    std::vector< std::vector< tk::real > >
    histOutput( const std::vector< HistData >& h,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const tk::Fields& U,
      const tk::Fields& P ) const
    { return self->histOutput( h, inpoel, coord, U, P ); }

    //! Public interface to returning analytic solution
    tk::InitializeFn::result_type
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return self->analyticSolution( xi, yi, zi, t ); }

    //! Public interface to returning the analytic solution for conserved vars
    tk::InitializeFn::result_type
    solution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return self->solution( xi, yi, zi, t ); }

    //! Public interface to returning the specific total energy
    tk::real
    sp_totalenergy( std::size_t e, const tk::Fields& unk ) const
    { return self->sp_totalenergy( e, unk ); }

    //! Public interface to returning the relevant sound speed in each cell
    void
    soundspeed(
      std::size_t nielem,
      const tk::Fields& U,
      const tk::Fields& P,
      std::vector< tk::real >& ss ) const
    { return self->soundspeed( nielem, U, P, ss ); }

    //! Copy assignment
    FVPDE& operator=( const FVPDE& x )
    { FVPDE tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    FVPDE( const FVPDE& x ) : self( x.self->copy() ) {}
    //! Move assignment
    FVPDE& operator=( FVPDE&& ) noexcept = default;
    //! Move constructor
    FVPDE( FVPDE&& ) noexcept = default;

  private:
    //! \brief Concept is a pure virtual base class specifying the requirements
    //!   of polymorphic objects deriving from it
    struct Concept {
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual std::size_t nprim() const = 0;
      virtual std::size_t nmat() const = 0;
      virtual void IcBoxElems( const tk::Fields&,
        std::size_t,
        std::vector< std::unordered_set< std::size_t > >& ) const = 0;
      virtual void initialize(
        const tk::Fields&,
        const std::vector< std::size_t >&,
        const tk::UnsMesh::Coords&,
        const std::vector< std::unordered_set< std::size_t > >&,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&,
        tk::Fields&,
        tk::real,
        const std::size_t nielem ) const = 0;
      virtual void lhs( const tk::Fields&, tk::Fields& ) const = 0;
      virtual void updatePrimitives( const tk::Fields&,
                                     tk::Fields&,
                                     std::size_t ) const = 0;
      virtual void cleanTraceMaterial( tk::real,
                                       const tk::Fields&,
                                       tk::Fields&,
                                       tk::Fields&,
                                       std::size_t ) const = 0;
      virtual void reconstruct( const tk::Fields&,
                                const inciter::FaceData&,
                                const std::map< std::size_t,
                                  std::vector< std::size_t > >&,
                                const std::vector< std::size_t >&,
                                const tk::UnsMesh::Coords&,
                                tk::Fields&,
                                tk::Fields& ) const = 0;
      virtual void limit( const tk::Fields&,
                          const inciter::FaceData&,
                          const std::map< std::size_t,
                            std::vector< std::size_t > >&,
                          const std::vector< std::size_t >&,
                          const tk::UnsMesh::Coords&,
                          const std::vector< int >&,
                          tk::Fields&,
                          tk::Fields& ) const = 0;
      virtual void rhs( tk::real,
        const tk::Fields&,
        const tk::Fields&,
        const inciter::FaceData&,
        const std::vector< std::size_t >&,
        const tk::UnsMesh::Coords&,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&,
        const tk::Fields&,
        const tk::Fields&,
        tk::Fields&,
        std::vector< int >& ) const = 0;
      virtual tk::real dt( const inciter::FaceData&,
                           const tk::Fields&,
                           const tk::Fields&,
                           const tk::Fields&,
                           const tk::Fields&,
                           const std::size_t,
                           const std::vector< int >&,
                           std::vector< tk::real >& ) const = 0;
      virtual std::vector< std::string > analyticFieldNames() const = 0;
      virtual std::vector< std::string > surfNames() const = 0;
      virtual std::vector< std::string > histNames() const = 0;
      virtual std::vector< std::string > names() const = 0;
      virtual std::vector< std::vector< tk::real > > surfOutput(
        const inciter::FaceData&,
        const tk::Fields&,
        const tk::Fields& ) const = 0;
      virtual std::vector< std::vector< tk::real > > histOutput(
        const std::vector< HistData >&,
        const std::vector< std::size_t >&,
        const tk::UnsMesh::Coords&,
        const tk::Fields&,
        const tk::Fields& ) const = 0;
      virtual tk::InitializeFn::result_type analyticSolution(
        tk::real xi, tk::real yi, tk::real zi, tk::real t ) const = 0;
      virtual tk::InitializeFn::result_type solution(
        tk::real xi, tk::real yi, tk::real zi, tk::real t ) const = 0;
      virtual tk::real sp_totalenergy(
        std::size_t, const tk::Fields& ) const = 0;
      virtual void soundspeed(
        std::size_t,
        const tk::Fields&,
        const tk::Fields&,
        std::vector< tk::real >& ) const = 0;
    };

    //! \brief Model models the Concept above by deriving from it and overriding
    //!   the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      explicit Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      std::size_t nprim() const override
      { return data.nprim(); }
      std::size_t nmat() const override
      { return data.nmat(); }
      void IcBoxElems( const tk::Fields& geoElem,
        std::size_t nielem,
        std::vector< std::unordered_set< std::size_t > >& inbox )
      const override { data.IcBoxElems( geoElem, nielem, inbox ); }
      void initialize(
        const tk::Fields& L,
        const std::vector< std::size_t >& inpoel,
        const tk::UnsMesh::Coords& coord,
        const std::vector< std::unordered_set< std::size_t > >& inbox,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&
          elemblkid,
        tk::Fields& unk,
        tk::real t,
        const std::size_t nielem )
      const override { data.initialize( L, inpoel, coord, inbox, elemblkid, unk,
        t, nielem ); }
      void lhs( const tk::Fields& geoElem, tk::Fields& l ) const override
      { data.lhs( geoElem, l ); }
      void updatePrimitives( const tk::Fields& unk,
                             tk::Fields& prim,
                             std::size_t nielem )
      const override { data.updatePrimitives( unk, prim, nielem ); }
      void cleanTraceMaterial( tk::real t,
                               const tk::Fields& geoElem,
                               tk::Fields& unk,
                               tk::Fields& prim,
                               std::size_t nielem )
      const override { data.cleanTraceMaterial( t, geoElem, unk, prim, nielem ); }
      void reconstruct( const tk::Fields& geoElem,
                        const inciter::FaceData& fd,
                        const std::map< std::size_t,
                          std::vector< std::size_t > >& esup,
                        const std::vector< std::size_t >& inpoel,
                        const tk::UnsMesh::Coords& coord,
                        tk::Fields& U,
                        tk::Fields& P ) const override
      {
        data.reconstruct( geoElem, fd, esup, inpoel, coord, U, P );
      }
      void limit( const tk::Fields& geoFace,
                  const inciter::FaceData& fd,
                  const std::map< std::size_t, std::vector< std::size_t > >&
                    esup,
                  const std::vector< std::size_t >& inpoel,
                  const tk::UnsMesh::Coords& coord,
                  const std::vector< int >& srcFlag,
                  tk::Fields& U,
                  tk::Fields& P ) const override
      {
        data.limit( geoFace, fd, esup, inpoel, coord, srcFlag, U, P );
      }
      void rhs(
        tk::real t,
        const tk::Fields& geoFace,
        const tk::Fields& geoElem,
        const inciter::FaceData& fd,
        const std::vector< std::size_t >& inpoel,
        const tk::UnsMesh::Coords& coord,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&
          elemblkid,
        const tk::Fields& U,
        const tk::Fields& P,
        tk::Fields& R,
        std::vector< int >& srcFlag ) const override
      {
        data.rhs( t, geoFace, geoElem, fd, inpoel, coord, elemblkid, U, P, R,
          srcFlag );
      }
      tk::real dt( const inciter::FaceData& fd,
                   const tk::Fields& geoFace,
                   const tk::Fields& geoElem,
                   const tk::Fields& U,
                   const tk::Fields& P,
                   const std::size_t nielem,
                   const std::vector< int >& srcFlag,
                   std::vector< tk::real >& local_dte ) const override
      { return data.dt( fd, geoFace, geoElem, U, P, nielem, srcFlag,
          local_dte ); }
      std::vector< std::string > analyticFieldNames() const override
      { return data.analyticFieldNames(); }
      std::vector< std::string > surfNames() const override
      { return data.surfNames(); }
      std::vector< std::string > histNames() const override
      { return data.histNames(); }
      std::vector< std::string > names() const override
      { return data.names(); }
      std::vector< std::vector< tk::real > > surfOutput(
        const inciter::FaceData& fd,
        const tk::Fields& U,
        const tk::Fields& P ) const override
      { return data.surfOutput( fd, U, P ); }
      std::vector< std::vector< tk::real > > histOutput(
        const std::vector< HistData >& h,
        const std::vector< std::size_t >& inpoel,
        const tk::UnsMesh::Coords& coord,
        const tk::Fields& U,
        const tk::Fields& P ) const override
      { return data.histOutput( h, inpoel, coord, U, P ); }
      tk::InitializeFn::result_type
      analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t )
       const override { return data.analyticSolution( xi, yi, zi, t ); }
      tk::InitializeFn::result_type
      solution( tk::real xi, tk::real yi, tk::real zi, tk::real t )
       const override { return data.solution( xi, yi, zi, t ); }
      tk::real sp_totalenergy( std::size_t e, const tk::Fields& unk )
       const override { return data.sp_totalenergy( e, unk ); }
      void soundspeed(
        std::size_t nielem,
        const tk::Fields& U,
        const tk::Fields& P,
        std::vector< tk::real >& ss )
       const override { return data.soundspeed( nielem, U, P, ss ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // inciter::

#endif // FVPDE_h
