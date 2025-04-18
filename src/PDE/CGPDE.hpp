// *****************************************************************************
/*!
  \file      src/PDE/CGPDE.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Partial differential equation base for continuous Galerkin PDEs
  \details   This file defines a generic partial differential equation (PDE)
    class for PDEs that use continuous Galerkin spatial discretization.
    The class uses runtime polymorphism without client-side inheritance:
    inheritance is confined to the internals of the class, invisible to
    client-code. The class exclusively deals with ownership enabling client-side
    value semantics. Credit goes to Sean Parent at Adobe:
    https://github.com/sean-parent/sean-parent.github.com/wiki/
    Papers-and-Presentations.
*/
// *****************************************************************************
#ifndef CGPDE_h
#define CGPDE_h

#include <array>
#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <unordered_set>
#include <unordered_map>

#include "Types.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "FunctionPrototypes.hpp"
#include "Mesh/CommMap.hpp"
#include "History.hpp"
#include "Table.hpp"

namespace inciter {

namespace cg {

using ncomp_t = tk::ncomp_t;

//! \brief Evaluate the increment from t to t+dt of an analytical solution at
//!   (x,y,z) for all components
std::vector< tk::real >
solinc( tk::ncomp_t ncomp, const std::vector< EOS >&, tk::real x, tk::real y,
        tk::real z, tk::real t, tk::real dt, tk::InitializeFn solution );

//! Compute boundary point normals
std::unordered_map< int,
  std::unordered_map< std::size_t, std::array< tk::real, 4 > > >
bnorm( const std::map< int, std::vector< std::size_t > >& bface,
       const std::vector< std::size_t >& triinpoel,
       const std::array< std::vector< tk::real >, 3 >& coord,
       const std::vector< std::size_t >& gid,
       const std::unordered_map< int,
         std::unordered_set< std::size_t > >& bcnodes );

} // cg::

//! \brief Partial differential equation base for continuous Galerkin PDEs
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   invisible to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a CGPDE,
//!   see inciter::CompFlow.
class CGPDE {

  private:
    using ncomp_t = tk::ncomp_t;
    using real = tk::real;

  public:
    //! Default constructor taking no arguments for Charm++
    explicit CGPDE() = default;

    //! Constructor taking an object modeling Concept.
    //! \details The object of class T comes pre-constructed.
    //! \param[in] x Instantiated object of type T given by the template
    //!   argument.
    template< typename T > explicit CGPDE( T x ) :
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
    explicit CGPDE( std::function<T(Args...)> x, Args&&... args ) :
      self( std::make_unique< Model<T> >(
              std::move( x( std::forward<Args>(args)... ) ) ) ) {}

    //! Public interface to determining which nodes are in IC box
    void IcBoxNodes( const tk::UnsMesh::Coords& coord,
      const std::vector< std::size_t >& inpoel,
      const std::unordered_map< std::size_t, std::set< std::size_t > >& elemblkid,
      std::vector< std::unordered_set< std::size_t > >& inbox,
      std::unordered_map< std::size_t, std::set< std::size_t > >& nodeblkid,
      std::size_t& nuserblk )
    { self->IcBoxNodes( coord, inpoel, elemblkid, inbox, nodeblkid, nuserblk );
    }

    //! Public interface to setting the initial conditions for the diff eq
    void initialize(
      const std::array< std::vector< real >, 3 >& coord,
      tk::Fields& unk,
      real t,
      real V,
      const std::vector< std::unordered_set< std::size_t > >& inbox,
      const std::vector< tk::real >& blkvols,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&
        nodeblkid )
    { self->initialize( coord, unk, t, V, inbox, blkvols, nodeblkid ); }

    //! Public interface to querying a velocity
    void velocity( const tk::Fields& u, tk::UnsMesh::Coords& v ) const
    { self->velocity(u,v); }

    //! Public interface to querying a sound speed
    void soundspeed( const tk::Fields& u, std::vector< tk::real >& s ) const
    { self->soundspeed(u,s); }

    //! Public interface to computing the nodal gradients for ALECG
    void chBndGrad( const std::array< std::vector< real >, 3 >& coord,
      const std::vector< std::size_t >& inpoel,
      const std::vector< std::size_t >& bndel,
      const std::vector< std::size_t >& gid,
      const std::unordered_map< std::size_t, std::size_t >& bid,
      const tk::Fields& U,
      tk::Fields& G ) const
    { self->chBndGrad( coord, inpoel, bndel, gid, bid, U, G ); }

    //! Public interface to computing the right-hand side vector for ALECG
    void rhs(
      real t,
      const std::array< std::vector< real >, 3 >& coord,
      const std::vector< std::size_t >& inpoel,
      const std::vector< std::size_t >& triinpoel,
      const std::vector< std::size_t >& gid,
      const std::unordered_map< std::size_t, std::size_t >& bid,
      const std::unordered_map< std::size_t, std::size_t >& lid,
      const std::vector< real >& dfn,
      const std::pair< std::vector< std::size_t >,
                       std::vector< std::size_t > >& psup,
      const std::pair< std::vector< std::size_t >,
                       std::vector< std::size_t > >& esup,
      const std::vector< int >& symbctri,
      const std::vector< real >& vol,
      const std::vector< std::size_t >& edgenode,
      const std::vector< std::size_t >& edgeid,
      const std::vector< std::unordered_set< std::size_t > >& boxnodes,
      const tk::Fields& G,
      const tk::Fields& U,
      const tk::Fields& W,
      const std::vector< real >& tp,
      real V,
      tk::Fields& R ) const
    { self->rhs( t, coord, inpoel, triinpoel, gid, bid, lid, dfn, psup,
        esup, symbctri, vol, edgenode, edgeid,
        boxnodes, G, U, W, tp, V, R ); }

    //! Public interface to compute boundary surface integrals of pressure
    void bndPressureInt(
      const std::array< std::vector< real >, 3 >& coord,
      const std::vector< std::size_t >& triinpoel,
      const std::vector< int >& symbctri,
      const tk::Fields& U,
      const std::array< tk::real, 3 >& CM,
      std::vector< real >& F ) const
    { self->bndPressureInt( coord, triinpoel, symbctri, U, CM, F ); }

    //! Public interface for computing the minimum time step size
    real dt( const std::array< std::vector< real >, 3 >& coord,
             const std::vector< std::size_t >& inpoel,
             tk::real t,
             tk::real dtn,
             const tk::Fields& U,
             const std::vector< tk::real >& vol,
             const std::vector< tk::real >& voln ) const
    { return self->dt( coord, inpoel, t, dtn, U, vol, voln ); }

    //! Public interface for computing a time step size for each mesh node
    void dt( uint64_t it,
             const std::vector< real >& vol,
             const tk::Fields& U,
             std::vector< real >& dtp ) const
    { self->dt( it, vol, U, dtp ); }

    //! \brief Public interface for querying Dirichlet boundary condition values
    //!  set by the user on a given side set for all components in a PDE system
    std::map< std::size_t, std::vector< std::pair<bool,real> > >
    dirbc( real t,
           real deltat,
           const std::vector< real >& tp,
           const std::vector< real >& dtp,
           const std::pair< const int, std::vector< std::size_t > >& sides,
           const std::array< std::vector< real >, 3 >& coord,
           bool increment ) const
    { return self->dirbc( t, deltat, tp, dtp, sides, coord, increment ); }

    //! Public interface to set symmetry boundary conditions at nodes
    void
    symbc( tk::Fields& U,
           const std::array< std::vector< real >, 3 >& coord,
           const std::unordered_map< int,
                   std::unordered_map< std::size_t,
                     std::array< real, 4 > > >& bnorm,
           const std::unordered_set< std::size_t >& nodes ) const
    { self->symbc( U, coord, bnorm, nodes ); }

    //! Public interface to set farfield boundary conditions at nodes
    void
    farfieldbc( tk::Fields& U,
                const std::array< std::vector< real >, 3 >& coord,
                const std::unordered_map< int,
                        std::unordered_map< std::size_t,
                          std::array< real, 4 > > >& bnorm,
                const std::unordered_set< std::size_t >& nodes ) const
    { self->farfieldbc( U, coord, bnorm, nodes ); }

    //! Public interface to set slip wall boundary conditions at nodes
    void
    slipwallbc( tk::Fields& U,
           const tk::Fields& W,
           const std::array< std::vector< real >, 3 >& coord,
           const std::unordered_map< int,
                   std::unordered_map< std::size_t,
                     std::array< real, 4 > > >& bnorm,
           const std::unordered_set< std::size_t >& nodes ) const
    { self->slipwallbc( U, W, coord, bnorm, nodes ); }

    //! Public interface to applying time dependent boundary conditions at nodes
    void
    timedepbc( tk::real t,
      tk::Fields& U,
      const std::vector< std::unordered_set< std::size_t > >& nodes,
      const std::vector< tk::Table<5> >& timedepfn ) const
    { self->timedepbc( t, U, nodes, timedepfn ); }

    //! Public interface to returning maps of output var functions
    std::map< std::string, tk::GetVarFn > OutVarFn() const
    { return self->OutVarFn(); }

    //! Public interface to returning analytic field output labels
    std::vector< std::string > analyticFieldNames() const
    { return self->analyticFieldNames(); }

    //! Public interface to returning surface field output labels
    std::vector< std::string > surfNames() const { return self->surfNames(); }

    //! Public interface to returning time history field output labels
    std::vector< std::string > histNames() const { return self->histNames(); }

    //! Public interface to returning variable names
    std::vector< std::string > names() const { return self->names(); }

    //! Public interface to returning nodal surface field output
    std::vector< std::vector< real > >
    surfOutput( const std::map< int, std::vector< std::size_t > >& bnd,
                const tk::Fields& U ) const
    { return self->surfOutput( bnd, U ); }

    //! Public interface to returning elemental surface field output
    std::vector< std::vector< real > >
    elemSurfOutput( const std::map< int, std::vector< std::size_t > >& bface,
      const std::vector< std::size_t >& triinpoel,
      const tk::Fields& U ) const
    { return self->elemSurfOutput( bface, triinpoel, U ); }

    //! Public interface to returning time history output
    std::vector< std::vector< real > >
    histOutput( const std::vector< HistData >& h,
                const std::vector< std::size_t >& inpoel,
                const tk::Fields& U ) const
    { return self->histOutput( h, inpoel, U ); }

    //! Public interface to returning analytic solution
    tk::InitializeFn::result_type
    analyticSolution( real xi, real yi, real zi, real t ) const
    { return self->analyticSolution( xi, yi, zi, t ); }

    //! Public interface to returning the analytic solution for conserved vars
    tk::InitializeFn::result_type
    solution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return self->solution( xi, yi, zi, t ); }

    //! Copy assignment
    CGPDE& operator=( const CGPDE& x )
    { CGPDE tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    CGPDE( const CGPDE& x ) : self( x.self->copy() ) {}
    //! Move assignment
    CGPDE& operator=( CGPDE&& ) noexcept = default;
    //! Move constructor
    CGPDE( CGPDE&& ) noexcept = default;

  private:
    //! \brief Concept is a pure virtual base class specifying the requirements
    //!   of polymorphic objects deriving from it
    struct Concept {
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void IcBoxNodes( const tk::UnsMesh::Coords&,
        const std::vector< std::size_t >&,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&,
        std::vector< std::unordered_set< std::size_t > >&,
        std::unordered_map< std::size_t, std::set< std::size_t > >&,
        std::size_t& ) = 0;
      virtual void initialize(
        const std::array< std::vector< real >, 3 >&,
        tk::Fields&,
        real,
        real,
        const std::vector< std::unordered_set< std::size_t > >&,
        const std::vector< tk::real >&,
        const std::unordered_map< std::size_t, std::set< std::size_t > >& ) = 0;
      virtual void velocity( const tk::Fields&, tk::UnsMesh::Coords& )
        const = 0;
      virtual void soundspeed( const tk::Fields&, std::vector< tk::real >& )
        const = 0;
      virtual void chBndGrad( const std::array< std::vector< real >, 3 >&,
        const std::vector< std::size_t >&,
        const std::vector< std::size_t >&,
        const std::vector< std::size_t >&,
        const std::unordered_map< std::size_t, std::size_t >&,
        const tk::Fields&,
        tk::Fields& ) const = 0;
      virtual void rhs(
        real,
        const std::array< std::vector< real >, 3 >&,
        const std::vector< std::size_t >&,
        const std::vector< std::size_t >&,
        const std::vector< std::size_t >&,
        const std::unordered_map< std::size_t, std::size_t >&,
        const std::unordered_map< std::size_t, std::size_t >&,
        const std::vector< real >&,
        const std::pair< std::vector< std::size_t >,
                         std::vector< std::size_t > >&,
        const std::pair< std::vector< std::size_t >,
                         std::vector< std::size_t > >&,
        const std::vector< int >&,
        const std::vector< real >&,
        const std::vector< std::size_t >&,
        const std::vector< std::size_t >&,
        const std::vector< std::unordered_set< std::size_t > >&,
        const tk::Fields&,
        const tk::Fields&,
        const tk::Fields&,
        const std::vector< real >&,
        real,
        tk::Fields& ) const = 0;
      virtual void bndPressureInt(
        const std::array< std::vector< real >, 3 >&,
        const std::vector< std::size_t >&,
        const std::vector< int >&,
        const tk::Fields&,
        const std::array< tk::real, 3 >&,
        std::vector< real >& ) const = 0;
      virtual real dt( const std::array< std::vector< real >, 3 >&,
                       const std::vector< std::size_t >&,
                       tk::real,
                       tk::real,
                       const tk::Fields&,
                       const std::vector< tk::real >& ,
                       const std::vector< tk::real >& ) const = 0;
      virtual void dt( uint64_t,
                       const std::vector< real > &,
                       const tk::Fields&,
                       std::vector< real >& ) const = 0;
      virtual std::map< std::size_t, std::vector< std::pair<bool,real> > >
      dirbc( real,
             real,
             const std::vector< real >&,
             const std::vector< real >&,
             const std::pair< const int, std::vector< std::size_t > >&,
             const std::array< std::vector< real >, 3 >&,
             bool ) const = 0;
      virtual void symbc(
        tk::Fields& U,
        const std::array< std::vector< real >, 3 >&,
        const std::unordered_map< int,
                std::unordered_map< std::size_t,
                  std::array< real, 4 > > >&,
        const std::unordered_set< std::size_t >& ) const = 0;
      virtual void farfieldbc(
        tk::Fields&,
        const std::array< std::vector< real >, 3 >&,
        const std::unordered_map< int,
                std::unordered_map< std::size_t,
                  std::array< real, 4 > > >&,
        const std::unordered_set< std::size_t >& ) const = 0;
      virtual void slipwallbc(
        tk::Fields& U,
        const tk::Fields& W,
        const std::array< std::vector< real >, 3 >&,
        const std::unordered_map< int,
                std::unordered_map< std::size_t,
                  std::array< real, 4 > > >&,
        const std::unordered_set< std::size_t >& ) const = 0;
      virtual void timedepbc(
        tk::real,
        tk::Fields&,
        const std::vector< std::unordered_set< std::size_t > >&,
        const std::vector< tk::Table<5> >& ) const = 0;
      virtual std::map< std::string, tk::GetVarFn > OutVarFn() const = 0;
      virtual std::vector< std::string > analyticFieldNames() const = 0;
      virtual std::vector< std::string > surfNames() const = 0;
      virtual std::vector< std::string > histNames() const = 0;
      virtual std::vector< std::string > names() const = 0;
      virtual std::vector< std::vector< real > > surfOutput(
        const std::map< int, std::vector< std::size_t > >&,
        const tk::Fields& ) const = 0;
      virtual std::vector< std::vector< real > > elemSurfOutput(
        const std::map< int, std::vector< std::size_t > >&,
        const std::vector< std::size_t >&,
        const tk::Fields& ) const = 0;
      virtual std::vector< std::vector< real > > histOutput(
        const std::vector< HistData >&,
        const std::vector< std::size_t >&,
        const tk::Fields& ) const = 0;
      virtual tk::InitializeFn::result_type analyticSolution(
        real xi, real yi, real zi, real t ) const = 0;
      virtual tk::InitializeFn::result_type solution(
        tk::real xi, tk::real yi, tk::real zi, tk::real t ) const = 0;
    };

    //! \brief Model models the Concept above by deriving from it and overriding
    //!   the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      explicit Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      void IcBoxNodes( const tk::UnsMesh::Coords& coord,
        const std::vector< std::size_t >& inpoel,
        const std::unordered_map< std::size_t, std::set< std::size_t > >& elemblkid,
        std::vector< std::unordered_set< std::size_t > >& inbox,
        std::unordered_map< std::size_t, std::set< std::size_t > >& nodeblkid,
        std::size_t& nuserblk )
      override { data.IcBoxNodes( coord, inpoel, elemblkid, inbox, nodeblkid,
        nuserblk ); }
      void initialize(
        const std::array< std::vector< real >, 3 >& coord,
        tk::Fields& unk,
        real t,
        real V,
        const std::vector< std::unordered_set< std::size_t > >& inbox,
        const std::vector< tk::real >& blkvols,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&
          nodeblkid)
      override { data.initialize( coord, unk, t, V, inbox, blkvols, nodeblkid ); }
      void velocity( const tk::Fields& u, tk::UnsMesh::Coords& v ) const
      override { data.velocity(u,v); }
      void soundspeed( const tk::Fields& u, std::vector< tk::real >& s ) const
      override { data.soundspeed(u,s); }
      void chBndGrad( const std::array< std::vector< real >, 3 >& coord,
        const std::vector< std::size_t >& inpoel,
        const std::vector< std::size_t >& bndel,
        const std::vector< std::size_t >& gid,
        const std::unordered_map< std::size_t, std::size_t >& bid,
        const tk::Fields& U,
        tk::Fields& G ) const override
      { data.chBndGrad( coord, inpoel, bndel, gid, bid, U, G ); }
      void rhs(
        real t,
        const std::array< std::vector< real >, 3 >& coord,
        const std::vector< std::size_t >& inpoel,
        const std::vector< std::size_t >& triinpoel,
        const std::vector< std::size_t >& gid,
        const std::unordered_map< std::size_t, std::size_t >& bid,
        const std::unordered_map< std::size_t, std::size_t >& lid,
        const std::vector< real >& dfn,
        const std::pair< std::vector< std::size_t >,
                         std::vector< std::size_t > >& psup,
        const std::pair< std::vector< std::size_t >,
                         std::vector< std::size_t > >& esup,
        const std::vector< int >& symbctri,
        const std::vector< real >& vol,
        const std::vector< std::size_t >& edgenode,
        const std::vector< std::size_t >& edgeid,
        const std::vector< std::unordered_set< std::size_t > >& boxnodes,
        const tk::Fields& G,
        const tk::Fields& U,
        const tk::Fields& W,
        const std::vector< real >& tp,
        real V,
        tk::Fields& R ) const override
      { data.rhs( t, coord, inpoel, triinpoel, gid, bid, lid, dfn, psup,
                  esup, symbctri, vol, edgenode,
                  edgeid, boxnodes, G, U, W, tp, V, R ); }
      void bndPressureInt(
        const std::array< std::vector< real >, 3 >& coord,
        const std::vector< std::size_t >& triinpoel,
        const std::vector< int >& symbctri,
        const tk::Fields& U,
        const std::array< tk::real, 3 >& CM,
        std::vector< real >& F ) const override
      { data.bndPressureInt( coord, triinpoel, symbctri, U, CM, F ); }
      real dt( const std::array< std::vector< real >, 3 >& coord,
               const std::vector< std::size_t >& inpoel,
               tk::real t,
               tk::real dtn,
               const tk::Fields& U,
               const std::vector< tk::real >& vol,
               const std::vector< tk::real >& voln ) const override
      { return data.dt( coord, inpoel, t, dtn, U, vol, voln ); }
      void dt( uint64_t it,
               const std::vector< real > & vol,
               const tk::Fields& U,
               std::vector< real >& dtp ) const override
      { data.dt( it, vol, U, dtp ); }
      std::map< std::size_t, std::vector< std::pair<bool,real> > >
      dirbc( real t,
             real deltat,
             const std::vector< real >& tp,
             const std::vector< real >& dtp,
             const std::pair< const int, std::vector< std::size_t > >& sides,
             const std::array< std::vector< real >, 3 >& coord,
             bool increment ) const
        override { return data.dirbc( t, deltat, tp, dtp, sides, coord,
                                      increment ); }
      void symbc(
        tk::Fields& U,
        const std::array< std::vector< real >, 3 >& coord,
        const std::unordered_map< int,
                std::unordered_map< std::size_t,
                  std::array< real, 4 > > >& bnorm,
        const std::unordered_set< std::size_t >& nodes ) const override
      { data.symbc( U, coord, bnorm, nodes ); }
      void farfieldbc(
        tk::Fields& U,
        const std::array< std::vector< real >, 3 >& coord,
        const std::unordered_map< int,
                std::unordered_map< std::size_t,
                  std::array< real, 4 > > >& bnorm,
        const std::unordered_set< std::size_t >& nodes ) const override
      { data.farfieldbc( U, coord, bnorm, nodes ); }
      void slipwallbc(
        tk::Fields& U,
        const tk::Fields& W,
        const std::array< std::vector< real >, 3 >& coord,
        const std::unordered_map< int,
                std::unordered_map< std::size_t,
                  std::array< real, 4 > > >& bnorm,
        const std::unordered_set< std::size_t >& nodes ) const override
      { data.slipwallbc( U, W, coord, bnorm, nodes ); }
      void
      timedepbc(
        tk::real t,
        tk::Fields& U,
        const std::vector< std::unordered_set< std::size_t > >& nodes,
        const std::vector< tk::Table<5> >& timedepfn ) const override
      { data.timedepbc( t, U, nodes, timedepfn ); }
      std::map< std::string, tk::GetVarFn > OutVarFn() const override
      { return data.OutVarFn(); }
      std::vector< std::string > analyticFieldNames() const override
      { return data.analyticFieldNames(); }
      std::vector< std::string > surfNames() const override
      { return data.surfNames(); }
      std::vector< std::string > histNames() const override
      { return data.histNames(); }
      std::vector< std::string > names() const override
      { return data.names(); }
      std::vector< std::vector< real > > surfOutput(
        const std::map< int, std::vector< std::size_t > >& bnd,
        const tk::Fields& U ) const override
      { return data.surfOutput( bnd, U ); }
      std::vector< std::vector< real > > elemSurfOutput(
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::vector< std::size_t >& triinpoel,
        const tk::Fields& U ) const override
      { return data.elemSurfOutput( bface, triinpoel, U ); }
      std::vector< std::vector< real > > histOutput(
        const std::vector< HistData >& h,
        const std::vector< std::size_t >& inpoel,
        const tk::Fields& U ) const override
      { return data.histOutput( h, inpoel, U ); }
      tk::InitializeFn::result_type
      analyticSolution( real xi, real yi, real zi, real t )
       const override { return data.analyticSolution( xi, yi, zi, t ); }
      tk::InitializeFn::result_type
      solution( real xi, real yi, real zi, real t )
       const override { return data.solution( xi, yi, zi, t ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // inciter::

#endif // CGPDE_h
