// *****************************************************************************
/*!
  \file      src/PDE/Transport/DGTransport.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Scalar transport using disccontinous Galerkin discretization
  \details   This file implements the physics operators governing transported
     scalars using disccontinuous Galerkin discretization.
*/
// *****************************************************************************
#ifndef DGTransport_h
#define DGTransport_h

#include <vector>
#include <array>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <map>

#include "Macro.h"
#include "Exception.h"
#include "Vector.h"
#include "Inciter/Options/BC.h"
#include "UnsMesh.h"
#include "Integrate/Quadrature.h"
#include "Integrate/Initialize.h"
#include "Integrate/Mass.h"
#include "Integrate/Surface.h"
#include "Integrate/Boundary.h"
#include "Integrate/Volume.h"
#include "Integrate/Riemann/Upwind.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief Transport equation used polymorphically with tk::DGPDE
//! \details The template argument(s) specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/Transport/Physics.h
//!   - Problem - problem configuration, see PDE/Transport/Problem.h
//! \note The default physics is DGAdvection, set in
//!    inciter::deck::check_transport()
template< class Physics, class Problem >
class Transport {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;
    using bcconf_t = kw::sideset::info::expect::type;
    using eq = tag::transport;

    //! Extract BC configuration ignoring if BC not specified
    //! \param[in] c Equation system index (among multiple systems configured)
    //! \return Vector of BC config of type bcconf_t used to apply BCs for all
    //!   scalar components this Transport eq system is configured for
    //! \note A more preferable way of catching errors such as this function
    //!   hides is during parsing, so that we don't even get here if BCs are not
    //!   correctly specified. For now we simply ignore if BCs are not
    //!   specified by allowing empty BC vectors from the user input.
    template< typename bctag >
    std::vector< bcconf_t >
    config( ncomp_t c ) {
      std::vector< bcconf_t > bc;
      const auto& v = g_inputdeck.get< tag::param, eq, bctag >();
      if (v.size() > c) bc = v[c];
      return bc;
    }

  public:
    //! Constructor
    //! \param[in] c Equation system index (among multiple systems configured)
    explicit Transport( ncomp_t c ) :
      m_system( c ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< eq >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< eq >(c) ),
      m_bcextrapolate( config< tag::bcextrapolate >( c ) ),
      m_bcinlet( config< tag::bcinlet >( c ) ),
      m_bcoutlet( config< tag::bcoutlet >( c ) ),
      m_bcdir( config< tag::bcdir >( c ) )
    {
      Problem::errchk( m_system, m_ncomp );
    }

    //! Initalize the transport equations for DG
    //! \param[in] L Element mass matrix
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    void initialize( const tk::Fields& L,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord,
                     tk::Fields& unk,
                     tk::real t ) const
    {
      tk::initialize( m_system, m_ncomp, m_offset, L, inpoel, coord,
                      Problem::solution, unk, t );
    }

    //! Compute the left hand side mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const {
      tk::mass( m_ncomp, m_offset, geoElem, l );
    }

    //! Compute right hand side
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] U Solution vector at recent time step
    //! \param[in] limFunc Limiter function for higher-order solution dofs
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              const tk::Fields& limFunc,
              tk::Fields& R ) const
    {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();

      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nprop() == ndof*m_ncomp && R.nprop() == ndof*m_ncomp,
              "Number of components in solution and right-hand side vector " 
              "must equal "+ std::to_string(ndof*m_ncomp) );
      Assert( inpoel.size()/4 == U.nunk(), "Connectivity inpoel has incorrect "
              "size" );
      Assert( fd.Inpofa().size()/3 == fd.Esuf().size()/2,
              "Mismatch in inpofa size" );

      // set rhs to zero
      R.fill(0.0);

      // supported boundary condition types and associated state functions
      std::vector< std::pair< std::vector< bcconf_t >, tk::StateFn > > bctypes{{
        { m_bcextrapolate, Extrapolate },
        { m_bcinlet, Inlet },
        { m_bcoutlet, Outlet },
        { m_bcdir, Dirichlet } }};

      // compute internal surface flux integrals
      tk::surfInt( m_system, m_ncomp, m_offset, inpoel, coord, fd, geoFace,
                   Upwind::flux, Problem::prescribedVelocity, U, limFunc, R );

      // compute boundary surface flux integrals
      for (const auto& b : bctypes)
        tk::sidesetInt( m_system, m_ncomp, m_offset, b.first, fd, geoFace,
          inpoel, coord, t, Upwind::flux, Problem::prescribedVelocity,
          b.second, U, limFunc, R );

      switch(ndof)
      {
        case 1:           //DG(P0)
          break;
 
        case 4:          //DG(P1)
          // compute volume integrals
          tk::volIntP1( m_system, m_ncomp, m_offset, inpoel, coord, geoElem,
                        flux, Problem::prescribedVelocity, U, limFunc, R );
          break;

        case 10:        //DG(P2)
          // compute volume integrals
          tk::volIntP2( m_system, m_ncomp, m_offset, inpoel, coord, geoElem,
                        flux, Problem::prescribedVelocity, U, R );
          break;

        default:
          Throw( "dg::Transport::rhs() not defined for NDOF=" +
                 std::to_string(ndof) );
      }
    }

    //! Compute the minimum time step size
//     //! \param[in] U Solution vector at recent time step
//     //! \param[in] coord Mesh node coordinates
//     //! \param[in] inpoel Mesh element connectivity
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& /*coord*/,
                 const std::vector< std::size_t >& /*inpoel*/,
                 const inciter::FaceData& /*fd*/,
                 const tk::Fields& /*geoFace*/,
                 const tk::Fields& /*geoElem*/,
                 const tk::Fields& /*limFunc*/,
                 const tk::Fields& /*U*/ ) const
    {
      tk::real mindt = std::numeric_limits< tk::real >::max();
      return mindt;
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    void side( std::unordered_set< int >& conf ) const
    { Problem::side( conf ); }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    //! \details This functions should be written in conjunction with
    //!   fieldOutput(), which provides the vector of fields to be output
    std::vector< std::string > fieldNames() const {
      std::vector< std::string > n;
      const auto& depvar =
      g_inputdeck.get< tag::param, eq, tag::depvar >().at(m_system);
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_numerical" );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_analytic" );
      // will output error for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_error" );
      return n;
    }

    //!
    std::vector< std::vector< tk::real > >
    avgElemToNode( const std::vector< std::size_t >& /*inpoel*/,
                   const tk::UnsMesh::Coords& /*coord*/,
                   const tk::Fields& /*geoElem*/,
                   const tk::Fields& /*limFunc*/,
                   const tk::Fields& /*U*/ ) const
    {
      std::vector< std::vector< tk::real > > out;
      return out;
    }

    //! Return field output going to file
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] t Physical time
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    //! \details This functions should be written in conjunction with names(),
    //!   which provides the vector of field names
    //! \note U is overwritten
    std::vector< std::vector< tk::real > >
    fieldOutput( const tk::Fields&,
                 const std::vector< std::size_t >& inpoel,
                 const tk::UnsMesh::Coords& coord,
                 tk::real t,
                 const tk::Fields& geoElem,
                 tk::Fields& U ) const
    {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      Assert( geoElem.nunk() == U.nunk(), "Size mismatch" );
      std::vector< std::vector< tk::real > > out;
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c*ndof, m_offset ) );
      // evaluate analytic solution at time t
      auto E = U;
      tk::initializeP0( m_system, m_ncomp, m_offset, inpoel, coord,
                        Problem::solution, E, t );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( E.extract( c*ndof, m_offset ) );
      // will output error for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) {
        auto mark = c*ndof;
        auto u = U.extract( mark, m_offset );
        auto e = E.extract( mark, m_offset );
        for (std::size_t i=0; i<u.size(); ++i)
          e[i] = std::pow( e[i] - u[i], 2.0 ) * geoElem(i,0,0);
        out.push_back( e );
      }
      return out;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const {
      std::vector< std::string > n;
      const auto& depvar =
      g_inputdeck.get< tag::param, eq, tag::depvar >().at(m_system);
      // construct the name of the numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) );
      return n;
    }

    //! Return analytic solution (if defined by Problem) at xi, yi, zi, t
    //! \param[in] xi X-coordinate at which to evaluate the analytic solution
    //! \param[in] yi Y-coordinate at which to evaluate the analytic solution
    //! \param[in] zi Z-coordinate at which to evaluate the analytic solution
    //! \param[in] t Physical time at which to evaluate the analytic solution
    //! \return Vector of analytic solution at given spatial location and time
    std::vector< tk::real >
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return Problem::solution( m_system, m_ncomp, xi, yi, zi, t ); }

  private:
    const ncomp_t m_system;             //!< Equation system index
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    const ncomp_t m_offset;             //!< Offset this PDE operates from
    //! Extrapolation BC configuration
    const std::vector< bcconf_t > m_bcextrapolate;
    //! Inlet BC configuration
    const std::vector< bcconf_t > m_bcinlet;
    //! Outlet BC configuration
    const std::vector< bcconf_t > m_bcoutlet;
    //! Dirichlet BC configuration
    const std::vector< bcconf_t > m_bcdir;

    //! Evaluate physical flux function for this PDE system
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ugp Numerical solution at the Gauss point at which to
    //!   evaluate the flux
    //! \param[in] v Prescribed velocity evaluated at the Gauss point at which
    //!   to evaluate the flux
    //! \return Flux vectors for all components in this PDE system
    //! \note The function signature must follow tk::FluxFn
    static tk::FluxFn::result_type
    flux( ncomp_t,
          ncomp_t ncomp,
          const std::vector< tk::real >& ugp,
          const std::vector< std::array< tk::real, 3 > >& v )
    {
      Assert( ugp.size() == ncomp, "Size mismatch" );
      Assert( v.size() == ncomp, "Size mismatch" );

      std::vector< std::array< tk::real, 3 > > fl( ugp.size() );

      for (ncomp_t c=0; c<ncomp; ++c)
        fl[c] = {{ v[c][0] * ugp[c], v[c][1] * ugp[c], v[c][2] * ugp[c] }};

      return fl;
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at extrapolation boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    Extrapolate( ncomp_t, ncomp_t, const std::vector< tk::real >& ul,
                 tk::real, tk::real, tk::real, tk::real,
                 const std::array< tk::real, 3 >& )
    {
      return {{ ul, ul }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at extrapolation boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    Inlet( ncomp_t, ncomp_t, const std::vector< tk::real >& ul,
           tk::real, tk::real, tk::real, tk::real,
           const std::array< tk::real, 3 >& )
    {
      auto ur = ul;
      std::fill( begin(ur), end(ur), 0.0 );
      return {{ ul, std::move(ur) }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at outlet boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    Outlet( ncomp_t, ncomp_t, const std::vector< tk::real >& ul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& )
    {
      return {{ ul, ul }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at Dirichlet boundaries
    //! \param[in] system Equation system index
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ul Left (domain-internal) state
    //! \param[in] x X-coordinate at which to compute the states
    //! \param[in] y Y-coordinate at which to compute the states
    //! \param[in] z Z-coordinate at which to compute the states
    //! \param[in] t Physical time
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    Dirichlet( ncomp_t system, ncomp_t ncomp, const std::vector< tk::real >& ul,
               tk::real x, tk::real y, tk::real z, tk::real t,
               const std::array< tk::real, 3 >& )
    {
      return {{ ul, Problem::solution( system, ncomp, x, y, z, t ) }};
    }
};

} // dg::
} // inciter::

#endif // DGTransport_h
