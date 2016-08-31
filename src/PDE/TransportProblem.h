// *****************************************************************************
/*!
  \file      src/PDE/TransportProblem.h
  \author    J. Bakosi
  \date      Wed 31 Aug 2016 08:05:59 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for transport equations
  \details   This file defines policy classes for a transport equation, defined
    in PDE/Transport.h.

    General requirements on transport equation problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::ProblemType type() noexcept {
          return ctr::ProblemType::SHEAR_DIFF;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef TransportProblem_h
#define TransportProblem_h

#include <cmath>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! Transport PDE problem: diffusion of a shear layer
class TransportProblemShearDiff {

  public:
    //! Do error checking on PDE parameters
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this transport equation
    template< class eq >
    static void errchk( tk::ctr::ncomp_type e, tk::ctr::ncomp_type ncomp )
    {
      const auto& u0 = g_inputdeck.get< tag::param, eq, tag::u0 >()[e];
      ErrChk( ncomp == u0.size(),
        "Wrong number of advection-diffusion PDE parameters 'u0'" );
      const auto& lambda = g_inputdeck.get< tag::param, eq, tag::lambda >()[e];
      ErrChk( ncomp == lambda.size(),
        "Wrong number of advection-diffusion PDE parameters 'lambda'" );
      const auto& diff =
        g_inputdeck.get< tag::param, eq, tag::diffusivity >()[e];
      ErrChk( ncomp == diff.size(),
        "Wrong number of advection-diffusion PDE parameters 'diffusivity'" );
    }

    //! Set initial conditions for dispersion in simple shear flow
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this transport equation
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    template< class eq >
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::MeshNodes& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type ncomp,
                      tk::ctr::ncomp_type offset,
                      tk::real t )
    {
      const auto& u0 = g_inputdeck.get< tag::param, eq, tag::u0 >()[e];
      const auto& lambda = g_inputdeck.get< tag::param, eq, tag::lambda >()[e];
      const auto& diff =
        g_inputdeck.get< tag::param, eq, tag::diffusivity >()[e];
      const tk::real X0 = 7200.0;        // x position of source
      const tk::real t0 = g_inputdeck.get< tag::discr, tag::t0 >(); // initial t
      const auto& x = coord[0];
      const auto& y = coord[1];
      for (ncomp_t c=0; c<ncomp; ++c) {
        const auto b = 1.0 + lambda[c]*lambda[c]*t*t/12.0;
        const auto M =
          4.0*M_PI*t0*std::sqrt( 1.0 + lambda[c]*lambda[c]*t0*t0/12.0 );
        for (ncomp_t i=0; i<x.size(); ++i) {
          tk::real a = x[i] - X0 - u0[c]*t - 0.5*lambda[c]*y[i]*t;
          unk( i, c, offset ) =
            M * exp( -a*a/(4.0*M_PI*diff[c]*t*b) - y[i]*y[i]/(4.0*diff[c]*t) )
              / ( 4.0*M_PI*t*std::sqrt(b) );
        }
      }
    }

    //! Assign prescribed shear velocity to nodes of tetrahedron element
    //! \param[in] N Element node indices
    //! \param[in] coord Mesh node coordinates
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this transport equation
    //! \return Velocity assigned to all vertices of a tetrehedron, size:
    //!   ncomp * ndim * nnode = [ncomp][3][4]
    template< class eq >
    static std::vector< std::array< std::array< tk::real, 4 >, 3 > >
    prescribedVelocity( const std::array< std::size_t, 4 >& N,
                        const std::array< std::vector< tk::real >, 3 >& coord,
                        tk::ctr::ncomp_type e,
                        tk::ctr::ncomp_type ncomp )
    {
      const auto& u0 = g_inputdeck.get< tag::param, eq, tag::u0 >()[e];
      const auto& lambda = g_inputdeck.get< tag::param, eq, tag::lambda >()[e];
      const auto& y = coord[1];
      std::vector< std::array< std::array< tk::real, 4 >, 3 > > vel( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c) {
        std::array< std::array< tk::real, 4 >, 3 > v;
        v[0][0] = u0[c] + lambda[c]*y[N[0]];  v[1][0] = 0.0;  v[2][0] = 0.0;
        v[0][1] = u0[c] + lambda[c]*y[N[1]];  v[1][1] = 0.0;  v[2][1] = 0.0;
        v[0][2] = u0[c] + lambda[c]*y[N[2]];  v[1][2] = 0.0;  v[2][2] = 0.0;
        v[0][3] = u0[c] + lambda[c]*y[N[3]];  v[1][3] = 0.0;  v[2][3] = 0.0;
        vel[c] = std::move(v);
      }
      return vel;
    }

    //! Return the velocity field at cell nodes
    //! \return Array of the four values of the three velocity coordinates
    static std::vector< std::array< tk::real, 4 > >
    velocity( const tk::MeshNodes&,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& )
    { return { {{0.0, 0.0, 0.0, 0.0}},
               {{0.0, 0.0, 0.0, 0.0}},
               {{0.0, 0.0, 0.0, 0.0}} }; }

    //! Problem type enum accessor
    //! \return Problem type as enum
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SHEAR_DIFF; }
};

//! Transport PDE problem: rotation of Zalesak's slotted cylinder
//! \see Zalesak, S. (1978) Fully Multidimensional Flux-Corrected Transport
//!   Algorithms for Fluids. Journal of Computational Physics 31, 335-362.
//! \see Leveque R.J. (1996) High-Resolution Conservative Algorithms for
//!   Advection in Incompressible Flow. SIAM Journal on Numerical Analyis 33,
//!   627-665.
class TransportProblemSlotCyl {

  public:
    //! Do error checking on PDE parameters
    template< class eq >
    static void errchk( tk::ctr::ncomp_type, tk::ctr::ncomp_type ) {}

    //! Set initial conditions for Zalesak's slotted cylinder case
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] ncomp Number of components in this transport equation
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    template< class eq >
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::MeshNodes& unk,
                      tk::ctr::ncomp_type,
                      tk::ctr::ncomp_type ncomp,
                      tk::ctr::ncomp_type offset,
                      tk::real t )
    {
      tk::real R0 = 0.15;
      // position of the center of the cone
      tk::real x0 = 0.5;
      tk::real y0 = 0.25;
      tk::real r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
      tk::real p0 = std::asin( (x0-0.5)/r );
      tk::real kx = 0.5 + r*std::sin( p0 + t );
      tk::real ky = 0.5 - r*std::cos( p0 + t );
      // position of the center of the hump
      x0 = 0.25;
      y0 = 0.5;
      r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
      p0 = std::asin( (x0-0.5)/r );
      tk::real hx = 0.5 + r*std::sin( p0+t );
      tk::real hy = 0.5 - r*std::cos( p0+t );
      const auto& x = coord[0];
      const auto& y = coord[1];
      for (ncomp_t c=0; c<ncomp; ++c) unk.fill( c, offset, 0.0 );
      for (ncomp_t i=0; i<x.size(); ++i) {
          // cone
          r = std::sqrt((x[i]-kx)*(x[i]-kx) + (y[i]-ky)*(y[i]-ky)) / R0;
          if (r<1.0)
            for (ncomp_t c=0; c<ncomp; ++c)
              unk( i, c, offset ) = 0.6*(1.0-r);
          // hump
          r = std::sqrt((x[i]-hx)*(x[i]-hx) + (y[i]-hy)*(y[i]-hy)) / R0;
          if (r<1.0)
            for (ncomp_t c=0; c<ncomp; ++c)
              unk( i, c, offset ) = 0.2*(1.0+std::cos(M_PI*std::min(r,1.0)));
        }
    }

    //! Assign prescribed velocity to nodes of tetrahedron element
    //! \param[in] N Element node indices
    //! \param[in] coord Mesh node coordinates
    //! \param[in] ncomp Number of components in this transport equation
    //! \return Velocity assigned to all vertices of a tetrehedron, size:
    //!   ncomp * ndim * nnode = [ncomp][3][4]
    template< class eq >
    static std::vector< std::array< std::array< tk::real, 4 >, 3 > >
    prescribedVelocity( const std::array< std::size_t, 4 >& N,
                        const std::array< std::vector< tk::real >, 3 >& coord,
                        tk::ctr::ncomp_type,
                        tk::ctr::ncomp_type ncomp )
    {
      const auto& x = coord[0];
      const auto& y = coord[1];
      std::vector< std::array< std::array< tk::real, 4 >, 3 > > vel( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c) {
        std::array< std::array< tk::real, 4 >, 3 > v;
        v[0][0] = 0.5-y[N[0]];  v[1][0] = x[N[0]]-0.5;  v[2][0] = 0.0;
        v[0][1] = 0.5-y[N[1]];  v[1][1] = x[N[1]]-0.5;  v[2][1] = 0.0;
        v[0][2] = 0.5-y[N[2]];  v[1][2] = x[N[2]]-0.5;  v[2][2] = 0.0;
        v[0][3] = 0.5-y[N[3]];  v[1][3] = x[N[3]]-0.5;  v[2][3] = 0.0;
        vel[c] = std::move(v);
      }
      return vel;
    }

    //! Return the velocity field at cell nodes
    //! \return Array of the four values of the three velocity coordinates
    static std::vector< std::array< tk::real, 4 > >
    velocity( const tk::MeshNodes&,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& )
    { return { {{0.0, 0.0, 0.0, 0.0}},
               {{0.0, 0.0, 0.0, 0.0}},
               {{0.0, 0.0, 0.0, 0.0}} }; }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SLOT_CYL; }
};

//! List of all transport equation problem policies
using TransportProblems = boost::mpl::vector< TransportProblemShearDiff
                                            , TransportProblemSlotCyl 
                                            >;

} // inciter::

#endif // TransportProblem_h
