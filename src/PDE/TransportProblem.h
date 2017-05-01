// *****************************************************************************
/*!
  \file      src/PDE/TransportProblem.h
  \author    J. Bakosi
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

#include <vector>
#include <array>

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
                      tk::Fields& unk,
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
    static std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields&,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& )
    { return {{ {{0.0, 0.0, 0.0, 0.0}},
                {{0.0, 0.0, 0.0, 0.0}},
                {{0.0, 0.0, 0.0, 0.0}} }}; }

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
    //! \param[in] time Physical time
    template< class eq >
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type,
                      tk::ctr::ncomp_type ncomp,
                      tk::ctr::ncomp_type offset,
                      tk::real time )
    {
      for (ncomp_t c=0; c<ncomp; ++c) {
        auto t = time + 2.0*M_PI/ncomp * c;
        const tk::real R0 = 0.15;
        // center of the cone
        tk::real x0 = 0.5;
        tk::real y0 = 0.25;
        tk::real r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real kx = 0.5 + r*std::sin( t );
        tk::real ky = 0.5 - r*std::cos( t );
        // center of the hump
        x0 = 0.25;
        y0 = 0.5;
        r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real hx = 0.5 + r*std::sin( t-M_PI/2.0 ),
                 hy = 0.5 - r*std::cos( t-M_PI/2.0 );
        // center of the slotted cylinder
        x0 = 0.5;
        y0 = 0.75;
        r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real cx = 0.5 + r*std::sin( t+M_PI ),
                 cy = 0.5 - r*std::cos( t+M_PI );
        // end points of the cylinder slot
        tk::real i1x = 0.525, i1y = cy - r*std::cos( std::asin(0.025/r) ),
                 i2x = 0.525, i2y = 0.8,
                 i3x = 0.475, i3y = 0.8;
        // rotate end points of cylinder slot
        tk::real ri1x = 0.5 + std::cos(t)*(i1x-0.5) - std::sin(t)*(i1y-0.5),
                 ri1y = 0.5 + std::sin(t)*(i1x-0.5) + std::cos(t)*(i1y-0.5),
                 ri2x = 0.5 + std::cos(t)*(i2x-0.5) - std::sin(t)*(i2y-0.5),
                 ri2y = 0.5 + std::sin(t)*(i2x-0.5) + std::cos(t)*(i2y-0.5),
                 ri3x = 0.5 + std::cos(t)*(i3x-0.5) - std::sin(t)*(i3y-0.5),
                 ri3y = 0.5 + std::sin(t)*(i3x-0.5) + std::cos(t)*(i3y-0.5);
        // direction of slot sides
        tk::real v1x = ri2x-ri1x, v1y = ri2y-ri1y,
                 v2x = ri3x-ri2x, v2y = ri3y-ri2y;
        // lengths of direction of slot sides vectors
        tk::real v1 = std::sqrt(v1x*v1x + v1y*v1y),
                 v2 = std::sqrt(v2x*v2x + v2y*v2y);
        const auto& x = coord[0];
        const auto& y = coord[1];
        unk.fill( c, offset, 0.0 );
        for (ncomp_t i=0; i<x.size(); ++i) {
          // cone
          r = std::sqrt((x[i]-kx)*(x[i]-kx) + (y[i]-ky)*(y[i]-ky)) / R0;
          if (r<1.0) unk( i, c, offset ) = 0.6*(1.0-r);
          // hump
          r = std::sqrt((x[i]-hx)*(x[i]-hx) + (y[i]-hy)*(y[i]-hy)) / R0;
          if (r<1.0)
              unk( i, c, offset ) = 0.2*(1.0+std::cos(M_PI*std::min(r,1.0)));
          // cylinder
          r = std::sqrt((x[i]-cx)*(x[i]-cx) + (y[i]-cy)*(y[i]-cy)) / R0;
          const std::array< tk::real, 2 > r1{{ v1x, v1y }},
                                          r2{{ x[i]-ri1x, y[i]-ri1y }};
          const auto d1 = (r1[0]*r2[1] - r2[0]*r1[1]) / v1;
          const std::array< tk::real, 2 > r3{{ v2x, v2y }},
                                          r4{{ x[i]-ri2x, y[i]-ri2y }};
          const auto d2 = (r3[0]*r4[1] - r4[0]*r3[1]) / v2;
          if (r<1.0 && (d1>0.05 || d1<0.0 || d2<0.0)) unk( i, c, offset ) = 0.6;
        }
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
    static std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields&,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::array< std::size_t, 4 >& N )
    {
      const auto& x = coord[0];
      const auto& y = coord[1];
      return {{ {{ 0.5-y[N[0]], 0.5-y[N[1]], 0.5-y[N[2]], 0.5-y[N[3]] }},
                {{ x[N[0]]-0.5, x[N[1]]-0.5, x[N[2]]-0.5, x[N[3]]-0.5 }},
                {{ 0.0,         0.0,         0.0,         0.0 }} }};
    }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SLOT_CYL; }
};

//! List of all transport equation problem policies
using TransportProblems = boost::mpl::vector< TransportProblemShearDiff
                                            , TransportProblemSlotCyl 
                                            >;

} // inciter::

#endif // TransportProblem_h
