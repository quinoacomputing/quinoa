// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/SlotCyl.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Problem configuration for transport equations
  \details   This file defines a Problem policy class for the scalar transport
    equations, defined in PDE/Transport/Transport.h. See
    PDE/Transport/Problems.h for general requirements on Problem policy classes
    for Transport.
*/
// *****************************************************************************
#ifndef TransportProblemSlotCyl_h
#define TransportProblemSlotCyl_h

#include <vector>
#include <array>

#include <cmath>

#include "Types.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! Transport PDE problem: rotation of Zalesak's slotted cylinder
//! \see Zalesak, S. (1978) Fully Multidimensional Flux-Corrected Transport
//!   Algorithms for Fluids. Journal of Computational Physics 31, 335-362.
//! \see Leveque R.J. (1996) High-Resolution Conservative Algorithms for
//!   Advection in Incompressible Flow. SIAM Journal on Numerical Analyis 33,
//!   627-665.
class TransportProblemSlotCyl {

  public:
    //! Evaluate analytical solution at (x,y,t) for all components
    //! \param[in] ncomp Number of components in this transport equation system
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] t Time where to evaluate the solution
    //! \return Values of all components evaluated at (x,y,t)
    static std::vector< tk::real >
    solution( ncomp_t, ncomp_t ncomp,
              tk::real x, tk::real y, tk::real, tk::real t )
    {
      using std::sin; using std::cos;
      std::vector< tk::real > s( ncomp, 0.0 );
      for (ncomp_t c=0; c<ncomp; ++c) {
        auto T = t + 2.0*M_PI/ncomp * c;
        const tk::real R0 = 0.15;
        // center of the cone
        tk::real x0 = 0.5;
        tk::real y0 = 0.25;
        tk::real r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real kx = 0.5 + r*sin( T );
        tk::real ky = 0.5 - r*cos( T );
        // center of the hump
        x0 = 0.25;
        y0 = 0.5;
        r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real hx = 0.5 + r*sin( T-M_PI/2.0 ),
                 hy = 0.5 - r*cos( T-M_PI/2.0 );
        // center of the slotted cylinder
        x0 = 0.5;
        y0 = 0.75;
        r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real cx = 0.5 + r*sin( T+M_PI ),
                 cy = 0.5 - r*cos( T+M_PI );
        // end points of the cylinder slot
        tk::real i1x = 0.525, i1y = cy - r*cos( std::asin(0.025/r) ),
                 i2x = 0.525, i2y = 0.8,
                 i3x = 0.475, i3y = 0.8;
        // rotate end points of cylinder slot
        tk::real ri1x = 0.5 + cos(T)*(i1x-0.5) - sin(T)*(i1y-0.5),
                 ri1y = 0.5 + sin(T)*(i1x-0.5) + cos(T)*(i1y-0.5),
                 ri2x = 0.5 + cos(T)*(i2x-0.5) - sin(T)*(i2y-0.5),
                 ri2y = 0.5 + sin(T)*(i2x-0.5) + cos(T)*(i2y-0.5),
                 ri3x = 0.5 + cos(T)*(i3x-0.5) - sin(T)*(i3y-0.5),
                 ri3y = 0.5 + sin(T)*(i3x-0.5) + cos(T)*(i3y-0.5);
        // direction of slot sides
        tk::real v1x = ri2x-ri1x, v1y = ri2y-ri1y,
                 v2x = ri3x-ri2x, v2y = ri3y-ri2y;
        // lengths of direction of slot sides vectors
        tk::real v1 = std::sqrt(v1x*v1x + v1y*v1y),
                 v2 = std::sqrt(v2x*v2x + v2y*v2y);
        // cone
        r = std::sqrt((x-kx)*(x-kx) + (y-ky)*(y-ky)) / R0;
        if (r<1.0) s[c] = 0.6*(1.0-r);
        // hump
        r = std::sqrt((x-hx)*(x-hx) + (y-hy)*(y-hy)) / R0;
        if (r<1.0) s[c] = 0.2*(1.0+cos(M_PI*std::min(r,1.0)));
        // cylinder
        r = std::sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)) / R0;
        const std::array< tk::real, 2 > r1{{ v1x, v1y }},
                                        r2{{ x-ri1x, y-ri1y }};
        const auto d1 = (r1[0]*r2[1] - r2[0]*r1[1]) / v1;
        const std::array< tk::real, 2 > r3{{ v2x, v2y }},
                                        r4{{ x-ri2x, y-ri2y }};
        const auto d2 = (r3[0]*r4[1] - r4[0]*r3[1]) / v2;
        if (r<1.0 && (d1>0.05 || d1<0.0 || d2<0.0)) s[c] = 0.6;
      }
      return s;
    }

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    //! \param[in] ncomp Number of components in this transport equation system
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] t Time where to evaluate the solution increment starting from
    //! \param[in] dt Time increment at which evaluate the solution increment to
    //! \return Increment in values of all components evaluated at (x,y,t+dt)
    static std::vector< tk::real >
    solinc( ncomp_t, ncomp_t ncomp, tk::real x, tk::real y, tk::real,
            tk::real t, tk::real dt )
    {
      auto st1 = solution( 0, ncomp, x, y, 0.0, t );
      auto st2 = solution( 0, ncomp, x, y, 0.0, t+dt );
      std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                      []( tk::real s, tk::real& d ){ return d -= s; } );
      return st2;
    }

    //! Do error checking on PDE parameters
    static void errchk( ncomp_t, ncomp_t ) {}

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::transport; using tag::bcdir;
      for (const auto& s : g_inputdeck.get< param, transport, bcdir >())
        for (const auto& i : s)
          conf.insert( std::stoi(i) );
    }

    //! Assign prescribed velocity at a point
    //! \param[in] x X coordinate at which to assign velocity
    //! \param[in] y y coordinate at which to assign velocity
    //! \param[in] ncomp Number of components in this transport equation
    //! \return Velocity assigned to all vertices of a tetrehedron, size:
    //!   ncomp * ndim = [ncomp][3]
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( tk::real x,
                        tk::real y,
                        tk::real,
                        ncomp_t,
                        ncomp_t ncomp )
    {
      std::vector< std::array< tk::real, 3 > > vel( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c)
        vel[c] = {{ 0.5-y, x-0.5, 0.0 }};
      return vel;
    }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SLOT_CYL; }
};

} // inciter::

#endif // TransportProblemSlotCyl_h
