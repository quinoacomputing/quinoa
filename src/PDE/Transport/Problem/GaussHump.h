// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/GaussHump.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Problem configuration for transport equations
  \details   This file defines a Problem policy class for the scalar transport
    equation, defined in PDE/Transport/Transport.h. See
    PDE/Transport/Problems.h for general requirements on Problem policy classes
    for Transport.
*/
// *****************************************************************************
#ifndef TransportProblemGaussHump_h
#define TransportProblemGaussHump_h

#include <vector>
#include <array>

#include <cmath>

#include "Types.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! Transport PDE problem: advection of two-dimensional Gaussian hump
class TransportProblemGaussHump {

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
      const auto vel = prescribedVelocity( x, y, 0.0, ncomp, ncomp );

      std::vector< tk::real > s( ncomp, 0.0 );
      for (ncomp_t c=0; c<ncomp; ++c) 
      {
        // center of the hump
        auto x0 = 0.25 + vel[c][0]*t;
        auto y0 = 0.25 + vel[c][1]*t;

        // hump
        s[c] = 1.0 * exp( -((x-x0)*(x-x0)
                     + (y-y0)*(y-y0))/(2.0 * 0.005) );
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
    //! \param[in] ncomp Number of components in this transport equation
    //! \return Velocity assigned to all vertices of a tetrehedron, size:
    //!   ncomp * ndim = [ncomp][3]
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( tk::real,
                        tk::real,
                        tk::real,
                        ncomp_t,
                        ncomp_t ncomp )
    {
      std::vector< std::array< tk::real, 3 > > vel( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c)
        vel[c] = {{ 0.1, 0.1, 0.0 }};
      return vel;
    }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::GAUSS_HUMP; }
};

} // inciter::

#endif // TransportProblemGaussHump_h
