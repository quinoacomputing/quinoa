// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/GaussHump.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for transport equations
  \details   This file declares a Problem policy class for the transport
    equations, defined in PDE/Transport/CGTransport.h implementing
    node-centered continuous Galerkin (CG) and PDE/Transport/DGTransport.h
    implementing cell-centered discontinuous Galerkin (DG) discretizations.
    See PDE/Transport/Problem.h for general requirements on Problem policy
    classes for cg::Transport and dg::Transport.
*/
// *****************************************************************************
#ifndef TransportProblemGaussHump_h
#define TransportProblemGaussHump_h

#include <vector>
#include <array>

#include "Types.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"
#include "EoS/EoS_Base.hpp"

namespace inciter {

//! Transport PDE problem: advection of two-dimensional Gaussian hump
class TransportProblemGaussHump {
  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::transport;

  public:
    //! Initialize numerical solution
    static std::vector< tk::real >
    initialize( ncomp_t system, ncomp_t ncomp,
                const std::vector< EoS_Base* >& mat_blk, tk::real x, tk::real y,
                tk::real, tk::real t );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static std::vector< tk::real >
    analyticSolution( ncomp_t system, ncomp_t ncomp,
                      const std::vector< EoS_Base* >& mat_blk, tk::real x,
                      tk::real y, tk::real z, tk::real t )
    { return initialize( system, ncomp, mat_blk, x, y, z, t ); }

    //! Do error checking on PDE parameters
    void errchk( ncomp_t, ncomp_t ) const {}

    //! Assign prescribed velocity at a point
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( ncomp_t,
                        ncomp_t ncomp,
                        tk::real,
                        tk::real,
                        tk::real,
                        tk::real );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::GAUSS_HUMP; }
};

} // inciter::

#endif // TransportProblemGaussHump_h
