// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/CylVortex.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for transport equations
  \details   This file declares a Problem policy class for the transport
    equations, defined in PDE/Transport/CGTransport.hpp implementing
    node-centered continuous Galerkin (CG) and PDE/Transport/DGTransport.hpp
    implementing cell-centered discontinuous Galerkin (DG) discretizations.
    See PDE/Transport/Problem.hpp for general requirements on Problem policy
    classes for cg::Transport and dg::Transport.
*/
// *****************************************************************************
#ifndef TransportProblemCylVortex_h
#define TransportProblemCylVortex_h

#include <vector>
#include <array>

#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/Options/Problem.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

//! Transport PDE problem: deformation of cylinder in a vortex
class TransportProblemCylVortex {
  private:
    using ncomp_t = tk::ncomp_t;
    using eq = tag::transport;

  public:
    //! Evaluate analytical solution at (x,y,t) for all components
    static std::vector< tk::real >
    initialize( ncomp_t ncomp,
                const std::vector< EOS >& mat_blk, tk::real x, tk::real y,
                tk::real, tk::real t );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static std::vector< tk::real >
    analyticSolution( ncomp_t ncomp,
                      const std::vector< EOS >& mat_blk, tk::real x,
                      tk::real y, tk::real z, tk::real t )
    { return initialize( ncomp, mat_blk, x, y, z, t ); }

    //! Do error checking on PDE parameters
    void errchk( ncomp_t ) const {}

    //! Assign prescribed velocity at a point
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( ncomp_t ncomp, tk::real x, tk::real y,
      tk::real,
      tk::real t );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::CYL_VORTEX; }
};

} // inciter::

#endif // TransportProblemCylVortex_h
