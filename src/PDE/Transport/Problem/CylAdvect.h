// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/CylAdvect.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#ifndef TransportProblemCylAdvect_h
#define TransportProblemCylAdvect_h

#include <vector>
#include <array>

#include "Types.h"
#include "SystemComponents.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! Transport PDE problem: advection of cylinder
class TransportProblemCylAdvect {
  private:
    using ncomp_t = tk::ctr::ncomp_type;
    using eq = tag::transport;

  public:
    //! Evaluate analytical solution at (x,y,t) for all components
    static std::vector< tk::real >
    solution( ncomp_t system, ncomp_t ncomp,
              tk::real x, tk::real y, tk::real, tk::real t );

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    std::vector< tk::real >
    solinc( ncomp_t, ncomp_t ncomp, tk::real x, tk::real y, tk::real,
            tk::real t, tk::real dt ) const;

    //! Do error checking on PDE parameters
    void errchk( ncomp_t, ncomp_t ) const {}

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    void side( std::unordered_set< int >& conf ) const;

    //! Assign prescribed velocity at a point
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( ncomp_t, ncomp_t ncomp, tk::real, tk::real, tk::real );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::CYL_ADVECT; }
};

} // inciter::

#endif // TransportProblemCylAdvect_h
