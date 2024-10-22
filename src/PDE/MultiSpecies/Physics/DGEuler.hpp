// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Physics/DGEuler.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for the equations governing multi-species flow
    using a DG method
  \details   This file defines a Physics policy class for the compressible
    flow equations class dg::MultiSpecies, defined in
    PDE/MultiSpecies/DGMultiSpecies.h. This specific algorithm solves the
    equations of multi-species flow and uses a DG discretization scheme. See
    PDE/MultiSpecies/Physics/DG.h for general requirements on Physics policy
    classes for dg::MultiSpecies.
*/
// *****************************************************************************
#ifndef MultiSpeciesPhysicsDGEuler_h
#define MultiSpeciesPhysicsDGEuler_h

#include "Types.hpp"
#include "Exception.hpp"
#include "Inciter/Options/Physics.hpp"

namespace inciter {

namespace dg {

//! MultiSpecies system of PDEs problem: Euler (inviscid)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in MultiSpecies the Euler equations
//!   governing multi-material compressible flow.
class MultiSpeciesPhysicsEuler {

  public:
    //! Compute the time step size restriction based on this physics
    //! \return A large time step size, i.e., ignore
    tk::real dtRestriction(
      const tk::Fields&,
      std::size_t,
      const std::vector< int >& ) const
    { return std::numeric_limits< tk::real >::max(); }

    //! Compute sources corresponding to this physics
    void physSrc(
      std::size_t,
      tk::real,
      const tk::Fields&,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&,
      tk::Fields&,
      std::vector< int >& ) const {}

    //! Return enum denoting physics policy
    //! \return Enum denoting physics policy.
    static ctr::PhysicsType type() noexcept { return ctr::PhysicsType::EULER; }
};

} // dg::

} // inciter::

#endif // CompFlowPhysicsDGMultiSpeciesEuler_h
