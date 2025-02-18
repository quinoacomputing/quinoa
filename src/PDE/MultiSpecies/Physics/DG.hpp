// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Physics/DG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics configurations for multi-material compressible flow using
    discontinuous Galerkin finite element methods
  \details   This file configures all Physics policy classes for multi-material
    compressible flow implementied using discontinuous Galerkin finite element
    discretizations, defined in PDE/MultiSpecies/DGMultiSpecies.h.

    General requirements on MultiSpecies Physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::EULER;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for Physics policies.
*/
// *****************************************************************************
#ifndef MultiSpeciesPhysicsDG_h
#define MultiSpeciesPhysicsDG_h

#include <brigand/sequences/list.hpp>

#include "DGEuler.hpp"

namespace inciter {
namespace dg {

//! MultiSpecies Physics policies implemented using discontinuous Galerkin
using MultiSpeciesPhysics = brigand::list< MultiSpeciesPhysicsEuler
                                         >;

} // dg::
} // inciter::

#endif // MultiSpeciesPhysicsDG_h
