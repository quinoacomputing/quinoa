// *****************************************************************************
/*!
  \file      src/PDE/PoissonPhysics.h
  \author    J. Bakosi
  \date      Fri 22 Jul 2016 02:49:02 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Physics configurations for the Poisson equation
  \details   This file defines policy classes for the Poisson equation, defined
    in PDE/Poisson.h.

    General requirements on Poisson equations physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::BASE;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef PoissonPhysics_h
#define PoissonPhysics_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

//! Poisson system of PDEs problem: basic (advection + diffusion)
class PoissonPhysicsBase {
  public:


    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::BASE; }
};

//! List of all Poisson problem policies
using PoissonPhysics = boost::mpl::vector< PoissonPhysicsBase >;

} // inciter::

#endif // AdvDifPoisson_h
