// *****************************************************************************
/*!
  \file      src/PDE/PoissonPhysics.h
  \author    J. Bakosi
  \date      Mon 29 Aug 2016 01:12:25 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Physics configurations for the Poisson equation
  \details   This file defines policy classes for the Poisson equation, defined
    in PDE/Poisson.h.

    General requirements on Poisson equations physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::LAPLACE;
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

//! Poisson system of PDEs problem: Laplace
class PoissonPhysicsLaplace {
  public:


    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::LAPLACE; }
};

//! List of all Poisson problem policies
using PoissonPhysics = boost::mpl::vector< PoissonPhysicsLaplace >;

} // inciter::

#endif // AdvDifPoisson_h
