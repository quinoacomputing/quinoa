// *****************************************************************************
/*!
  \file      src/PDE/AdvDiffPhysics.h
  \author    J. Bakosi
  \date      Fri 22 Jul 2016 02:49:02 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Physics configurations for the advection-diffusion equations
  \details   This file defines policy classes for advection-diffusion equations,
    defined in PDE/AdvDiff.h.

    General requirements on advection-diffusion equation physics policy classes:

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
#ifndef AdvDiffPhysics_h
#define AdvDiffPhysics_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

//! AdvDiff system of PDEs problem: basic (advection + diffusion)
class AdvDiffPhysicsBase {
  public:


    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::BASE; }
};

//! List of all AdvDiff problem policies
using AdvDiffPhysics = boost::mpl::vector< AdvDiffPhysicsBase >;

} // inciter::

#endif // AdvDifAdvDiff_h
