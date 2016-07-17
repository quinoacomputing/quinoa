// *****************************************************************************
/*!
  \file      src/PDE/CompNSProblem.h
  \author    J. Bakosi
  \date      Wed 03 Feb 2016 04:16:54 PM MST
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for the compressible Navier-Stokes equation
  \details   This file defines policy classes for the compressible Navier-Stokes
    equations, defined in PDE/CompNS.h.

    General requirements on Navier-Stokes equations problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::ProblemType type() noexcept {
          return ctr::ProblemType::USER_DEFINED;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef CompNSProblem_h
#define CompNSProblem_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompNS system of PDEs problem: user defined
class CompNSProblemUserDefined {
  public:
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::USER_DEFINED; }
};

//! List of all CompNS problem policies
using CompNSProblems = boost::mpl::vector< CompNSProblemUserDefined
                                        >;

} // inciter::

#endif // CompNSProblem_h
