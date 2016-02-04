//******************************************************************************
/*!
  \file      src/PDE/EulerProblem.h
  \author    J. Bakosi
  \date      Wed 03 Feb 2016 04:16:54 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Problem configurations for the Euler equation
  \details   This file defines policy classes for the Euler system of artial
    differential equations, defined in PDE/Euler.h.

    General requirements on Euler equations problem policy classes:

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
//******************************************************************************
#ifndef EulerProblem_h
#define EulerProblem_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! Euler system of PDEs problem: user defined
class EulerProblemUserDefined {
  public:
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::USER_DEFINED; }
};

//! List of all Euler system of PDE's problem policies
using EulerProblems = boost::mpl::vector< EulerProblemUserDefined
                                        >;

} // inciter::

#endif // EulerProblem_h
