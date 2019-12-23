// *****************************************************************************
/*!
  \file      src/PDE/ProblemCommon.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Common functionality for Problem classes in inciter
  \details   This file declares some common functionality for Problem
    definitions in inciter, generally used across all types of PDE definitions.
*/
// *****************************************************************************
#ifndef ProblemCommon_h
#define ProblemCommon_h

#include "Types.hpp"
#include "FunctionPrototypes.hpp"

namespace inciter {

//! \brief Evaluate the increment from t to t+dt of an analytical solution at
//!   (x,y,z) for all components
std::vector< tk::real >
solinc( tk::ncomp_t system, tk::ncomp_t ncomp, tk::real x, tk::real y,
        tk::real z, tk::real t, tk::real dt, tk::SolutionFn solution );

} // inciter::

#endif // ProblemCommon.h
