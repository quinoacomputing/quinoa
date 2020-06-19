// *****************************************************************************
/*!
  \file      src/PDE/CGPDE.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions common to ALECG
  \details   Functions common to ALECG.
*/
// *****************************************************************************

#include <array>
#include <vector>
#include <unordered_map>

#include "Vector.hpp"
#include "DerivedData.hpp"
#include "Exception.hpp"
#include "Around.hpp"
#include "Fields.hpp"
#include "CGPDE.hpp"
#include "FunctionPrototypes.hpp"

namespace inciter {
namespace cg {

std::vector< tk::real >
solinc( tk::ncomp_t system, tk::ncomp_t ncomp, tk::real x, tk::real y,
        tk::real z, tk::real t, tk::real dt, tk::SolutionFn solution )
// *****************************************************************************
// Evaluate the increment from t to t+dt of the analytical solution at (x,y,z)
// for all components
//! \param[in] system Equation system index, i.e., which equation system we
//!   operate on among the systems of PDEs configured by the user
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution increment starting from
//! \param[in] dt Time increment at which evaluate the solution increment to
//! \param[in] solution Function used to evaluate the solution
//! \return Increment in values of all components evaluated at (x,y,z,t+dt)
// *****************************************************************************
{
  int inbox = 0;
  auto st1 = solution( system, ncomp, x, y, z, t, inbox );
  auto st2 = solution( system, ncomp, x, y, z, t+dt, inbox );

  std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                  []( tk::real s, tk::real& d ){ return d -= s; } );

  return st2;
}

} // cg::
} // inciter::
