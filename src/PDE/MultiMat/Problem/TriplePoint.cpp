// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/TriplePoint.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the multi-material flow equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/MultiMat.hpp. See
    PDE/MultiMat/Problem.hpp for general requirements on Problem policy classes
    for MultiMat.
*/
// *****************************************************************************

#include "TriplePoint.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemTriplePoint;

tk::SolutionFn::result_type
MultiMatProblemTriplePoint::solution( ncomp_t system,
                                            ncomp_t ncomp,
                                            tk::real x,
                                            tk::real y,
                                            tk::real,
                                            tk::real,
                                            int& )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which multi-material
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::SolutionFn
//! \details This function only initializes the triple point problem,
//!   but does not actually give the analytical solution at time greater than 0.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 12, "Number of scalar components must be 12" );

  auto nmat =
    g_inputdeck.get< tag::param, eq, tag::nmat >()[system];

  std::vector< tk::real > s(ncomp, 0.0), r(nmat, 0.0);
  tk::real p, u, v, w, temp;
  auto alphamin = 1.0e-12;

  // velocity
  u = 0.0;
  v = 0.0;
  w = 0.0;

  if (x<1.0) {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = 1.0-2.0*alphamin;
    s[volfracIdx(nmat, 1)] = alphamin;
    s[volfracIdx(nmat, 2)] = alphamin;
    // pressure
    p = 1.0;
    // temperature
    temp = 4.3554007e-4;
  }
  else {
    if (y<1.5) {
      // volume-fraction
      s[volfracIdx(nmat, 0)] = alphamin;
      s[volfracIdx(nmat, 1)] = 1.0-2.0*alphamin;
      s[volfracIdx(nmat, 2)] = alphamin;
      // pressure
      p = 0.1;
      // temperature
      temp = 3.4843206e-4;
    }
    else {
      // volume-fraction
      s[volfracIdx(nmat, 0)] = alphamin;
      s[volfracIdx(nmat, 1)] = alphamin;
      s[volfracIdx(nmat, 2)] = 1.0-2.0*alphamin;
      // pressure
      p = 0.1;
      // temperature
      temp = 3.4843206e-4;
    }
  }

  for (std::size_t k=0; k<nmat; ++k)
  {
    // density
    r[k] = eos_density< eq >( system, p, temp, k );
    // partial density
    s[densityIdx(nmat, k)] = s[volfracIdx(nmat, k)]*r[k];
    // total specific energy
    s[energyIdx(nmat, k)] = s[volfracIdx(nmat, k)]*
      eos_totalenergy< eq >( system, r[k], u, v, w, p, k );
  }

  return s;
}

tk::SrcFn::result_type
MultiMatProblemTriplePoint::src( ncomp_t, ncomp_t ncomp, tk::real,
                                       tk::real, tk::real, tk::real )
// *****************************************************************************
//  Compute and return source term for manufactured solution
//! \param[in] ncomp Number of scalar components in this PDE system
//! \return Array of reals containing the source for all components
//! \note The function signature must follow tk::SrcFn
// *****************************************************************************
{
  std::vector< tk::real > s( ncomp, 0.0 );

  return s;
}

std::vector< std::string >
MultiMatProblemTriplePoint::fieldNames( ncomp_t )
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  auto nmat =
    g_inputdeck.get< tag::param, eq, tag::nmat >()[0];

  return MultiMatFieldNames(nmat);
}

std::vector< std::vector< tk::real > >
MultiMatProblemTriplePoint::fieldOutput(
  ncomp_t system,
  ncomp_t,
  ncomp_t offset,
  std::size_t nunk,
  tk::real,
  tk::real,
  const std::vector< tk::real >&,
  const std::array< std::vector< tk::real >, 3 >&,
  tk::Fields& U,
  const tk::Fields& P )
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  // number of degrees of freedom
  const std::size_t rdof =
    g_inputdeck.get< tag::discr, tag::rdof >();

  // number of materials
  auto nmat =
    g_inputdeck.get< tag::param, eq, tag::nmat >()[system];

  return MultiMatFieldOutput(system, nmat, offset, nunk, rdof, U, P);
}

std::vector< std::string >
MultiMatProblemTriplePoint::names( ncomp_t )
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
