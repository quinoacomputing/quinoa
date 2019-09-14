// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/InterfaceAdvection.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the multi-material compressible flow
    equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined under PDE/MultiMat. See
    PDE/MultiMat/Problem.hpp for general requirements on Problem policy classes
    for MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatProblemInterfaceAdvection_h
#define MultiMatProblemInterfaceAdvection_h

#include <string>
#include <unordered_set>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! MultiMat system of PDEs problem: interface advection
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class MultiMatProblemInterfaceAdvection {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::multimat;

  public:
    //! Evaluate analytical solution at (x,y,z,t) for all components
    static tk::SolutionFn::result_type
    solution( ncomp_t system,
              ncomp_t ncomp,
              tk::real x,
              tk::real y,
              tk::real /*z*/,
              tk::real t );

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    static std::vector< tk::real >
    solinc( ncomp_t system,
            ncomp_t ncomp,
            tk::real x,
            tk::real y,
            tk::real z,
            tk::real t,
            tk::real dt );

    //! Compute and return source term for interface advection
    static tk::SrcFn::result_type
    src( ncomp_t, ncomp_t ncomp, tk::real, tk::real, tk::real, tk::real );

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    static void side( std::unordered_set< int >& conf );

    //! Return field names to be output to file
    static std::vector< std::string > fieldNames( ncomp_t );

    //! Return field output going to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( ncomp_t system,
                 ncomp_t /*ncomp*/,
                 ncomp_t offset,
                 tk::real t,
                 tk::real,
                 const std::vector< tk::real >&,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 tk::Fields& U );

    //! Return names of integral variables to be output to diagnostics file
    static std::vector< std::string > names( ncomp_t );

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::INTERFACE_ADVECTION; }
};

} // inciter::

#endif // MultiMatProblemInterfaceAdvection_h
