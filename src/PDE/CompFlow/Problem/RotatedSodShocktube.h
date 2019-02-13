// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/RotatedSodShocktube.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Problem configuration for rotated Sod's shock-tube
  \details   This file defines a policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problems.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemRotatedSodShocktube_h
#define CompFlowProblemRotatedSodShocktube_h

#include <string>
#include <unordered_set>
#include <cmath>

#include "Types.h"
#include "Vector.h"
#include "Inciter/Options/Problem.h"
#include "SodShocktube.h"

namespace inciter {

//! CompFlow system of PDEs problem: rotated Sod shock-tube
//! \see G. A. Sod. A Survey of Several Finite Difference Methods for Systems of
//!   Nonlinear Hyperbolic Conservation Laws. J. Comput. Phys., 27:1–31, 1978.
class CompFlowProblemRotatedSodShocktube : public CompFlowProblemSodShocktube {

  private:
    //! Rotate vector about X axis by -45 degress
    //! \param[in] v Vector to rotate
    //! \return Rotated vector
    static std::array< tk::real, 3 >
    rotateX( const std::array< tk::real, 3 >& v ) {
      using std::cos;  using std::sin;
      tk::real angle = -45.0*M_PI/180.0;
      std::array< std::array< tk::real, 3 >, 3 >
        R{{ {{ 1.0,         0.0,          0.0 }},
            {{ 0.0,   cos(angle), -sin(angle) }},
            {{ 0.0,   sin(angle),  cos(angle) }} }};
      return {{ tk::dot(R[0],v), tk::dot(R[1],v), tk::dot(R[2],v) }};
    }

    //! Rotate vector about Y axis by -45 degress
    //! \param[in] v Vector to rotate
    //! \return Rotated vector
    static std::array< tk::real, 3 >
    rotateY( const std::array< tk::real, 3 >& v ) {
      using std::cos;  using std::sin;
      tk::real angle = -45.0*M_PI/180.0;
      std::array< std::array< tk::real, 3 >, 3 >
        R{{ {{ cos(angle),  0.0, sin(angle) }},
            {{ 0.0,         1.0,        0.0 }},
            {{ -sin(angle), 0.0, cos(angle) }} }};
      return {{ tk::dot(R[0],v), tk::dot(R[1],v), tk::dot(R[2],v) }};
    }

    //! Rotate vector about Z axis by -45 degress
    //! \param[in] v Vector to rotate
    //! \return Rotated vector
    static std::array< tk::real, 3 >
    rotateZ( const std::array< tk::real, 3 >& v ) {
      using std::cos;  using std::sin;
      tk::real angle = -45.0*M_PI/180.0;
      std::array< std::array< tk::real, 3 >, 3 >
        R{{ {{ cos(angle), -sin(angle), 0.0 }},
            {{ sin(angle),  cos(angle), 0.0 }},
            {{ 0.0,         0.0,        1.0 }} }};
      return {{ tk::dot(R[0],v), tk::dot(R[1],v), tk::dot(R[2],v) }};
    }

  public:
    //! Evaluate analytical solution at (x,y,0) for all components
    //! \param[in] system Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Time at which to evaluate the solution
    //! \return Values of all components evaluated at (x,y,0)
    //! \note The function signature must follow tk::SolutionFn
    static tk::SolutionFn::result_type
    solution( ncomp_t system, ncomp_t ncomp, tk::real x, tk::real y, tk::real z,
              tk::real t ) {
      // Assume the domain is rotated by 45 degrees about the X, Y, and then Z
      // axis compared to the original tube with largest dimension in X
      auto c = rotateX( rotateY( rotateZ( {{ x, y, z }} ) ) );
      return CompFlowProblemSodShocktube::solution( system, ncomp,
                                                    c[0], c[1], c[2], t );
    }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::ROTATED_SOD_SHOCKTUBE; }
};

} // inciter::

#endif // CompFlowProblemRotatedSodShocktube_h
