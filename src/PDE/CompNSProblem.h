// *****************************************************************************
/*!
  \file      src/PDE/CompNSProblem.h
  \author    J. Bakosi
  \date      Wed 20 Jul 2016 08:49:21 AM MDT
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

    //! Set initial conditions
    //! \param[in] deck Input deck
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which compressible
    //!   Navier-Stokes equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    template< class eq >
    static void init( const ctr::InputDeck& deck,
                      const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::MeshNodes& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset )
    {
      IGNORE(deck);
      IGNORE(e);
      const auto& x = coord[0];
      for (ncomp_t i=0; i<x.size(); ++i) {
         unk( i, 0, offset ) = 1.0;     // density
         unk( i, 1, offset ) = 0.0;     // density * velocity
         unk( i, 2, offset ) = 0.0;
         unk( i, 3, offset ) = -1.0;
         unk( i, 4, offset ) = 1.0;     // energy
      }
    }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::USER_DEFINED; }
};

//! List of all CompNS problem policies
using CompNSProblems = boost::mpl::vector< CompNSProblemUserDefined >;

} // inciter::

#endif // CompNSProblem_h
