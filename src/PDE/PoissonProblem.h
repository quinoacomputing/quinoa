// *****************************************************************************
/*!
  \file      src/PDE/PoissonProblem.h
  \author    J. Bakosi
  \date      Mon 18 Jul 2016 03:46:24 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for the Poisson equation
  \details   This file defines policy classes for the Poisson
    partial differential equation, defined in PDE/Poisson.h.

    General requirements on Poisson partial differential equation problem policy
    classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::ProblemType type() noexcept {
          return ctr::ProblemType::POISSON;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef PoissonProblem_h
#define PoissonProblem_h

#include <boost/mpl/vector.hpp>

#include "Macro.h"
#include "Types.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! \brief Poisson PDE problem: inhomogeneous Poisson equation with Dirichlet
//! and Neumann boundary conditions
class PoissonProblemDirNeu {

  public:
    //! Do error checking on PDE parameters
    //! \param[in] deck Input deck
    //! \param[in] e Equation system index, i.e., which Poisson equation
    //!   system we operate on among the systems of PDEs
    template< class eq >
    static void errchk( const ctr::InputDeck& deck,
                        tk::ctr::ncomp_type e,
                        tk::ctr::ncomp_type ncomp )
    {
      IGNORE(deck);
      IGNORE(e);
      IGNORE(ncomp);
    }

    //! Set initial conditions
    //! \param[in] deck Input deck
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which Poisson equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    template< class eq >
    static void init( const ctr::InputDeck& deck,
                      const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::MeshNodes& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type ncomp,
                      tk::ctr::ncomp_type offset,
                      tk::real t )
    {
      IGNORE(deck);
      IGNORE(e);
      IGNORE(t);
      const auto& x = coord[0];
      for (ncomp_t c=0; c<ncomp; ++c)
        for (ncomp_t i=0; i<x.size(); ++i)
          unk( i, c, offset ) = 0.0;
    }

    //! Problem type enum accessor
    //! \return Problem type as enum
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::DIR_NEU; }
};

//! List of all Poisson PDE's problem policies
using PoissonProblems = boost::mpl::vector< PoissonProblemDirNeu >;

} // inciter::

#endif // PoissonProblem_h
