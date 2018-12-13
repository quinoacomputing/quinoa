// *****************************************************************************
/*!
  \file      src/PDE/MultiMatCompFlow/Problem/VorticalFlow.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Problem configuration for the multi-material compressible flow
    equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined under PDE/MultiMatCompFlow. See
    PDE/MultiMatCompFlow/Problem.h for general requirements on Problem policy
    classes for MultiMatCompFlow.
*/
// *****************************************************************************
#ifndef MultiMatCompFlowProblemVorticalFlow_h
#define MultiMatCompFlowProblemVorticalFlow_h

#include <string>
#include <unordered_set>

#include "Types.h"
#include "FunctionPrototypes.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! MultiMatCompFlow system of PDEs problem: vortical flow
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class MultiMatCompFlowProblemVorticalFlow {

  private:
    using ncomp_t = tk::ctr::ncomp_type;
    using eq = tag::multimat_compflow;

  public:
    //! Evaluate analytical solution at (x,y,z) for all components
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \return Values of all components evaluated at (x,y,z)
    //! \note The function signature must follow tk::SolutionFn
    static tk::SolutionFn::result_type
    solution( ncomp_t, ncomp_t ncomp, tk::real, tk::real, tk::real, tk::real )
    {
      return std::vector< tk::real >( ncomp, 0.0 );
    }

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \return Increment in values of all components: all zero for this problem
    static std::vector< tk::real >
    solinc( ncomp_t, ncomp_t ncomp, tk::real, tk::real, tk::real, tk::real,
            tk::real )
    {
      return std::vector< tk::real >( ncomp, 0.0 );
    }

    //! Compute and return source term for vortical flow manufactured solution
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \return Array of reals containing the source for all components
    //! \note The function signature must follow tk::SrcFn
    static tk::SrcFn::result_type
    src( ncomp_t, ncomp_t ncomp, tk::real, tk::real, tk::real, tk::real ) {
      return std::vector< tk::real >( ncomp, 0.0 );
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::bcdir;
      for (const auto& s : g_inputdeck.get< param, eq, bcdir >())
        for (const auto& i : s)
          conf.insert( std::stoi(i) );
    }

    //! Return field names to be output to file
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > fieldNames( ncomp_t ncomp ) {
      return std::vector< std::string >( ncomp, "fieldvar" );
    }

    //! Return field output going to file
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( ncomp_t,
                 ncomp_t ncomp,
                 ncomp_t,
                 tk::real,
                 tk::real,
                 const std::vector< tk::real >&,
                 const std::array< std::vector< tk::real >, 3 >&,
                 tk::Fields& U )
    {
      return std::vector< std::vector< tk::real > >
                        ( ncomp, std::vector< tk::real >( U.nunk(), 0.0 ) );
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names( ncomp_t ncomp ) {
      return std::vector< std::string >( ncomp, "diagvar" );
    }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::VORTICAL_FLOW; }
};

} // inciter::

#endif // MultiMatCompFlowProblemVorticalFlow_h
