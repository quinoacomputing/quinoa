// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/RiemannFactory.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register available Riemann solvers for multimaterial compressible
             hydrodynamics into a factory
  \details   Register available Riemann solvers for multimaterial compressible
             hydrodynamics into a factory.
*/
// *****************************************************************************
#ifndef RiemannSolverFactory_h
#define RiemannSolverFactory_h

#include <map>
#include <functional>

#include "NoWarning/value_factory.hpp"

#include "Riemann/RiemannSolver.hpp"
#include "Inciter/Options/Flux.hpp"

namespace inciter {

//! Factory for Riemann solvers
//! \details This factory is used to store the constructors as a
//!   std::function of specific Riemann solvers that can be invoked at a
//!   later time compared to the point where the map is populated. The key
//!   is an enum, uniquely idenfitying a specific Riemann solver. The value
//!   is std::function storing a constructor to be invoked. The type of
//!   object stored in std::function is a generic (base) class constructor,
//!   which provides a polymorphyic interface (overridable functions) that
//!   specific (child) Riemann solvers override, yielding runtime polymorphism.
using MultiMatRiemannFactory =
  std::map< ctr::FluxType, std::function< RiemannSolver() > >;

//! Functor to register a Riemann solver into the Riemann solver factory
struct registerRiemannSolver {
  //! Factory to which to register the Riemann solver
  MultiMatRiemannFactory& factory;
  //! Constructor
  //! \param[in] f Factory
  explicit registerRiemannSolver( MultiMatRiemannFactory& f ) : factory( f ) {}
  //! \brief Function call operator templated on the type that implements
  //!   a specific Riemann solver
  template< typename U > void operator()( brigand::type_<U> ) {
     // Function object holding the (default) constructor to be called later
     // without bound arguments, since all specific Riemann solvers'
     // constructors are compiler-generated (default) constructors, and thus
     // taking no arguments.
     std::function< U() > c = boost::value_factory< U >();
     // Associate constructor function object to flux type in factory
     factory.emplace( U::type(),
       std::bind(boost::value_factory< RiemannSolver >(), std::move(c)) );
  }
};

//! Register available Riemann solvers into a factory
MultiMatRiemannFactory multimatRiemannSolvers();

} // inciter::

#endif // RiemannSolverFactory_h
