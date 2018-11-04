// *****************************************************************************
/*!
  \file      src/PDE/PDEFactory.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Differential equations factory
  \details   This file declares the type for a differential equations factory.
*/
// *****************************************************************************
#ifndef PDEFactory_h
#define PDEFactory_h

#include <map>
#include <set>

#include <brigand/algorithms/for_each.hpp>

#include "CGPDE.h"
#include "DGPDE.h"
#include "Factory.h"
#include "SystemComponents.h"
#include "Inciter/Options/PDE.h"

namespace inciter {

//! \brief Factory for PDEs using continuous Galerkin discretization storing
//!   keys associated to their constructors
using CGFactory =
  std::map< ctr::PDEKey, std::function< CGPDE(const tk::ctr::ncomp_type&) > >;

//! \brief Factory for PDEs using discontinuous Galerkin discretization storing
//!   keys associated to their constructors
using DGFactory =
  std::map< ctr::PDEKey, std::function< DGPDE(const tk::ctr::ncomp_type&) > >;

//! \brief Function object for registering a partial differential equation
//!   into the partial differential equation factory
//! \details This functor is repeatedly called by MPL's cartesian_product,
//!   sweeping all combinations of the partial differential equation
//!   policies. The purpose of template template is to simplify client code
//!   as that will not have to specify the template arguments of the
//!   template argument (the policies of Eq), since we can figure it out
//!   here. See also http://stackoverflow.com/a/214900. The template
//!   argument Eq specifies a "child" class that is used polymorphically
//!   with a "base" class modeling a concept defined in the base. The base
//!   is given by the template argument PDE. The template argument Factory
//!   specifies which factory to store the registered and configured child
//1   PDE.
template< template< class, class > class Eq, class Factory, class PDE >
struct registerPDE {
  //! Need to store the reference to factory we are registering into
  Factory& factory;
  //! Need to store which differential equation we are registering
  const ctr::PDEType type;
  //! Constructor, also count number of unique equation types registered
  //! \param[in,out] f Factory into which to register PDE
  //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
  //! \param[in] eqTypes Equation type counters
  explicit registerPDE( Factory& f,
                        ctr::PDEType t,
                        std::set< ctr::PDEType >& eqTypes ) :
    factory( f ), type( t ) { eqTypes.insert( t ); }
  //! \brief Function call operator called with tk::cartesian_product for
  //!   each unique sequence of policy combinations
  template< typename U > void operator()( brigand::type_<U> ) {
    // Get problem policy: first type of brigand::list U
    using Physics = typename brigand::front< U >;
    // Get problem policy: last type of brigand::list U
    using Problem = typename brigand::back< U >;
    // Build differential equation key
    ctr::PDEKey key{ type, Physics::type(), Problem::type() };
    // Register equation (with policies given by brigand::list U) into factory
    tk::recordModelLate< PDE, Eq< Physics, Problem > >
                       ( factory, key, static_cast<tk::ctr::ncomp_type>(0) );
  }
};

//! Wrapper of registerPDE specialized for registering CG PDEs
//! \details The sole reason for this functor is to simplify client-code
//!   calling registerPDE specialized to CG PDEs
template< template< class, class > class Eq >
struct registerCG : registerPDE< Eq, CGFactory, CGPDE > {
  //! Delegate constructor to base and specialize to CG
  //! \param[in] f Factory to register to
  //! \param[in] eqtypes Counters for equation types in factory
  //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
  explicit registerCG( CGFactory& f,
                       std::set< ctr::PDEType >& eqtypes,
                       ctr::PDEType t ) :
    registerPDE< Eq, CGFactory, CGPDE >( f, t, eqtypes ) {}
};

//! Wrapper of registerPDE specialized for registering DG PDEs
//! \details The sole reason for this functor is to simplify client-code
//!   calling registerPDE specialized to DG PDEs
template< template< class, class > class Eq >
struct registerDG : registerPDE< Eq, DGFactory, DGPDE > {
  //! Delegate constructor to base and specialize to CG
  //! \param[in] f Factory to register to
  //! \param[in] eqtypes Counters for equation types in factory
  //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
  explicit registerDG( DGFactory& f,
                       std::set< ctr::PDEType >& eqtypes,
                       ctr::PDEType t ) :
    registerPDE< Eq, DGFactory, DGPDE >( f, t, eqtypes ) {}
};

//! \brief Convert and return values from vector as string
//! \param[in] v Vector whose components to return as a string
//! \return Concatenated string of values read from a vector
template< typename V >
std::string parameters( const V& v ) {
  std::stringstream s;
  s << "{ ";
  for (auto p : v) s << p << ' ';
  s << "}";
  return s.str();
}

} // inciter::

#endif // PDEFactory_h
