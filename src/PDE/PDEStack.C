//******************************************************************************
/*!
  \file      src/PDE/PDEStack.C
  \author    J. Bakosi
  \date      Sun 03 Apr 2016 02:26:31 PM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Stack of partial differential equations
  \details   This file defines class PDEStack, which implements various
    functionality related to registering and instantiating partial differential
    equation types. Registration and instantiation use a partial differential
    equation factory, which is a std::map (an associative container),
    associating unique partial differential equation keys to their constructor
    calls. For more details, see the in-code documentation of the constructor.
*/
//******************************************************************************

#include <boost/mpl/cartesian_product.hpp>

#include "PDEStack.h"
#include "Tags.h"
#include "SystemComponents.h"
#include "Inciter/Options/Problem.h"

#include "AdvDiff.h"
#include "Euler.h"

#include "AdvDiffProblem.h"
#include "EulerProblem.h"

using inciter::PDEStack;

PDEStack::PDEStack()
//******************************************************************************
//  Constructor: register all partial differential equations into factory
//! \details This constructor consists of several blocks, each registering a
//!   potentially large number of entries in the partial differential equation
//!   factory, m_factory, which is of type inciter::PDEFactory, a std::map. At
//!   this time, each type of partial differential equation can be configured to
//!   use a unique _problem policy_. (More types of policies will most likely
//!   come in the future.) Policy classes are template arguments to the partial
//!   differential equation classes and influence their behavior in a different
//!   way, abstracting away certain functions, e.g., how to set problem-specific
//!   initial and/or boundary conditions and how to update their coefficients
//!   during time integration. For more information on policy-based design, see
//!   http://en.wikipedia.org/wiki/Policy-based_design. This abstraction allows
//!   [separation of concerns](http://en.wikipedia.org/wiki/Separation_of_concerns).
//!
//!   Since the functionality of the policies are orthogonal to each other,
//!   i.e., they do not depend on each other or their host (the partial
//!   differential equation class), a Cartesian product of combinations are
//!   possible, depending on which policies are selected. _This constructor
//!   registers all possible combinations of policies for all available
//!   differential equations._ By _register_, we mean, an entry is recorded in
//!   an associative container, a std::map, that associates a lightweight key of
//!   type inciter::ctr::PDEKey, consisting of only an enum for each policy
//!   type, to an std::function object that holds the constructor bound to its
//!   arguments corresponding to a particular partial differential equation +
//!   policies combination. Note that registering these entries in the map does
//!   not invoke the constructors. The mapped value simply stores how the
//!   constructors should be invoked at a later time. At some point later,
//!   based on user input, we then instantiate only the partial differential
//!   equations (and only those configurations) that are requested by the user.
//!
//!   Since all partial differential equation types (registered in the factory)
//!   "inherit" from a common "base", client-code is unform and generic, and
//!   thus immune to changes in the inner workings of the particular partial
//!   differential equations as long as they fullfill certain concepts, i.e.,
//!   implement certain member functinos, enforced by the _common base_, PDE.
//!   The words "inherit and "base" are quoted here, because the common base
//!   does not use inheritance in the normal OOP sense and does not use
//!   reference semantics, i.e., pointers, visible to client-code either. The
//!   relationship is more of a _models a_-type, which simplifies client-code
//!   and allows for the benfits of runtime inheritance with value-semantics
//!   which is less error prone and easier to read. See more about the
//!   _models-a_ relationship and its implementation in PDE/PDE.h.
//!
//!   The design discussed above allows the registration, instantiation, and
//!   use of the partial differential equations to be generic, which eliminates
//!   a lot of boiler-plate code and makes client-code uniform.
//!
//!   _Details of registration using mpl::cartesian_product:_
//!
//!   The template argument to mpl::cartesian_product requires a sequence of
//!   sequences of types. We use vector of vectors of types, listing all
//!   possible policies. The constructor argument to mpl::cartesian_product is a
//!   functor that is to be applied to all combinations. mpl::cartesian_product
//!   will then create all possible combinations of these types and call the
//!   user-supplied functor with each type of the created sequence as a template
//!   parameter. The user-supplied functor here is registerPDE, which, i.e.,
//!   its constructor call, needs a single template argument, a class templated
//!   on policy classes. This is the partial differential equation class to be
//!   configured by selecting policies and to be registered. The arguments to
//!   registerPDE's constructor are the factory, the enum denoting the
//!   differential equation type, and a reference to a variable of type
//!   std::set< ctr::PDEType >, which is only used internally to PDEStack
//!   for counting up the number of unique differential equation types
//!   registered, used for diagnostics purposes.
//! \author J. Bakosi
//******************************************************************************
{
  namespace mpl = boost::mpl;

  // Advection-diffusion PDE
  // Construct vector of vectors for all possible policies for PDE
  using AdvDiffPolicies = mpl::vector< AdvDiffProblems >;
  // Register PDE for all combinations of policies
  mpl::cartesian_product< AdvDiffPolicies >(
    registerPDE< AdvDiff >( m_factory, ctr::PDEType::ADV_DIFF, m_eqTypes ) );

  // Euler system of PDEs
  // Construct vector of vectors for all possible policies for PDE
  using EulerPolicies = mpl::vector< EulerProblems >;
  // Register PDE for all combinations of policies
  mpl::cartesian_product< EulerPolicies >(
    registerPDE< Euler >( m_factory, ctr::PDEType::EULER, m_eqTypes ) );
}

std::vector< inciter::PDE >
PDEStack::selected() const
//******************************************************************************
//  Instantiate all selected partial differential equations
//! \return std::vector of instantiated partial differential equation objects
//! \author J. Bakosi
//******************************************************************************
{
  std::map< ctr::PDEType, ncomp_t > cnt;    // count PDEs per type
  std::vector< PDE > pdes;                      // will store instantiated PDEs

  for (const auto& d : g_inputdeck.get< tag::selected, tag::pde >()) {
    if (d == ctr::PDEType::ADV_DIFF)
      pdes.push_back( createPDE< tag::advdiff >( d, cnt ) );
    else if (d == ctr::PDEType::EULER)
      pdes.push_back( createPDE< tag::euler >( d, cnt ) );
    else Throw( "Can't find selected PDE" );
  }

  return pdes;
}

std::vector< std::vector< std::pair< std::string, std::string > > >
PDEStack::info() const
//******************************************************************************
//  Return information on all selected partial differential equations
//! \return A vector of vector of pair of strings, containing the configuration
//!   for each selected partial differential equation
//! \author J. Bakosi
//******************************************************************************
{
  std::map< ctr::PDEType, ncomp_t > cnt; // count PDEs per type
  // will store info on all differential equations selected
  std::vector< std::vector< std::pair< std::string, std::string > > > info;

  for (const auto& d : g_inputdeck.get< tag::selected, tag::pde >()) {
    if (d == ctr::PDEType::ADV_DIFF)
      info.emplace_back( infoAdvDiff( cnt ) );
    else if (d == ctr::PDEType::EULER)
      info.emplace_back( infoEuler( cnt ) );
    else Throw( "Can't find selected PDE" );
  }

  return info;
}

std::vector< std::pair< std::string, std::string > >
PDEStack::infoAdvDiff( std::map< ctr::PDEType, ncomp_t >& cnt ) const
//******************************************************************************
//  Return information on the advection-diffusion PDE
//! \param[inout] cnt std::map of counters for all partial differential equation
//!   types
//! \return vector of string pairs describing the PDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::PDEType::ADV_DIFF ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::PDE().name( ctr::PDEType::ADV_DIFF ), "" );
  info.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, tag::advdiff, tag::problem >()[c] ) );
  info.emplace_back( "start offset in unknowns array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::advdiff >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::advdiff >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "coeff diffusivity [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::advdiff, tag::diffusivity >().at(c) ) );
  info.emplace_back( "coeff u0 [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::advdiff, tag::u0 >().at(c) ) );
  info.emplace_back( "coeff lambda [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::advdiff, tag::lambda >().at(c) ) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
PDEStack::infoEuler( std::map< ctr::PDEType, ncomp_t >& cnt ) const
//******************************************************************************
//  Return information on the Euler system of PDEs
//! \param[inout] cnt std::map of counters for all partial differential equation
//!   types
//! \return vector of string pairs describing the PDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::PDEType::EULER ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::PDE().name( ctr::PDEType::EULER ), "" );
  info.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, tag::euler, tag::problem >()[c] ) );
  info.emplace_back( "start offset in unknowns array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::euler >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::euler >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );

  return info;
}
