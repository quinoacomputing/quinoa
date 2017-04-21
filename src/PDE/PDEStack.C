// *****************************************************************************
/*!
  \file      src/PDE/PDEStack.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Stack of partial differential equations
  \details   This file defines class PDEStack, which implements various
    functionality related to registering and instantiating partial differential
    equation types. Registration and instantiation use a partial differential
    equation factory, which is a std::map (an associative container),
    associating unique partial differential equation keys to their constructor
    calls. For more details, see the in-code documentation of the constructor.
*/
// *****************************************************************************

#include "NoWarning/cartesian_product.h"

#include "PDEStack.h"
#include "Tags.h"
#include "SystemComponents.h"
#include "Inciter/Options/Problem.h"

#include "Transport.h"
#include "CompFlow.h"

#include "TransportPhysics.h"
#include "CompFlowPhysics.h"

#include "TransportProblem.h"
#include "CompFlowProblem.h"

using inciter::PDEStack;

PDEStack::PDEStack() : m_factory(), m_eqTypes()
// *****************************************************************************
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
// *****************************************************************************
{
  namespace mpl = boost::mpl;

  // Transport PDE
  // Construct vector of vectors for all possible policies for PDE
  using TransportPolicies = mpl::vector< TransportPhysics, TransportProblems >;
  // Register PDE for all combinations of policies
  mpl::cartesian_product< TransportPolicies >(
    registerPDE< Transport >( m_factory, ctr::PDEType::TRANSPORT, m_eqTypes ) );

  // Compressible flow system of PDEs
  // Construct vector of vectors for all possible policies for PDE
  using CompFlowPolicies = mpl::vector< CompFlowPhysics, CompFlowProblems >;
  // Register PDE for all combinations of policies
  mpl::cartesian_product< CompFlowPolicies >(
    registerPDE< CompFlow >( m_factory, ctr::PDEType::COMPFLOW, m_eqTypes ) );
}

std::vector< inciter::PDE >
PDEStack::selected() const
// *****************************************************************************
//  Instantiate all selected partial differential equations
//! \return std::vector of instantiated partial differential equation objects
//! \author J. Bakosi
// *****************************************************************************
{
  std::map< ctr::PDEType, ncomp_t > cnt;    // count PDEs per type
  std::vector< PDE > pdes;                      // will store instantiated PDEs

  for (const auto& d : g_inputdeck.get< tag::selected, tag::pde >()) {
    if (d == ctr::PDEType::TRANSPORT)
      pdes.push_back( createPDE< tag::transport >( d, cnt ) );
    else if (d == ctr::PDEType::COMPFLOW)
      pdes.push_back( createPDE< tag::compflow >( d, cnt ) );
    else Throw( "Can't find selected PDE" );
  }

  return pdes;
}

std::vector< std::vector< std::pair< std::string, std::string > > >
PDEStack::info() const
// *****************************************************************************
//  Return information on all selected partial differential equations
//! \return A vector of vector of pair of strings, containing the configuration
//!   for each selected partial differential equation
//! \author J. Bakosi
// *****************************************************************************
{
  std::map< ctr::PDEType, ncomp_t > cnt; // count PDEs per type
  // will store info on all differential equations selected
  std::vector< std::vector< std::pair< std::string, std::string > > > nfo;

  for (const auto& d : g_inputdeck.get< tag::selected, tag::pde >()) {
    if (d == ctr::PDEType::TRANSPORT)
      nfo.emplace_back( infoTransport( cnt ) );
    else if (d == ctr::PDEType::COMPFLOW)
      nfo.emplace_back( infoCompFlow( cnt ) );
    else Throw( "Can't find selected PDE" );
  }

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
PDEStack::infoTransport( std::map< ctr::PDEType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the transport PDE
//! \param[inout] cnt std::map of counters for all partial differential equation
//!   types
//! \return vector of string pairs describing the PDE configuration
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::PDEType::TRANSPORT ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::TRANSPORT ), "" );
  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, tag::transport, tag::problem >()[c] ) );
  nfo.emplace_back( "start offset in unknowns array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::transport >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::transport >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
  const auto& diff =
     g_inputdeck.get< tag::param, tag::transport, tag::diffusivity >();
  if (diff.size() > c)
    nfo.emplace_back( "coeff diffusivity [" + std::to_string( ncomp ) + "]",
                       parameters( diff[c] ) );
  const auto& u0 = g_inputdeck.get< tag::param, tag::transport, tag::u0 >();
  if (u0.size() > c)
    nfo.emplace_back( "coeff u0 [" + std::to_string( ncomp ) + "]",
                       parameters( u0[c] ) );
  const auto& lambda =
    g_inputdeck.get< tag::param, tag::transport, tag::lambda >();
  if (lambda.size() > c)
    nfo.emplace_back( "coeff lambda [" + std::to_string( ncomp ) + "]",
      parameters( lambda[c] ) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
PDEStack::infoCompFlow( std::map< ctr::PDEType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the compressible flow system of PDEs
//! \param[inout] cnt std::map of counters for all partial differential equation
//!   types
//! \return vector of string pairs describing the PDE configuration
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::PDEType::COMPFLOW ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::COMPFLOW ), "" );
  nfo.emplace_back( "physics", ctr::Physics().name(
    g_inputdeck.get< tag::param, tag::compflow, tag::physics >()[c] ) );
  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, tag::compflow, tag::problem >()[c] ) );
  nfo.emplace_back( "start offset in unknowns array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::compflow >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::compflow >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
  nfo.emplace_back( "material id", parameters(
    g_inputdeck.get< tag::param, tag::compflow, tag::id >() ) );
  nfo.emplace_back( "ratio of specific heats", parameters(
    g_inputdeck.get< tag::param, tag::compflow, tag::gamma >() ) );

  const auto& mu = g_inputdeck.get< tag::param, tag::compflow, tag::mu >();
  if (!mu.empty())
    nfo.emplace_back( "dynamic viscosity", parameters( mu ) );

  const auto& cv = g_inputdeck.get< tag::param, tag::compflow, tag::cv >();
  if (!cv.empty())
    nfo.emplace_back( "specific heat at const. volume", parameters( cv ) );

  const auto& k = g_inputdeck.get< tag::param, tag::compflow, tag::k >();
  if (!k.empty())
    nfo.emplace_back( "heat conductivity", parameters( k ) );

  const auto& npar = g_inputdeck.get< tag::param, tag::compflow, tag::npar >();
  if (!npar.empty())
    nfo.emplace_back( "number of tracker particles", parameters( npar ) );

  return nfo;
}
