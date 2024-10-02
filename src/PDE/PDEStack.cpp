// *****************************************************************************
/*!
  \file      src/PDE/PDEStack.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Stack of partial differential equations
  \details   This file defines class PDEStack, which implements various
    functionality related to registering and instantiating partial differential
    equation types. Registration and instantiation use a partial differential
    equation factory, which is a std::map (an associative container),
    associating unique partial differential equation keys to their constructor
    calls. For more details, see the in-code documentation of the constructor.
*/
// *****************************************************************************

#include "PDEStack.hpp"

#include "ConfigureTransport.hpp"
#include "ConfigureCompFlow.hpp"
#include "ConfigureMultiMat.hpp"
#include "ConfigureMultiSpecies.hpp"

using inciter::PDEStack;

PDEStack::PDEStack() : m_cgfactory(), m_dgfactory(), m_fvfactory(),
                       m_cgEqTypes(), m_dgEqTypes(), m_fvEqTypes()
// *****************************************************************************
//  Constructor: register all partial differential equations into factory
//! \details This constructor consists of several blocks, each registering a
//!   potentially large number of entries in a partial differential equation
//!   factory, a standard associative container. At this time, each type of
//!   partial differential equation can be configured to use a unique _physics
//!   policy_ and a unique _problem policy_. (More types of policies might be
//!   introduced in the future.) Policy classes are template arguments to the
//!   partial differential equation classes and influence their behavior in a
//!   different way, abstracting away certain functions, e.g., how to set
//!   problem-specific initial and/or boundary conditions and how to update
//!   their coefficients during time integration. For more information on
//!   policy-based design, see http://en.wikipedia.org/wiki/Policy-based_design.
//!   This abstraction allows [separation of concerns]
//!   (http://en.wikipedia.org/wiki/Separation_of_concerns).
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
//!   _models-a_ relationship and its implementation in, e.g., PDE/CGPDE.h.
//!
//!   The design discussed above allows the registration, instantiation, and
//!   use of the partial differential equations to be generic, which eliminates
//!   a lot of boiler-plate code and makes client-code uniform.
//!
//!   _Details of registration using brigand::for_each and
//!   tk::cartesian_product:_
//!
//!   The template argument to brigand::for_each, as used below, requires a
//!   list of list of types. We use brigand::list of brigand::list of types,
//!   listing all possible policies, where the inner list must have exactly two
//!   types, as the list of lists is constructed from two lists using the
//!   cartesian product, and the length of the outer list (the list of lists) is
//!   arbitrary. The constructor argument to brigand::for_each is a functor that
//!   is to be applied to all members of the outer list. tk::cartesian_product
//!   will create all possible combinations of these types and call the functor
//!   with each type of the created sequence as a template parameter. The
//!   functor here inherits from registerPDE, which, i.e., its constructor call,
//!   needs a single template argument, a class templated on policy classes.
//!   This is the partial differential equation class to be configured by
//!   selecting policies and to be registered. The arguments to registerPDE's
//!   constructor are the factory, the enum denoting the differential equation
//!   type, and a reference to a variable of type std::set< ctr::PDEType >,
//!   which is only used internally to PDEStack for counting up the number of
//!   unique differential equation types registered, used for diagnostics
//!   purposes.
// *****************************************************************************
{
  registerTransport( m_cgfactory, m_dgfactory, m_cgEqTypes, m_dgEqTypes );
  registerCompFlow( m_cgfactory, m_dgfactory, m_cgEqTypes, m_dgEqTypes );
  registerMultiMat( m_dgfactory, m_fvfactory, m_dgEqTypes, m_fvEqTypes );
}

std::vector< inciter::CGPDE >
PDEStack::selectedCG() const
// *****************************************************************************
//  Instantiate all selected PDEs using continuous Galerkin discretization
//! \return std::vector of instantiated partial differential equation objects
// *****************************************************************************
{
  std::map< ctr::PDEType, ncomp_t > cnt;    // count PDEs per type
  std::vector< CGPDE > pdes;                // will store instantiated PDEs

  const auto sch = g_inputdeck.get< tag::scheme >();
  if (sch == ctr::SchemeType::ALECG || sch == ctr::SchemeType::OversetFE) {

    const auto& d = g_inputdeck.get< tag::pde >();
    for (std::size_t i=0; i<g_inputdeck.get< tag::mesh >().size(); ++i) {
      if (d == ctr::PDEType::TRANSPORT)
        pdes.push_back( createCG< tag::transport >( d, cnt ) );
      else if (d == ctr::PDEType::COMPFLOW)
        pdes.push_back( createCG< tag::compflow >( d, cnt ) );
      else Throw( "Can't find selected CGPDE" );
    }

  }

  return pdes;
}

std::vector< inciter::DGPDE >
PDEStack::selectedDG() const
// *****************************************************************************
//  Instantiate all selected PDEs using discontinuous Galerkin discretization
//! \return std::vector of instantiated partial differential equation objects
// *****************************************************************************
{
  std::map< ctr::PDEType, ncomp_t > cnt;    // count PDEs per type
  std::vector< DGPDE > pdes;                // will store instantiated PDEs

  auto sch = g_inputdeck.get< tag::scheme >();
  if (sch == ctr::SchemeType::DG ||
      sch == ctr::SchemeType::P0P1 || sch == ctr::SchemeType::DGP1 ||
      sch == ctr::SchemeType::DGP2 || sch == ctr::SchemeType::PDG ||
      sch == ctr::SchemeType::FV) {

    const auto& d = g_inputdeck.get< tag::pde >();
    for (std::size_t i=0; i<g_inputdeck.get< tag::mesh >().size(); ++i) {
      if (d == ctr::PDEType::TRANSPORT)
        pdes.push_back( createDG< tag::transport >( d, cnt ) );
      else if (d == ctr::PDEType::COMPFLOW)
        pdes.push_back( createDG< tag::compflow >( d, cnt ) );
      else if (d == ctr::PDEType::MULTIMAT)
        pdes.push_back( createDG< tag::multimat >( d, cnt ) );
      else if (d == ctr::PDEType::MULTISPECIES)
        pdes.push_back( createDG< tag::multispecies >( d, cnt ) );
      else Throw( "Can't find selected DGPDE" );
    }

  }

  return pdes;
}

std::vector< inciter::FVPDE >
PDEStack::selectedFV() const
// *****************************************************************************
//  Instantiate all selected PDEs using finite volume discretization
//! \return std::vector of instantiated partial differential equation objects
// *****************************************************************************
{
  std::map< ctr::PDEType, ncomp_t > cnt;    // count PDEs per type
  std::vector< FVPDE > pdes;                // will store instantiated PDEs

  auto sch = g_inputdeck.get< tag::scheme >();
  if (sch == ctr::SchemeType::FV) {

    const auto& d = g_inputdeck.get< tag::pde >();
    for (std::size_t i=0; i<g_inputdeck.get< tag::mesh >().size(); ++i) {
      if (d == ctr::PDEType::MULTIMAT)
        pdes.push_back( createFV< tag::multimat >( d, cnt ) );
      else Throw( "Can't find selected FVPDE" );
    }

  }

  return pdes;
}

std::vector< std::vector< std::pair< std::string, std::string > > >
PDEStack::info() const
// *****************************************************************************
//  Return information on all selected partial differential equations
//! \return A vector of vector of pair of strings, containing the configuration
//!   for each selected partial differential equation
// *****************************************************************************
{
  std::map< ctr::PDEType, ncomp_t > cnt; // count PDEs per type
  // will store info on all differential equations selected
  std::vector< std::vector< std::pair< std::string, std::string > > > nfo;

  const auto& d = g_inputdeck.get< tag::pde >();
  if (d == ctr::PDEType::TRANSPORT)
    nfo.emplace_back( infoTransport( cnt ) );
  else if (d == ctr::PDEType::COMPFLOW)
    nfo.emplace_back( infoCompFlow( cnt ) );
  else if (d == ctr::PDEType::MULTIMAT)
    nfo.emplace_back( infoMultiMat( cnt ) );
  else if (d == ctr::PDEType::MULTISPECIES)
    nfo.emplace_back( infoMultiSpecies( cnt ) );
  else Throw( "Can't find selected PDE" );

  return nfo;
}
