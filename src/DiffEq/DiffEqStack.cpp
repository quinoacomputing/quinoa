// *****************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Stack of differential equations
  \details   This file defines class DiffEqStack, which implements various
    functionality related to registering and instantiating differential equation
    types. Registration and instantiation use a differential equation factory,
    which is a std::map (an associative container), associating unique
    differential equation keys to their constructor calls. For more details, see
    the in-code documentation of the constructor.
*/
// *****************************************************************************

#include "DiffEqStack.hpp"
#include "Tags.hpp"

#include "ConfigureDirichlet.hpp"
#include "ConfigureMixDirichlet.hpp"
#include "ConfigureGeneralizedDirichlet.hpp"
#include "ConfigureWrightFisher.hpp"
#include "ConfigureOrnsteinUhlenbeck.hpp"
#include "ConfigureDiagOrnsteinUhlenbeck.hpp"
#include "ConfigureBeta.hpp"
#include "ConfigureNumberFractionBeta.hpp"
#include "ConfigureMassFractionBeta.hpp"
#include "ConfigureMixNumberFractionBeta.hpp"
#include "ConfigureMixMassFractionBeta.hpp"
#include "ConfigureGamma.hpp"
#include "ConfigureSkewNormal.hpp"
#include "ConfigureVelocity.hpp"
#include "ConfigurePosition.hpp"
#include "ConfigureDissipation.hpp"

using walker::DiffEqStack;

DiffEqStack::DiffEqStack() : m_factory(), m_eqTypes()
// *****************************************************************************
//  Constructor: register all differential equations into factory
//! \details This constructor consists of several blocks, each registering a
//!   potentially large number of entries in the differential equation factory,
//!   m_factory, which is of type walker::DiffEqFactory, a std::map. At this
//!   time, each type of differential equation can be configured to use a unique
//!   _initialization policy_ and a unique _coefficients policy_. (More types
//!   of policies will most likely come in the future.) The policy classes are
//!   template arguments to the differential equation classes and influence
//!   their behavior in a different way, abstracting away certain functions,
//!   e.g., how to set initial conditions and how to update their coefficients
//!   during time integration. For more information on policy-based design, see
//!   http://en.wikipedia.org/wiki/Policy-based_design. This abstraction allows
//!   [separation of concerns]
//!   (http://en.wikipedia.org/wiki/Separation_of_concerns).
//!
//!   Since the functionality of the policies are orthogonal to each other,
//!   i.e., they do not depend on each other or their host (the differential
//!   equation class), a Cartesian product of combinations are possible,
//!   depending on which policies are selected. _This constructor registers all
//!   possible combinations of policies for all available differential
//!   equations._ By _register_, we mean, an entry is recorded in an associative
//!   container, a std::map, that associates a lightweight key of type
//!   walker::ctr::DiffEqKey, consisting of only an enum for each policy type,
//!   to an std::function object that holds the constructor bound to its
//!   arguments corresponding to a particular differential equation + policies
//!   combination. Note that registering these entries in the map does not
//!   invoke the constructors. The mapped value simply stores how the
//!   constructors should be invoked at a later time. At some point later,
//!   based on user input, we then instantiate only the differential equations
//!   (and only those configurations) that are requested by the user.
//!
//!   Since all differential equation types (registered in the factory)
//!   "inherit" from a common "base", client-code is unform and generic, and
//!   thus immune to changes in the inner workings of the particular
//!   differential equations as long as they fullfill certain concepts, i.e.,
//!   implement certain member functinos, enforced by the _common base_, DiffEq.
//!   The words "inherit and "base" are quoted here, because the common base
//!   does not use inheritance in the normal OOP sense and does not use
//!   reference semantics, i.e., pointers, visible to client-code either. The
//!   relationship is more of a _models a_-type, which simplifies client-code
//!   and allows for the benfits of runtime inheritance with value-semantics
//!   which is less error prone and easier to read. See more about the
//!   _models-a_ relationship and its implementation in DiffEq/DiffEq.h.
//!
//!   The design discussed above allows the registration, instantiation, and
//!   use of the differential equations to be generic, which eliminates a lot of
//!   boiler-plate code and makes client-code uniform.
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
//!   functor here inherits from registerDiffEq, which, i.e., its constructor
//!   call, needs a single template argument, a class templated on policy
//!   classes. This is the differential equation class to be configured by
//!   selecting policies and to be registered. The arguments to
//!   registerDiffEq's constructor are the factory, the enum denoting the
//!   differential equation type, and a reference to a variable of type
//!   std::set< ctr::DiffEqType >, which is only used internally to DiffEqStack
//!   for counting up the number of unique differential equation types
//!   registered, used for diagnostics purposes.
// *****************************************************************************
{
  registerDirichlet( m_factory, m_eqTypes );
  registerMixDirichlet( m_factory, m_eqTypes );
  registerGenDir( m_factory, m_eqTypes );
  registerWrightFisher( m_factory, m_eqTypes );
  registerOrnsteinUhlenbeck( m_factory, m_eqTypes );
  registerDiagOrnsteinUhlenbeck( m_factory, m_eqTypes );
  registerBeta( m_factory, m_eqTypes );
  registerNumberFractionBeta( m_factory, m_eqTypes );
  registerMassFractionBeta( m_factory, m_eqTypes );
  registerMixNumberFractionBeta( m_factory, m_eqTypes );
  registerMixMassFractionBeta( m_factory, m_eqTypes );
  registerGamma( m_factory, m_eqTypes );
  registerSkewNormal( m_factory, m_eqTypes );
  registerVelocity( m_factory, m_eqTypes );
  registerPosition( m_factory, m_eqTypes );
  registerDissipation( m_factory, m_eqTypes );
}

std::vector< walker::DiffEq >
DiffEqStack::selected() const
// *****************************************************************************
//  Instantiate all selected differential equations
//! \return std::vector of instantiated differential equation objects
// *****************************************************************************
{
  std::map< ctr::DiffEqType, ncomp_t > cnt; // count DiffEqs per type
  std::vector< DiffEq > diffeqs;            // will store instantiated DiffEqs

  for (const auto& d : g_inputdeck.get< tag::selected, tag::diffeq >()) {
    if (d == ctr::DiffEqType::DIRICHLET)
      diffeqs.push_back( createDiffEq< tag::dirichlet >( d, cnt ) );
    else if (d == ctr::DiffEqType::MIXDIRICHLET)
      diffeqs.push_back( createDiffEq< tag::mixdirichlet >( d, cnt ) );
    else if (d == ctr::DiffEqType::GENDIR)
      diffeqs.push_back( createDiffEq< tag::gendir >( d, cnt ) );
    else if (d == ctr::DiffEqType::WRIGHTFISHER)
      diffeqs.push_back( createDiffEq< tag::wrightfisher >( d, cnt ) );
    else if (d == ctr::DiffEqType::OU)
      diffeqs.push_back( createDiffEq< tag::ou >( d, cnt ) );
    else if (d == ctr::DiffEqType::DIAG_OU)
      diffeqs.push_back( createDiffEq< tag::diagou >( d, cnt ) );
    else if (d == ctr::DiffEqType::BETA)
      diffeqs.push_back( createDiffEq< tag::beta >( d, cnt ) );
    else if (d == ctr::DiffEqType::NUMFRACBETA)
      diffeqs.push_back( createDiffEq< tag::numfracbeta >( d, cnt ) );
    else if (d == ctr::DiffEqType::MASSFRACBETA)
      diffeqs.push_back( createDiffEq< tag::massfracbeta >( d, cnt ) );
    else if (d == ctr::DiffEqType::MIXNUMFRACBETA)
      diffeqs.push_back( createDiffEq< tag::mixnumfracbeta >( d, cnt ) );
    else if (d == ctr::DiffEqType::MIXMASSFRACBETA)
      diffeqs.push_back( createDiffEq< tag::mixmassfracbeta >( d, cnt ) );
    else if (d == ctr::DiffEqType::SKEWNORMAL)
      diffeqs.push_back( createDiffEq< tag::skewnormal >( d, cnt ) );
    else if (d == ctr::DiffEqType::GAMMA)
      diffeqs.push_back( createDiffEq< tag::gamma >( d, cnt ) );
    else if (d == ctr::DiffEqType::VELOCITY)
      diffeqs.push_back( createDiffEq< tag::velocity >( d, cnt ) );
    else if (d == ctr::DiffEqType::POSITION)
      diffeqs.push_back( createDiffEq< tag::position >( d, cnt ) );
    else if (d == ctr::DiffEqType::DISSIPATION)
      diffeqs.push_back( createDiffEq< tag::dissipation >( d, cnt ) );
    else Throw( "Can't find selected DiffEq" );
  }

  return diffeqs;
}

std::pair< std::vector< std::string >, std::vector< tk::Table<1> > >
DiffEqStack::tables() const
// *****************************************************************************
//  Instantiate tables from which extra statistics data to be output sampled for
//  all selected differential equations
//! \return Vector of names and tables to sample from during time stepping
// *****************************************************************************
{
  std::map< ctr::DiffEqType, ncomp_t > cnt;     // count DiffEqs per type
  std::vector< std::string > nam;               // names of instantiated tables
  std::vector< tk::Table<1> > tab;              // instantiated tables

  for (const auto& d : g_inputdeck.get< tag::selected, tag::diffeq >()) {
    std::pair< std::vector< std::string >, std::vector< tk::Table<1> > > t;

    if (d == ctr::DiffEqType::MIXMASSFRACBETA)
      t = createTables< tag::mixmassfracbeta >( d, cnt );
    else if (d == ctr::DiffEqType::VELOCITY)
      t = createTables< tag::velocity >( d, cnt );

    nam.insert( end(nam), begin(t.first), end(t.first) );
    tab.insert( end(tab), begin(t.second), end(t.second) );
  }

  return { nam, tab };
}

std::vector< std::vector< std::pair< std::string, std::string > > >
DiffEqStack::info() const
// *****************************************************************************
//  Return information on all selected differential equations
//! \return A vector of vector of pair of strings, containing the configuration
//!   for each selected differential equation
// *****************************************************************************
{
  std::map< ctr::DiffEqType, ncomp_t > cnt; // count DiffEqs per type
  // will store info on all differential equations selected
  std::vector< std::vector< std::pair< std::string, std::string > > > nfo;

  for (const auto& d : g_inputdeck.get< tag::selected, tag::diffeq >()) {
    if (d == ctr::DiffEqType::DIRICHLET)
      nfo.emplace_back( infoDirichlet( cnt ) );
    else if (d == ctr::DiffEqType::MIXDIRICHLET)
      nfo.emplace_back( infoMixDirichlet( cnt ) );
    else if (d == ctr::DiffEqType::GENDIR)
      nfo.emplace_back( infoGenDir( cnt ) );
    else if (d == ctr::DiffEqType::WRIGHTFISHER)
      nfo.emplace_back( infoWrightFisher( cnt ) );
    else if (d == ctr::DiffEqType::OU)
      nfo.emplace_back( infoOrnsteinUhlenbeck( cnt ) );
    else if (d == ctr::DiffEqType::DIAG_OU)
      nfo.emplace_back( infoDiagOrnsteinUhlenbeck( cnt ) );
    else if (d == ctr::DiffEqType::BETA)
      nfo.emplace_back( infoBeta( cnt ) );
    else if (d == ctr::DiffEqType::NUMFRACBETA)
      nfo.emplace_back( infoNumberFractionBeta( cnt ) );
    else if (d == ctr::DiffEqType::MASSFRACBETA)
      nfo.emplace_back( infoMassFractionBeta( cnt ) );
    else if (d == ctr::DiffEqType::MIXNUMFRACBETA)
      nfo.emplace_back( infoMixNumberFractionBeta( cnt ) );
    else if (d == ctr::DiffEqType::MIXMASSFRACBETA)
      nfo.emplace_back( infoMixMassFractionBeta( cnt ) );
    else if (d == ctr::DiffEqType::SKEWNORMAL)
      nfo.emplace_back( infoSkewNormal( cnt ) );
    else if (d == ctr::DiffEqType::GAMMA)
      nfo.emplace_back( infoGamma( cnt ) );
    else if (d == ctr::DiffEqType::VELOCITY)
      nfo.emplace_back( infoVelocity( cnt ) );
    else if (d == ctr::DiffEqType::POSITION)
      nfo.emplace_back( infoPosition( cnt ) );
    else if (d == ctr::DiffEqType::DISSIPATION)
      nfo.emplace_back( infoDissipation( cnt ) );
    else Throw( "Can't find selected DiffEq" );
  }

  return nfo;
}
