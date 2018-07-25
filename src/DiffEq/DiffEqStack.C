// *****************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Stack of differential equations
  \details   This file defines class DiffEqStack, which implements various
    functionality related to registering and instantiating differential equation
    types. Registration and instantiation use a differential equation factory,
    which is a std::map (an associative container), associating unique
    differential equation keys to their constructor calls. For more details, see
    the in-code documentation of the constructor.
*/
// *****************************************************************************

#include <brigand/algorithms/for_each.hpp>

#include "DiffEqStack.h"
#include "Tags.h"
#include "SystemComponents.h"
#include "Options/RNG.h"
#include "Walker/Options/CoeffPolicy.h"
#include "Walker/Options/InitPolicy.h"
#include "Walker/Options/HydroTimeScales.h"
#include "Walker/Options/HydroProductions.h"
#include "CartesianProduct.h"

#include "Beta.h"
#include "DiagOrnsteinUhlenbeck.h"
#include "Dirichlet.h"
#include "Gamma.h"
#include "GeneralizedDirichlet.h"
#include "MassFractionBeta.h"
#include "MixMassFractionBeta.h"
#include "MixNumberFractionBeta.h"
#include "NumberFractionBeta.h"
#include "OrnsteinUhlenbeck.h"
#include "SkewNormal.h"
#include "WrightFisher.h"
#include "Velocity.h"
#include "Position.h"
#include "Dissipation.h"

#include "BetaCoeffPolicy.h"
#include "DiagOrnsteinUhlenbeckCoeffPolicy.h"
#include "DirichletCoeffPolicy.h"
#include "GammaCoeffPolicy.h"
#include "GeneralizedDirichletCoeffPolicy.h"
#include "MassFractionBetaCoeffPolicy.h"
#include "MixMassFractionBetaCoeffPolicy.h"
#include "MixNumberFractionBetaCoeffPolicy.h"
#include "NumberFractionBetaCoeffPolicy.h"
#include "OrnsteinUhlenbeckCoeffPolicy.h"
#include "SkewNormalCoeffPolicy.h"
#include "WrightFisherCoeffPolicy.h"
#include "VelocityCoeffPolicy.h"
#include "PositionCoeffPolicy.h"
#include "DissipationCoeffPolicy.h"

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
//!   [separation of concerns](http://en.wikipedia.org/wiki/Separation_of_concerns).
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
  namespace mpl = boost::mpl;

  // Dirichlet SDE
  // Construct vector of vectors for all possible policies for SDE
  using DirPolicies =
    tk::cartesian_product< InitPolicies, DirichletCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< DirPolicies >(
    registerDiffEq< Dirichlet >
                  ( m_factory, ctr::DiffEqType::DIRICHLET, m_eqTypes ) );

  // Lochner's generalized Dirichlet SDE
  // Construct vector of vectors for all possible policies for SDE
  using GenDirPolicies =
    tk::cartesian_product< InitPolicies, GeneralizedDirichletCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< GenDirPolicies >(
    registerDiffEq< GeneralizedDirichlet >
                  ( m_factory, ctr::DiffEqType::GENDIR, m_eqTypes ) );

  // Wright-Fisher SDE
  // Construct vector of vectors for all possible policies for SDE
  using WrightFisherPolicies =
    tk::cartesian_product< InitPolicies, WrightFisherCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< WrightFisherPolicies >(
    registerDiffEq< WrightFisher >
                  ( m_factory, ctr::DiffEqType::WRIGHTFISHER, m_eqTypes ) );

  // Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using OrnsteinUhlenbeckPolicies =
    tk::cartesian_product< InitPolicies, OrnsteinUhlenbeckCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< OrnsteinUhlenbeckPolicies >(
    registerDiffEq< OrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::OU, m_eqTypes ) );

  // Diagonal Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using DiagOrnsteinUhlenbeckPolicies =
    tk::cartesian_product< InitPolicies, DiagOrnsteinUhlenbeckCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< DiagOrnsteinUhlenbeckPolicies >(
    registerDiffEq< DiagOrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::DIAG_OU, m_eqTypes ) );

  // beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using BetaPolicies =
    tk::cartesian_product< InitPolicies, BetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< BetaPolicies >(
    registerDiffEq< Beta >( m_factory, ctr::DiffEqType::BETA, m_eqTypes ) );

  // Number-fraction beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using NumberFractionBetaPolicies =
    tk::cartesian_product< InitPolicies, NumberFractionBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< NumberFractionBetaPolicies >(
    registerDiffEq< NumberFractionBeta >
                  ( m_factory, ctr::DiffEqType::NUMFRACBETA, m_eqTypes ) );

  // Mass-fraction beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using MassFractionBetaPolicies =
    tk::cartesian_product< InitPolicies, MassFractionBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< MassFractionBetaPolicies >(
    registerDiffEq< MassFractionBeta >
                  ( m_factory, ctr::DiffEqType::MASSFRACBETA, m_eqTypes ) );

  // Mix number-fraction beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using MixNumFracBetaPolicies =
    tk::cartesian_product< InitPolicies, MixNumFracBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< MixNumFracBetaPolicies >(
    registerDiffEq< MixNumberFractionBeta >
                  ( m_factory, ctr::DiffEqType::MIXNUMFRACBETA, m_eqTypes ) );

  // Mix mass-fraction beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using MixMassFracBetaPolicies =
    tk::cartesian_product< InitPolicies, MixMassFracBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< MixMassFracBetaPolicies >(
    registerDiffEq< MixMassFractionBeta >
                  ( m_factory, ctr::DiffEqType::MIXMASSFRACBETA, m_eqTypes ) );

  // Skew-normal SDE
  // Construct vector of vectors for all possible policies for SDE
  using SkewNormalPolicies =
    tk::cartesian_product< InitPolicies, SkewNormalCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< SkewNormalPolicies >(
    registerDiffEq< SkewNormal >
                  ( m_factory, ctr::DiffEqType::SKEWNORMAL, m_eqTypes ) );

  // Gamma SDE
  // Construct vector of vectors for all possible policies for SDE
  using GammaPolicies =
    tk::cartesian_product< InitPolicies, GammaCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< GammaPolicies >(
    registerDiffEq< Gamma >( m_factory, ctr::DiffEqType::GAMMA, m_eqTypes ) );

  // Velocity SDE
  // Construct vector of vectors for all possible policies for SDE
  using VelocityPolicies =
    tk::cartesian_product< InitPolicies, VelocityCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< VelocityPolicies >(
    registerDiffEq< Velocity >
                  ( m_factory, ctr::DiffEqType::VELOCITY, m_eqTypes ) );

  // Position equation
  // Construct vector of vectors for all possible policies for equation
  using PositionPolicies =
    tk::cartesian_product< InitPolicies, PositionCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< PositionPolicies >(
    registerDiffEq< Position >
                  ( m_factory, ctr::DiffEqType::POSITION, m_eqTypes ) );

  // Dissipation equation
  // Construct vector of vectors for all possible policies for equation
  using DissipationPolicies =
    tk::cartesian_product< InitPolicies, DissipationCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< DissipationPolicies >(
    registerDiffEq< Dissipation >
                  ( m_factory, ctr::DiffEqType::DISSIPATION, m_eqTypes ) );
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

std::pair< std::vector< std::string >, std::vector< tk::Table > >
DiffEqStack::tables() const
// *****************************************************************************
//  Instantiate tables from which extra statistics data to be output sampled for
//  all selected differential equations
//! \return Vector of names and tables to sample from during time stepping
// *****************************************************************************
{
  std::map< ctr::DiffEqType, ncomp_t > cnt;     // count DiffEqs per type
  std::vector< std::string > nam;               // names of instantiated tables
  std::vector< tk::Table > tab;                 // instantiated tables

  for (const auto& d : g_inputdeck.get< tag::selected, tag::diffeq >()) {
    std::pair< std::vector< std::string >, std::vector< tk::Table > > t;

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
    else if (d == ctr::DiffEqType::GENDIR)
      nfo.emplace_back( infoGenDir( cnt ) );
    else if (d == ctr::DiffEqType::WRIGHTFISHER)
      nfo.emplace_back( infoWrightFisher( cnt ) );
    else if (d == ctr::DiffEqType::OU)
      nfo.emplace_back( infoOU( cnt ) );
    else if (d == ctr::DiffEqType::DIAG_OU)
      nfo.emplace_back( infoDiagOU( cnt ) );
    else if (d == ctr::DiffEqType::BETA)
      nfo.emplace_back( infoBeta( cnt ) );
    else if (d == ctr::DiffEqType::NUMFRACBETA)
      nfo.emplace_back( infoNumberFractionBeta( cnt ) );
    else if (d == ctr::DiffEqType::MASSFRACBETA)
      nfo.emplace_back( infoMassFractionBeta( cnt ) );
    else if (d == ctr::DiffEqType::MIXNUMFRACBETA)
      nfo.emplace_back( infoMixNumFracBeta( cnt ) );
    else if (d == ctr::DiffEqType::MIXMASSFRACBETA)
      nfo.emplace_back( infoMixMassFracBeta( cnt ) );
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

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoDirichlet( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the Dirichlet SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DIRICHLET ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DIRICHLET ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::dirichlet >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::dirichlet >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::dirichlet, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::dirichlet, tag::b >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::dirichlet, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::dirichlet, tag::kappa >().at(c) )
  );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoGenDir( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on Lochner's generalized Dirichlet SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::GENDIR ];  // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::GENDIR ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::gendir >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::gendir >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::gendir, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::gendir, tag::b >().at(c) ) );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::gendir, tag::S >().at(c) ) );
  nfo.emplace_back(
    "coeff kappa [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::gendir, tag::kappa >().at(c) )
  );
  nfo.emplace_back(
    "coeff c [" + std::to_string( ncomp*(ncomp-1)/2 ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::gendir, tag::c >().at(c) ) );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::gendir, tag::spike >().at(c) );
  betapdfs( nfo,
            g_inputdeck.get< tag::param, tag::gendir, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoWrightFisher( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the Wright-Fisher SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::WRIGHTFISHER ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::WRIGHTFISHER ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::wrightfisher >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::wrightfisher >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff omega [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::wrightfisher, tag::omega >().at(c) ) );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::wrightfisher, tag::spike >().at(c)
  );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoOU( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the Ornstein-Uhlenbeck SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::OU ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::OU ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::ou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::ou >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::ou, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::ou, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp*(ncomp+1)/2 ) + ", upper tri]",
    parameters( g_inputdeck.get< tag::param, tag::ou, tag::sigmasq >().at(c) )
  );
  nfo.emplace_back( "coeff theta [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::ou, tag::theta >().at(c) ) );
  nfo.emplace_back( "coeff mu [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::ou, tag::mu >().at(c) ) );
  spikes( nfo, g_inputdeck.get< tag::param, tag::ou, tag::spike >().at(c) );
  betapdfs( nfo,
            g_inputdeck.get< tag::param, tag::ou, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoDiagOU( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the diagonal Ornstein-Uhlenbeck SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DIAG_OU ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DIAG_OU ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::diagou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::diagou >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::diagou, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::diagou, tag::sigmasq >().at(c) ) );
  nfo.emplace_back( "coeff theta [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::diagou, tag::theta >().at(c) )
  );
  nfo.emplace_back( "coeff mu [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::diagou, tag::mu >().at(c) ) );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::diagou, tag::spike >().at(c) );
  betapdfs( nfo,
            g_inputdeck.get< tag::param, tag::diagou, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::BETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::BETA ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::beta >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::beta >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::beta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::beta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::beta, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::beta, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::beta, tag::b >().at(c) ) );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::beta, tag::S >().at(c) ) );
  nfo.emplace_back(
    "coeff kappa [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::beta, tag::kappa >().at(c) )
  );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::beta, tag::spike >().at(c) );
  betapdfs( nfo,
            g_inputdeck.get< tag::param, tag::beta, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoNumberFractionBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt )
const
// *****************************************************************************
//  Return information on the number-fraction beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::NUMFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::NUMFRACBETA ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::numfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::numfracbeta >()[c] / 3;
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::coeffpolicy >()[c] ) );

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::numfracbeta, tag::b >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::numfracbeta, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::numfracbeta, tag::kappa >().at(c) ) );
  nfo.emplace_back(
    "coeff rho2 [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::numfracbeta, tag::rho2 >().at(c) ) );
  nfo.emplace_back(
    "coeff rcomma [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::numfracbeta, tag::rcomma >().at(c) ) );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::numfracbeta, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoMassFractionBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt )
const
// *****************************************************************************
//  Return information on the mass-fraction beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::MASSFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::MASSFRACBETA ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::massfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::massfracbeta >()[c] / 3;
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::coeffpolicy >()[c] ) );

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b [" + std::to_string( ncomp ) + "]",
    parameters(g_inputdeck.get< tag::param, tag::massfracbeta, tag::b >().at(c))
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::massfracbeta, tag::S >().at(c) ) );
  nfo.emplace_back(
    "coeff kappa [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::massfracbeta, tag::kappa >().at(c) ) );
  nfo.emplace_back(
    "coeff rho2 [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::massfracbeta, tag::rho2 >().at(c) ) );
  nfo.emplace_back(
    "coeff r [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::massfracbeta, tag::r >().at(c) ) );
  spikes(
    nfo,
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoMixNumFracBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt )
const
// *****************************************************************************
//  Return information on the mix number-fraction beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::MIXNUMFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back(
    ctr::DiffEq().name( ctr::DiffEqType::MIXNUMFRACBETA ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::mixnumfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::mixnumfracbeta >()[c] / 3;
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::coeffpolicy >()[c] )
  );

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b' [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::bprime >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa' [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::kappaprime >().at(c)
    )
  );
  nfo.emplace_back(
    "coeff rho2 [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::rho2 >().at(c) ) );
  nfo.emplace_back(
    "coeff rcomma [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::rcomma >().at(c) )
  );
  spikes( nfo,
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoMixMassFracBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt )
const
// *****************************************************************************
//  Return information on the mix mass-fraction beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::MIXMASSFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back(
    ctr::DiffEq().name( ctr::DiffEqType::MIXMASSFRACBETA ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::mixmassfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::mixmassfracbeta >()[c];
  nfo.emplace_back( "number of components",
    std::to_string( ncomp ) + " (" + std::to_string(ncomp/4) + "*4) " );

  // Optional coupled
  const auto& coupled_velocity =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::velocity >();
  const auto& coupled_velocity_id =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::velocity_id >();
  Assert( coupled_velocity.size() == coupled_velocity_id.size(),
          "Size mismatch" );

  if (coupled_velocity.size() > c) {
    nfo.emplace_back( "coupled velocity depvar", std::string( 1,
      coupled_velocity[c] ) );
    nfo.emplace_back( "coupled velocity depvar offset", std::to_string(
      coupled_velocity_id[c] ) );
  }

  // Optional coupled
  const auto& coupled_dissipation =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::dissipation >();
  const auto& coupled_dissipation_id =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::dissipation_id >();
  Assert( coupled_dissipation.size() == coupled_dissipation_id.size(),
          "Size mismatch" );

  if (coupled_dissipation.size() > c) {
    nfo.emplace_back( "coupled dissipation depvar", std::string( 1,
      coupled_dissipation[c] ) );
    nfo.emplace_back( "coupled dissipation depvar ofs", std::to_string(
    coupled_dissipation_id[c] ) );
  }

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::initpolicy >()[c] )
  );
  auto cp =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::coeffpolicy >()[c];
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name( cp ) );
  if (cp == ctr::CoeffPolicyType::HYDROTIMESCALE) {
    nfo.emplace_back(
      "inverse hydro time scales [" + std::to_string( ncomp/4 ) + "]",
      options( ctr::HydroTimeScales(),
               g_inputdeck.get< tag::param,
                                tag::mixmassfracbeta,
                                tag::hydrotimescales >().at(c) ) );
    nfo.emplace_back(
      "production/dissipation [" + std::to_string( ncomp/4 ) + "]",
      options( ctr::HydroProductions(),
               g_inputdeck.get< tag::param,
                                tag::mixmassfracbeta,
                                tag::hydroproductions >().at(c) ) );
  }

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b' [" + std::to_string( ncomp/4 ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::bprime >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp/4 ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa' [" + std::to_string( ncomp/4 ) + "]",
    parameters( g_inputdeck.get< tag::param,
                                 tag::mixmassfracbeta,
                                 tag::kappaprime >().at(c) ) );
  nfo.emplace_back(
    "coeff rho2 [" + std::to_string( ncomp/4 ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::rho2 >().at(c) ) );
  nfo.emplace_back(
    "coeff r [" + std::to_string( ncomp/4 ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::r >().at(c) )
  );
  spikes( nfo,
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoSkewNormal( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the skew-normal SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::SKEWNORMAL ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::SKEWNORMAL ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::skewnormal >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::skewnormal >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::skewnormal, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff T [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::skewnormal, tag::timescale >().at(c) )
  );
  nfo.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::skewnormal, tag::sigmasq >().at(c) ) );
  nfo.emplace_back(
    "coeff lambda [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::skewnormal, tag::lambda >().at(c) ) );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::skewnormal, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::skewnormal, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoGamma( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the gamma SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::GAMMA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::GAMMA ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::gamma >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::gamma >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::gamma, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::gamma, tag::b >().at(c) ) );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::gamma, tag::S >().at(c) ) );
  nfo.emplace_back(
    "coeff kappa [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::gamma, tag::kappa >().at(c) )
  );
  spikes( nfo, g_inputdeck.get< tag::param, tag::gamma, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::gamma, tag::betapdf >().at(c) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoVelocity( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the velocity SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::VELOCITY ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::VELOCITY ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::velocity >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::velocity >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  // Required coupled
  nfo.emplace_back( "coupled position depvar", std::string( 1,
    g_inputdeck.get< tag::param, tag::velocity, tag::position >()[c] ) );
  nfo.emplace_back( "coupled position depvar offset", std::to_string(
    g_inputdeck.get< tag::param, tag::velocity, tag::position_id >()[c] ) );
  nfo.emplace_back( "coupled dissipation depvar", std::string( 1,
    g_inputdeck.get< tag::param, tag::velocity, tag::dissipation >()[c] ) );
  nfo.emplace_back( "coupled dissipation depv offs", std::to_string(
    g_inputdeck.get< tag::param, tag::velocity, tag::dissipation_id >()[c] ) );

  // Optional coupled
  const auto& coupled_mixmassfracbeta =
    g_inputdeck.get< tag::param, tag::velocity, tag::mixmassfracbeta >();
  const auto& coupled_mixmassfracbeta_id =
    g_inputdeck.get< tag::param, tag::velocity, tag::mixmassfracbeta_id >();
  Assert( coupled_mixmassfracbeta.size() == coupled_mixmassfracbeta_id.size(),
          "Size mismatch" );
  if (coupled_mixmassfracbeta.size() > c) {
    nfo.emplace_back( "coupled mixmassfracbeta depvar", std::string( 1,
      coupled_mixmassfracbeta[c] ) );
    nfo.emplace_back( "coupled mixmassfracbeta dv of", std::to_string(
      coupled_mixmassfracbeta_id[c] ) );
  }

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::velocity, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::velocity, tag::initpolicy >()[c] ) );
  auto cp = g_inputdeck.get< tag::param, tag::velocity, tag::coeffpolicy >()[c];
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name( cp ) );

  auto solve = g_inputdeck.get< tag::param, tag::velocity, tag::solve >()[c];
  auto depvar = ctr::Depvar();
  nfo.emplace_back( depvar.group(), depvar.name(solve) );

  auto variant =
    g_inputdeck.get< tag::param, tag::velocity, tag::variant >()[c];
  auto velocity = ctr::VelocityVariant();
  nfo.emplace_back( velocity.group(), velocity.name(variant) );

  if (cp == ctr::CoeffPolicyType::HYDROTIMESCALE) {
    nfo.emplace_back(
      "inverse hydro time scale",
      options( ctr::HydroTimeScales(),
               g_inputdeck.get< tag::param,
                                tag::velocity,
                                tag::hydrotimescales >().at(c) ) );
    nfo.emplace_back(
      "production/dissipation",
      options( ctr::HydroProductions(),
               g_inputdeck.get< tag::param,
                                tag::velocity,
                                tag::hydroproductions >().at(c) ) );
  }

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::velocity, tag::rng >()[c] ) );
  nfo.emplace_back( "coeff C0", std::to_string(
    g_inputdeck.get< tag::param, tag::velocity, tag::c0 >().at(c) ) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoPosition( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the position eq
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::POSITION ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::POSITION ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::position >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::position >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  // Required coupled
  nfo.emplace_back( "coupled velocity depvar", std::string( 1,
    g_inputdeck.get< tag::param, tag::position, tag::velocity >()[c] ) );
  nfo.emplace_back( "coupled velocity depvar offset", std::to_string(
    g_inputdeck.get< tag::param, tag::position, tag::velocity_id >()[c] ) );

  nfo.emplace_back( "kind", "deterministic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::position, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::position, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::position, tag::coeffpolicy >()[c] ) );
  auto solve = g_inputdeck.get< tag::param, tag::position, tag::solve >()[c];
  nfo.emplace_back( "solve for", ctr::Depvar().name( solve ) );

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::position, tag::rng >()[c] ) );

  return nfo;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoDissipation( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
// *****************************************************************************
//  Return information on the dissipation eq
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DISSIPATION ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DISSIPATION ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::dissipation >(c) ) );
  nfo.emplace_back( "number of components", std::to_string(
    g_inputdeck.get< tag::component >().get< tag::dissipation >()[c] ) );

  // Required coupled
  nfo.emplace_back( "coupled velocity depvar", std::string( 1,
    g_inputdeck.get< tag::param, tag::dissipation, tag::velocity >()[c] ) );
  nfo.emplace_back( "coupled velocity depvar offset", std::to_string(
    g_inputdeck.get< tag::param, tag::dissipation, tag::velocity_id >()[c] ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::dissipation, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::dissipation, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::dissipation, tag::coeffpolicy >()[c] ) );

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::dissipation, tag::rng >()[c] ) );
  nfo.emplace_back( "coeff C3", std::to_string(
    g_inputdeck.get< tag::param, tag::dissipation, tag::c3 >().at(c) ) );
  nfo.emplace_back( "coeff C4", std::to_string(
    g_inputdeck.get< tag::param, tag::dissipation, tag::c4 >().at(c) ) );
  nfo.emplace_back( "coeff COM1", std::to_string(
    g_inputdeck.get< tag::param, tag::dissipation, tag::com1 >().at(c) ) );
  nfo.emplace_back( "coeff COM2", std::to_string(
    g_inputdeck.get< tag::param, tag::dissipation, tag::com2 >().at(c) ) );

  return nfo;
}
