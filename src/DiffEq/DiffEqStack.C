// *****************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Stack of differential equations
  \details   This file defines class DiffEqStack, which implements various
    functionality related to registering and instantiating differential equation
    types. Registration and instantiation use a differential equation factory,
    which is a std::map (an associative container), associating unique
    differential equation keys to their constructor calls. For more details, see
    the in-code documentation of the constructor.
*/
// *****************************************************************************

#include "NoWarning/cartesian_product.h"

#include "DiffEqStack.h"
#include "Tags.h"
#include "SystemComponents.h"
#include "Options/RNG.h"
#include "Walker/Options/CoeffPolicy.h"
#include "Walker/Options/InitPolicy.h"
#include "Walker/Options/HydroTimeScales.h"
#include "Walker/Options/HydroProductions.h"

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
//!   _Details of registration using mpl::cartesian_product:_
//!
//!   The template argument to mpl::cartesian_product requires a sequence of
//!   sequences of types. We use vector of vectors of types, listing all
//!   possible policies. The constructor argument to mpl::cartesian_product is a
//!   functor that is to be applied to all combinations. mpl::cartesian_product
//!   will then create all possible combinations of these types and call the
//!   user-supplied functor with each type of the created sequence as a template
//!   parameter. The user-supplied functor here is registerDiffEq, which, i.e.,
//!   its constructor call, needs a single template argument, a class templated
//!   on policy classes. This is the differential equation class to be
//!   configured by selecting policies and to be registered. The arguments to
//!   registerDiffEq's constructor are the factory, the enum denoting the
//!   differential equation type, and a reference to a variable of type
//!   std::set< ctr::DiffEqType >, which is only used internally to DiffEqStack
//!   for counting up the number of unique differential equation types
//!   registered, used for diagnostics purposes.
//! \author J. Bakosi
// *****************************************************************************
{
  namespace mpl = boost::mpl;

  // Dirichlet SDE
  // Construct vector of vectors for all possible policies for SDE
  using DirPolicies = mpl::vector< InitPolicies, DirichletCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< DirPolicies >(
    registerDiffEq< Dirichlet >
                  ( m_factory, ctr::DiffEqType::DIRICHLET, m_eqTypes ) );

  // Lochner's generalized Dirichlet SDE
  // Construct vector of vectors for all possible policies for SDE
  using GenDirPolicies = mpl::vector< InitPolicies,
                                      GeneralizedDirichletCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< GenDirPolicies >(
    registerDiffEq< GeneralizedDirichlet >
                  ( m_factory, ctr::DiffEqType::GENDIR, m_eqTypes ) );

  // Wright-Fisher SDE
  // Construct vector of vectors for all possible policies for SDE
  using WrightFisherPolicies =
    mpl::vector< InitPolicies, WrightFisherCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< WrightFisherPolicies >(
    registerDiffEq< WrightFisher >
                  ( m_factory, ctr::DiffEqType::WRIGHTFISHER, m_eqTypes ) );

  // Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using OrnsteinUhlenbeckPolicies =
    mpl::vector< InitPolicies, OrnsteinUhlenbeckCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< OrnsteinUhlenbeckPolicies >(
    registerDiffEq< OrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::OU, m_eqTypes ) );

  // Diagonal Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using DiagOrnsteinUhlenbeckPolicies =
    mpl::vector< InitPolicies, DiagOrnsteinUhlenbeckCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< DiagOrnsteinUhlenbeckPolicies >(
    registerDiffEq< DiagOrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::DIAG_OU, m_eqTypes ) );

  // beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using BetaPolicies = mpl::vector< InitPolicies, BetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< BetaPolicies >(
    registerDiffEq< Beta >( m_factory, ctr::DiffEqType::BETA, m_eqTypes ) );

  // Number-fraction beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using NumberFractionBetaPolicies =
    mpl::vector< InitPolicies, NumberFractionBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< NumberFractionBetaPolicies >(
    registerDiffEq< NumberFractionBeta >
                  ( m_factory, ctr::DiffEqType::NUMFRACBETA, m_eqTypes ) );

  // Mass-fraction beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using MassFractionBetaPolicies =
    mpl::vector< InitPolicies, MassFractionBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< MassFractionBetaPolicies >(
    registerDiffEq< MassFractionBeta >
                  ( m_factory, ctr::DiffEqType::MASSFRACBETA, m_eqTypes ) );

  // Mix number-fraction beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using MixNumFracBetaPolicies =
    mpl::vector< InitPolicies, MixNumFracBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< MixNumFracBetaPolicies >(
    registerDiffEq< MixNumberFractionBeta >
                  ( m_factory, ctr::DiffEqType::MIXNUMFRACBETA, m_eqTypes ) );

  // Mix mass-fraction beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using MixMassFracBetaPolicies =
    mpl::vector< InitPolicies, MixMassFracBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< MixMassFracBetaPolicies >(
    registerDiffEq< MixMassFractionBeta >
                  ( m_factory, ctr::DiffEqType::MIXMASSFRACBETA, m_eqTypes ) );

  // Skew-normal SDE
  // Construct vector of vectors for all possible policies for SDE
  using SkewNormalPolicies =
    mpl::vector< InitPolicies, SkewNormalCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< SkewNormalPolicies >(
    registerDiffEq< SkewNormal >
                  ( m_factory, ctr::DiffEqType::SKEWNORMAL, m_eqTypes ) );

  // Gamma SDE
  // Construct vector of vectors for all possible policies for SDE
  using GammaPolicies = mpl::vector< InitPolicies, GammaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< GammaPolicies >(
    registerDiffEq< Gamma >( m_factory, ctr::DiffEqType::GAMMA, m_eqTypes ) );
}

std::vector< walker::DiffEq >
DiffEqStack::selected() const
// *****************************************************************************
//  Instantiate all selected differential equations
//! \return std::vector of instantiated differential equation objects
//! \author J. Bakosi
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
//! \author J. Bakosi
// *****************************************************************************
{
  std::map< ctr::DiffEqType, ncomp_t > cnt;     // count DiffEqs per type
  std::vector< std::string > nam;               // names of instantiated tables
  std::vector< tk::Table > tab;                 // instantiated tables

  for (const auto& d : g_inputdeck.get< tag::selected, tag::diffeq >()) {
    std::pair< std::vector< std::string >, std::vector< tk::Table > > t;

    if (d == ctr::DiffEqType::DIRICHLET)
      t = createTables< tag::dirichlet >( d, cnt );
    else if (d == ctr::DiffEqType::GENDIR)
      t = createTables< tag::gendir >( d, cnt );
    else if (d == ctr::DiffEqType::WRIGHTFISHER)
      t = createTables< tag::wrightfisher >( d, cnt );
    else if (d == ctr::DiffEqType::OU)
      t = createTables< tag::ou >( d, cnt );
    else if (d == ctr::DiffEqType::DIAG_OU)
      t = createTables< tag::diagou >( d, cnt );
    else if (d == ctr::DiffEqType::BETA)
      t = createTables< tag::beta >( d, cnt );
    else if (d == ctr::DiffEqType::NUMFRACBETA)
      t = createTables< tag::numfracbeta >( d, cnt );
    else if (d == ctr::DiffEqType::MASSFRACBETA)
      t = createTables< tag::massfracbeta >( d, cnt );
    else if (d == ctr::DiffEqType::MIXNUMFRACBETA)
      t = createTables< tag::mixnumfracbeta >( d, cnt );
    else if (d == ctr::DiffEqType::MIXMASSFRACBETA)
      t = createTables< tag::mixmassfracbeta >( d, cnt );
    else if (d == ctr::DiffEqType::SKEWNORMAL)
      t = createTables< tag::skewnormal >( d, cnt );
    else if (d == ctr::DiffEqType::GAMMA)
      t = createTables< tag::gamma >( d, cnt );
    else Throw( "Can't find selected DiffEq" );

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
//! \author J. Bakosi
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DIRICHLET ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DIRICHLET ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::dirichlet, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::dirichlet >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::dirichlet >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::GENDIR ];  // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::GENDIR ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::gendir, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::gendir >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::gendir >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::WRIGHTFISHER ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::WRIGHTFISHER ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::wrightfisher >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::wrightfisher >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::OU ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::OU ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::ou, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::ou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::ou >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DIAG_OU ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DIAG_OU ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::diagou, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::diagou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::diagou >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::BETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::BETA ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::beta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::beta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::beta, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::beta >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::beta >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::NUMFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::NUMFRACBETA ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::numfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::numfracbeta >()[c] / 3;
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::MASSFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::MASSFRACBETA ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::massfracbeta, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::massfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::massfracbeta >()[c] / 3;
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::MIXNUMFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back(
    ctr::DiffEq().name( ctr::DiffEqType::MIXNUMFRACBETA ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::coeffpolicy >()[c] )
  );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::mixnumfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::mixnumfracbeta >()[c] / 3;
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::MIXMASSFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back(
    ctr::DiffEq().name( ctr::DiffEqType::MIXMASSFRACBETA ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::initpolicy >()[c] )
  );
  auto cp =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::coeffpolicy >()[c];
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name( cp ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::mixmassfracbeta >()[c] / 4;
  if (cp == ctr::CoeffPolicyType::HYDROTIMESCALE_HOMOGENEOUS_DECAY) {
    nfo.emplace_back(
      "inverse hydro time scales [" + std::to_string( ncomp ) + "]",
      options( ctr::HydroTimeScales(),
               g_inputdeck.get< tag::param,
                                tag::mixmassfracbeta,
                                tag::hydrotimescales >().at(c) ) );
    nfo.emplace_back(
      "production/dissipation [" + std::to_string( ncomp ) + "]",
      options( ctr::HydroProductions(),
               g_inputdeck.get< tag::param,
                                tag::mixmassfracbeta,
                                tag::hydroproductions >().at(c) ) );
  }
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::mixmassfracbeta >(c) ) );
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b' [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::bprime >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa' [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param,
                                 tag::mixmassfracbeta,
                                 tag::kappaprime >().at(c) ) );
  nfo.emplace_back(
    "coeff rho2 [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::rho2 >().at(c) ) );
  nfo.emplace_back(
    "coeff r [" + std::to_string( ncomp ) + "]",
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::SKEWNORMAL ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::SKEWNORMAL ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::skewnormal, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::skewnormal >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::skewnormal >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
//! \author J. Bakosi
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::GAMMA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::GAMMA ), "" );
  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::gamma, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::gamma >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::gamma >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );
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
