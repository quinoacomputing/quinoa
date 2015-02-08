//******************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.C
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 05:49:37 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Stack of differential equations
  \details   This file defines class DiffEqStack, which implements various
    functionality related to registering and instantiating differential equation
    types. Registration and instantiation use a differential equation factory,
    which is a std::map (an associative container), associating unique
    differential equation keys to their constructor calls. For more details, see
    the in-code documentation of the constructor.
*/
//******************************************************************************

#include <boost/mpl/cartesian_product.hpp>

#include <DiffEqStack.h>
#include <OrnsteinUhlenbeck.h>
#include <DiagOrnsteinUhlenbeck.h>
#include <Dirichlet.h>
#include <GeneralizedDirichlet.h>
#include <WrightFisher.h>
#include <Beta.h>
#include <SkewNormal.h>
#include <Gamma.h>
#include <Factory.h>

using walker::DiffEqStack;

DiffEqStack::DiffEqStack()
//******************************************************************************
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
//!   to an std::function object that holds the a constructor bound to its
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
//******************************************************************************
{
  namespace mpl = boost::mpl;

  // Dirichlet SDE
  // Construct vector of vectors for all possible policies for SDE
  using DirPolicies = mpl::vector< tk::InitPolicies, DirichletCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< DirPolicies >(
    registerDiffEq< Dirichlet >
                  ( m_factory, ctr::DiffEqType::DIRICHLET, m_eqTypes ) );

  // Lochner's generalized Dirichlet SDE
  // Construct vector of vectors for all possible policies for SDE
  using GenDirPolicies =
    mpl::vector< tk::InitPolicies, GeneralizedDirichletCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< GenDirPolicies >(
    registerDiffEq< GeneralizedDirichlet >
                  ( m_factory, ctr::DiffEqType::GENDIR, m_eqTypes ) );

  // Wright-Fisher SDE
  // Construct vector of vectors for all possible policies for SDE
  using WrightFisherPolicies =
    mpl::vector< tk::InitPolicies, WrightFisherCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< WrightFisherPolicies >(
    registerDiffEq< WrightFisher >
                  ( m_factory, ctr::DiffEqType::WRIGHTFISHER, m_eqTypes ) );

  // Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using OrnsteinUhlenbeckPolicies =
    mpl::vector< tk::InitPolicies, OrnsteinUhlenbeckCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< OrnsteinUhlenbeckPolicies >(
    registerDiffEq< OrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::OU, m_eqTypes ) );

  // Diagonal Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using DiagOrnsteinUhlenbeckPolicies =
    mpl::vector< tk::InitPolicies, DiagOrnsteinUhlenbeckCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< DiagOrnsteinUhlenbeckPolicies >(
    registerDiffEq< DiagOrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::DIAG_OU, m_eqTypes ) );

  // Beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using BetaPolicies = mpl::vector< tk::InitPolicies, BetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< BetaPolicies >(
    registerDiffEq< Beta >
                  ( m_factory, ctr::DiffEqType::BETA, m_eqTypes ) );

  // Skew-normal SDE
  // Construct vector of vectors for all possible policies for SDE
  using SkewNormalPolicies =
    mpl::vector< tk::InitPolicies, SkewNormalCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< SkewNormalPolicies >(
    registerDiffEq< SkewNormal >
                  ( m_factory, ctr::DiffEqType::SKEWNORMAL, m_eqTypes ) );

  // Gamma SDE
  // Construct vector of vectors for all possible policies for SDE
  using GammaPolicies = mpl::vector< tk::InitPolicies, GammaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< GammaPolicies >(
    registerDiffEq< Gamma >
                  ( m_factory, ctr::DiffEqType::GAMMA, m_eqTypes ) );

}

std::vector< walker::DiffEq >
DiffEqStack::selected() const
//******************************************************************************
//  Instantiate all selected differential equations
//! \return std::vector of instantiated differential equation objects
//! \author J. Bakosi
//******************************************************************************
{
  std::map< ctr::DiffEqType, int > cnt; // count DiffEqs per type
  std::vector< DiffEq > diffeqs;        // will store instantiated DiffEqs

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
    else if (d == ctr::DiffEqType::SKEWNORMAL)
      diffeqs.push_back( createDiffEq< tag::skewnormal >( d, cnt ) );
    else if (d == ctr::DiffEqType::GAMMA)
      diffeqs.push_back( createDiffEq< tag::gamma >( d, cnt ) );
    else Throw( "Can't find selected DiffEq" );
  }

  return diffeqs;
}

std::vector< std::vector< std::pair< std::string, std::string > > >
DiffEqStack::info() const
//******************************************************************************
//  Return information on all selected differential equations
//! \return 
//! \author J. Bakosi
//******************************************************************************
{
  std::map< ctr::DiffEqType, int > cnt;         // count DiffEqs per type
  // will store info on all differential equations selected
  std::vector< std::vector< std::pair< std::string, std::string > > > info;

  for (const auto& d : g_inputdeck.get< tag::selected, tag::diffeq >()) {
    if (d == ctr::DiffEqType::DIRICHLET)
      info.emplace_back( infoDirichlet( cnt ) );
    else if (d == ctr::DiffEqType::GENDIR)
      info.emplace_back( infoGenDir( cnt ) );
    else if (d == ctr::DiffEqType::WRIGHTFISHER)
      info.emplace_back( infoWrightFisher( cnt ) );
    else if (d == ctr::DiffEqType::OU)
      info.emplace_back( infoOU( cnt ) );
    else if (d == ctr::DiffEqType::DIAG_OU)
      info.emplace_back( infoDiagOU( cnt ) );
    else if (d == ctr::DiffEqType::BETA)
      info.emplace_back( infoBeta( cnt ) );
    else if (d == ctr::DiffEqType::SKEWNORMAL)
      info.emplace_back( infoSkewNormal( cnt ) );
    else if (d == ctr::DiffEqType::GAMMA)
      info.emplace_back( infoGamma( cnt ) );
    else Throw( "Can't find selected DiffEq" );
  }

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoDirichlet( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the Dirichlet SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DIRICHLET ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DIRICHLET ), "" );
  info.emplace_back( "kind", "stochastic" );
  info.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::dirichlet, tag::depvar >()[c] ) );
  info.emplace_back( "initialization policy", tk::ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", tk::ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::dirichlet >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::dirichlet >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::rng >()[c] ) );
  info.emplace_back( "coeff b [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::dirichlet, tag::b >(c) );
  info.emplace_back( "coeff S [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::dirichlet, tag::S >(c) );
  info.emplace_back( "coeff kappa [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::dirichlet, tag::kappa >(c) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoGenDir( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on Lochner's generalized Dirichlet SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::GENDIR ];  // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::GENDIR ), "" );
  info.emplace_back( "kind", "stochastic" );
  info.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::gendir, tag::depvar >()[c] ) );
  info.emplace_back( "initialization policy", tk::ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", tk::ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::gendir >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::gendir >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::rng >()[c] ) );
  info.emplace_back( "coeff b [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gendir, tag::b >(c) );
  info.emplace_back( "coeff S [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gendir, tag::S >(c) );
  info.emplace_back( "coeff kappa [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gendir, tag::kappa >(c) );
  info.emplace_back( "coeff c [" + std::to_string( ncomp*(ncomp-1)/2 ) + "]",
                     parameters< tag::param, tag::gendir, tag::c >( c ) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoWrightFisher( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the Wright-Fisher SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::WRIGHTFISHER ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::WRIGHTFISHER ), "" );
  info.emplace_back( "kind", "stochastic" );
  info.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::depvar >()[c] ) );
  info.emplace_back( "initialization policy", tk::ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", tk::ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::wrightfisher >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::wrightfisher >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::rng >()[c] ) );
  info.emplace_back( "coeff omega [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::wrightfisher, tag::omega >(c) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoOU( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the Ornstein-Uhlenbeck SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::OU ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::OU ), "" );
  info.emplace_back( "kind", "stochastic" );
  info.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::ou, tag::depvar >()[c] ) );
  info.emplace_back( "initialization policy", tk::ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", tk::ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::ou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::ou >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::ou, tag::rng >()[c] ) );
  info.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp*(ncomp+1)/2 ) + ", upper tri]",
    parameters< tag::param, tag::ou, tag::sigmasq >(c) );
  info.emplace_back( "coeff theta [" + std::to_string( ncomp ) + "]",
    parameters< tag::param, tag::ou, tag::theta >(c) );
  info.emplace_back( "coeff mu [" + std::to_string( ncomp ) + "]",
    parameters< tag::param, tag::ou, tag::mu >(c) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoDiagOU( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the diagonal Ornstein-Uhlenbeck SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DIAG_OU ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DIAG_OU ), "" );
  info.emplace_back( "kind", "stochastic" );
  info.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::diagou, tag::depvar >()[c] ) );
  info.emplace_back( "initialization policy", tk::ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", tk::ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::diagou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::diagou >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::rng >()[c] ) );
  info.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp ) + "]",
    parameters< tag::param, tag::diagou, tag::sigmasq >(c) );
  info.emplace_back( "coeff theta [" + std::to_string( ncomp ) + "]",
    parameters< tag::param, tag::diagou, tag::theta >(c) );
  info.emplace_back( "coeff mu [" + std::to_string( ncomp ) + "]",
    parameters< tag::param, tag::diagou, tag::mu >(c) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoBeta( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::BETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::BETA ), "" );
  info.emplace_back( "kind", "stochastic" );
  info.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::beta, tag::depvar >()[c] ) );
  info.emplace_back( "initialization policy", tk::ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::beta, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", tk::ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::beta, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::beta >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::beta >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::beta, tag::rng >()[c] ) );
  info.emplace_back( "coeff b [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::beta, tag::b >(c) );
  info.emplace_back( "coeff S [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::beta, tag::S >(c) );
  info.emplace_back( "coeff kappa [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::beta, tag::kappa >(c) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoSkewNormal( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the skew-normal SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::SKEWNORMAL ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::SKEWNORMAL ), "" );
  info.emplace_back( "kind", "stochastic" );
  info.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::skewnormal, tag::depvar >()[c] ) );
  info.emplace_back( "initialization policy", tk::ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", tk::ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::skewnormal >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::skewnormal >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::rng >()[c] ) );
  info.emplace_back( "coeff T [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::skewnormal, tag::timescale >(c) );
  info.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp ) + "]",
    parameters< tag::param, tag::skewnormal, tag::sigmasq >(c) );
  info.emplace_back( "coeff lambda [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::skewnormal, tag::lambda >(c) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoGamma( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the gamma SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
//! \author J. Bakosi
//******************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::GAMMA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > info;

  info.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::GAMMA ), "" );
  info.emplace_back( "kind", "stochastic" );
  info.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::gamma, tag::depvar >()[c] ) );
  info.emplace_back( "initialization policy", tk::ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", tk::ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::gamma >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::gamma >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::rng >()[c] ) );
  info.emplace_back( "coeff b [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gamma, tag::b >(c) );
  info.emplace_back( "coeff S [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gamma, tag::S >(c) );
  info.emplace_back( "coeff kappa [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gamma, tag::kappa >(c) );

  return info;
}
