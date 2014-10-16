//******************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.C
  \author    J. Bakosi
  \date      Fri 10 Oct 2014 03:04:44 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Stack of differential equations
  \details   Stack of differential equations
*/
//******************************************************************************

#include <boost/mpl/cartesian_product.hpp>

#include <DiffEqStack.h>
#include <OrnsteinUhlenbeck.h>
#include <Dirichlet.h>
#include <GenDirichlet.h>
#include <Factory.h>

using quinoa::DiffEqStack;

DiffEqStack::DiffEqStack()
//******************************************************************************
//  Constructor: register all differential equations into factory
//! \author J. Bakosi
//******************************************************************************
{
  namespace mpl = boost::mpl;

  // Dirichlet SDE
  // Construct vector of vectors for all possible policies for SDE
  using DirPolicies = mpl::vector< InitPolicies, DirCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< DirPolicies >(
    registerDiffEq< Dirichlet >
                  ( m_factory, ctr::DiffEqType::DIRICHLET, m_eqTypes ) );

  // Lochner's generalized Dirichlet SDE
  // Construct vector of vectors for all possible policies for SDE
  using GenDirPolicies = mpl::vector< InitPolicies, GenDirCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< GenDirPolicies >(
    registerDiffEq< GenDirichlet >
                  ( m_factory, ctr::DiffEqType::GENDIR, m_eqTypes ) );

  // Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using OUPolicies = mpl::vector< InitPolicies, OUCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< OUPolicies >(
    registerDiffEq< OrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::OU, m_eqTypes ) );
}

std::vector< quinoa::DiffEq >
DiffEqStack::selected() const
//******************************************************************************
//  Instantiate all selected differential equations
//! \author J. Bakosi
//******************************************************************************
{
  std::map< ctr::DiffEqType, int > cnt; //!< Count DiffEqs per type
  std::vector< DiffEq > diffeqs;

  for (const auto& d : g_inputdeck.get< tag::selected, tag::diffeq >()) {
    if (d == ctr::DiffEqType::DIRICHLET)
      diffeqs.push_back( createDiffEq< tag::dirichlet >( m_factory, d, cnt ) );
    else if (d == ctr::DiffEqType::GENDIR)
      diffeqs.push_back( createDiffEq< tag::gendir >( m_factory, d, cnt ) );
    else if (d == ctr::DiffEqType::OU)
      diffeqs.push_back( createDiffEq< tag::ou >( m_factory, d, cnt ) );
    else Throw( "Can't find selected DiffEq" );
  }

  return diffeqs;
}

std::vector< std::vector< std::pair< std::string, std::string > > >
DiffEqStack::info() const
//******************************************************************************
//  Return information on all selected differential equations
//! \author J. Bakosi
//******************************************************************************
{
  std::map< ctr::DiffEqType, int > cnt; //!< Count DiffEqs per type
  std::vector< std::vector< std::pair< std::string, std::string > > > info;

  for (const auto& d : g_inputdeck.get< tag::selected, tag::diffeq >()) {
    if (d == ctr::DiffEqType::DIRICHLET)
      info.emplace_back( infoDirichlet( cnt ) );
    else if (d == ctr::DiffEqType::GENDIR)
      info.emplace_back( infoGenDir( cnt ) );
    else if (d == ctr::DiffEqType::OU)
      info.emplace_back( infoOU( cnt ) );
    else Throw( "Can't find selected DiffEq" );
  }

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoDirichlet( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the Dirichlet SDE
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
  info.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::dirichlet >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::dirichlet >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tk::tag::rng >()[c] ) );
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
  info.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::gendir, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::gendir >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::gendir >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::gendir, tk::tag::rng >()[c] ) );
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
DiffEqStack::infoOU( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the Ornstein-Uhlenbeck SDE
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
  info.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::ou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::ou >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::ou, tk::tag::rng >()[c] ) );
  info.emplace_back( "coeff sigma [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::ou, tag::sigma >(c) );
  info.emplace_back( "coeff theta [" + std::to_string( ncomp ) + "]",
    parameters< tag::param, tag::ou, tag::theta >(c) );
  info.emplace_back( "coeff mu [" + std::to_string( ncomp ) + "]",
    parameters< tag::param, tag::ou, tag::mu >(c) );

  return info;
}
