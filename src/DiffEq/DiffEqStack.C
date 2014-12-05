//******************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.C
  \author    J. Bakosi
  \date      Fri 05 Dec 2014 03:13:24 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Stack of differential equations
  \details   Stack of differential equations
*/
//******************************************************************************

#include <boost/mpl/cartesian_product.hpp>

#include <DiffEqStack.h>
#include <OrnsteinUhlenbeck.h>
#include <DiagOrnsteinUhlenbeck.h>
#include <Dirichlet.h>
#include <GenDirichlet.h>
#include <WrightFisher.h>
#include <Beta.h>
#include <SkewNormal.h>
#include <Gamma.h>
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

  // Wright-Fisher SDE
  // Construct vector of vectors for all possible policies for SDE
  using WFPolicies = mpl::vector< InitPolicies, WFCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< WFPolicies >(
    registerDiffEq< WrightFisher >
                  ( m_factory, ctr::DiffEqType::WRIGHTFISHER, m_eqTypes ) );

  // Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using OUPolicies = mpl::vector< InitPolicies, OUCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< OUPolicies >(
    registerDiffEq< OrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::OU, m_eqTypes ) );

  // Diagonal Ornstein-Uhlenbeck SDE
  // Construct vector of vectors for all possible policies for SDE
  using DiagOUPolicies = mpl::vector< InitPolicies, DiagOUCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< DiagOUPolicies >(
    registerDiffEq< DiagOrnsteinUhlenbeck >
                  ( m_factory, ctr::DiffEqType::DIAG_OU, m_eqTypes ) );

  // Beta SDE
  // Construct vector of vectors for all possible policies for SDE
  using BetaPolicies = mpl::vector< InitPolicies, BetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< BetaPolicies >(
    registerDiffEq< Beta >
                  ( m_factory, ctr::DiffEqType::BETA, m_eqTypes ) );

  // Skew-normal SDE
  // Construct vector of vectors for all possible policies for SDE
  using SkewNormalPolicies = mpl::vector< InitPolicies, SkewNormalCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< SkewNormalPolicies >(
    registerDiffEq< SkewNormal >
                  ( m_factory, ctr::DiffEqType::SKEWNORMAL, m_eqTypes ) );

  // Gamma SDE
  // Construct vector of vectors for all possible policies for SDE
  using GammaPolicies = mpl::vector< InitPolicies, GammaCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< GammaPolicies >(
    registerDiffEq< Gamma >
                  ( m_factory, ctr::DiffEqType::GAMMA, m_eqTypes ) );

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
    else if (d == ctr::DiffEqType::WRIGHTFISHER)
      diffeqs.push_back( createDiffEq< tag::wrightfisher >( m_factory, d, cnt ) );
    else if (d == ctr::DiffEqType::OU)
      diffeqs.push_back( createDiffEq< tag::ou >( m_factory, d, cnt ) );
    else if (d == ctr::DiffEqType::DIAG_OU)
      diffeqs.push_back( createDiffEq< tag::diagou >( m_factory, d, cnt ) );
    else if (d == ctr::DiffEqType::BETA)
      diffeqs.push_back( createDiffEq< tag::beta >( m_factory, d, cnt ) );
    else if (d == ctr::DiffEqType::SKEWNORMAL)
      diffeqs.push_back( createDiffEq< tag::skewnormal >( m_factory, d, cnt ) );
    else if (d == ctr::DiffEqType::GAMMA)
      diffeqs.push_back( createDiffEq< tag::gamma >( m_factory, d, cnt ) );
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
DiffEqStack::infoWrightFisher( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the Wright-Fisher SDE
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
  info.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::wrightfisher >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::wrightfisher >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tk::tag::rng >()[c] ) );
  info.emplace_back( "coeff omega [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::wrightfisher, tag::omega >(c) );

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
  info.emplace_back( "coeff sigma [" + std::to_string( ncomp*(ncomp+1)/2 )
                     + ", upper tri]",
                     parameters< tag::param, tag::ou, tag::sigma >(c) );
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
  info.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::diagou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::diagou >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::diagou, tk::tag::rng >()[c] ) );
  info.emplace_back( "coeff sigma [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::diagou, tag::sigma >(c) );
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
  info.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::beta, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::beta, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::beta >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::beta >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::beta, tk::tag::rng >()[c] ) );
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
  info.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::skewnormal >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::skewnormal >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tk::tag::rng >()[c] ) );
  info.emplace_back( "coeff T [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::skewnormal, tag::timescale >(c) );
  info.emplace_back( "coeff sigma [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::skewnormal, tag::sigma >(c) );
  info.emplace_back( "coeff lambda [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::skewnormal, tag::lambda >(c) );

  return info;
}

std::vector< std::pair< std::string, std::string > >
DiffEqStack::infoGamma( std::map< ctr::DiffEqType, int >& cnt ) const
//******************************************************************************
//  Return information on the gamma SDE
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
  info.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::initpolicy >()[c] ) );
  info.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::gamma, tag::coeffpolicy >()[c] ) );
  info.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::gamma >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::gamma >()[c];
  info.emplace_back( "number of components", std::to_string( ncomp ) );
  info.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::gamma, tk::tag::rng >()[c] ) );
  info.emplace_back( "coeff b [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gamma, tag::b >(c) );
  info.emplace_back( "coeff S [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gamma, tag::S >(c) );
  info.emplace_back( "coeff kappa [" + std::to_string( ncomp ) + "]",
                     parameters< tag::param, tag::gamma, tag::kappa >(c) );

  return info;
}
