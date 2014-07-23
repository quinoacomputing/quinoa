//******************************************************************************
/*!
  \file      src/Main/QuinoaPrint.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:32:40 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     QuinoaPrint
  \details   QuinoaPrint
*/
//******************************************************************************

#include <QuinoaPrint.h>
#include <Quinoa/Options/SDE.h>

using quinoa::QuinoaPrint;

void
QuinoaPrint::RequestedStats(const std::string& msg) const
//******************************************************************************
//  Echo requested statistics if differs from default
//! \author J. Bakosi
//******************************************************************************
{
  if (m_ctr.get< tag::stat >() != ctr::InputDeckDefaults.get< tag::stat >()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for ( auto& v : m_ctr.get< tag::stat >() ) {
      m_stream <<= v;
    }
    m_stream << '\n';
  }
}

void
QuinoaPrint::EstimatedStats(const std::string& msg) const
//******************************************************************************
//  Echo estimated statistics if differs from default
//! \author J. Bakosi
//******************************************************************************
{
  if (m_ctr.get<tag::stat>() != ctr::InputDeckDefaults.get<tag::stat>()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for (auto& v : m_ctr.get<tag::stat>()) {
       m_stream << v;
    }
    m_stream << '\n';
  }
}

std::vector< std::string >
QuinoaPrint::SDEPolicyNames( const ctr::SDEKey& key ) const
//******************************************************************************
//  Return SDE policies names
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< std::string > names;
  names.push_back( ctr::InitPolicy().name( key.get<tag::initpolicy>() ) );
  names.push_back( ctr::CoeffPolicy().name( key.get<tag::coeffpolicy>() ) );
  return names;
}

void
QuinoaPrint::printModel( const quinoa::Model& model ) const
//******************************************************************************
//  Echo configuration of a model
//! \author J. Bakosi
//******************************************************************************
{
  // Echo dependent variable
  m_stream << m_item_name_value_fmt % m_item_indent
                                    % "Dependent variable"
                                    % model.depvar();

  // Echo equation type and RNG if model is stochastic
  if (model.stochastic()) {
    m_stream << m_item_name_value_fmt % m_item_indent
                                      % "Equation"
                                      % "stochastic";
    tk::Option< tk::ctr::RNG > rng;
    m_stream << m_item_name_value_fmt % m_item_indent
                                      % rng.group()
                                      % rng.name( model.rng() );
  } else {
    // Only echo equation type if model is deterministic
    m_stream << m_item_name_value_fmt % m_item_indent
                                      % "Equation"
                                      % "deterministic";
  }

  // Echo initialization policy
  m_stream << m_item_name_value_fmt % m_item_indent
                                    % "Init policy"
                                    % model.initPolicy();
  // Echo coefficients policy
  m_stream << m_item_name_value_fmt % m_item_indent
                                    % "Coefficients policy"
                                    % model.coeffPolicy();
  // Echo number of components
  m_stream << m_item_name_value_fmt % m_item_indent
                                    % "Number of components"
                                    % model.ncomp();
  m_stream << '\n';
}
