//******************************************************************************
/*!
  \file      src/Base/Print.C
  \author    J. Bakosi
  \date      Tue 28 Jan 2014 03:56:15 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Print
  \details   Print
*/
//******************************************************************************

#include <Print.h>

using tk::Print;

void
Print::echoMKLParams( const ctr::MKLRNGParam& p ) const
//******************************************************************************
//  Echo information on MKL random number generator
//! \author J. Bakosi
//******************************************************************************
{
  Option< ctr::MKLUniformMethod > um;
  Option< ctr::MKLGaussianMethod > gm;

  m_stream << m_item_name_value_fmt
              % m_item_indent
              % "seed"
              % p.get< tag::seed >();

  m_stream << m_item_name_value_fmt
              % m_item_indent
              % um.group()
              % um.name( p.get< tag::uniform_method >() );

  m_stream << m_item_name_value_fmt
              % m_item_indent
              % gm.group()
              % gm.name( p.get< tag::gaussian_method >() );
}

void
Print::echoRNGSSEParams( const ctr::RNGSSEParam& p,
                         const ctr::RNG& rng,
                         const ctr::RNGType& r ) const
//******************************************************************************
//  Echo information on RNGSSE random number generator
//! \author J. Bakosi
//******************************************************************************
{
  m_stream << m_item_name_value_fmt
              % m_item_indent
              % "seed"
              % p.get< tag::seed >();

  if ( rng.supportsSeq(r) ) {
    Option< ctr::RNGSSESeqLen > seq;
    m_stream << m_item_name_value_fmt
                % m_item_indent
                % seq.group()
                % seq.name( p.get< tag::seqlen >() );
  }
}


void
Print::MKLParams( const std::vector< ctr::RNGType >& vec,
                  const ctr::MKLRNGParameters& map ) const
//******************************************************************************
//  Print all fields of MKL RNG parameters
//! \author J. Bakosi
//******************************************************************************
{
  ctr::RNG rng;

  for (auto& r : vec) {
    if (rng.lib(r) == ctr::RNGLibType::MKL) {
      subsection( rng.name(r) );
      const auto& m = map.find(r);
      if (m == map.end()) {   // no parameter map entry, print defaults
        echoMKLParams( ctr::MKLRNGParam() );
      } else {
        echoMKLParams( m->second );
      }
    }
  }
}

void
Print::RNGSSEParams( const std::vector< ctr::RNGType >& vec,
                     const ctr::RNGSSEParameters& map ) const
//******************************************************************************
//  Print all fields of RNGSSE RNG parameters
//! \author J. Bakosi
//******************************************************************************
{
  ctr::RNG rng;

  for (auto& r : vec) {
    if (rng.lib(r) == ctr::RNGLibType::RNGSSE) {
      subsection( rng.name(r) );
      const auto& m = map.find(r);
      if (m == map.end()) {   // no parameter map entry, print defaults
        echoRNGSSEParams( ctr::RNGSSEParam(), rng, r );
      } else {
        echoRNGSSEParams( m->second, rng, r );
      }
    }
  }
}
