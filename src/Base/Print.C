//******************************************************************************
/*!
  \file      src/Base/Print.C
  \author    J. Bakosi
  \date      Fri 07 Feb 2014 10:07:17 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Print
  \details   Print
*/
//******************************************************************************

#include <Print.h>

using tk::Print;

void
Print::header(const std::string& title) const
//******************************************************************************
//  Print header
//! \author J. Bakosi
//******************************************************************************
{
  m_stream << m_header_fmt % boost::io::group(std::setfill('='), "");
  m_stream << std::endl;
  m_stream << m_header_fmt % title;
  m_stream << std::endl;
  m_stream << m_header_fmt % boost::io::group(std::setfill('='), "");
}

void
Print::part(const std::string& title) const
//******************************************************************************
//  Print part header: title
//! \author J. Bakosi
//******************************************************************************
{
  using std::operator+;
  std::string::size_type half_length = title.size()/2 + 1;
  std::string s(half_length, '-');
  std::string underline(s + " o " + s);
  std::string upper(title);
  std::transform(title.begin(), title.end(), upper.begin(), ::toupper);
  upper = "< " + upper + " >";
  m_stream << m_part_fmt % upper;
  m_stream << m_part_underline_fmt % underline;
}

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
