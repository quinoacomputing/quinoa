//******************************************************************************
/*!
  \file      src/Main/RNGPrint.C
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 12:19:57 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Pretty printer base for pretty printers supporting RNGs
  \details   Pretty printer base for pretty printers supporting RNGs.
*/
//******************************************************************************

#include <RNGPrint.h>

using tk::RNGPrint;

#ifdef HAS_MKL
void
RNGPrint::echoMKLParams( const ctr::RNGMKLParam& p ) const
//******************************************************************************
//  Echo information on MKL random number generator
//! \param[in] p MKL RNG parameters
//! \author J. Bakosi
//******************************************************************************
{
  ctr::MKLUniformMethod um;
  ctr::MKLGaussianMethod gm;

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
#endif

void
RNGPrint::echoRNGSSEParams( const ctr::RNGSSEParam& p,
                            const ctr::RNG& rng,
                            const ctr::RNGType& r ) const
//******************************************************************************
//  Echo information on RNGSSE random number generator
//! \param[in] p RNGSSE RNG parameters
//! \param[in] rng RNG options object
//! \param[in] r RNG type enum
//! \author J. Bakosi
//******************************************************************************
{
  m_stream << m_item_name_value_fmt
              % m_item_indent
              % "seed"
              % p.get< tag::seed >();

  if ( rng.supportsSeq(r) ) {
    ctr::RNGSSESeqLen seq;
    m_stream << m_item_name_value_fmt
                % m_item_indent
                % seq.group()
                % seq.name( p.get< tag::seqlen >() );
  }
}

#ifdef HAS_MKL
void
RNGPrint::MKLParams( const std::vector< ctr::RNGType >& vec,
                     const ctr::RNGMKLParameters& map ) const
//******************************************************************************
//  Print all fields of MKL RNG parameters
//! \param[in] vec Vector of RNG type enums to print
//! \param[in] map MKL RNG parameters map
//! \author J. Bakosi
//******************************************************************************
{
  ctr::RNG rng;

  for (auto& r : vec) {
    if (rng.lib(r) == ctr::RNGLibType::MKL) {
      subsection( rng.name(r) );
      const auto& m = map.find(r);
      if (m == map.end()) {   // no parameter map entry, print defaults
        echoMKLParams( ctr::RNGMKLParam() );
      } else {
        echoMKLParams( m->second );
      }
    }
  }
}
#endif

void
RNGPrint::RNGSSEParams( const std::vector< ctr::RNGType >& vec,
                        const ctr::RNGSSEParameters& map ) const
//******************************************************************************
//  Print all fields of RNGSSE RNG parameters
//! \param[in] vec Vector of RNG type enums to print
//! \param[in] map RNGSSE RNG parameters map
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
