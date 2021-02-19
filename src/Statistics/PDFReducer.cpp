// *****************************************************************************
/*!
  \file      src/Statistics/PDFReducer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging PDFs across PEs
  \details   Custom Charm++ reducer for merging PDFs across PEs.
*/
// *****************************************************************************

#include <memory>

#include "PDFReducer.hpp"

namespace tk {

std::pair< int, std::unique_ptr<char[]> >
serialize( std::size_t meshid, const std::vector< tk::UniPDF >& u )
// *****************************************************************************
// Serialize univariate PDFs to raw memory stream
//! \param[in] meshid Mesh ID
//! \param[in] u Univariate PDFs
//! \return Pair of the length and the raw stream containing the serialized PDFs
// *****************************************************************************
{
  // Prepare for serializing PDF to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | meshid;
  sizer | const_cast< std::vector< tk::UniPDF >& >( u );

  // Create raw character stream to store the serialized PDF
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize PDF, the message will contain a univariate PDF
  PUP::toMem packer( flatData.get() );
  packer | meshid;
  packer | const_cast< std::vector< tk::UniPDF >& >( u );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeUniPDFs( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging a univariate PDFs during reduction across
// PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized PDF
//! \return Aggregated PDF built for further aggregation if needed
// *****************************************************************************
{
  // Will store deserialized univariate PDFs
  std::size_t meshid;
  std::vector< tk::UniPDF > updf;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize PDFs from raw stream
  creator | meshid;
  creator | updf;

  for (int m=1; m<nmsg; ++m) {
    // Unpack PDF
    std::size_t mid;
    std::vector< tk::UniPDF > u;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | mid;
    curCreator | u;
    // Merge PDFs
    meshid = mid;
    std::size_t i = 0;
    for (const auto& p : u) updf[i++].addPDF( p );
  }

  // Serialize vector of merged PDF to raw stream
  auto stream = tk::serialize( meshid, updf );

  // Forward serialized PDFs
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< tk::UniPDF >& u,
           const std::vector< tk::BiPDF >& b,
           const std::vector< tk::TriPDF >& t )
// *****************************************************************************
// Serialize vectors of PDFs to raw memory stream
//! \param[in] u Vector of univariate PDFs
//! \param[in] b Vector of bivariate PDFs
//! \param[in] t Vector of trivariate PDFs
//! \return Pair of the length and the raw stream containing the serialized PDFs
// *****************************************************************************
{
  // Prepare for serializing PDFs to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::vector< tk::UniPDF >& >( u );
  sizer | const_cast< std::vector< tk::BiPDF >& >( b );
  sizer | const_cast< std::vector< tk::TriPDF >& >( t );

  // Create raw character stream to store the serialized PDFs
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize PDFs, each message will contain a vector of PDFs
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::vector< tk::UniPDF >& >( u );
  packer | const_cast< std::vector< tk::BiPDF >& >( b );
  packer | const_cast< std::vector< tk::TriPDF >& >( t );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergePDF( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging PDFs during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized PDFs
//! \return Aggregated PDFs built for further aggregation if needed
// *****************************************************************************
{
  // Will store deserialized uni-, bi-, and tri-variate PDFs
  std::vector< tk::UniPDF > updf;
  std::vector< tk::BiPDF > bpdf;
  std::vector< tk::TriPDF > tpdf;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize PDFs from raw stream
  creator | updf;
  creator | bpdf;
  creator | tpdf;

  for (int m=1; m<nmsg; ++m) {
    // Unpack PDFs
    std::vector< tk::UniPDF > u;
    std::vector< tk::BiPDF > b;
    std::vector< tk::TriPDF > t;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | u;
    curCreator | b;
    curCreator | t;
    // Merge PDFs
    std::size_t i = 0;
    for (const auto& p : u) updf[i++].addPDF( p );
    i = 0;
    for (const auto& p : b) bpdf[i++].addPDF( p );
    i = 0;
    for (const auto& p : t) tpdf[i++].addPDF( p );
  }

  // Serialize vector of merged PDFs to raw stream
  auto stream = tk::serialize( updf, bpdf, tpdf );

  // Forward serialized PDFs
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::
