//******************************************************************************
/*!
  \file      src/Statistics/PDFUtil.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:52:59 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     PDF utilities
  \brief     PDF utilities.
*/
//******************************************************************************

#include "PDFUtil.h"
#include "Make_unique.h"

namespace tk {

std::pair< int, std::unique_ptr<char[]> >
serialize( const std::tuple< std::vector< tk::UniPDF >,
                             std::vector< tk::BiPDF >,
                             std::vector< tk::TriPDF > >& pdf )
//******************************************************************************
// Serialize vectors of PDFs to raw memory stream
//! \param[in] pdf Tuple of vectors of uni-, bi-, and tri-variate PDFs
//! \return Pair of the length and the raw stream containing the serialized PDFs
//! \author J. Bakosi
//******************************************************************************
{
  // Prepare for serializing PDFs to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::vector<tk::UniPDF>& >( std::get<0>(pdf) );
  sizer | const_cast< std::vector<tk::BiPDF>& >( std::get<1>(pdf) );
  sizer | const_cast< std::vector<tk::TriPDF>& >( std::get<2>(pdf) );

  // Create raw character stream to store the serialized PDFs
  std::unique_ptr<char[]> flatData = tk::make_unique<char[]>( sizer.size() );

  // Serialize PDFs, each message will contain a vector of PDFs
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::vector<tk::UniPDF>& >( std::get<0>(pdf) );
  packer | const_cast< std::vector<tk::BiPDF>& >( std::get<1>(pdf) );
  packer | const_cast< std::vector<tk::TriPDF>& >( std::get<2>(pdf) );

  // Return size of and raw stream
  return { static_cast< int >( sizer.size() ), std::move(flatData) };
}

std::tuple< std::vector< tk::UniPDF >,
            std::vector< tk::BiPDF >,
            std::vector< tk::TriPDF > >
merge( CkReductionMsg* msg )
//******************************************************************************
// Deserialize and merge vectors of PDFs from Charm's CkReductionMsg
//! \param[in] msg Charm++ reduction message containing the serialized PDFs
//! \return Vector of merged PDFs
//! \author J. Bakosi
//******************************************************************************
{
  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msg->getData() );

  // Will store deserialized uni-, bi-, and tri-variate PDFs
  std::vector< tk::UniPDF > updf;
  std::vector< tk::BiPDF > bpdf;
  std::vector< tk::TriPDF > tpdf;

  // Deserialize PDFs from raw stream
  creator | updf;
  creator | bpdf;
  creator | tpdf;

  // Create tuple to hold all PDFs
  std::tuple< std::vector< tk::UniPDF >,
              std::vector< tk::BiPDF >,
              std::vector< tk::TriPDF > >
    res( std::vector< tk::UniPDF >( updf.size() ),
         std::vector< tk::BiPDF >( bpdf.size() ),
         std::vector< tk::TriPDF >( tpdf.size() ) );

  // Merge PDFs into tuple
  std::size_t i = 0;
  for (const auto& p : updf) std::get<0>(res)[i++].addPDF(p);
  i = 0;
  for (const auto& p : bpdf) std::get<1>(res)[i++].addPDF(p);
  i = 0;
  for (const auto& p : tpdf) std::get<2>(res)[i++].addPDF(p);

  // Return merged PDFs
  return res;
}

CkReductionMsg*
mergePDF( int nmsg, CkReductionMsg **msgs )
//******************************************************************************
// Charm++ custom reducer for merging PDFs during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized PDFs
//! \return Aggregated PDFs built for further aggregation if needed
//! \author J. Bakosi
//******************************************************************************
{
  // Will store merged PDFs in deserialized form
  std::tuple< std::vector< tk::UniPDF >,
              std::vector< tk::BiPDF >,
              std::vector< tk::TriPDF > > pdf;

  // Deserialize and merge vector of PDFs with partial sums
  for (int i=0; i<nmsg; ++i) pdf = tk::merge( msgs[i] );

  // Serialize vector of merged PDFs to raw stream
  auto stream = tk::serialize( pdf );

  // Forward serialized PDFs
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::
