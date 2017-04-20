// *****************************************************************************
/*!
  \file      src/Statistics/PDFReducer.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Custom Charm++ reducer for merging PDFs across PEs
  \details   Custom Charm++ reducer for merging PDFs across PEs.
*/
// *****************************************************************************
#ifndef PDFReducer_h
#define PDFReducer_h

#include <tuple>
#include <vector>

#include "NoWarning/charm++.h"

#include "UniPDF.h"
#include "BiPDF.h"
#include "TriPDF.h"

namespace tk {

//! Serialize vectors of PDFs to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< tk::UniPDF >& u,
           const std::vector< tk::BiPDF >& b,
           const std::vector< tk::TriPDF >& t );

//! Charm++ custom reducer for merging PDFs during reduction across PEs
CkReductionMsg*
mergePDF( int nmsg, CkReductionMsg **msgs );

} // tk::

#endif // PDFReducer_h
