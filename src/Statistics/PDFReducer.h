// *****************************************************************************
/*!
  \file      src/Statistics/PDFReducer.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
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

//! Serialize univariate PDF to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< tk::UniPDF >& u );

//! \brief Charm++ custom reducer for merging a univariate PDF during reduction
//!    across PEs
CkReductionMsg*
mergeUniPDFs( int nmsg, CkReductionMsg **msgs );

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
