//******************************************************************************
/*!
  \file      src/Statistics/PDFUtil.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:10:29 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     PDF utilities
  \brief     PDF utilities.
*/
//******************************************************************************
#ifndef PDFUtil_h
#define PDFUtil_h

#include <tuple>
#include <vector>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <charm++.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include "UniPDF.h"
#include "BiPDF.h"
#include "TriPDF.h"

namespace tk {

//! Serialize vectors of PDFs to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::tuple< std::vector< tk::UniPDF >,
                             std::vector< tk::BiPDF >,
                             std::vector< tk::TriPDF > >& pdf );

//! Deserialize and merge vectors of PDFs from Charm's CkReductionMsg
std::tuple< std::vector< tk::UniPDF >,
            std::vector< tk::BiPDF >,
            std::vector< tk::TriPDF > >
merge( CkReductionMsg* msg );

//! Charm++ custom reducer for merging PDFs during reduction across PEs
CkReductionMsg*
mergePDF( int nmsg, CkReductionMsg **msgs );

} // tk::

#endif // PDFUtil_h
