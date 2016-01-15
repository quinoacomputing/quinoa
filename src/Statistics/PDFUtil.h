//******************************************************************************
/*!
  \file      src/Statistics/PDFUtil.h
  \author    J. Bakosi
  \date      Thu 14 Jan 2016 02:38:42 PM MST
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
serialize( const std::vector< tk::UniPDF >& u,
           const std::vector< tk::BiPDF >& b,
           const std::vector< tk::TriPDF >& t );

//! Charm++ custom reducer for merging PDFs during reduction across PEs
CkReductionMsg*
mergePDF( int nmsg, CkReductionMsg **msgs );

} // tk::

#endif // PDFUtil_h
