//******************************************************************************
/*!
  \file      src/Control/PDFInfo.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 09:35:00 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     PDF information
  \details   PDF information
*/
//******************************************************************************
#ifndef PDFInfo_h
#define PDFInfo_h

#include <ControlTypes.h>

namespace tk {
namespace ctr {

//! PDF information bundle
struct PDFInfo {
  const std::string& name;                  //!< identifier
  const std::vector< tk::real >& exts;      //!< extents
  std::vector< std::string > vars;          //!< dependent variables
};

//! Return sample space variables for PDF
template< std::size_t d > std::vector< std::string >
vars( const std::vector< Probability >& pdfs, long int idx )
{
  long int n = -1;
  std::vector< std::string > v;
  for (const auto& probability : pdfs) {
    if (probability.size() == d) ++n;
    if (n == idx) {
      for (const auto& term : probability)
        v.push_back( term.var + std::to_string(term.field+1) );
      return v;
    }
  }
  Throw( "Cannot find PDF." );
}

//! Find PDF information given the sample space dimension and its index
//! \param[in]  idx  Index of the PDF with given sample space dimension
template< std::size_t d >
PDFInfo pdf( const std::vector< std::vector< tk::real > >& binsizes,
             const std::vector< std::string >& names,
             const std::vector< std::vector< tk::real > >& exts,
             const std::vector< Probability >& pdfs,
             long int idx )
{
  Assert( binsizes.size() == names.size(),
          "Number of binsizes vector and the number of PDF names must "
          "equal in InputDeck::pdfname()." );
  long int n = -1;
  long int i = 0;
  for (const auto& bs : binsizes) {
    if (bs.size() == d) ++n;
    if (n == idx) return { names[i], exts[i], vars<d>(pdfs,idx) };
    ++i;
  }
  Throw( "Cannot find PDF name." );
}

} // ctr::
} // tk::

#endif // PDFInfo_h
