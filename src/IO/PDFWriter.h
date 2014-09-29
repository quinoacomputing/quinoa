//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Sun 28 Sep 2014 09:52:40 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     PDF writer
  \details   PDF writer
*/
//******************************************************************************
#ifndef PDFWriter_h
#define PDFWriter_h

#include <string>

#include <Writer.h>
#include <UniPDF.h>
#include <BiPDF.h>
#include <Quinoa/Options/PDFCentering.h>

namespace quinoa {

//! PDFWriter : Writer
class PDFWriter : public tk::Writer {

  public:
    //! Constructor
    explicit PDFWriter( const std::string& filename ) : Writer( filename ) {}

    //! Write univariate PDF to text file
    void writeTxt( const UniPDF& pdf, const std::vector< tk::real >& uext )
    const;

    //! Write bivariate PDF to text file
    void writeTxt( const BiPDF& pdf, const std::vector< tk::real >& uext )
    const;

    //! Write bivariate PDF to gmsh (text) file format
    void writeGmsh( const BiPDF& pdf,
                    const std::string& pdfname,
                    ctr::PDFCenteringType centering,
                    const std::vector< tk::real >& uext ) const;

  private:
    //! Query extents and other metadata of univariate PDF sample space
    void extents( const UniPDF& pdf,
                  const std::vector< tk::real >& uext,
                  std::size_t& nbi,
                  tk::real& min,
                  tk::real& max,
                  tk::real& binsize,
                  std::pair< UniPDF::key_type, UniPDF::key_type >& ext,
                  std::vector< tk::real >& outpdf ) const;

    //! Query extents and other metadata of bivariate PDF sample space
    void extents( const BiPDF& pdf,
                  const std::vector< tk::real >& uext,
                  std::size_t& nbix,
                  std::size_t& nbiy,
                  tk::real& xmin,
                  tk::real& xmax,
                  tk::real& ymin,
                  tk::real& ymax,
                  std::array< tk::real, 2 >& binsize,
                  std::pair< BiPDF::key_type, BiPDF::key_type >& ext,
                  std::vector< tk::real >& outpdf ) const;
};

} // quinoa::

#endif // PDFWriter_h
