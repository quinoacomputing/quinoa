//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Wed 08 Oct 2014 08:31:36 AM MDT
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
#include <TriPDF.h>
#include <Quinoa/Options/PDFCentering.h>
#include <Quinoa/InputDeck/InputDeck.h>

namespace quinoa {

//! PDFWriter : Writer
class PDFWriter : public tk::Writer {

  public:
    //! Constructor
    explicit PDFWriter( const std::string& filename,
                        ctr::TxtFloatFormatType format =
                          ctr::TxtFloatFormatType::DEFAULT,
                        std::streamsize precision = std::cout.precision() );

    //! Write univariate PDF to text file
    void writeTxt( const UniPDF& pdf, const ctr::InputDeck::PDFInfo& info )
    const;

    //! Write bivariate PDF to text file
    void writeTxt( const BiPDF& pdf, const ctr::InputDeck::PDFInfo& info )
    const;

    //! Write trivariate PDF to text file
    void writeTxt( const TriPDF& pdf, const ctr::InputDeck::PDFInfo& info )
    const;

    //! Write bivariate PDF to gmsh (text) file format
    void writeGmshTxt( const BiPDF& pdf, const ctr::InputDeck::PDFInfo& info,
                       ctr::PDFCenteringType centering ) const;

    //! Write trivariate PDF to gmsh (text) file format
    void writeGmshTxt( const TriPDF& pdf, const ctr::InputDeck::PDFInfo& info,
                       ctr::PDFCenteringType centering ) const;

    //! Write bivariate PDF to gmsh (binary) file format
    void writeGmshBin( const BiPDF& pdf, const ctr::InputDeck::PDFInfo& info,
                       ctr::PDFCenteringType centering ) const;

    //! Write trivariate PDF to gmsh (binary) file format
    void writeGmshBin( const TriPDF& pdf, const ctr::InputDeck::PDFInfo& info,
                       ctr::PDFCenteringType centering ) const;

  private:
    //! Assert the number of sample space dimensions given
    template< std::size_t size, class Container >
    void assertSampleSpaceDimensions( const Container& c ) const {
      Assert( c.size() == size,
              "Number of sample space variables must equal " +
              std::to_string( size ) + " in PDF writer." );
    }

    //! Assert the number of sample space extents given
    template< std::size_t size, class Container >
    void assertSampleSpaceExtents( const Container& c ) const {
      if (!c.empty())
        Assert( c.size() == size*2,
                "PDF user-specified sample space extents must be defined by " +
                std::to_string( size*2 ) +" real numbers: minx, maxx, ..." );
    }

    //! Query extents and other metadata of univariate PDF sample space
    void extents( const UniPDF& pdf,
                  const std::vector< tk::real >& uext,
                  std::size_t& nbi,
                  tk::real& min,
                  tk::real& max,
                  tk::real& binsize,
                  std::array< long, 2*UniPDF::dim >& ext,
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
                  std::array< tk::real, BiPDF::dim >& binsize,
                  std::array< long, 2*BiPDF::dim >& ext,
                  std::vector< tk::real >& outpdf ) const;

    //! Query extents and other metadata of trivariate PDF sample space
    void extents( const TriPDF& pdf,
                  const std::vector< tk::real >& uext,
                  std::size_t& nbix,
                  std::size_t& nbiy,
                  std::size_t& nbiz,
                  tk::real& xmin,
                  tk::real& xmax,
                  tk::real& ymin,
                  tk::real& ymax,
                  tk::real& zmin,
                  tk::real& zmax,
                  std::array< tk::real, TriPDF::dim >& binsize,
                  std::array< long, 2*TriPDF::dim >& ext,
                  std::vector< tk::real >& outpdf ) const;
};

} // quinoa::

#endif // PDFWriter_h
