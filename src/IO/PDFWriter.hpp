// *****************************************************************************
/*!
  \file      src/IO/PDFWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     PDF writer class declaration
  \details   This file declares a PDF writer class that facilitates outputing
    probability density functions (PDFs) into files in various formats using
    various configurations.
*/
// *****************************************************************************
#ifndef PDFWriter_h
#define PDFWriter_h

#include <string>
#include <iostream>

#include "Macro.hpp"
#include "Writer.hpp"
#include "UniPDF.hpp"
#include "BiPDF.hpp"
#include "TriPDF.hpp"
#include "StatCtr.hpp"
#include "Options/PDFCentering.hpp"
#include "Options/TxtFloatFormat.hpp"

namespace tk {

//! PDFWriter : Writer
class PDFWriter : public tk::Writer {

  public:
    //! Constructor
    explicit PDFWriter(
      const std::string& filename,
      tk::ctr::TxtFloatFormatType format = tk::ctr::TxtFloatFormatType::DEFAULT,
      std::streamsize precision = std::cout.precision() );

    //! Write univariate PDF to text file
    void writeTxt( const UniPDF& pdf, const tk::ctr::PDFInfo& info ) const;

    //! Write bivariate PDF to text file
    void writeTxt( const BiPDF& pdf, const tk::ctr::PDFInfo& info ) const;

    //! Write trivariate PDF to text file
    void writeTxt( const TriPDF& pdf, const tk::ctr::PDFInfo& info ) const;

    //! Write bivariate PDF to gmsh (text) file format
    void writeGmshTxt( const BiPDF& pdf, const tk::ctr::PDFInfo& info,
                       ctr::PDFCenteringType centering ) const;

    //! Write trivariate PDF to gmsh (text) file format
    void writeGmshTxt( const TriPDF& pdf, const tk::ctr::PDFInfo& info,
                       ctr::PDFCenteringType centering ) const;

    //! Write bivariate PDF to gmsh (binary) file format
    void writeGmshBin( const BiPDF& pdf, const tk::ctr::PDFInfo& info,
                       ctr::PDFCenteringType centering ) const;

    //! Write trivariate PDF to gmsh (binary) file format
    void writeGmshBin( const TriPDF& pdf, const tk::ctr::PDFInfo& info,
                       ctr::PDFCenteringType centering ) const;

    //! Write bivariate PDF to Exodus II file format
    void writeExodusII( const BiPDF& pdf, const tk::ctr::PDFInfo& info,
                        ctr::PDFCenteringType centering ) const;

    //! Write trivariate PDF to Exodus II file format
    void writeExodusII( const TriPDF& pdf, const tk::ctr::PDFInfo& info,
                        ctr::PDFCenteringType centering ) const;

  private:
    //! Assert the number of sample space dimensions given
    template< std::size_t size, class Container >
    void assertSampleSpaceDimensions( [[maybe_unused]] const Container& c )
    const {
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

    // Create Exodus II file
    int createExFile() const;

    // Write Exodus II file header
    void writeExHdr( int outFileId, int nnode, int nelem ) const;

    // Output probability density function as Exodus II results field
    void writeExVar( int exoFile,
                     ctr::PDFCenteringType centering,
                     const std::vector< tk::real >& probability ) const;

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
                  std::vector< tk::real >& outpdf,
                  ctr::PDFCenteringType centering ) const;

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
                  std::vector< tk::real >& outpdf,
                  ctr::PDFCenteringType centering ) const;
};

} // tk::

#endif // PDFWriter_h
