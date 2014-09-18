//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Thu 18 Sep 2014 11:51:46 AM MDT
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
    void writeTxt( const UniPDF& pdf ) const;

    //! Write bivariate PDF to text file
    void writeTxt( const BiPDF& pdf ) const;

    //! Write bivariate PDF to gmsh (text) file format
    void writeGmsh( const BiPDF& pdf,
                    const std::string& pdfname,
                    ctr::PDFCenteringType centering ) const;
};

} // quinoa::

#endif // PDFWriter_h
