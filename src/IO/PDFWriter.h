//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Wed 10 Sep 2014 03:47:16 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     PDF writer
  \details   PDF writer
*/
//******************************************************************************
#ifndef PDFWriter_h
#define PDFWriter_h

#include <string>

#include <Writer.h>
#include <PDF.h>
#include <JPDF.h>

namespace quinoa {

//! PDFWriter : Writer
class PDFWriter : public tk::Writer {

  public:
    //! Constructor
    explicit PDFWriter( const std::string& filename ) : Writer( filename ) {}

    //! Write PDF to file
    void writeTxt( const PDF& pdf ) const;

    //! Write joint PDF to text file
    void writeTxt( const JPDF& jpdf ) const;

    //! Write joint PDF to gmsh (text) file format
    void writeGmsh( const JPDF& jpdf ) const;
};

} // quinoa::

#endif // PDFWriter_h
