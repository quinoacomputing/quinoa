//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Wed Apr 23 11:17:19 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
    explicit PDFWriter(const std::string& filename) :
      Writer(filename) {}

    //! Destructor: Release PDF file handle
    ~PDFWriter() noexcept override = default;

    //! Write PDF to file
    void writeTxt(const tk::PDF& pdf);

    //! Write joint PDF to text file
    void writeTxt(const tk::JPDF& jpdf);

    //! Write joint PDF to gmsh (text) file format
    void writeGmsh(const tk::JPDF& jpdf);

  private:
    //! Don't permit copy constructor
    PDFWriter(const PDFWriter&) = delete;
    //! Don't permit copy assigment
    PDFWriter& operator=(const PDFWriter&) = delete;
    //! Don't permit move constructor
    PDFWriter(PDFWriter&&) = delete;
    //! Don't permit move assigment
    PDFWriter& operator=(PDFWriter&&) = delete;
};

} // quinoa::

#endif // PDFWriter_h
