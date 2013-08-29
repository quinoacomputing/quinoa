//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:07:21 2013
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

class Memory;

namespace quinoa {

//! PDFWriter : Writer
class PDFWriter : public Writer {

  public:
    //! Constructor
    explicit PDFWriter(const std::string& filename) :
      Writer(filename) {}

    //! Destructor: Release PDF file handle
    virtual ~PDFWriter() noexcept = default;

    //! Write PDF to file
    void write(const PDF* pdf);

    //! Write joint PDF to text file
    void writeTxt(const JPDF* jpdf);

    //! Write joint PDF to gmsh (text) file format
    void writeGmsh(const JPDF* jpdf);

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

} // namespace quinoa

#endif // PDFWriter_h
