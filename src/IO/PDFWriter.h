//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Thu 25 Oct 2012 06:13:46 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     PDF writer
  \details   PDF writer
*/
//******************************************************************************
#ifndef PDFWriter_h
#define PDFWriter_h

#include <string>
#include <fstream>

#include <PDF.h>

using namespace std;

namespace Quinoa {

//! PDFWriter base class
class PDFWriter {

  public:
    //! Constructor: Acquire PDF file handle
    PDFWriter(const string filename);

    //! Destructor: Release PDF file handle
    ~PDFWriter();

    //! Write PDF to file
    void write(const PDF* pdf);

    //! PDF file name
    const string m_filename;

    //! PDF file output stream
    ofstream m_outPDF;

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

} // namespace Quinoa

#endif // PDFWriter_h
