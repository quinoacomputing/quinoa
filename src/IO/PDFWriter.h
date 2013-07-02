//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Tue Jul  2 15:35:29 2013
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
#include <JPDF.h>

class Memory;

namespace Quinoa {

//! PDFWriter base class
class PDFWriter {

  public:
    //! Constructor: Acquire PDF file handle
    explicit PDFWriter(const std::string filename);

    //! Destructor: Release PDF file handle
    ~PDFWriter() noexcept;

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

    const std::string m_filename;            //!< PDF file name
    std::ofstream m_outPDF;                  //!< PDF file output stream
};

} // namespace Quinoa

#endif // PDFWriter_h
