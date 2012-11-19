//******************************************************************************
/*!
  \file      src/IO/PDFWriter.h
  \author    J. Bakosi
  \date      Sun 18 Nov 2012 06:44:27 PM MST
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

using namespace std;

class Memory;

namespace Quinoa {

//! PDFWriter base class
class PDFWriter {

  public:
    //! Constructor: Acquire PDF file handle
    PDFWriter(Memory* memory, const string filename);

    //! Destructor: Release PDF file handle
    ~PDFWriter();

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

    Memory* m_memory;                   //!< Memory object pointer 
    const string m_filename;            //!< PDF file name
    ofstream m_outPDF;                  //!< PDF file output stream
};

} // namespace Quinoa

#endif // PDFWriter_h
