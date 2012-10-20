//******************************************************************************
/*!
  \file      src/IO/PDFWriter.C
  \author    J. Bakosi
  \date      Sat 20 Oct 2012 10:43:15 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     PDF writer
  \details   PDF writer
*/
//******************************************************************************

#include <string>
#include <fstream>

#include <PDFWriter.h>
#include <IOException.h>

using namespace Quinoa;

PDFWriter::PDFWriter(const string filename) : m_filename(filename)
//******************************************************************************
//  Constructor: Acquire PDF file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outPDF.open(m_filename, ofstream::out);
  if (!m_outPDF.good()) throw IOException(FATAL, FAILED_OPEN, m_filename);
}

PDFWriter::~PDFWriter()
//******************************************************************************
//  Destructor: Release PDF file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outPDF.close();
  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
  if (m_outPDF.fail())
    cerr << "WARNING: Failed to close file: " << m_filename << endl;
}

void
PDFWriter::write(const PDF* pdf)
//******************************************************************************
//  Write out standardized PDF to file
//! \author  J. Bakosi
//******************************************************************************
{
  const Pdf* f = pdf->getPDF();
  const real binsize = pdf->getBinsize();
  const real sp = pdf->getNsample()*binsize;
  for (auto& p : *f) {
    m_outPDF << p.first*binsize << "\t" << p.second/sp << endl;
  }
}
