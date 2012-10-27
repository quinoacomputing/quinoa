//******************************************************************************
/*!
  \file      src/IO/PDFWriter.C
  \author    J. Bakosi
  \date      Sat 27 Oct 2012 03:25:25 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Univariate PDF writer
  \details   Univariate PDF writer
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
//! \param[in]  filename  File name to open for writing
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
  // a warning, thus we avoid terminate if an exception is propagating through
  if (m_outPDF.fail())
    cerr << "WARNING: Failed to close file: " << m_filename << endl;
}

void
PDFWriter::write(const PDF* pdf)
//******************************************************************************
//  Write out standardized PDF to file
//! \param[in]  pdf  Object pointer to univariate PDF
//! \author  J. Bakosi
//******************************************************************************
{
  auto f = pdf->getMap();
  const real binsize = pdf->getBinsize();
  const real sp = pdf->getNsample()*binsize;
  for (auto& p : *f) m_outPDF << p.first*binsize << "\t" << p.second/sp << endl;
}

void
PDFWriter::write(const JPDF* jpdf)
//******************************************************************************
//  Write out standardized joint PDF to file
//! \param[in]  jpdf  Object pointer to joint PDF
//! \author  J. Bakosi
//******************************************************************************
{
  auto f = jpdf->getMap();
  real binsize = jpdf->getBinsize();
  const real sp = jpdf->getNsample()*binsize*binsize;
  for (auto& p : *f) {
    m_outPDF << p.first[0]*binsize << " " << p.first[1]*binsize
             << " " << p.second/sp << endl;
  }
}
