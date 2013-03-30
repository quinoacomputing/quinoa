//******************************************************************************
/*!
  \file      src/IO/PDFWriter.C
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 06:26:11 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Univariate PDF writer
  \details   Univariate PDF writer
*/
//******************************************************************************

#include <string>
#include <fstream>
#include <map>

#include <Memory.h>
#include <PDFWriter.h>
#include <IOException.h>

using namespace Quinoa;

PDFWriter::PDFWriter(const string filename) :
  m_filename(filename)
//******************************************************************************
//  Constructor: Acquire PDF file handle
//! \param[in]  filename  File name to open for writing
//! \author J. Bakosi
//******************************************************************************
{
  m_outPDF.open(m_filename, ofstream::out);
  Assert(m_outPDF.good(), IOException,FATAL,IO_FAILED_OPEN,m_filename);
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
    cout << "WARNING: Failed to close file: " << m_filename << endl;
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
PDFWriter::writeTxt(const JPDF* jpdf)
//******************************************************************************
//  Write out standardized joint PDF to text file (only for 2D)
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

void
PDFWriter::writeGmsh(const JPDF* jpdf)
//******************************************************************************
//  Write out standardized joint PDF to Gmsh (text) format (only for 2D)
//! \param[in]  jpdf  Object pointer to joint PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Output mesh header: mesh version, file type, data size
  m_outPDF << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  Assert(!m_outPDF.bad(), IOException,FATAL,IO_FAILED_WRITE,m_filename);

  auto f = jpdf->getMap();
  real binsize = jpdf->getBinsize();
  const real sp = jpdf->getNsample()*binsize*binsize;

  // Find sample space extents
  real xmin = numeric_limits<real>::max();
  real xmax = numeric_limits<real>::min();
  real ymin = numeric_limits<real>::max();
  real ymax = numeric_limits<real>::min();
  for (auto& p : *f) {
    if (binsize*p.first[0] < xmin) xmin = binsize*p.first[0];
    if (binsize*p.first[0] > xmax) xmax = binsize*p.first[0];
    if (binsize*p.first[1] < ymin) ymin = binsize*p.first[1];
    if (binsize*p.first[1] > ymax) ymax = binsize*p.first[1];
  }
  int nbix = static_cast<int>((xmax - xmin)/binsize + 1);
  int nbiy = static_cast<int>((ymax - ymin)/binsize + 1);

  // Output points of discretized sample space
  m_outPDF << "$Nodes\n" << nbix*nbiy << endl;
  int k=0;
  for (int i=0; i<nbix; i++ ) {
    real x = xmin + i*binsize;
    for (int j=0; j<nbiy; j++ ) {
      real y = ymin + j*binsize;
      m_outPDF << k++ << " " << x << " " << y << " 0\n";
    }
  }
  m_outPDF << "$EndNodes\n";

  // Output elements of discretized sample space
  --nbix;  --nbiy;
  m_outPDF << "$Elements\n" << nbix*nbiy << "\n";
  for (int i=0; i<nbix*nbiy; ++i) {
    m_outPDF << i << " 3 2 1 1 " << i+i/nbiy << " " << i+1+i/nbiy << " "
             << i+2+nbiy+i/nbiy << " " << i+1+nbiy+i/nbiy << "\n";
  }
  m_outPDF << "$EndElements\n";

  // Output function values
  ++nbiy;
  m_outPDF << "$NodeData\n1\n\"Computed\"\n1\n0.0\n3\n0\n1\n"
           << f->size() << "\n";
  for (auto& p : *f) {
    int bin = static_cast<int>(
                (binsize*p.first[0] - xmin)/(xmax - xmin)*nbix*nbiy +
                (binsize*p.first[1] - ymin)/(ymax - ymin)*nbiy );
    m_outPDF << bin << " " << p.second/sp << endl;
  }

  m_outPDF << "$NodetData\n";
}
