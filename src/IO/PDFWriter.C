//******************************************************************************
/*!
  \file      src/IO/PDFWriter.C
  \author    J. Bakosi
  \date      Sun 18 Nov 2012 08:37:09 PM MST
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

PDFWriter::PDFWriter(Memory* memory, const string filename) :
  m_memory(memory), m_filename(filename)
//******************************************************************************
//  Constructor: Acquire PDF file handle
//! \param[in]  filename  File name to open for writing
//! \param[in]  memory    MEmory object pointer
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
  int xmin = numeric_limits<int>::max();
  int xmax = numeric_limits<int>::min();
  int ymin = numeric_limits<int>::max();
  int ymax = numeric_limits<int>::min();
  for (auto& p : *f) {
    if (p.first[0] < xmin) xmin = p.first[0];
    if (p.first[0] > xmax) xmax = p.first[0];
    if (p.first[1] < ymin) ymin = p.first[1];
    if (p.first[1] > ymax) ymax = p.first[1];
  }
  int nbix = xmax-xmin+1;
  int nbiy = ymax-ymin+1;

//   // Interpolate function values of joint PDF to points
//   JPDF::pdf ijpdf;
//   for (auto& p : *f) {
//     vector<int> s(2,0);
//   }

//   for (int i=0; i<nbix*nb; ++i) {
//     sjpdf[i+i/nbiy] += 
//   }

  // Output points
  m_outPDF << "$Nodes\n" << nbix*nbiy << endl;
  int k=0;
  for (int i=0; i<nbix; i++ )
    for (int j=0; j<nbiy; j++ ) {
      real x = xmin + i*binsize;
      real y = ymin + j*binsize;
      m_outPDF << k++ << " " << x << " " << y << " 0\n";
    }
  m_outPDF << "$EndNodes\n";

  --nbix;  --nbiy;
  m_outPDF << "$Elements\n" << nbix*nbiy << "\n";
  for (int i=0; i<nbix*nbiy; ++i) {
    m_outPDF << i << " 3 2 1 1 " << i+i/nbiy << " " << i+1+i/nbiy << " "
             << i+2+nbiy+i/nbiy << " " << i+1+nbiy+i/nbiy << "\n";
  }
  m_outPDF << "$EndElements\n";

  // Output function values
  m_outPDF << "$ElementData\n1\n\"Computed\"\n1\n0.0\n3\n0\n1\n" << f->size()
           << "\n";
  for (auto& p : *f) {
    m_outPDF << p.first[0]*nbiy + p.first[1] << " " << p.second/sp << endl;
  }
  m_outPDF << "$ElementData\n";
}
