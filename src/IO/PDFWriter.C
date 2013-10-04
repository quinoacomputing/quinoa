//******************************************************************************
/*!
  \file      src/IO/PDFWriter.C
  \author    J. Bakosi
  \date      Thu 03 Oct 2013 08:48:52 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Univariate PDF writer
  \details   Univariate PDF writer
*/
//******************************************************************************

#include <PDFWriter.h>
#include <Exception.h>

using namespace quinoa;

void
PDFWriter::write(const PDF& pdf)
//******************************************************************************
//  Write out standardized PDF to file
//! \param[in]  pdf  Univariate PDF
//! \author  J. Bakosi
//******************************************************************************
{
  auto f = pdf.getMap();
  const real binsize = pdf.getBinsize();
  const real sp = pdf.getNsample()*binsize;
  for (auto& p : *f) {
    m_outFile << p.first*binsize << "\t" << p.second/sp << std::endl;
  }
}

void
PDFWriter::writeTxt(const JPDF& jpdf)
//******************************************************************************
//  Write out standardized joint PDF to text file (only for 2D)
//! \param[in]  jpdf  Joint PDF
//! \author  J. Bakosi
//******************************************************************************
{
  auto f = jpdf.getMap();
  real binsize = jpdf.getBinsize();
  const real sp = jpdf.getNsample()*binsize*binsize;
  for (auto& p : *f) {
    m_outFile << p.first[0]*binsize << " " << p.first[1]*binsize
              << " " << p.second/sp << std::endl;
  }
}

void
PDFWriter::writeGmsh(const JPDF& jpdf)
//******************************************************************************
//  Write out standardized joint PDF to Gmsh (text) format (only for 2D)
//! \param[in]  jpdf  Joint PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Output mesh header: mesh version, file type, data size
  m_outFile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  ErrChk(!m_outFile.bad(), ExceptType::FATAL,
         "Failed to write to file: " + m_filename);

  auto f = jpdf.getMap();
  real binsize = jpdf.getBinsize();
  const real sp = jpdf.getNsample()*binsize*binsize;

  // Find sample space extents
  real xmin = std::numeric_limits<real>::max();
  real xmax = std::numeric_limits<real>::min();
  real ymin = std::numeric_limits<real>::max();
  real ymax = std::numeric_limits<real>::min();
  for (auto& p : *f) {
    if (binsize*p.first[0] < xmin) xmin = binsize*p.first[0];
    if (binsize*p.first[0] > xmax) xmax = binsize*p.first[0];
    if (binsize*p.first[1] < ymin) ymin = binsize*p.first[1];
    if (binsize*p.first[1] > ymax) ymax = binsize*p.first[1];
  }
  int nbix = static_cast<int>((xmax - xmin)/binsize + 1);
  int nbiy = static_cast<int>((ymax - ymin)/binsize + 1);

  // Output points of discretized sample space
  m_outFile << "$Nodes\n" << nbix*nbiy << std::endl;
  int k=0;
  for (int i=0; i<nbix; i++ ) {
    real x = xmin + i*binsize;
    for (int j=0; j<nbiy; j++ ) {
      real y = ymin + j*binsize;
      m_outFile << k++ << " " << x << " " << y << " 0\n";
    }
  }
  m_outFile << "$EndNodes\n";

  // Output elements of discretized sample space
  --nbix;  --nbiy;
  m_outFile << "$Elements\n" << nbix*nbiy << "\n";
  for (int i=0; i<nbix*nbiy; ++i) {
    m_outFile << i << " 3 2 1 1 " << i+i/nbiy << " " << i+1+i/nbiy << " "
             << i+2+nbiy+i/nbiy << " " << i+1+nbiy+i/nbiy << "\n";
  }
  m_outFile << "$EndElements\n";

  // Output function values
  ++nbiy;
  m_outFile << "$NodeData\n1\n\"Computed\"\n1\n0.0\n3\n0\n1\n"
           << f->size() << "\n";
  for (auto& p : *f) {
    int bin = static_cast<int>(
                (binsize*p.first[0] - xmin)/(xmax - xmin)*nbix*nbiy +
                (binsize*p.first[1] - ymin)/(ymax - ymin)*nbiy );
    m_outFile << bin << " " << p.second/sp << std::endl;
  }

  m_outFile << "$NodetData\n";
  ErrChk(!m_outFile.bad(), ExceptType::FATAL,
         "Failed to write to file: " + m_filename);
}
