//******************************************************************************
/*!
  \file      src/IO/PDFWriter.C
  \author    J. Bakosi
  \date      Thu 18 Sep 2014 07:17:02 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Univariate PDF writer
  \details   Univariate PDF writer
*/
//******************************************************************************

#include <PDFWriter.h>
#include <Exception.h>

using quinoa::PDFWriter;

void
PDFWriter::writeTxt( const UniPDF& pdf ) const
//******************************************************************************
//  Write out standardized univariate PDF to file
//! \param[in]  pdf  Univariate PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Output header
  m_outFile << "# Univariate PDF. Columns: x, probability\n#" << std::endl;

  const auto binsize = pdf.binsize();
  const auto sp = binsize * pdf.nsample();
  for (const auto& p : pdf.map())
    m_outFile << binsize * p.first << "\t" << p.second / sp << std::endl;
}

void
PDFWriter::writeTxt( const BiPDF& pdf ) const
//******************************************************************************
//  Write out standardized bivariate PDF to text file
//! \param[in]  pdf  Bivariate PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Output header
  m_outFile << "# Joint bivariate PDF. Columns: x, y, probability\n#\n"
            << "# Example visualization with gnuplot:\n#\n"
            << "# gnuplot> set dgrid3d 50,50,1\n"
            << "# gnuplot> set cntrparam levels 20\n"
            << "# gnuplot> set contour\n"
            << "# gnuplot> splot <filename.txt> with lines\n#" << std::endl;

  const auto binsize = pdf.binsize();
  const auto sp = binsize[0] * binsize[1] * pdf.nsample();

  // Output data
  for (const auto& p : pdf.map())
    m_outFile << binsize[0] * p.first[0] << " " << binsize[1] * p.first[1]
              << " " << p.second / sp << std::endl;
}

void
PDFWriter::writeGmsh( const BiPDF& pdf,
                      const std::string& pdfname,
                      ctr::PDFCenteringType centering ) const
//******************************************************************************
//  Write out standardized bivariate PDF to Gmsh (text) format
//! \param[in]  pdf  Bivariate PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Output mesh header: mesh version, file type, data size
  m_outFile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );

  // Get access to PDF underlying map and binsizes
  auto f = pdf.map();
  const auto binsize = pdf.binsize();
  const tk::real sp = binsize[0] * binsize[1] * pdf.nsample();

  // Find extents of sample space
  const auto ext = pdf.extents();
  std::size_t nbix = ext.second[0] - ext.first[0] + 1;
  std::size_t nbiy = ext.second[1] - ext.first[1] + 1;

  const tk::real xmin = binsize[0] * ext.first[0];
  const tk::real xmax = binsize[0] * ext.second[0];
  const tk::real ymin = binsize[1] * ext.first[1];
  const tk::real ymax = binsize[1] * ext.second[1];

  // Output grid points of discretized sample space (2D Cartesian grid)
  m_outFile << "$Nodes\n" << (nbix+1)*(nbiy+1) << std::endl;
  int k=0;
  for (int i=0; i<=nbix; i++) {
    tk::real x = xmin + i*binsize[0];
    for (int j=0; j<=nbiy; j++) {
      tk::real y = ymin + j*binsize[1];
      m_outFile << ++k << " " << x << " " << y << " 0\n";
    }
  }
  m_outFile << "$EndNodes\n";

  // Output elements of discretized sample space (2D Cartesian grid)
  m_outFile << "$Elements\n" << nbix*nbiy << "\n";
  for (int i=0; i<nbix*nbiy; ++i) {
    m_outFile << i+1 << " 3 2 1 1 " << i+i/nbiy+1 << " "
              << i+2+i/nbiy << " " << i+3+nbiy+i/nbiy << " "
              << i+2+nbiy+i/nbiy << std::endl;
  }
  m_outFile << "$EndElements\n";

  // Output PDF function values
  std::string c( "Element" );
  if (centering == ctr::PDFCenteringType::NODE) {
    ++nbix; ++nbiy;
    c = "Node";
  }
  std::vector< bool > out( nbix*nbiy, false ); // indicate bins filled
  m_outFile << "$" << c << "Data\n1\n\"" << pdfname << "\"\n1\n0.0\n3\n0\n1\n"
            << nbix*nbiy << "\n";
  for (const auto& p : f) {
    const auto bin = (p.first[0] - ext.first[0]) * nbiy +
                     (p.first[1] - ext.first[1]) % nbiy;
    Assert( bin < nbix*nbiy, "bin overflow in PDFWriter::writeGmsh()." );
    out[ bin ] = true;
    m_outFile << bin+1 << " " << p.second/sp << std::endl;
  }
  // Output bins nonexistent in PDF (looks better as zero than holes in gmsh)
  for (std::size_t i=0; i<out.size(); ++i) {
    if (!out[i]) m_outFile << i+1 << " 0" << std::endl;
  }
  m_outFile << "$End" << c << "Data\n";

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}
