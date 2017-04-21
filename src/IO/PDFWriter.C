// *****************************************************************************
/*!
  \file      src/IO/PDFWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Univariate PDF writer
  \brief     PDF writer class definition
  \details   This file defines a PDF writer class that facilitates outputing
    probability density functions (PDFs) into files in various formats using
    various configurations.
*/
// *****************************************************************************

#include <iostream>
#include <iomanip>

#include "NoWarning/exodusII.h"

#include "PDFWriter.h"
#include "Exception.h"

using tk::PDFWriter;

PDFWriter::PDFWriter( const std::string& filename,
                      ctr::TxtFloatFormatType format,
                      kw::precision::info::expect::type precision ) :
  Writer( filename )
// *****************************************************************************
//  Constructor
//! \param[in] filename Output filename to which output the PDF
//! \param[in] format Configure floating-point output format for ASCII output
//! \param[in] precision Configure precision for floating-point ASCII output
//! \author J. Bakosi
// *****************************************************************************
{
  // Set floating-point format for output file stream
  if (format == ctr::TxtFloatFormatType::DEFAULT)
    {} //m_outFile << std::defaultfloat;   GCC does not yet support this
  else if (format == ctr::TxtFloatFormatType::FIXED)
    m_outFile << std::fixed;
  else if (format == ctr::TxtFloatFormatType::SCIENTIFIC)
    m_outFile << std::scientific;
  else Throw( "Text floating-point format not recognized." );

  // Set numeric precision for output file stream if the input makes sense
  if (precision > 0 && precision < std::numeric_limits< tk::real >::digits10+2)
    m_outFile << std::setprecision( static_cast<int>(precision) );
}

void
PDFWriter::writeTxt( const UniPDF& pdf, const tk::ctr::PDFInfo& info ) const
// *****************************************************************************
//  Write out standardized univariate PDF to file
//! \param[in] pdf Univariate PDF
//! \param[in] info PDF metadata
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 1 >( vars );
  assertSampleSpaceExtents< 1 >( uext );

  // Query and optionally override number of bins and minimum of sample space if
  // user-specified extents were given and copy probabilities from pdf to an
  // array for output
  std::size_t nbi;
  tk::real min, max;
  std::vector< tk::real > outpdf;
  tk::real binsize;
  std::array< long, 2*UniPDF::dim > ext;
  extents( pdf, uext, nbi, min, max, binsize, ext, outpdf );

  // Output header
  m_outFile << "# vim: filetype=sh:\n#\n"
            << "# Univariate PDF: " << name << '(' << vars[0] << ')' << '\n'
            << "# -----------------------------------------------\n"
            << "# Numeric precision: " << m_outFile.precision() << '\n'
            << "# Bin size: " << binsize << '\n'
            << "# Number of bins estimated: " << ext[1] - ext[0] + 1
            << '\n'
            << "# Number of bins output: " << nbi << '\n'
            << "# Sample space extent: [" << min << " : " << max << "]\n#\n"
            << "# Example step-by-step visualization with gnuplot\n"
            << "# -----------------------------------------------\n"
            << "# gnuplot> set grid\n"
            << "# gnuplot> unset key\n"
            << "# gnuplot> set xlabel \"" << vars[0] << "\"\n"
            << "# gnuplot> set ylabel \"" << name << "(" << vars[0] << ")\"\n"
            << "# gnuplot> plot ";
  if (!uext.empty()) m_outFile << "[" << uext[0] << ':' << uext[1] << "] ";
  m_outFile << "\"" << m_filename << "\" with points\n#\n"
            << "# Gnuplot one-liner for quick copy-paste\n"
            << "# -----------------------------------------------\n"
            << "# set grid; unset key; set xlabel \"" << vars[0]
            << "\"; set ylabel \"" << name << "(" << vars[0]
            << ")\"; plot";
  if (!uext.empty()) m_outFile << " [" << uext[0] << ':' << uext[1] << "]";
  m_outFile << " \"" << m_filename << "\" w p\n#\n"
            << "# Data columns: " << vars[0] << ", " << name << "(" << vars[0]
            << ")\n# -----------------------------------------------\n";

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {
    for (const auto& p : pdf.map())
      m_outFile << binsize * static_cast<tk::real>(p.first) << '\t'
                << static_cast<tk::real>(p.second) / binsize /
                   static_cast<tk::real>(pdf.nsample())
                << std::endl;
  } else { // If user-specified sample space extents, output outpdf array
    std::size_t bin = 0;
    for (const auto& p : outpdf)
      m_outFile << binsize * static_cast<tk::real>(bin++) + uext[0] << '\t'
                << p << std::endl;
  }
}

void
PDFWriter::writeTxt( const BiPDF& pdf, const tk::ctr::PDFInfo& info ) const
// *****************************************************************************
//  Write out standardized bivariate PDF to text file
//! \param[in] pdf Bivariate PDF
//! \param[in] info PDF metadata
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 2 >( vars );
  assertSampleSpaceExtents< 2 >( uext );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 2D array for output
  std::size_t nbix, nbiy;
  tk::real xmin, xmax, ymin, ymax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 2 > binsize;
  std::array< long, 2*BiPDF::dim > ext;
  extents( pdf, uext, nbix, nbiy, xmin, xmax, ymin, ymax, binsize, ext, outpdf,
           ctr::PDFCenteringType::ELEM );

  // Output metadata
  m_outFile << "# vim: filetype=sh:\n#\n"
            << "# Joint bivariate PDF: " << name << '(' << vars[0] << ','
            << vars[1] << ")\n"
            << "# -----------------------------------------------\n"
            << "# Numeric precision: " << m_outFile.precision() << '\n'
            << "# Bin sizes: " << binsize[0] << ", " << binsize[1] << '\n'
            << "# Number of bins estimated: " << ext[1] - ext[0] + 1 << " x "
            << ext[3] - ext[2] + 1 << '\n'
            << "# Number of bins output: " << nbix << " x " << nbiy << '\n'
            << "# Sample space extents: [" << xmin << " : " << xmax
            << "], [" << ymin << " : " << ymax << "]\n#\n"
            << "# Example step-by-step visualization with gnuplot\n"
            << "# -----------------------------------------------\n"
            << "# gnuplot> set grid\n"
            << "# gnuplot> unset key\n"
            << "# gnuplot> set xlabel \"" << vars[0] << "\"\n"
            << "# gnuplot> set ylabel \"" << vars[1] << "\"\n"
            << "# gnuplot> set zlabel \"" << name << "(" << vars[0] << ","
            << vars[1] << ")\"\n"
            << "# gnuplot> set dgrid3d 50,50,1\n"
            << "# gnuplot> set cntrparam levels 20\n"
            << "# gnuplot> set contour\n";
  if (!uext.empty())
    m_outFile << "# gnuplot> set xrange [" << uext[0] << ':' << uext[1] << "]\n"
              << "# gnuplot> set yrange [" << uext[2] << ':' << uext[3] << "]\n";
         
  m_outFile << "# gnuplot> splot \"" << m_filename << "\" with lines\n#\n"
            << "# Gnuplot one-liner for quick copy-paste\n"
            << "# --------------------------------------\n"
            << "# set grid; unset key; set xlabel \"" << vars[0]
            << "\"; set ylabel \"" << vars[1] << "\"; set zlabel \"" << name
            << "(" << vars[0] << ',' << vars[1] << ")\"; set dgrid3d 50,50,1; "
               "set cntrparam levels 20; set contour; ";
  if (!uext.empty())
    m_outFile << "set xrange [" << uext[0] << ':' << uext[1] << "]; set yrange "
                 "[" << uext[2] << ':' << uext[3] << "]; ";
  m_outFile << "splot \"" << m_filename << "\" w l\n#\n"
            << "# Data columns: " << vars[0] << ", " << vars[1] << ", "
            << name << '(' << vars[0] << ',' << vars[1] << ")\n"
            << "# -----------------------------------------------\n";

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {
    for (const auto& p : pdf.map())
      m_outFile << binsize[0] * static_cast<tk::real>(p.first[0]) << '\t'
                << binsize[1] * static_cast<tk::real>(p.first[1]) << '\t'
                << p.second / binsize[0] / binsize[1] /
                   static_cast<tk::real>(pdf.nsample())
                << std::endl;
  } else { // If user-specified sample space extents, output outpdf array
    std::size_t bin = 0;
    for (const auto& p : outpdf) {
      m_outFile << binsize[0] * static_cast<tk::real>(bin % nbix) + uext[0]
                << '\t'
                << binsize[1] * static_cast<tk::real>(bin / nbix) + uext[2]
                << '\t'
                << p
                << std::endl;
      ++bin;
    }
  }

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}

void
PDFWriter::writeTxt( const TriPDF& pdf, const tk::ctr::PDFInfo& info ) const
// *****************************************************************************
//  Write out standardized trivariate PDF to text file
//! \param[in] pdf Trivariate PDF
//! \param[in] info PDF metadata
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 3 >( vars );
  assertSampleSpaceExtents< 3 >( uext );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 3D array for output
  std::size_t nbix, nbiy, nbiz;
  tk::real xmin, xmax, ymin, ymax, zmin, zmax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 3 > binsize;
  std::array< long, 2*TriPDF::dim > ext;
  extents( pdf, uext, nbix, nbiy, nbiz, xmin, xmax, ymin, ymax, zmin, zmax,
           binsize, ext, outpdf, ctr::PDFCenteringType::ELEM );

  // Output header
  m_outFile << "# vim: filetype=sh:\n#\n"
            << "# Joint trivariate PDF: " << name << '(' << vars[0] << ','
            << vars[1] << ',' << vars[2] << ")\n"
            << "# -----------------------------------------------\n"
            << "# Numeric precision: " << m_outFile.precision() << '\n'
            << "# Bin sizes: " << binsize[0] << ", " << binsize[1] << ", "
            << binsize[2] << '\n'
            << "# Number of bins estimated: " << ext[1] - ext[0] + 1 << " x "
            << ext[3] - ext[2] + 1 << " x " << ext[5] - ext[4] + 1 << '\n'
            << "# Number of bins output: " << nbix << " x " << nbiy << " x "
            << nbiz << '\n'
            << "# Sample space extents: [" << xmin << " : " << xmax << "], ["
            << ymin << " : " << ymax << "], [" << zmin << " : " << zmax
            << "]\n#\n"
            << "# Example step-by-step visualization with gnuplot\n"
            << "# -----------------------------------------------\n"
            << "# gnuplot> set grid\n"
            << "# gnuplot> set xlabel \"" << vars[0] << "\"\n"
            << "# gnuplot> set ylabel \"" << vars[1] << "\"\n"
            << "# gnuplot> set zlabel \"" << vars[2] << "\"\n";
  if (!uext.empty())
    m_outFile << "# gnuplot> set xrange [" << uext[0] << ':' << uext[1] << "]\n"
              << "# gnuplot> set yrange [" << uext[2] << ':' << uext[3] << "]\n"
              << "# gnuplot> set zrange [" << uext[4] << ':' << uext[5] << "]\n";
  m_outFile << "# gnuplot> splot \"" << m_filename << "\" pointtype 7 "
               "linecolor palette title \"" << name << '(' << vars[0] << ','
            << vars[1] << ',' << vars[2] << ")\"\n#\n"
            << "# Gnuplot one-liner for quick copy-paste\n"
            << "# --------------------------------------\n"
            << "# set grid; set xlabel \"" << vars[0] << "\"; set ylabel \""
            << vars[1] << "\"; set zlabel \"" << vars[2] << "\"; ";
  if (!uext.empty())
    m_outFile << "set xrange [" << uext[0] << ':' << uext[1] << "]; set yrange "
                 "[" << uext[2] << ':' << uext[3] << "]; set zrange ["
              << uext[4] << ':' << uext[5] << "]; ";
  m_outFile << "splot \"" << m_filename << "\" pt 7 linecolor palette title \""
            << name << '(' << vars[0] << ',' << vars[1] << ',' << vars[2] << ')'
            << "\"\n#\n"
            << "# Data columns: " << vars[0] << ", " << vars[1] << ", "
            << vars[2] << ", " << name << '(' << vars[0] << ',' << vars[1]
            << ',' << vars[2] << ")\n"
            << "# -----------------------------------------------\n";

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {
    for (const auto& p : pdf.map())
      m_outFile << binsize[0] * static_cast<tk::real>(p.first[0]) << '\t'
                << binsize[1] * static_cast<tk::real>(p.first[1]) << '\t'
                << binsize[2] * static_cast<tk::real>(p.first[2]) << '\t'
                << p.second / binsize[0] / binsize[1] / binsize[2]
                            / static_cast<tk::real>(pdf.nsample())
                << std::endl;
  } else { // If user-specified sample space extents, output outpdf array
    std::size_t bin = 0;
    const auto n = nbix*nbiy;
    for (const auto& p : outpdf) {
      m_outFile << binsize[0] * static_cast<tk::real>(bin % n % nbix) + uext[0]
                << '\t'
                << binsize[1] * static_cast<tk::real>(bin % n / nbix) + uext[2]
                << '\t'
                << binsize[2] * static_cast<tk::real>(bin / n) + uext[4] << '\t'
                << p
                << std::endl;
      ++bin;
    }
  }

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}

void
PDFWriter::writeGmshTxt( const BiPDF& pdf,
                         const tk::ctr::PDFInfo& info,
                         ctr::PDFCenteringType centering ) const
// *****************************************************************************
//  Write out standardized bivariate PDF to Gmsh (text) format
//! \param[in] pdf Bivariate PDF
//! \param[in] info PDF metadata
//! \param[in] centering Bin centering on sample space mesh
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 2 >( vars );
  assertSampleSpaceExtents< 2 >( uext );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 2D array for output
  std::size_t nbix, nbiy;
  tk::real xmin, xmax, ymin, ymax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 2 > binsize;
  std::array< long, 2*BiPDF::dim > ext;
  extents( pdf, uext, nbix, nbiy, xmin, xmax, ymin, ymax, binsize, ext, outpdf,
           centering );

  // Output metadata. The #s are unnecessary, but vi will color it differently.
  m_outFile << "$Comments\n"
            << "# vim: filetype=sh:\n"
            << "# Joint bivariate PDF: " << name << '(' << vars[0] << ','
            << vars[1] << ")\n"
            << "# -----------------------------------------------\n"
            << "# Numeric precision: " << m_outFile.precision() << '\n'
            << "# Bin sizes: " << binsize[0] << ", " << binsize[1] << '\n'
            << "# Number of bins estimated: " << ext[1] - ext[0] + 1 << " x "
            << ext[3] - ext[2] + 1 << '\n'
            << "# Number of bins output: " << nbix << " x " << nbiy << '\n'
            << "# Sample space extents: [" << xmin << " : " << xmax
            << "], [" << ymin << " : " << ymax << "]\n"
            << "$EndComments\n";

  // Output mesh header: mesh version, file type, data size
  m_outFile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );

  // Output grid points of discretized sample space (2D Cartesian grid)
  m_outFile << "$Nodes\n" << (nbix+1)*(nbiy+1) << std::endl;
  int k=0;
  for (std::size_t i=0; i<=nbiy; i++) {
    tk::real y = ymin + i*binsize[1];
    for (std::size_t j=0; j<=nbix; j++) {
      tk::real x = xmin + j*binsize[0];
      m_outFile << ++k << ' ' << x << ' ' << y << " 0\n";
    }
  }
  m_outFile << "$EndNodes\n";

  // Output elements of discretized sample space (2D Cartesian grid)
  m_outFile << "$Elements\n" << nbix*nbiy << "\n";
  for (std::size_t i=0; i<nbix*nbiy; ++i) {
    const auto y = i/nbix;
    m_outFile << i+1 << " 3 2 1 1 " << i+y+1 << ' ' << i+y+2 << ' '
              << i+y+nbix+3 << ' ' << i+y+nbix+2 << std::endl;
  }
  m_outFile << "$EndElements\n";

  // Output PDF function values in element or node centers
  std::string c( "Element" );
  if (centering == ctr::PDFCenteringType::NODE) {
    ++nbix; ++nbiy;
    c = "Node";
  }
  m_outFile << '$' << c << "Data\n1\n\"" << name << "\"\n1\n0.0\n3\n0\n1\n"
            << nbix*nbiy << "\n";

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {

    std::vector< bool > out( nbix*nbiy, false ); // indicate bins filled
    for (const auto& p : pdf.map()) {
      const auto bin = (p.first[1] - ext[2]) * static_cast<long>(nbix) +
                       (p.first[0] - ext[0]) % static_cast<long>(nbix);
      Assert( bin >= 0, "Bin underflow in PDFWriter::writeGmshTxt()." );
      Assert( static_cast<std::size_t>(bin) < nbix*nbiy,
              "Bin overflow in PDFWriter::writeGmshTxt()." );
      out[ static_cast<std::size_t>(bin) ] = true;
      m_outFile << bin+1 << '\t'
                << p.second / binsize[0] / binsize[1]
                            / static_cast<tk::real>(pdf.nsample())
                << std::endl;
    }
    // Output bins nonexistent in PDF (gmsh sometimes fails to plot the exiting
    // bins if holes exist in the data, it also looks better as zero than holes)
    for (std::size_t i=0; i<out.size(); ++i)
      if (!out[i]) m_outFile << i+1 << "\t0" << std::endl;

  } else { // If user-specified sample space extents, output outpdf array

    std::size_t bin = 0;
    for (const auto& p : outpdf) m_outFile << ++bin << ' ' << p << std::endl;

  }

  m_outFile << "$End" << c << "Data\n";

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}

void
PDFWriter::writeGmshTxt( const TriPDF& pdf,
                         const tk::ctr::PDFInfo& info,
                         ctr::PDFCenteringType centering ) const
// *****************************************************************************
//  Write out standardized trivariate PDF to Gmsh (text) format
//! \param[in] pdf Trivariate PDF
//! \param[in] info PDF metadata
//! \param[in] centering Bin centering on sample space mesh
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 3 >( vars );
  assertSampleSpaceExtents< 3 >( uext );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 3D array for output
  std::size_t nbix, nbiy, nbiz;
  tk::real xmin, xmax, ymin, ymax, zmin, zmax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 3 > binsize;
  std::array< long, 2*TriPDF::dim > ext;
  extents( pdf, uext, nbix, nbiy, nbiz, xmin, xmax, ymin, ymax, zmin, zmax,
           binsize, ext, outpdf, centering );

  // Output metadata. The #s are unnecessary, but vi will color it differently.
  m_outFile << "$Comments\n"
            << "# vim: filetype=sh:\n#\n"
            << "# Joint trivariate PDF: " << name << '(' << vars[0] << ','
            << vars[1] << ',' << vars[2] << ")\n"
            << "# -----------------------------------------------\n"
            << "# Numeric precision: " << m_outFile.precision() << '\n'
            << "# Bin sizes: " << binsize[0] << ", " << binsize[1] << ", "
            << binsize[2] << '\n'
            << "# Number of bins estimated: " << ext[1] - ext[0] + 1 << " x "
            << ext[3] - ext[2] + 1 << " x " << ext[5] - ext[4] + 1 << '\n'
            << "# Number of bins output: " << nbix << " x " << nbiy << " x "
            << nbiz << '\n'
            << "# Sample space extents: [" << xmin << " : " << xmax << "], ["
            << ymin << " : " << ymax << "], [" << zmin << " : " << zmax << "]\n"
            << "$EndComments\n";

  // Output mesh header: mesh version, file type, data size
  m_outFile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );

  // Output grid points of discretized sample space (3D Cartesian grid)
  m_outFile << "$Nodes\n" << (nbix+1)*(nbiy+1)*(nbiz+1) << std::endl;
  int l=0;
  for (std::size_t k=0; k<=nbiz; k++) {
    tk::real z = zmin + k*binsize[2];
    for (std::size_t j=0; j<=nbiy; j++) {
      tk::real y = ymin + j*binsize[1];
      for (std::size_t i=0; i<=nbix; i++) {
        tk::real x = xmin + i*binsize[0];
        m_outFile << ++l << ' ' << x << ' ' << y << ' ' << z << '\n';
      }
    }
  }
  m_outFile << "$EndNodes\n";

  // Output elements of discretized sample space (3D Cartesian grid)
  m_outFile << "$Elements\n" << nbix*nbiy*nbiz << "\n";
  const auto n = nbix*nbiy;
  const auto p = (nbix+1)*(nbiy+1);
  for (std::size_t i=0; i<nbix*nbiy*nbiz; ++i) {
    const auto y = i/nbix + i/n*(nbix+1);
    m_outFile << i+1 << " 5 2 1 1 " << i+y+1 << ' ' << i+y+2 << ' '
              << i+y+nbix+3 << ' ' << i+y+nbix+2 << ' '
              << i+y+p+1 << ' ' << i+y+p+2 << ' '
              << i+y+p+nbix+3 << ' ' << i+y+p+nbix+2 << ' '
              << std::endl;
  }
  m_outFile << "$EndElements\n";

  // Output PDF function values in element or node centers
  std::string c( "Element" );
  if (centering == ctr::PDFCenteringType::NODE) {
    ++nbix; ++nbiy; ++nbiz;
    c = "Node";
  }
  m_outFile << '$' << c << "Data\n1\n\"" << name << "\"\n1\n0.0\n3\n0\n1\n"
            << nbix*nbiy*nbiz << "\n";

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {

    std::vector< bool > out( nbix*nbiy*nbiz, false ); // indicate bins filled
    for (const auto& q : pdf.map()) {
      const auto bin = (q.first[2] - ext[4]) * static_cast<long>(nbix*nbiy) +
                       (q.first[1] - ext[2]) * static_cast<long>(nbix) +
                       (q.first[0] - ext[0]) % static_cast<long>(nbix);
      Assert( bin >= 0, "Bin underflow in PDFWriter::writeGmshTxt()." );
      Assert( static_cast<std::size_t>(bin) < nbix*nbiy*nbiz,
              "Bin overflow in PDFWriter::writeGmshTxt()." );
      out[ static_cast<std::size_t>(bin) ] = true;
      m_outFile << bin+1 << '\t'
                << q.second / binsize[0] / binsize[1] / binsize[2]
                            / static_cast<tk::real>(pdf.nsample())
                << std::endl;
    }
    // Output bins nonexistent in PDF (gmsh sometimes fails to plot the exiting
    // bins if holes exist in the data, it also looks better as zero than holes)
    for (std::size_t i=0; i<out.size(); ++i)
      if (!out[i]) m_outFile << i+1 << "\t0" << std::endl;

  } else { // If user-specified sample space extents, output outpdf array

    std::size_t bin = 0;
    for (const auto& q : outpdf) m_outFile << ++bin << ' ' << q << std::endl;

  }

  m_outFile << "$End" << c << "Data\n";

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}

void
PDFWriter::writeGmshBin( const BiPDF& pdf,
                         const tk::ctr::PDFInfo& info,
                         ctr::PDFCenteringType centering ) const
// *****************************************************************************
//  Write out standardized bivariate PDF to Gmsh (binary) format
//! \param[in] pdf Bivariate PDF
//! \param[in] info PDF metadata
//! \param[in] centering Bin centering on sample space mesh
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 2 >( vars );
  assertSampleSpaceExtents< 2 >( uext );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 2D array for output
  std::size_t nbix, nbiy;
  tk::real xmin, xmax, ymin, ymax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 2 > binsize;
  std::array< long, 2*BiPDF::dim > ext;
  extents( pdf, uext, nbix, nbiy, xmin, xmax, ymin, ymax, binsize, ext, outpdf,
           centering );

  // Output metadata. The #s are unnecessary, but vi will color it differently.
  m_outFile << "$Comments\n"
            << "# vim: filetype=sh:\n"
            << "# Joint bivariate PDF: " << name << '(' << vars[0] << ','
            << vars[1] << ")\n"
            << "# -----------------------------------------------\n"
            << "# Numeric precision: 64-bit binary\n"
            << "# Bin sizes: " << binsize[0] << ", " << binsize[1] << '\n'
            << "# Number of bins estimated: " << ext[1] - ext[0] + 1 << " x "
            << ext[3] - ext[2] + 1 << '\n'
            << "# Number of bins output: " << nbix << " x " << nbiy << '\n'
            << "# Sample space extents: [" << xmin << " : " << xmax
            << "], [" << ymin << " : " << ymax << "]\n"
            << "$EndComments\n";

  // Output mesh header: mesh version, file type, data size
  m_outFile << "$MeshFormat\n2.2 1 8\n";
  int one = 1;
  m_outFile.write( reinterpret_cast<char*>(&one), sizeof(int) );
  m_outFile << "\n$EndMeshFormat\n";
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );

  // Output grid points of discretized sample space (2D Cartesian grid)
  m_outFile << "$Nodes\n" << (nbix+1)*(nbiy+1) << std::endl;
  int k = 0;
  tk::real z = 0.0;
  for (std::size_t i=0; i<=nbiy; i++) {
    tk::real y = ymin + i*binsize[1];
    for (std::size_t j=0; j<=nbix; j++) {
      tk::real x = xmin + j*binsize[0];
      ++k;
      m_outFile.write( reinterpret_cast< char* >( &k ), sizeof(int) );
      m_outFile.write( reinterpret_cast< char* >( &x ), sizeof(tk::real) );
      m_outFile.write( reinterpret_cast< char* >( &y ), sizeof(tk::real) );
      m_outFile.write( reinterpret_cast< char* >( &z ), sizeof(tk::real) );
    }
  }
  m_outFile << "\n$EndNodes\n";

  // Output elements of discretized sample space (2D Cartesian grid)
  m_outFile << "$Elements\n" << nbix*nbiy << "\n";
  int type = 3;                 // gmsh elem type: 4-node quadrangle
  std::size_t n = nbix*nbiy;    // number of elements in (this single) block
  int ntags = 2;                // number of element tags
  m_outFile.write( reinterpret_cast< char* >( &type ), sizeof(int) );
  m_outFile.write( reinterpret_cast< char* >( &n ), sizeof(int) );
  m_outFile.write( reinterpret_cast< char* >( &ntags ), sizeof(int) );
  for (std::size_t i=0; i<n; ++i) {
    const auto y = i/nbix;
    auto id = i+1;
    int tag[2] = { 1, 1 };
    int con[4] = { static_cast< int >( i+y+1 ),
                   static_cast< int >( i+y+2 ),
                   static_cast< int >( i+y+nbix+3 ),
                   static_cast< int >( i+y+nbix+2 ) };
    m_outFile.write( reinterpret_cast< char* >( &id ), sizeof(int) );
    m_outFile.write( reinterpret_cast< char* >( tag ), 2*sizeof(int) );
    m_outFile.write( reinterpret_cast< char* >( con ), 4*sizeof(int) );
  }
  m_outFile << "\n$EndElements\n";

  // Output PDF function values in element or node centers
  std::string c( "Element" );
  if (centering == ctr::PDFCenteringType::NODE) {
    ++nbix; ++nbiy;
    c = "Node";
  }
  m_outFile << '$' << c << "Data\n1\n\"" << name << "\"\n1\n0.0\n3\n0\n1\n"
            << nbix*nbiy << "\n";

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {

    std::vector< bool > out( nbix*nbiy, false ); // indicate bins filled
    for (const auto& p : pdf.map()) {
      const auto bin = (p.first[1] - ext[2]) * static_cast<long>(nbix) +
                       (p.first[0] - ext[0]) % static_cast<long>(nbix);
      Assert( bin >= 0, "Bin underflow in PDFWriter::writeGmshBin()." );
      Assert( static_cast<std::size_t>(bin) < nbix*nbiy,
              "Bin overflow in PDFWriter::writeGmshBin()." );
      out[ static_cast<std::size_t>(bin) ] = true;
      auto id = static_cast<int>(bin+1);
      tk::real prob = p.second / binsize[0] / binsize[1]
                               / static_cast<tk::real>(pdf.nsample());
      m_outFile.write( reinterpret_cast< char* >( &id ), sizeof(int) );
      m_outFile.write( reinterpret_cast< char* >( &prob ), sizeof(tk::real) );
    }
    // Output bins nonexistent in PDF (gmsh sometimes fails to plot the exiting
    // bins if holes exist in the data, it also looks better as zero than holes)
    tk::real prob = 0.0;
    for (std::size_t i=0; i<out.size(); ++i)
      if (!out[i]) {
        auto id = static_cast<int>(i+1);
        m_outFile.write( reinterpret_cast< char* >( &id ), sizeof(int) );
        m_outFile.write( reinterpret_cast< char* >( &prob ), sizeof(tk::real) );
      }

  } else { // If user-specified sample space extents, output outpdf array

    int bin = 0;
    for (auto& p : outpdf) {
      ++bin;
      m_outFile.write( reinterpret_cast< char* >( &bin ), sizeof(int) );
      m_outFile.write( reinterpret_cast< char* >( &p ), sizeof(tk::real) );
    }

  }

  m_outFile << "$End" << c << "Data\n";

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}

void
PDFWriter::writeGmshBin( const TriPDF& pdf,
                         const tk::ctr::PDFInfo& info,
                         ctr::PDFCenteringType centering ) const
// *****************************************************************************
//  Write out standardized trivariate PDF to Gmsh (binary) format
//! \param[in] pdf Trivariate PDF
//! \param[in] info PDF metadata
//! \param[in] centering Bin centering on sample space mesh
//! \author  J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 3 >( vars );
  assertSampleSpaceExtents< 3 >( uext );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 3D array for output
  std::size_t nbix, nbiy, nbiz;
  tk::real xmin, xmax, ymin, ymax, zmin, zmax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 3 > binsize;
  std::array< long, 2*TriPDF::dim > ext;
  extents( pdf, uext, nbix, nbiy, nbiz, xmin, xmax, ymin, ymax, zmin, zmax,
           binsize, ext, outpdf, centering );

  // Output metadata. The #s are unnecessary, but vi will color it differently.
  m_outFile << "$Comments\n"
            << "# vim: filetype=sh:\n#\n"
            << "# Joint trivariate PDF: " << name << '(' << vars[0] << ','
            << vars[1] << ',' << vars[2] << ")\n"
            << "# -----------------------------------------------\n"
            << "# Numeric precision: 64-bit binary\n"
            << "# Bin sizes: " << binsize[0] << ", " << binsize[1] << ", "
            << binsize[2] << '\n'
            << "# Number of bins estimated: " << ext[1] - ext[0] + 1 << " x "
            << ext[3] - ext[2] + 1 << " x " << ext[5] - ext[4] + 1 << '\n'
            << "# Number of bins output: " << nbix << " x " << nbiy << " x "
            << nbiz << '\n'
            << "# Sample space extents: [" << xmin << " : " << xmax << "], ["
            << ymin << " : " << ymax << "], [" << zmin << " : " << zmax << "]\n"
            << "$EndComments\n";

  // Output mesh header: mesh version, file type, data size
  m_outFile << "$MeshFormat\n2.2 1 8\n";
  int one = 1;
  m_outFile.write( reinterpret_cast<char*>(&one), sizeof(int) );
  m_outFile << "\n$EndMeshFormat\n";
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );

  // Output grid points of discretized sample space (3D Cartesian grid)
  m_outFile << "$Nodes\n" << (nbix+1)*(nbiy+1)*(nbiz+1) << std::endl;
  int l=0;
  for (std::size_t k=0; k<=nbiz; k++) {
    tk::real z = zmin + k*binsize[2];
    for (std::size_t j=0; j<=nbiy; j++) {
      tk::real y = ymin + j*binsize[1];
      for (std::size_t i=0; i<=nbix; i++) {
        tk::real x = xmin + i*binsize[0];
        ++l;
        m_outFile.write( reinterpret_cast< char* >( &l ), sizeof(int) );
        m_outFile.write( reinterpret_cast< char* >( &x ), sizeof(tk::real) );
        m_outFile.write( reinterpret_cast< char* >( &y ), sizeof(tk::real) );
        m_outFile.write( reinterpret_cast< char* >( &z ), sizeof(tk::real) );
      }
    }
  }
  m_outFile << "\n$EndNodes\n";

  // Output elements of discretized sample space (3D Cartesian grid)
  m_outFile << "$Elements\n" << nbix*nbiy*nbiz << "\n";
  int type = 5;                       // gmsh elem type: 8-node hexahedron
  std::size_t nelem = nbix*nbiy*nbiz; // num of elements in (this single) block
  int ntags = 2;                      // number of element tags
  m_outFile.write( reinterpret_cast< char* >( &type ), sizeof(int) );
  m_outFile.write( reinterpret_cast< char* >( &nelem ), sizeof(int) );
  m_outFile.write( reinterpret_cast< char* >( &ntags ), sizeof(int) );
  const auto n = nbix*nbiy;
  const auto p = (nbix+1)*(nbiy+1);
  for (std::size_t i=0; i<nelem; ++i) {
    const auto y = i/nbix + i/n*(nbix+1);
    auto id = i+1;
    int tag[2] = { 1, 1 };
    int con[8] = { static_cast< int >( i+y+1 ),
                   static_cast< int >( i+y+2 ),
                   static_cast< int >( i+y+nbix+3 ),
                   static_cast< int >( i+y+nbix+2 ),
                   static_cast< int >( i+y+p+1 ),
                   static_cast< int >( i+y+p+2 ),
                   static_cast< int >( i+y+p+nbix+3 ),
                   static_cast< int >( i+y+p+nbix+2 ) };
    m_outFile.write( reinterpret_cast< char* >( &id ), sizeof(int) );
    m_outFile.write( reinterpret_cast< char* >( tag ), 2*sizeof(int) );
    m_outFile.write( reinterpret_cast< char* >( con ), 8*sizeof(int) );
  }
  m_outFile << "\n$EndElements\n";

  // Output PDF function values in element or node centers
  std::string c( "Element" );
  if (centering == ctr::PDFCenteringType::NODE) {
    ++nbix; ++nbiy; ++nbiz;
    c = "Node";
  }
  m_outFile << '$' << c << "Data\n1\n\"" << name << "\"\n1\n0.0\n3\n0\n1\n"
            << nbix*nbiy*nbiz << "\n";

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {

    std::vector< bool > out( nbix*nbiy*nbiz, false ); // indicate bins filled
    for (const auto& q : pdf.map()) {
      const auto bin = (q.first[2] - ext[4]) * static_cast<long>(nbix*nbiy) +
                       (q.first[1] - ext[2]) * static_cast<long>(nbix) +
                       (q.first[0] - ext[0]) % static_cast<long>(nbix);
      Assert( bin >= 0, "Bin underflow in PDFWriter::writeGmshBin()." );
      Assert( static_cast<std::size_t>(bin) < nbix*nbiy*nbiz,
              "Bin overflow in PDFWriter::writeGmshBin()." );
      out[ static_cast<std::size_t>(bin) ] = true;
      auto id = static_cast<int>(bin+1);
      tk::real prob = q.second / binsize[0] / binsize[1] / binsize[2]
                               / static_cast<tk::real>(pdf.nsample());
      m_outFile.write( reinterpret_cast< char* >( &id ), sizeof(int) );
      m_outFile.write( reinterpret_cast< char* >( &prob ), sizeof(tk::real) );
    }
    // Output bins nonexistent in PDF (gmsh sometimes fails to plot the exiting
    // bins if holes exist in the data, it also looks better as zero than holes)
    tk::real prob = 0.0;
    for (std::size_t i=0; i<out.size(); ++i)
      if (!out[i]) {
        auto id = static_cast<int>(i+1);
        m_outFile.write( reinterpret_cast< char* >( &id ), sizeof(int) );
        m_outFile.write( reinterpret_cast< char* >( &prob ), sizeof(tk::real) );
      }

  } else { // If user-specified sample space extents, output outpdf array

    int bin = 0;
    for (auto& q : outpdf) {
      ++bin;
      m_outFile.write( reinterpret_cast< char* >( &bin ), sizeof(int) );
      m_outFile.write( reinterpret_cast< char* >( &q ), sizeof(tk::real) );
    }

  }

  m_outFile << "$End" << c << "Data\n";

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}

void
PDFWriter::writeExodusII( const BiPDF& pdf,
                          const tk::ctr::PDFInfo& info,
                          std::size_t it,
                          ctr::PDFCenteringType centering ) const
// *****************************************************************************
//  Write out standardized bivariate PDF to Exodus II format
//! \param[in] pdf Bivariate PDF
//! \param[in] info PDF metadata
//! \param[in] it Iteration count
//! \param[in] centering Bin centering on sample space mesh
//! \author  J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 2 >( vars );
  assertSampleSpaceExtents< 2 >( uext );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 2D array for output
  std::size_t nbix, nbiy;
  tk::real xmin, xmax, ymin, ymax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 2 > binsize;
  std::array< long, 2*BiPDF::dim > ext;
  extents( pdf, uext, nbix, nbiy, xmin, xmax, ymin, ymax, binsize, ext, outpdf,
           centering );

  // Create ExodusII file
  int outFile = createExFile();

  // Compute number of nodes and number of elements in sample space mesh
  std::size_t nelem = nbix*nbiy;
  std::size_t nnode = (nbix+1)*(nbiy+1);

  // Write ExodusII header
  writeExHdr( outFile, static_cast<int>(nnode), static_cast<int>(nelem) );

  // Write sample space variables as coordinate names
  char* coordnames[] = { const_cast< char* >( vars[0].c_str() ),
                         const_cast< char* >( vars[1].c_str() ),
                         const_cast< char* >( "probability" ) };
  ErrChk( ex_put_coord_names( outFile, coordnames ) == 0,
          "Failed to write coordinate names to file: " + m_filename );

  // Output grid points of discretized sample space (2D Cartesian grid)
  std::vector< tk::real > x( nnode, 0.0 );
  std::vector< tk::real > y( nnode, 0.0 );
  std::vector< tk::real > z( nnode, 0.0 );
  std::size_t k = 0;
  for (std::size_t i=0; i<=nbiy; i++)
    for (std::size_t j=0; j<=nbix; j++) {
      x[k] = xmin + j*binsize[0];
      y[k] = ymin + i*binsize[1];
      ++k;
    }
  ErrChk( ex_put_coord( outFile, x.data(), y.data(), z.data() ) == 0,
          "Failed to write coordinates to file: " + m_filename );

  // Output elements of discretized sample space (2D Cartesian grid)
  // Write element block information
  ErrChk( ex_put_block( outFile,
                        EX_ELEM_BLOCK,
                        1,
                        "QUADRANGLES",
                        static_cast<int64_t>(nelem),
                        4,
                        0,
                        0,
                        0 ) == 0,
          "Failed to write QUDARANGLE element block to file: " + m_filename );
  // Write element connectivity
  for (std::size_t i=0; i<nelem; ++i) {
    auto ye = i/nbix;
    int con[4] = { static_cast< int >( i+ye+1 ),
                   static_cast< int >( i+ye+2 ),
                   static_cast< int >( i+ye+nbix+3 ),
                   static_cast< int >( i+ye+nbix+2 ) };
    ErrChk( ex_put_partial_conn( outFile, EX_ELEM_BLOCK, 1,
              static_cast<int64_t>(i+1), 1, con, nullptr, nullptr ) == 0,
      "Failed to write element connectivity to file: " + m_filename );
  }

  // Output PDF function values in element or node centers
  ex_entity_type c = EX_ELEM_BLOCK;
  if (centering == ctr::PDFCenteringType::NODE) {
    ++nbix; ++nbiy;
    c = EX_NODE_BLOCK;
  }

  // Write PDF function values metadata
  ErrChk( ex_put_variable_param( outFile, c, 1 ) == 0,
            "Failed to write results metadata to file: " + m_filename );
  char* probname[1];
  std::string pdfname( name + '(' + vars[0] + ',' + vars[1] + ')' );
  probname[0] = const_cast< char* >( pdfname.c_str() );
  ErrChk( ex_put_variable_names( outFile, c, 1, probname ) == 0,
            "Failed to write results metadata to file: " + m_filename );

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {

    // Output PDF function values in element centers
    std::vector< tk::real > prob( nbix*nbiy, 0.0 );
    for (const auto& p : pdf.map()) {
      const auto bin = (p.first[1] - ext[2]) * static_cast<long>(nbix) +
                       (p.first[0] - ext[0]) % static_cast<long>(nbix);
      Assert( bin >= 0, "Bin underflow in PDFWriter::writeExodusII()." );
      Assert( static_cast<std::size_t>(bin) < nbix*nbiy,
              "Bin overflow in PDFWriter::writeExodusII()." );
      prob[ static_cast<std::size_t>(bin) ] =
        p.second / binsize[0] / binsize[1]
                 / static_cast<tk::real>(pdf.nsample());
    }
    writeExVar( outFile, it+1, centering, prob );

  } else { // If user-specified sample space extents, output outpdf array

    writeExVar( outFile, it+1, centering, outpdf );

  }

  ErrChk( ex_close(outFile) == 0, "Failed to close file: " + m_filename );
}

void
PDFWriter::writeExodusII( const TriPDF& pdf,
                          const tk::ctr::PDFInfo& info,
                          std::size_t it,
                          ctr::PDFCenteringType centering ) const
// *****************************************************************************
//  Write out standardized trivariate PDF to Exodus II format
//! \param[in] pdf Trivariate PDF
//! \param[in] info PDF metadata
//! \param[in] it Iteration count
//! \param[in] centering Bin centering on sample space mesh
//! \author  J. Bakosi
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;

  assertSampleSpaceDimensions< 3 >( vars );
  assertSampleSpaceExtents< 3 >( uext );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 3D array for output
  std::size_t nbix, nbiy, nbiz;
  tk::real xmin, xmax, ymin, ymax, zmin, zmax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 3 > binsize;
  std::array< long, 2*TriPDF::dim > ext;
  extents( pdf, uext, nbix, nbiy, nbiz, xmin, xmax, ymin, ymax, zmin, zmax,
           binsize, ext, outpdf, centering );

  // Create ExodusII file
  int outFile = createExFile();

  // Compute number of nodes and number of elements in sample space mesh
  std::size_t nelem = nbix*nbiy*nbiz;
  std::size_t nnode = (nbix+1)*(nbiy+1)*(nbiz+1);

  // Write ExodusII header
  writeExHdr( outFile, static_cast<int>(nnode), static_cast<int>(nelem) );

  // Write sample space variables as coordinate names
  char* coordnames[] = { const_cast< char* >( vars[0].c_str() ),
                         const_cast< char* >( vars[1].c_str() ),
                         const_cast< char* >( vars[2].c_str() ) };
  ErrChk( ex_put_coord_names( outFile, coordnames ) == 0,
          "Failed to write coordinate names to file: " + m_filename );

  // Output grid points of discretized sample space (2D Cartesian grid)
  std::vector< tk::real > x( nnode, 0.0 );
  std::vector< tk::real > y( nnode, 0.0 );
  std::vector< tk::real > z( nnode, 0.0 );
  std::size_t l=0;
  for (std::size_t k=0; k<=nbiz; k++)
    for (std::size_t i=0; i<=nbiy; i++)
      for (std::size_t j=0; j<=nbix; j++) {
        x[l] = xmin + static_cast<tk::real>(j) * binsize[0];
        y[l] = ymin + static_cast<tk::real>(i) * binsize[1];
        z[l] = zmin + static_cast<tk::real>(k) * binsize[2];
        ++l;
      }
  ErrChk( ex_put_coord( outFile, x.data(), y.data(), z.data() ) == 0,
          "Failed to write coordinates to file: " + m_filename );

  // Output elements of discretized sample space (2D Cartesian grid)
  // Write element block information
  ErrChk( ex_put_block( outFile,
                        EX_ELEM_BLOCK,
                        1,
                        "HEXAHEDRA",
                        static_cast<int64_t>(nelem),
                        8,
                        0,
                        0,
                        0 ) == 0,
          "Failed to write HEXAHEDRA element block to file: " + m_filename );
  // Write element connectivity
  const auto n = nbix*nbiy;
  const auto p = (nbix+1)*(nbiy+1);
  for (std::size_t i=0; i<nelem; ++i) {
    const auto ye = i/nbix + i/n*(nbix+1);
    int con[8] = { static_cast< int >( i+ye+1 ),
                   static_cast< int >( i+ye+2 ),
                   static_cast< int >( i+ye+nbix+3 ),
                   static_cast< int >( i+ye+nbix+2 ),
                   static_cast< int >( i+ye+p+1 ),
                   static_cast< int >( i+ye+p+2 ),
                   static_cast< int >( i+ye+p+nbix+3 ),
                   static_cast< int >( i+ye+p+nbix+2 ) };
    ErrChk(
      ex_put_partial_conn( outFile, EX_ELEM_BLOCK, 1,
        static_cast<int64_t>(i+1), 1, con, nullptr, nullptr ) == 0,
      "Failed to write element connectivity to file: " + m_filename );
  }

  // Output PDF function values in element or node centers
  ex_entity_type c = EX_ELEM_BLOCK;
  if (centering == ctr::PDFCenteringType::NODE) {
    ++nbix; ++nbiy; ++nbiz;
    c = EX_NODE_BLOCK;
  }

  // Write PDF function values metadata
  ErrChk( ex_put_variable_param( outFile, c, 1 ) == 0,
            "Failed to write results metadata to file: " + m_filename );
  char* probname[1];
  std::string pdfname( name + '(' + vars[0] + ',' +
                         vars[1] + ',' + vars[2] + ')' );
  probname[0] = const_cast< char* >( pdfname.c_str() );
  ErrChk( ex_put_variable_names( outFile, c, 1, probname ) == 0,
            "Failed to write results metadata to file: " + m_filename );

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {

    // Output PDF function values in element centers
    std::vector< tk::real > prob( nbix*nbiy*nbiz, 0.0 );
    for (const auto& q : pdf.map()) {
      const auto bin = (q.first[2] - ext[4]) * static_cast<long>(nbix*nbiy) +
                       (q.first[1] - ext[2]) * static_cast<long>(nbix) +
                       (q.first[0] - ext[0]) % static_cast<long>(nbix);
      Assert( bin >= 0, "Bin underflow in PDFWriter::writeExodusII()." );
      Assert( static_cast<std::size_t>(bin) < nbix*nbiy*nbiz,
              "Bin overflow in PDFWriter::writeExodusII()." );
      prob[ static_cast<std::size_t>(bin) ] =
        q.second / binsize[0] / binsize[1] / binsize[2]
                 / static_cast<tk::real>(pdf.nsample());
    }
    writeExVar( outFile, it+1, centering, prob );

  } else { // If user-specified sample space extents, output outpdf array

    writeExVar( outFile, it+1, centering, outpdf );

  }

  ErrChk( ex_close(outFile) == 0, "Failed to close file: " + m_filename );
}

int
PDFWriter::createExFile() const
// *****************************************************************************
//  Create Exodus II file
//! \return ExodusII file handle
//! \author J. Bakosi
// *****************************************************************************
{
  int cpuwordsize = sizeof( double );
  int iowordsize = sizeof( double );
  int outFileId = ex_create( m_filename.c_str(),
                             EX_CLOBBER | EX_NORMAL_MODEL,
                             &cpuwordsize,
                             &iowordsize );
  ErrChk( outFileId > 0, "Failed to create file: " + m_filename );
  return outFileId;
}

void
PDFWriter::writeExHdr( int outFileId, int nnode, int nelem ) const
// *****************************************************************************
//  Write Exodus II file header
//! \param[in] outFileId Output file ExodusII Id
//! \param[in] nnode Number of nodes in mesh to write
//! \param[in] nelem Number of elements in mesh to write
//! \author J. Bakosi
// *****************************************************************************
{
  ErrChk( ex_put_init( outFileId,
                       "Written by Quinoa",
                       3,                     // number of dimensions
                       nnode,                 // number of nodes
                       nelem,                 // number of elements
                       1,                     // number of element blocks
                       0,                     // number of node sets
                       0 ) == 0,              // number of side sets
          "Failed to write header to file: " + m_filename );
}

void
PDFWriter::writeExVar( int exoFile,
                       std::size_t,
                       ctr::PDFCenteringType centering,
                       const std::vector< tk::real >& probability ) const
// *****************************************************************************
//  Output probability density function as Exodus II results field
//! \param[in] exoFile ExodusII file handle to write to
//! \param[in] centering Node-, or element-centering to use on sample space mesh
//! \param[in] probability Probabilities at each sample space location
//! \author J. Bakosi
// *****************************************************************************
{
  if (centering == ctr::PDFCenteringType::NODE)
    ErrChk( ex_put_var( exoFile,
                        1,
                        EX_NODE_BLOCK,
                        1,
                        1,
                        static_cast<int64_t>(probability.size()),
                        probability.data() ) == 0,
            "Failed to write node-centered bivariate PDF to file: " +
            m_filename );
  else
    ErrChk( ex_put_var( exoFile,
                        1,
                        EX_ELEM_BLOCK,
                        1,
                        1,
                        static_cast<int64_t>(probability.size()),
                        probability.data() ) == 0,
            "Failed to write elem-centered bivariate PDF to file: " +
            m_filename );
}

void
PDFWriter::extents( const UniPDF& pdf,
                    const std::vector< tk::real >& uext,
                    std::size_t& nbi,
                    tk::real& min,
                    tk::real& max,
                    tk::real& binsize,
                    std::array< long, 2*UniPDF::dim >& ext,
                    std::vector< tk::real >& outpdf ) const
// *****************************************************************************
//  Query extents and other metadata of univariate PDF sample space
//! \details Query and optionally override number of bins and minimum of sample
//!    space if user-specified extents were given and copy probabilities from
//!    pdf to an array for output for plotting univariate PDF.
//! \param[in] pdf Univariate PDF object
//! \param[in] uext User-specified extents of sample space
//! \param[inout] nbi Number of bins
//! \param[inout] min Minimum value of sample space
//! \param[inout] max Maximum value of sample space
//! \param[inout] binsize Bin size
//! \param[inout] ext Extents of sample space
//! \param[inout] outpdf PDF ready to be written out to file
//! \author J. Bakosi
// *****************************************************************************
{
  assertSampleSpaceExtents< 1 >( uext );

  // Query bin size and extents of sample space from PDF
  binsize = pdf.binsize();
  ext = pdf.extents();

  // Compute number of bins of sample space (min bins: 1)
  Assert( ext[1] >= ext[0], "Wrong extents in PDFWriter::extents" );
  nbi = static_cast< std::size_t >( ext[1] - ext[0] + 1 );

  // Compute minimum and maximum of sample space
  min = binsize * static_cast< tk::real >( ext[0] );
  max = binsize * static_cast< tk::real >( ext[1] );

  // Override number of bins and minimum if user-specified extents were given,
  // and copy probabilities from pdf to an array for output
  if (!uext.empty()) {
    // Override number of bins by that based on user-specified extents
    Assert( uext[1] >= uext[0],
            "Wrong user-defined extents in PDFWriter::extents" );
    nbi = static_cast< std::size_t >(
            std::lround( (uext[1] - uext[0]) / binsize ) );
    // Override extents
    min = uext[0];
    max = uext[1];

    // Size output pdf to user-requested dimensions to overridden nbi and
    // initialize output probabilities to zero
    outpdf = std::vector< tk::real >( nbi, 0.0 );

    // Fill requested region of pdf to be output from computed pdf
    for (const auto& p : pdf.map()) {
      // Compute (i.e., shift) bin indices relative to user-requested extents
      const auto bin = p.first - std::lround( uext[0] / binsize );
      // Only copy probability value if shifted bin indices fall within
      // user-requested extents (lower inclusive, upper exclusive)
      if (bin >= 0 && bin < std::lround( (uext[1] - uext[0]) / binsize )) {
        Assert( static_cast<std::size_t>(bin) < nbi,
                "Bin overflow in user-specified-extent-based bin "
                "calculation of univariate PDF extents." );
        // Copy normalized probability to output pdf
        outpdf[ static_cast<std::size_t>(bin) ] =
          p.second / binsize / static_cast<tk::real>(pdf.nsample());
      }
    }
  }
}

void
PDFWriter::extents( const BiPDF& pdf,
                    const std::vector< tk::real >& uext,
                    std::size_t& nbix,
                    std::size_t& nbiy,
                    tk::real& xmin,
                    tk::real& xmax,
                    tk::real& ymin,
                    tk::real& ymax,
                    std::array< tk::real, BiPDF::dim >& binsize,
                    std::array< long, 2*BiPDF::dim >& ext,
                    std::vector< tk::real >& outpdf,
                    ctr::PDFCenteringType centering ) const
// *****************************************************************************
//  Query extents and other metadata of bivariate PDF sample space
//! \details Query and optionally override number of bins and minima of sample
//!    space if user-specified extents were given and copy probabilities from
//!    pdf to a logically 2D array for output for plotting bivariate joint PDF.
//! \param[in] pdf Bivariate PDF object
//! \param[in] uext User-specified extents of sample space
//! \param[inout] nbix Number of bins in x dimension
//! \param[inout] nbiy Number of bins in y dimension
//! \param[inout] xmin Minimum x value of sample space
//! \param[inout] xmax Maximum x value of sample space
//! \param[inout] ymin Minimum y value of sample space
//! \param[inout] ymax Maximum y value of sample space
//! \param[inout] binsize Bin size
//! \param[inout] ext Extents of sample space
//! \param[inout] outpdf PDF ready to be written out to file
//! \param[in] centering Bin centering on sample space mesh
//! \author J. Bakosi
// *****************************************************************************
{
  assertSampleSpaceExtents< 2 >( uext );

  // Query bin sizes and extents of sample space from PDF
  binsize = pdf.binsize();
  ext = pdf.extents();

  // Compute number of bins in sample space directions (min bins: 1)
  Assert( ext[1] >= ext[0], "Wrong extents in PDFWriter::extents" );
  Assert( ext[3] >= ext[2], "Wrong extents in PDFWriter::extents" );
  nbix = static_cast< std::size_t >( ext[1] - ext[0] + 1 );
  nbiy = static_cast< std::size_t >( ext[3] - ext[2] + 1 );

  // Compute minima and maxima of sample space
  xmin = binsize[0] * static_cast< tk::real >( ext[0] );
  xmax = binsize[0] * static_cast< tk::real >( ext[1] );
  ymin = binsize[1] * static_cast< tk::real >( ext[2] );
  ymax = binsize[1] * static_cast< tk::real >( ext[3] );

  // Override number of bins and minima if user-specified extents were given,
  // and copy probabilities from pdf to a logically 2D array for output
  if (!uext.empty()) {
    // Override number of bins by that based on user-specified extents
    Assert( uext[1] >= uext[0],
            "Wrong user-defined extents in PDFWriter::extents" );
    Assert( uext[3] >= uext[2],
            "Wrong user-defined extents in PDFWriter::extents" );
    nbix = static_cast< std::size_t >(
             std::lround( (uext[1] - uext[0]) / binsize[0] ) );
    nbiy = static_cast< std::size_t >(
             std::lround( (uext[3] - uext[2]) / binsize[1] ) );
    // Override extents
    xmin = uext[0];
    xmax = uext[1];
    ymin = uext[2];
    ymax = uext[3];

    // Temporarily increase number of bins if node-centered output required
    if (centering == ctr::PDFCenteringType::NODE) { ++nbix; ++nbiy; }

    // Size output pdf to user-requested dimensions to overridden nbiy * nbix
    // and initialize output probabilities to zero
    outpdf = std::vector< tk::real >( nbix*nbiy, 0.0 );

    // Fill requested region of pdf to be output from computed pdf
    for (const auto& p : pdf.map()) {
      // Compute (i.e., shift) bin indices relative to user-requested extents
      const auto x = p.first[0] - std::lround( uext[0] / binsize[0] );
      const auto y = p.first[1] - std::lround( uext[2] / binsize[1] );
      // Only copy probability value if shifted bin indices fall within
      // user-requested extents (lower inclusive, upper exclusive)
      if (x >= 0 && x < std::lround( (uext[1] - uext[0]) / binsize[0] ) &&
          y >= 0 && y < std::lround( (uext[3] - uext[2]) / binsize[1] ))
      {
        const auto bin =
          static_cast<std::size_t>(y)*nbix + static_cast<std::size_t>(x);
        Assert( bin < nbix*nbiy, "Bin overflow in user-specified-extent-based "
                "bin calculation of bivariate PDF." );
        // Copy normalized probability to output pdf
        outpdf[ bin ] = p.second / binsize[0] / binsize[1]
                                 / static_cast<tk::real>(pdf.nsample());
      }
    }

    // Revert number of bins if node-centered output required
    if (centering == ctr::PDFCenteringType::NODE) { --nbix; --nbiy; }
  }
}

void
PDFWriter::extents( const TriPDF& pdf,
                    const std::vector< tk::real >& uext,
                    std::size_t& nbix,
                    std::size_t& nbiy,
                    std::size_t& nbiz,
                    tk::real& xmin,
                    tk::real& xmax,
                    tk::real& ymin,
                    tk::real& ymax,
                    tk::real& zmin,
                    tk::real& zmax,
                    std::array< tk::real, TriPDF::dim >& binsize,
                    std::array< long, 2*TriPDF::dim >& ext,
                    std::vector< tk::real >& outpdf,
                    ctr::PDFCenteringType centering ) const
// *****************************************************************************
//  Query extents and other metadata of trivariate PDF sample space
//! \details Query and optionally override number of bins and minima of sample
//!   space if user-specified extents were given and copy probabilities from
//!   pdf to a logically 3D array for output for plotting trivariate joint PDF.
//! \param[in] pdf Trivariate PDF object
//! \param[in] uext User-specified extents of sample space
//! \param[inout] nbix Number of bins in x dimension
//! \param[inout] nbiy Number of bins in y dimension
//! \param[inout] nbiz Number of bins in z dimension
//! \param[inout] xmin Minimum x value of sample space
//! \param[inout] xmax Maximum x value of sample space
//! \param[inout] ymin Minimum y value of sample space
//! \param[inout] ymax Maximum y value of sample space
//! \param[inout] zmin Minimum z value of sample space
//! \param[inout] zmax Maximum z value of sample space
//! \param[inout] binsize Bin size
//! \param[inout] ext Extents of sample space
//! \param[inout] outpdf PDF ready to be written out to file
//! \param[in] centering Bin centering on sample space mesh
//! \author J. Bakosi
// *****************************************************************************
{
  assertSampleSpaceExtents< 3 >( uext );

  // Query bin sizes and extents of sample space from PDF
  binsize = pdf.binsize();
  ext = pdf.extents();

  // Compute number of bins in sample space directions (min bins: 1)
  Assert( ext[1] >= ext[0], "Wrong extents in PDFWriter::extents" );
  Assert( ext[3] >= ext[2], "Wrong extents in PDFWriter::extents" );
  Assert( ext[5] >= ext[4], "Wrong extents in PDFWriter::extents" );
  nbix = static_cast< std::size_t >( ext[1] - ext[0] + 1 );
  nbiy = static_cast< std::size_t >( ext[3] - ext[2] + 1 );
  nbiz = static_cast< std::size_t >( ext[5] - ext[4] + 1 );

  // Compute minima and maxima of sample space
  xmin = binsize[0] * static_cast< tk::real >( ext[0] );
  xmax = binsize[0] * static_cast< tk::real >( ext[1] );
  ymin = binsize[1] * static_cast< tk::real >( ext[2] );
  ymax = binsize[1] * static_cast< tk::real >( ext[3] );
  zmin = binsize[2] * static_cast< tk::real >( ext[4] );
  zmax = binsize[2] * static_cast< tk::real >( ext[5] );

  // Override number of bins and minima if user-specified extents were given,
  // and copy probabilities from pdf to a logically 3D array for output
  if (!uext.empty()) {
    // Override number of bins by that based on user-specified extents
    Assert( uext[1] >= uext[0],
            "Wrong user-defined extents in PDFWriter::extents" );
    Assert( uext[3] >= uext[2],
            "Wrong user-defined extents in PDFWriter::extents" );
    Assert( uext[5] >= uext[4],
            "Wrong user-defined extents in PDFWriter::extents" );
    nbix = static_cast< std::size_t >(
             std::lround( (uext[1] - uext[0]) / binsize[0] ) );
    nbiy = static_cast< std::size_t >(
             std::lround( (uext[3] - uext[2]) / binsize[1] ) );
    nbiz = static_cast< std::size_t >(
             std::lround( (uext[5] - uext[4]) / binsize[2] ) );
    // Override extents
    xmin = uext[0];
    xmax = uext[1];
    ymin = uext[2];
    ymax = uext[3];
    zmin = uext[4];
    zmax = uext[5];

    // Temporarily increase number of bins if node-centered output required
    if (centering == ctr::PDFCenteringType::NODE) { ++nbix; ++nbiy; ++nbiz; }

    // Size output pdf to user-requested dimensions to overridden nbiz * nbiy *
    // nbix and initialize output probabilities to zero
    outpdf = std::vector< tk::real >( nbiz * nbiy * nbix, 0.0 );

    // Fill requested region of pdf to be output from computed pdf
    for (const auto& p : pdf.map()) {
      // Compute (i.e., shift) bin indices relative to user-requested extents
      const auto x = p.first[0] - std::lround( uext[0] / binsize[0] );
      const auto y = p.first[1] - std::lround( uext[2] / binsize[1] );
      const auto z = p.first[2] - std::lround( uext[4] / binsize[2] );
      // Only copy probability value if shifted bin indices fall within
      // user-requested extents (lower inclusive, upper exclusive)
      if (x >= 0 && x < std::lround( (uext[1] - uext[0]) / binsize[0] ) &&
          y >= 0 && y < std::lround( (uext[3] - uext[2]) / binsize[1] ) &&
          z >= 0 && z < std::lround( (uext[5] - uext[4]) / binsize[2] ))
      {
        const auto X = static_cast< std::size_t >( x );
        const auto Y = static_cast< std::size_t >( y );
        const auto Z = static_cast< std::size_t >( z );
        const auto bin = nbix*(Z*nbiy + Y) + X;
        Assert( bin < nbix*nbiy*nbiz, "Bin overflow in "
              "user-specified-extent-based bin calculation of bivariate PDF." );
        // Copy normalized probability to output pdf
        outpdf[ bin ] =
          p.second / binsize[0] / binsize[1] / binsize[2]
                   / static_cast<tk::real>(pdf.nsample());
      }
    }

    // Revert number of bins if node-centered output required
    if (centering == ctr::PDFCenteringType::NODE) { --nbix; --nbiy; --nbiz; }
  }
}
