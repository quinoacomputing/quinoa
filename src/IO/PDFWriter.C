//******************************************************************************
/*!
  \file      src/IO/PDFWriter.C
  \author    J. Bakosi
  \date      Sun 28 Sep 2014 09:52:14 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Univariate PDF writer
  \details   Univariate PDF writer
*/
//******************************************************************************

#include <PDFWriter.h>
#include <Exception.h>

using quinoa::PDFWriter;

void
PDFWriter::writeTxt( const UniPDF& pdf, const std::vector< tk::real >& uext )
const
//******************************************************************************
//  Write out standardized univariate PDF to file
//! \param[in]  pdf   Univariate PDF
//! \param[in]  uext  Optional user-specified extents of the sample space
//! \author  J. Bakosi
//******************************************************************************
{
  if (!uext.empty())
    Assert( uext.size() == 2, "Univariate PDF user-specified sample space "
            "extents must be defined by two real numbers: min, max" );

  // Query and optionally override number of bins and minimum of sample space if
  // user-specified extents were given and copy probabilities from pdf to an
  // array for output
  std::size_t nbi;
  tk::real min, max;
  std::vector< tk::real > outpdf;
  tk::real binsize;
  std::pair< UniPDF::key_type, UniPDF::key_type > ext;
  extents( pdf, uext, nbi, min, max, binsize, ext, outpdf );

  // Output header
  m_outFile << "# Univariate PDF\n"
            << "# --------------\n"
            << "# Bin size: " << binsize << '\n'
            << "# Number of bins estimated: " << ext.second - ext.first + 1
            << '\n'
            << "# Number of bins output: " << nbi << '\n'
            << "# Sample space extent: [" << min << " : " << max << "]\n#\n"
            << "# Example step-by-step visualization with gnuplot\n"
            << "# -----------------------------------------------\n"
            << "# gnuplot> set grid\n"
            << "# gnuplot> plot ";
  if (!uext.empty()) m_outFile << "[" << uext[0] << ':' << uext[1] << "] ";
  m_outFile << "\"" << m_filename << "\" with points\n#\n"
            << "# Gnuplot one-liner for quick copy-paste\n"
            << "# --------------------------------------\n"
            << "# set grid; plot";
  if (!uext.empty()) m_outFile << " [" << uext[0] << ':' << uext[1] << "]";
  m_outFile << " \"" << m_filename << "\" w p\n#\n"
            << "# -----------< data columns: x, probability >-----------\n";

  // If no user-specified sample space extents, output pdf map directly
  m_outFile << std::scientific << std::setprecision(12);
  if (uext.empty()) {
    for (const auto& p : pdf.map())
      m_outFile << binsize * p.first << '\t'
                << p.second / binsize / pdf.nsample() << std::endl;
  } else { // If user-specified sample space extents, output outpdf array
    std::size_t bin = 0;
    for (const auto& p : outpdf)
      m_outFile << binsize * bin++ + uext[0] << '\t' << p << std::endl;
  }
}

void
PDFWriter::writeTxt( const BiPDF& pdf, const std::vector< tk::real >& uext )
const
//******************************************************************************
//  Write out standardized bivariate PDF to text file
//! \param[in]  pdf   Bivariate PDF
//! \param[in]  uext  Optional user-specified extents of the sample space
//! \author  J. Bakosi
//******************************************************************************
{
  if (!uext.empty())
    Assert( uext.size() == 4, "Bivariate joint PDF user-specified sample space "
       "extents must be defined by four real numbers: minx, maxx, miny, maxy" );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 2D array for output
  std::size_t nbix, nbiy;
  tk::real xmin, xmax, ymin, ymax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 2 > binsize;
  std::pair< BiPDF::key_type, BiPDF::key_type > ext;
  extents( pdf, uext, nbix, nbiy, xmin, xmax, ymin, ymax, binsize, ext, outpdf );

  // Output header
  m_outFile << "# Joint bivariate PDF\n"
            << "# -------------------\n"
            << "# Bin sizes: " << binsize[0] << ", " << binsize[1] << '\n'
            << "# Number of bins estimated: " << ext.first[1] - ext.first[0] + 1
            << " x " << ext.second[1] - ext.second[0] + 1 << '\n'
            << "# Number of bins output: " << nbix << " x " << nbiy << '\n'
            << "# Sample space extents: [" << xmin << " : " << xmax
            << "], [" << ymin << " : " << ymax << "]\n#\n"
            << "# Example step-by-step visualization with gnuplot\n"
            << "# -----------------------------------------------\n"
            << "# gnuplot> set grid\n"
            << "# gnuplot> set dgrid3d 50,50,1\n"
            << "# gnuplot> set cntrparam levels 20\n"
            << "# gnuplot> set contour\n";
  if (!uext.empty())
    m_outFile << "# gnuplot> set xrange [" << uext[0] << ':' << uext[1] << "]\n"
              << "# gnuplot> set yrange [" << uext[2] << ':' << uext[3] << "]\n";
         
  m_outFile << "# gnuplot> splot \"" << m_filename << "\" with lines\n#\n"
            << "# Gnuplot one-liner for quick copy-paste\n"
            << "# --------------------------------------\n"
            << "# set grid; set dgrid3d 50,50,1; set cntrparam levels 20; set "
               "contour; ";
  if (!uext.empty())
    m_outFile << "set xrange [" << uext[0] << ':' << uext[1] << "]; set yrange "
                 "[" << uext[2] << ':' << uext[3] << "]; ";
  m_outFile << "splot \"" << m_filename << "\" w l\n#\n"
            << "# -----------< data columns: x, y, probability >-----------\n";

  // If no user-specified sample space extents, output pdf map directly
  m_outFile << std::scientific << std::setprecision(12);
  if (uext.empty()) {
    for (const auto& p : pdf.map())
      m_outFile << binsize[0] * p.first[0] << '\t'
                << binsize[1] * p.first[1] << '\t'
                << p.second / binsize[0] / binsize[1] / pdf.nsample()
                << std::endl;
  } else { // If user-specified sample space extents, output outpdf array
    std::size_t bin = 0;
    for (const auto& p : outpdf) {
      m_outFile << binsize[0] * (bin % nbix) + uext[0] << '\t'
                << binsize[1] * (bin / nbix) + uext[2] << '\t'
                << p
                << std::endl;
      ++bin;
    }
  }

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}

void
PDFWriter::writeGmsh( const BiPDF& pdf,
                      const std::string& pdfname,
                      ctr::PDFCenteringType centering,
                      const std::vector< tk::real >& uext ) const
//******************************************************************************
//  Write out standardized bivariate PDF to Gmsh (text) format
//! \param[in]  pdf   Bivariate PDF
//! \param[in]  uext  Optional user-specified extents of the sample space
//! \author  J. Bakosi
//******************************************************************************
{
  // Output mesh header: mesh version, file type, data size
  m_outFile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );

  // Query and optionally override number of bins and minima of sample space if
  // user-specified extents were given and copy probabilities from pdf to a
  // logically 2D array for output
  std::size_t nbix, nbiy;
  tk::real xmin, xmax, ymin, ymax;
  std::vector< tk::real > outpdf;
  std::array< tk::real, 2 > binsize;
  std::pair< BiPDF::key_type, BiPDF::key_type > ext;
  extents( pdf, uext, nbix, nbiy, xmin, xmax, ymin, ymax, binsize, ext, outpdf );

  // Output grid points of discretized sample space (2D Cartesian grid)
  m_outFile << "$Nodes\n" << (nbix+1)*(nbiy+1) << std::endl;
  int k=0;
  for (int i=0; i<=nbiy; i++) {
    tk::real y = ymin + i*binsize[1];
    for (int j=0; j<=nbix; j++) {
      tk::real x = xmin + j*binsize[0];
      m_outFile << ++k << " " << x << " " << y << " 0\n";
    }
  }
  m_outFile << "$EndNodes\n";

  // Output elements of discretized sample space (2D Cartesian grid)
  m_outFile << "$Elements\n" << nbix*nbiy << "\n";
  for (int i=0; i<nbix*nbiy; ++i) {
    m_outFile << i+1 << " 3 2 1 1 " << i+i/nbix+1 << " "
              << i+2+i/nbix << " " << i+3+nbix+i/nbix << " "
              << i+2+nbix+i/nbix << std::endl;
  }
  m_outFile << "$EndElements\n";

  // Output PDF function values in element or node centers
  std::string c( "Element" );
  if (centering == ctr::PDFCenteringType::NODE) {
    ++nbix; ++nbiy;
    c = "Node";
  }
  m_outFile << "$" << c << "Data\n1\n\"" << pdfname << "\"\n1\n0.0\n3\n0\n1\n"
            << nbix*nbiy << "\n";

  // If no user-specified sample space extents, output pdf map directly
  m_outFile << std::scientific << std::setprecision(12);
  if (uext.empty()) {

    std::vector< bool > out( nbix*nbiy, false ); // indicate bins filled
    const auto ext = pdf.extents();
    for (const auto& p : pdf.map()) {
      const auto bin = (p.first[1] - ext.second[0]) * nbix +
                       (p.first[0] - ext.first[0]) % nbix;
      Assert( bin < nbix*nbiy, "Bin overflow in PDFWriter::writeGmsh()." );
      out[ bin ] = true;
      m_outFile << bin+1 << '\t'
                << p.second / binsize[0] / binsize[1] / pdf.nsample()
                << std::endl;
    }
    // Output bins nonexistent in PDF (gmsh sometimes fails to plot the exiting
    // bins if holes exist in the data, it also looks better as zero than holes)
    for (std::size_t i=0; i<out.size(); ++i)
      if (!out[i]) m_outFile << i+1 << "\t0" << std::endl;

  } else { // If user-specified sample space extents, output outpdf array

    std::size_t bin = 0;
    for (const auto& p : outpdf) m_outFile << ++bin << " " << p << std::endl;

  }

  m_outFile << "$End" << c << "Data\n";

  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
}

void
PDFWriter::extents( const UniPDF& pdf,
                    const std::vector< tk::real >& uext,
                    std::size_t& nbi,
                    tk::real& min,
                    tk::real& max,
                    tk::real& binsize,
                    std::pair< UniPDF::key_type, UniPDF::key_type >& ext,
                    std::vector< tk::real >& outpdf ) const
//******************************************************************************
//  Query and optionally override number of bins and minimum of sample space if
//  user-specified extents were given and copy probabilities from pdf to an
//  array for output for plotting univariate PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Query bin size and extents of sample space from PDF
  binsize = pdf.binsize();
  ext = pdf.extents();

  // Compute number of bins of sample space (min bins: 1)
  nbi = ext.second - ext.first + 1;

  // Compute minimum and maximum of sample space
  min = binsize * ext.first;
  max = binsize * ext.second;

  // Override number of bins and minimum if user-specified extents were given,
  // and copy probabilities from pdf to an array for output
  if (!uext.empty()) {
    Assert( uext.size() == 2, "Bivariate joint PDF user-specified sample space "
       "extents must be defined by four real numbers: minx, maxx, miny, maxy" );

    // Override number of bins by that based on user-specified extents
    nbi = std::lround( (uext[1] - uext[0]) / binsize );
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
        Assert( bin < nbi, "Bin overflow in user-specified-extent-based bin "
                           "calculation of univariate PDF extents." );
        // Copy normalized probability to output pdf
        outpdf[ bin ] = p.second / binsize / pdf.nsample();
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
                    std::array< tk::real, 2 >& binsize,
                    std::pair< BiPDF::key_type, BiPDF::key_type >& ext,
                    std::vector< tk::real >& outpdf ) const
//******************************************************************************
//  Query and optionally override number of bins and minima of sample space if
//  user-specified extents were given and copy probabilities from pdf to a
//  logically 2D array for output for plotting bivariate joint PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Query bin sizes and extents of sample space from PDF
  binsize = pdf.binsize();
  ext = pdf.extents();

  // Compute number of bins in sample space directions (min bins: 1)
  nbix = ext.first[1] - ext.first[0] + 1;
  nbiy = ext.second[1] - ext.second[0] + 1;

  // Compute minima and maxima of sample space
  xmin = binsize[0] * ext.first[0];
  xmax = binsize[0] * ext.first[1];
  ymin = binsize[1] * ext.second[0];
  ymax = binsize[1] * ext.second[1];

  // Override number of bins and minima if user-specified extents were given,
  // and copy probabilities from pdf to a logically 2D array for output
  if (!uext.empty()) {
    Assert( uext.size() == 4, "Bivariate joint PDF user-specified sample space "
          "extents must be defined by four real numbers: minx,maxx,miny,maxy" );

    // Override number of bins by that based on user-specified extents
    nbix = std::lround( (uext[1] - uext[0]) / binsize[0] );
    nbiy = std::lround( (uext[3] - uext[2]) / binsize[1] );
    // Override extents
    xmin = uext[0];
    xmax = uext[1];
    ymin = uext[2];
    ymax = uext[3];

    // Size output pdf to user-requested dimensions to overridden nbiy * nbix
    // and initialize output probabilities to zero
    outpdf = std::vector< tk::real >( nbix*nbiy, 0.0 );

    // Fill requested region of pdf to be output from computed pdf
    for (const auto& p : pdf.map()) {
      // Compute (i.e., shift) bin indices relative to user-requested extents
      const auto binx = p.first[0] - std::lround( uext[0] / binsize[0] );
      const auto biny = p.first[1] - std::lround( uext[2] / binsize[1] );
      // Only copy probability value if shifted bin indices fall within
      // user-requested extents (lower inclusive, upper exclusive)
      if (binx >= 0 && binx < std::lround( (uext[1] - uext[0]) / binsize[0] ) &&
          biny >= 0 && biny < std::lround( (uext[3] - uext[2]) / binsize[1] ))
      {
        Assert( binx + biny*nbix < nbix*nbiy, "Bin overflow in user-specified-"
                "-extent-based bin calculation of bivariate PDF extents." );
        // Copy normalized probability to output pdf
        outpdf[ binx + biny*nbix ] =
          p.second / binsize[0] / binsize[1] / pdf.nsample();
      }
    }
  }
}
