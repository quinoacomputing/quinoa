// *****************************************************************************
/*!
  \file      src/IO/TxtStatWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Text statistics writer declaration
  \details   This file declares the ASCII statistics writer class that
     facilitates outputing statistics to text files.
*/
// *****************************************************************************

#include <iostream>
#include <iomanip>

#include "TxtStatWriter.h"

using tk::TxtStatWriter;

TxtStatWriter::TxtStatWriter( const std::string& filename,
                              ctr::TxtFloatFormatType format,
                              kw::precision::info::expect::type precision,
                              std::ios_base::openmode mode ) :
  Writer( filename, mode ),
  m_precision( static_cast<int>(precision) ),
  m_width( std::max( 16, m_precision+8 ) )
// *****************************************************************************
//  Constructor
//! \param[in] filename Output filename to which output the statistics
//! \param[in] format Configure floating-point output format ASCII output
//! \param[in] precision Configure precision for floating-point ASCII output
//! \param[in] mode Configure file open mode
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
TxtStatWriter::header( const std::vector< std::string >& nameOrd,
                       const std::vector< std::string >& nameCen,
                       const std::vector< std::string >& nameExt ) const
// *****************************************************************************
//  Write out statistics file header
//! \param[in] nameOrd Vector of strings with the names of ordinary moments
//! \param[in] nameCen Vector of strings with the names of central moments
//! \param[in] nameExt Vector of strings with the names of extra data
//! \author J. Bakosi
// *****************************************************************************
{
  m_outFile << "#" << std::setw(9) << "1:it";
  m_outFile << std::setw(m_width) << "2:t";
  std::stringstream out;

  // Output names of ordinary moments
  std::size_t column = 3;
  for (const auto& n : nameOrd) {
    out << column++ << ":<" << n << ">";
    m_outFile << std::setw(m_width) << out.str();
    out.str("");
  }

  // Output name of central moments
  for (const auto& c : nameCen) {
    out << column++ << ":<" << c << ">";
    m_outFile << std::setw(m_width) << out.str();
    out.str("");
  }

  // Output name of extra data
  for (const auto& c : nameExt) {
    out << column++ << ":<" << c << ">";
    m_outFile << std::setw(m_width) << out.str();
    out.str("");
  }

  m_outFile << std::endl;
}

std::size_t
TxtStatWriter::stat( uint64_t it,
                     tk::real t,
                     const std::vector< tk::real >& ordinary,
                     const std::vector< tk::real >& central,
                     const std::vector< tk::real >& extra )
// *****************************************************************************
//  Write out statistics
//! \param[in] it Iteration counter
//! \param[in] t Time
//! \param[in] ordinary Vector with the ordinary moment statistics
//! \param[in] central Vector with the central moment statistics
//! \param[in] extra Vector with extra data to be also written (besides stats)
//! \return The total number of statistics written to the output file
//! \author J. Bakosi
// *****************************************************************************
{
  m_outFile << std::setw(10) << it;
  m_outFile << std::setw(m_width) << t;

  // Output ordinary moments
  for (const auto& o : ordinary) m_outFile << std::setw(m_width) << o;

  // Output central moments
  for (const auto& c : central) m_outFile << std::setw(m_width) << c;

  // Output extra data
  for (const auto& c : extra) m_outFile << std::setw(m_width) << c;

  m_outFile << std::endl;

  return ordinary.size() + central.size() + extra.size();
}
