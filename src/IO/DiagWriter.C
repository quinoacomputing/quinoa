// *****************************************************************************
/*!
  \file      src/IO/DiagWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Text diagnostics writer declaration
  \details   This file declares the ASCII diagnostics writer class that
     facilitates outputing diagnostics to text files.
*/
// *****************************************************************************

#include <iostream>
#include <iomanip>

#include "DiagWriter.h"

using tk::DiagWriter;

DiagWriter::DiagWriter( const std::string& filename,
                              ctr::TxtFloatFormatType format,
                              kw::precision::info::expect::type precision,
                              std::ios_base::openmode mode ) :
  Writer( filename, mode ),
  m_precision( static_cast<int>(precision) ),
  m_width( std::max( 16, m_precision+8 ) )
// *****************************************************************************
//  Constructor
//! \param[in] filename Output filename to which output the diagnostics
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
DiagWriter::header( const std::vector< std::string >& name ) const
// *****************************************************************************
//  Write out diagnostics file header
//! \param[in] name Vector of strings with the names of diagnostics
//! \author J. Bakosi
// *****************************************************************************
{
  m_outFile << "#" << std::setw(9) << "1:it";
  m_outFile << std::setw(m_width) << "2:t";
  std::stringstream out;

  // Output names of diagnostics
  std::size_t column = 3;
  for (const auto& n : name) {
    out << column++ << ':' << n;
    m_outFile << std::setw(m_width) << out.str();
    out.str("");
  }

  m_outFile << std::endl;
}

std::size_t
DiagWriter::diag( uint64_t it,
                  tk::real t,
                  const std::vector< tk::real >& diagnostics )
// *****************************************************************************
//  Write out diagnostics
//! \param[in] it Iteration counter
//! \param[in] t Time
//! \param[in] diagnostics Vector with the diagnostics
//! \return The total number of diagnostics written to the output file
//! \author J. Bakosi
// *****************************************************************************
{
  m_outFile << std::setw(10) << it;
  m_outFile << std::setw(m_width) << t;

  // Output diagnostics
  for (const auto& d : diagnostics) m_outFile << std::setw(m_width) << d;

  m_outFile << std::endl;

  return diagnostics.size();
}
