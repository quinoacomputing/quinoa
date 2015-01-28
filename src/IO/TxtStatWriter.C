//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.C
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 10:16:38 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Text statistics writer declaration
  \details   This file declares the ASCII statistics writer class that
     facilitates outputing statistics to text files.
*/
//******************************************************************************

#include <iostream>
#include <iomanip>

#include <TxtStatWriter.h>

using tk::TxtStatWriter;

void
TxtStatWriter::header( const std::vector< std::string >& nameOrd,
                       const std::vector< std::string >& nameCen ) const
//******************************************************************************
//  Write out statistics file header
//! \param[in] nameOrd Vector of strings with the names of ordinary moments
//! \param[in] nameCen Vector of strings with the names of central moments
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile << "#     it             t";

  // Output names of ordinary moments
  for (const auto& n : nameOrd) m_outFile << std::setw(12) << '<' << n << ">";

  // Output name of central moments
  for (const auto& c : nameCen) m_outFile << std::setw(10) << '<' << c << ">";

  m_outFile << std::endl;
}

std::size_t
TxtStatWriter::stat( int it,
                     tk::real t,
                     const std::vector< tk::real >& ordinary,
                     const std::vector< tk::real >& central )
//******************************************************************************
//  Write out statistics
//! \param[in] it Iteration counter
//! \param[in] t Time
//! \param[in] ordinary Vector with the ordinary moment statistics
//! \param[in] central Vector with the central moment statistics
//! \return The total number of statistics written to the output file
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile << std::setfill(' ') << std::setw(8) << it << "  "
            << std::scientific << std::setprecision(6) << std::setw(12) << t
            << "  ";

  // Output ordinary moments
  for (const auto& o : ordinary) m_outFile << " " << o << "  ";

  // Output central moments
  for (const auto& c : central) m_outFile << " " << c << "  ";

  m_outFile << std::endl;

  return ordinary.size() + central.size();
}
