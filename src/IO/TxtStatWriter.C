//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.C
  \author    J. Bakosi
  \date      Wed 21 Jan 2015 03:57:31 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     TxtStat writer base class definition
  \details   TxtStat writer base class definition
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
//! \param[in]  it         Iteration counter
//! \param[in]  t          Time
//! \return Bool signaling whether anything was written into file
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
