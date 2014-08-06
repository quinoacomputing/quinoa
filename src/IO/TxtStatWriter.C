//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.C
  \author    J. Bakosi
  \date      Tue 05 Aug 2014 04:48:26 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     TxtStat writer base class definition
  \details   TxtStat writer base class definition
*/
//******************************************************************************

#include <iostream>
#include <iomanip>

#include <TxtStatWriter.h>

using quinoa::TxtStatWriter;

void
TxtStatWriter::header() const
//******************************************************************************
//  Write out statistics file header
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile << "#     it             t";

  for (int i=0; i<m_nord; ++i)
    if ( m_plotOrdinary[i] )
      m_outFile << std::setw(12) << '<' << m_nameOrdinary[i] << ">";

  for (int i=0; i<m_ncen; ++i)
    m_outFile << std::setw(10) << '<' << m_nameCentral[i] << ">";

  m_outFile << std::endl;
}

void
TxtStatWriter::writeStat(const int it, const tk::real t)
//******************************************************************************
//  Write out statistics
//! \param[in]  it         Iteration counter
//! \param[in]  t          Time
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile << std::setfill(' ') << std::setw(8) << it << "  "
            << std::scientific << std::setprecision(6) << std::setw(12) << t
            << "  ";

  // Output ordinary moments
  for (int i=0; i<m_nord; ++i) {
    if ( m_plotOrdinary[i] )
      m_outFile << (m_ordinary[i]>0 ? " " : "") << m_ordinary[i] << "  ";
  }

  // Output central moments
  for (int i=0; i<m_ncen; ++i) {
    m_outFile << (m_central[i]>0 ? " " : "") << m_central[i] << "  ";
  }

  m_outFile << std::endl;
}
