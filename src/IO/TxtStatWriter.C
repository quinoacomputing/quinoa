//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.C
  \author    J. Bakosi
  \date      Fri 05 Sep 2014 12:04:51 PM MDT
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
TxtStatWriter::header( const std::vector< bool >& plotOrd,
                       const std::vector< std::string >& nameOrd,
                       const std::vector< std::string >& nameCen ) const
//******************************************************************************
//  Write out statistics file header
//! \author J. Bakosi
//******************************************************************************
{
  Assert( plotOrd.size() == nameOrd.size(), "plotOrd.size()!=nameOrd.size()" );

  m_outFile << "#     it             t";

  const auto nord = nameOrd.size();
  for (auto i=decltype(nord){0}; i<nord; ++i)
    if (plotOrd[i]) m_outFile << std::setw(12) << '<' << nameOrd[i] << ">";

  for (const auto& c : nameCen) m_outFile << std::setw(10) << '<' << c << ">";

  m_outFile << std::endl;
}

void
TxtStatWriter::writeStat( int it,
                          tk::real t,
                          const std::vector< tk::real >& ordinary,
                          const std::vector< tk::real >& central,
                          const std::vector< bool >& plotOrd )
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
  const auto nord = ordinary.size();
  for (auto i=decltype(nord){0}; i<nord; ++i)
    if (plotOrd[i])
      //m_outFile << (ordinary[i] > 0 ? " " : "") << ordinary[i] << "  ";
      m_outFile << " " << ordinary[i] << "  ";

  // Output central moments
  for (const auto& c : central) m_outFile << " " << c << "  ";

  m_outFile << std::endl;
}
