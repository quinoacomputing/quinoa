//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.C
  \author    J. Bakosi
  \date      Thu Aug 29 15:32:48 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TxtStat writer base class definition
  \details   TxtStat writer base class definition
*/
//******************************************************************************

#include <iostream>
#include <iomanip>

#include <TxtStatWriter.h>
#include <Statistics.h>

using namespace quinoa;

TxtStatWriter::TxtStatWriter(const std::string& filename,
                             Statistics* const statistics) :
  Writer(filename),
  m_statistics(statistics),
  m_nord(statistics->nord()),
  m_ncen(statistics->ncen()),
  m_ordinary(statistics->ordinary()),
  m_central(statistics->central())
//******************************************************************************
//  Constructor
//! \param[in]  filename     Filename of txt statistics output
//! \param[in]  statistics   Statistics estimator
//! \author J. Bakosi
//******************************************************************************
{
}

void
TxtStatWriter::header()
//******************************************************************************
//  Write out statistics file header
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile << "#     it             t";

  for (int i=0; i<m_nord; ++i)
    if (m_statistics->plotOrdinary(i))
      m_outFile << std::setw(12) << '<' << m_statistics->nameOrdinary(i) << ">";

  for (int i=0; i<m_ncen; ++i)
    m_outFile << std::setw(10) << '<' << m_statistics->nameCentral(i) << ">";

  m_outFile << std::endl;
}

void
TxtStatWriter::write(const int it, const real t)
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
    if (m_statistics->plotOrdinary(i))
      m_outFile << (m_ordinary[i]>0 ? " " : "") << m_ordinary[i] << "  ";
  }

  // Output central moments
  for (int i=0; i<m_ncen; ++i) {
    m_outFile << (m_central[i]>0 ? " " : "") << m_central[i] << "  ";
  }

  m_outFile << std::endl;
}
