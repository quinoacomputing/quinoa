//******************************************************************************
/*!
  \file      src/IO/TxtPlotWriter.C
  \author    J. Bakosi
  \date      Tue Jul  2 16:30:01 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TxtPlot writer base class definition
  \details   TxtPlot writer base class definition
*/
//******************************************************************************

#include <iostream>
#include <iomanip>

#include <TxtPlotWriter.h>
#include <Statistics.h>

using namespace Quinoa;

TxtPlotWriter::TxtPlotWriter(const std::string& filename,
                             Statistics* const statistics) :
  PlotWriter(filename),
  m_statistics(statistics),
  m_nord(statistics->nord()),
  m_ncen(statistics->ncen()),
  m_ordinary(statistics->ordinary()),
  m_central(statistics->central())
//******************************************************************************
//  Constructor
//! \param[in]  filename     Filename of txt plot output
//! \param[in]  statistics   Statistics estimator
//! \author J. Bakosi
//******************************************************************************
{
}

void
TxtPlotWriter::header()
//******************************************************************************
//  Write out plot header
//! \author J. Bakosi
//******************************************************************************
{
  m_outPlot << "#     it             t";

  for (int i=0; i<m_nord; ++i)
    if (m_statistics->plotOrdinary(i))
      m_outPlot << std::setw(12) << '<' << m_statistics->nameOrdinary(i) << ">";

  for (int i=0; i<m_ncen; ++i)
    m_outPlot << std::setw(10) << '<' << m_statistics->nameCentral(i) << ">";

  m_outPlot << std::endl;
}

void
TxtPlotWriter::write(const int it, const real t)
//******************************************************************************
//  Write out plot
//! \param[in]  it         Iteration counter
//! \param[in]  t          Time
//! \author J. Bakosi
//******************************************************************************
{
  m_outPlot << std::setfill(' ') << std::setw(8) << it << "  "
            << std::scientific << std::setprecision(6) << std::setw(12) << t
            << "  ";

  // Plot ordinary moments
  for (int i=0; i<m_nord; ++i) {
    if (m_statistics->plotOrdinary(i))
      m_outPlot << (m_ordinary[i]>0 ? " " : "") << m_ordinary[i] << "  ";
  }

  // Plot central moments
  for (int i=0; i<m_ncen; ++i) {
    m_outPlot << (m_central[i]>0 ? " " : "") << m_central[i] << "  ";
  }

  m_outPlot << std::endl;
}
