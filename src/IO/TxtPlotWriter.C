//******************************************************************************
/*!
  \file      src/IO/TxtPlotWriter.C
  \author    J. Bakosi
  \date      Tue 12 Mar 2013 10:51:12 PM MDT
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

TxtPlotWriter::TxtPlotWriter(const string& filename,
                             Statistics* const statistics) :
  PlotWriter(filename),
  m_statistics(statistics),
  m_nord(statistics->nord()),
  m_ordinary(statistics->ordinary())
//******************************************************************************
//  Constructor
//! \param[in]  filename     Filename of txt plot output
//! \param[in]  statistics   Statistics estimator
//! \author J. Bakosi
//******************************************************************************
{
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
  m_outPlot << setfill(' ') << setw(8) << it << "  " << scientific
            << setprecision(6) << setw(12) << t << "  ";

  for (int i=0; i<m_nord; ++i) {
    m_outPlot << m_ordinary[i] << "  ";
  }

  m_outPlot << endl;
}
