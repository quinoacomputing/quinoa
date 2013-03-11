//******************************************************************************
/*!
  \file      src/IO/TxtPlotWriter.C
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 08:50:55 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TxtPlot writer base class definition
  \details   TxtPlot writer base class definition
*/
//******************************************************************************

#include <iostream>
#include <iomanip>

#include <TxtPlotWriter.h>

using namespace Quinoa;

TxtPlotWriter::TxtPlotWriter(string filename) : PlotWriter(filename)
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
}

void
TxtPlotWriter::write()
//******************************************************************************
//  Write out plot
//! \author J. Bakosi
//******************************************************************************
{
//   m_outPlot << setfill(' ') << setw(8) << it << "  " << scientific
//             << setprecision(6) << setw(12) << t << "  ";
// 
//   for (int i=0; i<nord; ++i) {
//     m_outPlot << ordinary[i] << "  ";
//   }
// 
//   m_outPlot << endl;
}
