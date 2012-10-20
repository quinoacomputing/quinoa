//******************************************************************************
/*!
  \file      src/Statistics/PDF.C
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 09:48:06 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     PDF estimator base class
  \details   PDF estimator base class
*/
//******************************************************************************

#include <string>
#include <fstream>

#include <PDF.h>
#include <MemoryException.h>

using namespace Quinoa;

void
PDF::save(const string& filename)
//******************************************************************************
//  Insert new value into PDF
//! \param[in]   value    Value to insert
//! \author  J. Bakosi
//******************************************************************************
{
  ofstream pdf(filename, ofstream::out);

  //const real sx = m_nbin;///(m_max-m_min);
  const real sp = m_sample*m_binsize;
  for (pair<int,real> p : m_pdf) {
    pdf << p.first*m_binsize << "\t" << p.second/sp << endl;
  }

  pdf.close();
}
