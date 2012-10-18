//******************************************************************************
/*!
  \file      src/Statistics/PDF.C
  \author    J. Bakosi
  \date      Wed 17 Oct 2012 09:30:52 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     PDF estimator base class
  \details   PDF estimator base class
*/
//******************************************************************************

#include <PDF.h>
#include <MemoryException.h>

using namespace Quinoa;

PDF::PDF(const Real min, const Real max, const Real nbin) :
  m_min(min), m_max(max), m_nbin(nbin), m_binSize((max-min)/nbin),
  m_halfBinSize(m_binSize/2)
//******************************************************************************
//  Constructor: Initialize PDF container
//! \param[in]   min     Most negative value of the sample space
//! \param[in]   max     Most positive value of the sample space
//! \param[in]   nbin    Number of bins of the sample space
//! \author  J. Bakosi
//******************************************************************************
{
  // Initialize PDF container with the center-point of the bins as keys and
  // zeros as mapped values
  for (Int i=0; i<m_nbin; ++i) {
    pair<Pdf::iterator,Bool> n =
      m_pdf.emplace(m_min + m_binSize*i + m_halfBinSize, 0.0);
    if (!n.second)
      throw MemoryException(FATAL, BAD_INSERT);
  }
}

PDF::~PDF()
//******************************************************************************
//  Destructor: Clear PDF container
//! \author  J. Bakosi
//******************************************************************************
{
  m_pdf.clear();
}
