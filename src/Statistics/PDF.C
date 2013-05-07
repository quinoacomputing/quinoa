//******************************************************************************
/*!
  \file      src/Statistics/PDF.C
  \author    J. Bakosi
  \date      Tue May  7 12:16:51 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Univariate PDF estimator
  \details   Univariate PDF estimator
*/
//******************************************************************************

#include <iostream>

#include <PDF.h>

using namespace Quinoa;

PDF::PDF(const real& binsize) :
  m_binsize(binsize),
  m_pdf()
//******************************************************************************
//  Constructor: Initialize joint PDF container
//! \param[in]   binsize    Sample space bin size
//! \author J. Bakosi
//******************************************************************************
{
}

PDF::~PDF() noexcept
//******************************************************************************
//  Destructor: Clear joint PDF container
//! \author J. Bakosi
//******************************************************************************
{
  m_pdf.clear();
}

void
PDF::insert(const real& sample)
//******************************************************************************
//  Insert new sample into joint PDF
//! \param[in]   sample    Value to insert
//! \author  J. Bakosi
//******************************************************************************
{
  // Increase number of samples in joint PDF
  ++m_nsample;

  // Add sample to joint PDF
  ++m_pdf[floor(sample/m_binsize+0.5)];
}
