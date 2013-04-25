//******************************************************************************
/*!
  \file      src/Statistics/PDF.C
  \author    J. Bakosi
  \date      Wed 24 Apr 2013 11:22:15 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Univariate PDF estimator
  \details   Univariate PDF estimator
*/
//******************************************************************************

#include <iostream>

#include <PDF.h>
#include <StatException.h>

using namespace Quinoa;

PDF::PDF(const real& binsize) :
  m_binsize(binsize)
//******************************************************************************
//  Constructor: Initialize joint PDF container
//! \param[in]   binsize    Sample space bin size
//! \author J. Bakosi
//******************************************************************************
{
}

PDF::~PDF()
//******************************************************************************
//  Destructor: Clear joint PDF container
//! \author J. Bakosi
//******************************************************************************
{
  m_pdf.clear();
}

void
PDF::insert(const real& value)
//******************************************************************************
//  Insert new value into joint PDF
//! \param[in]   value    Value to insert
//! \author  J. Bakosi
//******************************************************************************
{
  // Increase number of samples in joint PDF
  ++m_nsample;

  // Add sample to joint PDF
  ++m_pdf[floor(value/m_binsize+0.5)];
}
