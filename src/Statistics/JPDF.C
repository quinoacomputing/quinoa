//******************************************************************************
/*!
  \file      src/Statistics/JPDF.C
  \author    J. Bakosi
  \date      Fri 03 May 2013 07:10:59 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Joint PDF estimator
  \details   Joint PDF estimator
*/
//******************************************************************************

#include <iostream>
#include <cassert>

#include <JPDF.h>

using namespace Quinoa;

JPDF::JPDF(const int dim, const real binsize) :
  m_binsize(binsize)
//******************************************************************************
//  Constructor: Initialize joint PDF container
//! \param[in]   dim        Dimension of sample space
//! \param[in]   binsize    Sample space bin size
//! \author J. Bakosi
//******************************************************************************
{
  // Grow temporary key vector size in advance
  m_key.resize(dim);
}

JPDF::~JPDF() noexcept
//******************************************************************************
//  Destructor: Clear joint PDF container
//! \author J. Bakosi
//******************************************************************************
{
  m_pdf.clear();
}

void
JPDF::insert(const vector<real>& sample)
//******************************************************************************
//  Insert new sample into joint PDF
//! \param[in]   sample    Sample to insert
//! \author  J. Bakosi
//******************************************************************************
{
  // Make sure sample has the same dimension as the joint PDF
  assert(sample.size() == m_key.size());

  // Find bin ids in all dimensions
  transform(sample.begin(), sample.end(), m_key.begin(),
            [&](const real& val)->int {
               return static_cast<int>(floor(val/m_binsize+0.5)); } );

  // Increase number of samples in joint PDF
  ++m_nsample;

  // Add sample to joint PDF
  ++m_pdf[m_key];
}
