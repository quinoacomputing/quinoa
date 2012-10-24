//******************************************************************************
/*!
  \file      src/Statistics/JPD.C
  \author    J. Bakosi
  \date      Wed Oct 24 08:11:24 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Joint PDF estimator
  \details   Joint PDF estimator
*/
//******************************************************************************

#include <JPDF.h>

using namespace Quinoa;

JPDF::JPDF(const real dim, const real binsize) :
      m_dim(dim), m_binsize(binsize), m_nsample(0)
//******************************************************************************
//  Constructor: Initialize joint PDF container
//! \param[in]   dim        Dimension of sample space
//! \param[in]   binsize    Sample space bin size
//! \author J. Bakosi
//******************************************************************************
{
  // Grow temporary key vector size in advance
  m_key.reserve(dim);
}

JPDF::~JPDF()
//******************************************************************************
//  Destructor: Clear joint PDF container
//! \author J. Bakosi
//******************************************************************************
{
  m_pdf.clear();
}

void
JPDF::insert(const vector<real>& value)
//******************************************************************************
//  Insert new value into joint PDF
//! \param[in]   value    Value to insert
//! \author  J. Bakosi
//******************************************************************************
{
  // Make sure sample has the same dimension as the joint PDF
  assert(value.size() == m_key.size());
  // Find bin ids in all dimensions
  transform(value.begin(), value.end(), m_key.begin(),
            [&](const int& v)->int { return floor(v/m_binsize+0.5); } );
  // Increase number of samples in joint PDF
  ++m_nsample;
  // Add sample to joint PDF
  ++m_pdf[m_key];
}
