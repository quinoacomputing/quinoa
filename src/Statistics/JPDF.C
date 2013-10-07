//******************************************************************************
/*!
  \file      src/Statistics/JPDF.C
  \author    J. Bakosi
  \date      Mon Oct  7 08:45:47 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Joint PDF estimator
  \details   Joint PDF estimator
*/
//******************************************************************************

#include <iostream>

#include <JPDF.h>
#include <Exception.h>

using namespace tk;

void
JPDF::insert(const std::vector<real>& sample)
//******************************************************************************
//  Insert new sample into joint PDF
//! \param[in]   sample    Sample to insert
//! \author  J. Bakosi
//******************************************************************************
{
  // Make sure sample has the same dimension as the joint PDF
  Assert(sample.size() == m_key.size(), ExceptType::FATAL,
         "Sample to be inserted and joint PDF have different dimensions");

  // Find bin ids in all dimensions
  std::transform(sample.begin(), sample.end(), m_key.begin(),
                 [&](const real& val)->int {
                   return static_cast<int>(floor(val/m_binsize+0.5)); } );

  // Increase number of samples in joint PDF
  ++m_nsample;

  // Add sample to joint PDF
  ++m_pdf[m_key];
}
