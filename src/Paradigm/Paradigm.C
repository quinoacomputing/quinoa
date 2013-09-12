//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.C
  \author    J. Bakosi
  \date      Wed 11 Sep 2013 10:22:02 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parallel programming paradigms
  \details   Parallel programming paradigms
*/
//******************************************************************************

#include <iostream>

#include <Paradigm.h>

using namespace quinoa;

void
Paradigm::echo() const
//******************************************************************************
//  Echo paradigm and configuration
//! \author  J. Bakosi
//******************************************************************************
{
  m_print.section("Compute environment");

  // OpenMP
  if (m_omp.available()) {
    m_print.item("OpenMP", "found");
    if (m_omp.used()) {
      m_print.item("Using threads", "yes");
      m_print.item("Number of threads", m_omp.nthread());
    } else {
      m_print.item("Using threads", "no");
    }
  } else {
    m_print.item("OpenMP", "not found");
  }
}
