//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.C
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 06:21:56 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parallel programming paradigms
  \details   Parallel programming paradigms
*/
//******************************************************************************

#include <iostream>

#include <Paradigm.h>

using tk::Paradigm;

Paradigm::Paradigm(const Print& print)
//******************************************************************************
//  Constructor
//! \param[in] print     Simple pretty printer
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Compute environment");

  // OpenMP
  if (m_omp.available()) {
    print.item("OpenMP", "found");
    if (m_omp.used()) {
      print.item("Using threads", "yes");
      print.item("Number of threads", m_omp.nthreads());
    } else {
      print.item("Using threads", "no");
    }
  } else {
    print.item("OpenMP", "not found");
  }
}
