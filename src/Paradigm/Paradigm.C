//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.C
  \author    J. Bakosi
  \date      Sun 25 May 2014 06:05:49 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parallel programming paradigms
  \details   Parallel programming paradigms
*/
//******************************************************************************

#include <Paradigm.h>

using tk::Paradigm;

void
Paradigm::info( const Print& print )
//******************************************************************************
//  Output info on compute environment
//! \param[in] print     Simple pretty printer
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Compute environment");

  // OpenMP
  OpenMP omp;
  if (omp.available()) {
    print.item("OpenMP", "found");
    if (omp.used()) {
      print.item("Using OpenMP threads", "yes");
      print.item("Number of OpenMP threads", omp.nthreads());
    } else {
      print.item("Using OpenMP threads", "no");
    }
  } else {
    print.item("OpenMP", "not found");
  }
}
