//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.C
  \author    J. Bakosi
  \date      Wed Sep 11 15:46:50 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parallel programming paradigms
  \details   Parallel programming paradigms
*/
//******************************************************************************

#include <iostream>

#include <Paradigm.h>

using namespace std;
using namespace quinoa;

void
Paradigm::echo(const QuinoaPrinter& print) const
//******************************************************************************
//  Echo paradigm and configuration
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Compute environment");

  // OpenMP
  if (m_omp.available()) {
    print.item("OpenMP", "found");
    if (m_omp.used()) {
      print.item("Using threads", "yes");
      print.item("Number of threads", m_omp.nthread());
    } else {
      print.item("Using threads", "no");
    }
  } else {
    print.item("OpenMP", "not found");
  }

  cout << endl;
}
