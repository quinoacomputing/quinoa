//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.C
  \author    J. Bakosi
  \date      Thu May 23 10:30:32 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parallel programming paradigms
  \details   Parallel programming paradigms
*/
//******************************************************************************

#include <iostream>

#include <Paradigm.h>

using namespace std;
using namespace Quinoa;

void
Paradigm::echo() const
//******************************************************************************
//  Echo paradigm and configuration
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Compute environment:"
        "\n--------------------" << endl;

  // OpenMP
  cout << " * OpenMP: ";
  if (m_omp.available()) {
    cout << "found";
    if (m_omp.used()) {
      cout << ", using ";
    } else {
      cout << ", not using ";
    }
    cout << m_omp.nthread() << " threads" << endl;
  } else {
    cout << "not found" << endl;
  }

  cout << endl;
}
