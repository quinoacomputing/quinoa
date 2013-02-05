//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.C
  \author    J. Bakosi
  \date      Mon 04 Feb 2013 08:55:40 PM MST
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
Paradigm::echo()
//******************************************************************************
//  Echo paradigm and configuration
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Parallel environment:"
        "\n---------------------" << endl;

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
