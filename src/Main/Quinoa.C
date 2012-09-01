//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Sat 01 Sep 2012 03:21:38 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>
#include <limits>

#include <QuinoaTypes.h>
#include <Memory.h>

using namespace std;
using namespace Quinoa;

int main(int argc, char* argv[]) {

  Memory<Int> intStore;
  Memory<Real> realStore;

  Int e1 = intStore.newEntry(10, SCALAR_VAR, "scalars");
  Int e2 = intStore.newEntry(10000, VECTOR_VAR, "vectors");
  Int e3 = realStore.newEntry(11234, TENSOR_VAR, "tensors");

  cout << intStore.getBytes() + realStore.getBytes() << endl;

  Int* e1ptr = intStore.getPtr(e1);
  cout << e1ptr << endl;

  intStore.freeAll();
  intStore.freeAll();
}
