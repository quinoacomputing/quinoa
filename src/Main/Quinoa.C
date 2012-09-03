//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Sun 02 Sep 2012 06:41:04 PM MDT
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

  Memory store;

  MemoryEntry* e1 = store.newEntry(1000, INT_VAL, SCALAR_VAR, "scalars");
  MemoryEntry* e2 = store.newEntry(1000, REAL_VAL, VECTOR_VAR, "vectors");
  MemoryEntry* e3 = store.newEntry(1000000, REAL_VAL, TENSOR_VAR, "tensors");

  cout << store.getBytes() << endl;

  Int* e1ptr = store.getPtr<Int>(e1);
  Real* e2ptr = store.getPtr<Real>(e2);
  Real* e3ptr = store.getPtr<Real>(e3);
  cout << e1ptr << endl;
  cout << e2ptr << endl;
  cout << e3ptr << endl;

  store.freeEntry(e1);
  store.freeEntry(e2);
  store.freeEntry(e3);

  store.freeAllEntries();
}
