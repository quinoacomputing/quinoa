//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Thu Sep  6 16:40:58 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>
#include <limits>

#include <QuinoaTypes.h>
#include <Memory.h>
#include <Driver.h>

using namespace std;
using namespace Quinoa;

int main(int argc, char* argv[]) {

  Memory memStore(1);
  Driver driver(&memStore);

  ErrorCode error = NO_ERROR;
  try {
    MemoryEntry* e1 = memStore.newEntry(100, INT_VAL, SCALAR_VAR, "scalars");
    MemoryEntry* e2 = memStore.newEntry(1000, REAL_VAL, VECTOR_VAR, "vectors");
    MemoryEntry* e3 = memStore.newEntry(1000000, REAL_VAL, TENSOR_VAR, "tensors");

    cout << "Memory usage: " << memStore.getBytes() << " Bytes" << endl;

    Int* e1ptr = memStore.getPtr<Int>(e1);
    Real* e2ptr = memStore.getPtr<Real>(e2);
    Real* e3ptr = memStore.getPtr<Real>(e3);
    cout << e1ptr << endl;
    cout << e2ptr << endl;
    cout << e3ptr << endl;

    cout << "---" << endl;
    cout << memStore.getPtr<Int>(memStore.getID("scalars")) << endl;

    memStore.freeEntry(e1);
    //memStore.freeEntry(e2);
    //memStore.freeEntry(e3);

  } catch (MemoryException& m) { error = m.handleException(&driver); }

  if (error != FATAL_ERROR) {
    cout << "still running..." << endl;
  }
}
