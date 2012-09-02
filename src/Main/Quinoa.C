//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Sat 01 Sep 2012 11:46:43 PM MDT
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

  Memory Store;

  Int e1 = Store.newEntry(10, INT_VAL, SCALAR_VAR, "scalars");
  Int e2 = Store.newEntry(10000, BOOL_VAL, VECTOR_VAR, "vectors");
  Int e3 = Store.newEntry(11234, REAL_VAL, TENSOR_VAR, "tensors");

  cout << Store.getBytes() << endl;

  Int* e1ptr = Store.getPtr<Int>(e1);
  Bool* e2ptr = Store.getPtr<Bool>(e2);
  Real* e3ptr = Store.getPtr<Real>(e3);
  cout << e1ptr << endl;
  cout << e2ptr << endl;
  cout << e3ptr << endl;

  Store.freeAll();
}
