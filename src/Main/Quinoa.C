//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 03:32:03 AM KST
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
#include <GmshReader.h>
#include <MeshException.h>

using namespace std;
using namespace Quinoa;

int main(int argc, char* argv[]) {

  Memory memStore(1);   // arg: nthreads
  Driver driver(&memStore);

  ErrorCode error = NO_ERROR;
  try {

    MemoryEntry* e = memStore.newEntry(10, INT, SCALAR, "scalars");

    UnsMesh mesh(&memStore);
    GmshReader gmsh("../../tmp/cylinder.msh", &mesh, &memStore);
    gmsh.read();

    memStore.echoAllEntries();
    memStore.freeEntry(e);

  } catch (MemoryException& e) { error = e.handleException(&driver); }
    catch (MeshException& e) { error = e.handleException(&driver); }

  if (error != FATAL_ERROR) {
    cout << "still running..." << endl;
  }
}
