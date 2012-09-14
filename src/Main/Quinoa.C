//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Fri Sep 14 16:20:37 2012
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
#include <GmshMeshWriter.h>
#include <MeshException.h>
#include <IOException.h>

using namespace std;
using namespace Quinoa;

int main(int argc, char* argv[]) {

  Memory memStore(1);   // arg: nthreads
  Driver driver(&memStore);

  //cout << sizeof(ValType) << endl;

  ErrorCode error = NO_ERROR;
  try {

    MemoryEntry* e =
      memStore.newEntry(10, ValType::INT, VarType::SCALAR, "scalars");

    UnsMesh mesh(&memStore);
    GmshReader inMesh("../../tmp/cylinder.msh", &mesh, &memStore);
    inMesh.read();
    GmshMeshWriter outMesh("../../tmp/cylinder_out.msh", &mesh, &memStore);
    outMesh.write();

    memStore.echoAllEntries();
    cout << "Allocated memory: " << memStore.getBytes() << endl;
    memStore.freeEntry(e);

  } // catch different exception types
    catch (MemoryException& e) { error = e.handleException(&driver); }
    catch (MeshException& e)   { error = e.handleException(&driver); }
    catch (IOException& e)     { error = e.handleException(&driver); }
    // catch uncaught exceptions
    catch (...) { Exception e(UNCAUGHT); error = e.handleException(&driver); }

  if (error != FATAL_ERROR) {
    cout << "still running..." << endl;
  }
}
