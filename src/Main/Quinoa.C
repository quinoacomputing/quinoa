//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Sat 15 Sep 2012 02:18:39 PM MDT
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

  ErrCode error = ErrCode::NO_ERROR;
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
    catch (...) {
      Exception e(ExceptType::UNCAUGHT);
      error = e.handleException(&driver);
    }

  if (error != ErrCode::FATAL) {
    cout << "still running..." << endl;
  }
}
