//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Sat 06 Oct 2012 08:08:51 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#include <QuinoaTypes.h>
#include <Memory.h>
#include <Driver.h>
#include <GmshTxtMeshReader2D.h>
//#include <GmshTxtMeshWriter.h>
#include <MeshException.h>
#include <IOException.h>

using namespace std;
using namespace Quinoa;

int main(int argc, char* argv[]) {

  Memory memStore(1);   // arg: nthreads
  Driver driver(&memStore);

  ErrCode error = ErrCode::NO_ERROR;
  try {

    //MemoryEntry* e =
    //  memStore.newEntry(10, ValType::INT, VarType::SCALAR, "scalars");

    UnsMesh mesh(&memStore);
    GmshTxtMeshReader2D inMesh("../../tmp/cylinder_stripped.msh", &mesh, &memStore);
    inMesh.read();
    //GmshTxtMeshWriter outMesh("../../tmp/cylinder_out.msh", &mesh, &memStore);
    //outMesh.write();

    memStore.echoAllEntries(MemoryEntryField::NAME);
    cout << "Allocated memory: " << memStore.getBytes() << endl;
    //memStore.freeEntry(e);

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
