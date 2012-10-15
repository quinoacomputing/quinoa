//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Mon Oct 15 08:41:42 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <QuinoaTypes.h>
#include <Memory.h>
#include <Driver.h>
#include <GmshTxtMeshReader.h>
#include <GmshTxtMeshWriter.h>
#include <MKLRandom.h>
#include <MeshException.h>
#include <IOException.h>
#include <MKLException.h>

using namespace std;
using namespace Quinoa;

int main(int argc, char* argv[]) {

  // query number of threads available
  #ifdef _OPENMP
  Int nthreads = 1;//omp_get_max_threads();
  #else
  Int nthreads = 1;
  #endif
  cout << "* Using number of OpenMP threads: " << nthreads << endl;

  Memory memStore(nthreads);   // arg: nthreads
  Driver driver(&memStore);

  ErrCode error = NO_ERROR;
  try {

    UnsMesh mesh(&memStore);
    GmshTxtMeshReader inMesh("../../tmp/cylinder.msh", &mesh, &memStore);
    inMesh.read();
    GmshTxtMeshWriter outMesh("../../tmp/cylinder_out.msh", &mesh, &memStore);
    outMesh.write();

    MKLRandom random(nthreads,1,&memStore);      // nthreads, seed
    random.addTable(UNIFORM,1000000,"Gauss");
    random.addTable(UNIFORM,1000000,"Uniform");
    random.regenTables();

    memStore.echoAllEntries(MemoryEntryField::NAME);
    cout << "Allocated memory: " << memStore.getBytes() << endl;

  } // catch different exception types
    catch (MemoryException& e) { error = e.handleException(&driver); }
    catch (MeshException& e)   { error = e.handleException(&driver); }
    catch (IOException& e)     { error = e.handleException(&driver); }
    catch (MKLException& e)    { error = e.handleException(&driver); }
    // catch uncaught exceptions
    catch (...) {
      Exception e(UNCAUGHT);
      error = e.handleException(&driver);
    }

  if (error != FATAL_ERROR) {
    cout << "still running..." << endl;
  }
}
