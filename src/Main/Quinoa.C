//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Wed 24 Oct 2012 05:37:48 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <mkl_vsl.h>

#include <QuinoaTypes.h>
#include <Memory.h>
#include <Driver.h>
#include <GmshTxtMeshReader.h>
#include <GmshTxtMeshWriter.h>
#include <MKLRandom.h>
#include <MeshException.h>
#include <IOException.h>
#include <MKLException.h>
#include <JPDF.h>
#include <PDFWriter.h>

using namespace std;
using namespace Quinoa;

int main(int argc, char* argv[]) {

  // query number of threads available
  #ifdef _OPENMP
  int nthread = 1;//omp_get_max_threads();
  #else
  int nthread = 1;
  #endif
  cout << "* Using number of OpenMP threads: " << nthread << endl;

  Memory memStore(nthread);   // arg: nthread
  Driver driver(&memStore);

  ErrCode error = NO_ERROR;
  try {

    // Memory
    //MemoryEntry* a = memStore.newEntry(10, INT, SCALAR, "_ja");

    // Mesh
    UnsMesh mesh(&memStore);
    GmshTxtMeshReader inMesh("../../tmp/cylinder.msh", &mesh, &memStore);
    inMesh.read();
    GmshTxtMeshWriter outMesh("../../tmp/cylinder_out.msh", &mesh, &memStore);
    outMesh.write();

    // Random
    int num = 100000;
    MKLRandom random(nthread, &memStore);
    MKLRndTable* tg = random.addTable(VSL_BRNG_MCG59, GAUSSIAN,
                                      VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                                      1, num, "Gaussian");
    MKLRndTable* tu = random.addTable(VSL_BRNG_MCG59, UNIFORM,
                                      VSL_RNG_METHOD_UNIFORM_STD,
                                      1, num, "Uniform");
    //MKLRndTable* t = random.addTable(VSL_BRNG_MCG59, GAMMA,
    //                                 VSL_RNG_METHOD_GAMMA_GNORM,
    //                                 1, num, "Gamma");
    random.regenTables();

    // PDF
    const real* rndg = random.getRnd(tg);
    const real* rndu = random.getRnd(tu);
    JPDF pdf(2,0.1);
    for (int i=0; i<num; ++i) {
      vector<real> v;
      v.push_back(rndg[i]);
      v.push_back(rndu[i]);
      v.push_back(rndu[i]);
      pdf.insert(v);
    }
    PDFWriter pw("pdf");
    //pw.write(&pdf);

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
