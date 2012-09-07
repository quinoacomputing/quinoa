//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 04:18:15 PM MDT
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
#include <GmshException.h>

using namespace std;
using namespace Quinoa;

int main(int argc, char* argv[]) {

  Memory memStore(1);   // arg: nthreads
  Driver driver(&memStore);

  ErrorCode error = NO_ERROR;
  try {

    UnsMesh mesh;
    GmshReader gmsh("../../tmp/cylinder.msh");
    gmsh.open();
    gmsh.read(&mesh);
    gmsh.close();

  } catch (MemoryException& m) { error = m.handleException(&driver); }
    catch (GmshException& g) { error = g.handleException(&driver); }

  if (error != FATAL_ERROR) {
    cout << "still running..." << endl;
  }
}
