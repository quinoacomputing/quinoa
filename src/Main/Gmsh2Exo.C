//******************************************************************************
/*!
  \file      src/Utilities/Gmsh2Exo.C
  \author    J. Bakosi
  \date      Wed Mar 19 08:08:48 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh to Exodus II mesh file converter
  \details   Gmsh to Exodus II mesh file converter
*/
//******************************************************************************

#include <Init.h>
#include <InitGmsh2Exo.h>
#include <Config.h>
#include <Gmsh2ExoDriver.h>

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Gmsh to Exodus II mesh file converter
//! \author  J. Bakosi
//******************************************************************************
{
  return tk::Main< gmsh2exo::Gmsh2ExoDriver, gmsh2exo::echoTPL >
                 ( argc,
                   argv,
                   "Quinoa: Gmsh to Exodus II mesh file converter",
                   GMSH2EXO_EXECUTABLE );
}
