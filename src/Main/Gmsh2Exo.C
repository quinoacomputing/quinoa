//******************************************************************************
/*!
  \file      src/Utilities/Gmsh2Exo.C
  \author    J. Bakosi
  \date      Thu Mar 27 16:16:39 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh to Exodus II mesh file converter
  \details   Gmsh to Exodus II mesh file converter
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <Gmsh2ExoDriver.h>
#include <TPLInfo/ExodusII.h>

namespace gmsh2exo {

void echoTPL(const tk::Print& print)
//******************************************************************************
//  Echo TPL version informaion for libs specific to Gmsh2Exo
//! \author  J. Bakosi
//******************************************************************************
{
  print.raw("\n");
  echoExodusII(print, "ExodusII library");
}

} // gmsh2exo::

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Gmsh to Exodus II mesh file converter
//! \author  J. Bakosi
//******************************************************************************
{
  return tk::Main< gmsh2exo::Gmsh2ExoDriver >
                 ( argc,
                   argv,
                   "Quinoa: Gmsh to Exodus II mesh file converter",
                   GMSH2EXO_EXECUTABLE,
                   gmsh2exo::echoTPL );
}
